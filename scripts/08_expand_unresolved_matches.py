#!/usr/bin/env python3
"""
Phase 8: Unresolved-entry expansion.

Phase 6 is the validated conservative baseline. Phase 8 is a SEPARATE expansion
layer that tries to recover additional biologically plausible links from the
~19k unresolved entries, while preserving reliability. It never overwrites
Phase 5/6/7 outputs; all Phase 8 results land in new files, plus merged
`*_v2.csv` files that union Phase 6 + Phase 8 recovered rows with provenance.

Recovery strategy (per unresolved row, in priority order):
  A. Pre-filter unrecoverable entries (Bucket E, numeric/no-alpha, sequence
     motifs) -> skipped_unrecoverable.
  F. Manual overrides (optional data/manual_name_overrides.csv) applied first.
  D. Cross-reference: id-like entities -> UniProt accession / xref columns and
     Phytozome gene IDs (sequence-confirmed is strongest).
  C. Description matching (Bucket D) via a TEMPORARY SQLite FTS5 index over
     UniProt descriptions, ranked by BM25 then re-scored with token overlap.
  E. Unknown species-prefix entities -> global UniProt gene-symbol lookup.

Realistic expectations (PI guidance):
  - Bucket E: generally unrecoverable; skipped outright.
  - Bucket D: main target; ~15-30% recovery via description matching.
  - Bucket C: limited; <10% via cross-reference fallback.
  - Bucket B unknown-prefix: ~10-20% via global symbol lookup.
  - Total realistic recovery: ~2,500-4,500 additional matches, mostly coverage.

Why FTS5 + token overlap: FTS5 gives fast BM25 candidate retrieval over millions
of UniProt descriptions without scanning the 21M-row index per entry. Token
overlap on top of BM25 is a transparent, tunable PRECISION filter so we can
separate strong (-> sometimes confidence) from partial (-> coverage) matches.

Why pre-filter Bucket E / generics: terms like "complex/domain/subunit/family"
with no specific name produce many spurious high-BM25 hits. Forcing them into
matches would inflate false positives and pollute the conservative baseline.

Plant-organism filter: this is a plant-protein dataset, but the FTS index spans
all of UniProt, so descriptive names readily match fungi/algae/dinoflagellates.
Fuzzy (description + unknown-prefix) matches are therefore restricted to an
allowlist of organism mnemonics already accepted by the validated Phase 6
datasets (plus the crosswalk). This bounds Phase 8 to baseline-vouched taxa
without needing external taxonomy. Disable with --no-plant-filter.

NOTE on the "<3 informative tokens" guidance: applied literally it would skip
legitimate 2-token names ("larreatricin hydroxylase") and, given the data
(only 5,129 of 14,528 D entries have >=3 informative tokens), make the
2,500-4,500 recovery target unreachable. The minimum is therefore exposed as
--min-informative-tokens (default 2). A 2-token name still needs BOTH specific
tokens present to reach the high-overlap threshold, so this stays conservative.

Standard library only. Does NOT train MTGNN, generate embeddings, or build edges.
Does NOT modify the read-only SQLite indexes or any Phase 1-7 output.
"""

import argparse
import csv
import os
import random
import re
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

# ── Paths ────────────────────────────────────────────────────────────────────
UNRESOLVED_CSV  = Path("outputs/unresolved_entries.csv")
RAW_CSV         = Path("outputs/raw_matches.csv")
P6_CONFIDENCE   = Path("outputs/plant_protein_dataset_confidence.csv")
P6_COVERAGE     = Path("outputs/plant_protein_dataset_coverage.csv")
CROSSWALK_CSV   = Path("data/species_crosswalk.csv")
UNIPROT_DB      = Path("indexes/uniprot_index.db")
PHYTOZOME_DB    = Path("indexes/phytozome_index.db")
MANUAL_OVERRIDE = Path("data/manual_name_overrides.csv")

OUT_EXPANDED    = Path("outputs/phase8_expanded_matches.csv")
OUT_REC_CONF    = Path("outputs/phase8_recovered_confidence.csv")
OUT_REC_COV     = Path("outputs/phase8_recovered_coverage.csv")
OUT_STILL       = Path("outputs/phase8_still_unresolved.csv")
OUT_AUDIT       = Path("outputs/phase8_audit_sample.csv")
OUT_V2_CONF     = Path("outputs/final_confidence_dataset_v2.csv")
OUT_V2_COV      = Path("outputs/final_coverage_dataset_v2.csv")
OUT_SUMMARY     = Path("outputs/phase8_summary.md")

ALL_OUTPUTS = [OUT_EXPANDED, OUT_REC_CONF, OUT_REC_COV, OUT_STILL,
               OUT_AUDIT, OUT_V2_CONF, OUT_V2_COV, OUT_SUMMARY]

# ── Vocabulary ───────────────────────────────────────────────────────────────
# Stopwords removed for description token-overlap scoring.
STOPWORDS = {
    'protein', 'gene', 'enzyme', 'type', 'like', 'related', 'putative',
    'probable', 'predicted', 'hypothetical', 'family', 'domain',
}
# UniProt description label tokens (RecName: Full=...; AltName: ...; etc.).
DESC_LABELS = {
    'recname', 'altname', 'subname', 'full', 'short', 'ec', 'flags', 'includes',
    'contains', 'allergen', 'biotech', 'cdantigen', 'innname',
}
# Generic terms that, alone, do not constitute a specific name.
GENERIC_TERMS = {
    'complex', 'family', 'domain', 'transport', 'chain', 'chains', 'subunit',
    'component', 'fragment', 'hypothetical', 'unknown', 'unnamed', 'predicted',
    'putative', 'probable', 'isoform', 'partial',
}
GREEK = {
    'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta', 'ε': 'epsilon',
    'ζ': 'zeta', 'η': 'eta', 'θ': 'theta', 'ι': 'iota', 'κ': 'kappa',
    'λ': 'lambda', 'μ': 'mu', 'ν': 'nu', 'ξ': 'xi', 'ο': 'omicron',
    'π': 'pi', 'ρ': 'rho', 'σ': 'sigma', 'ς': 'sigma', 'τ': 'tau',
    'υ': 'upsilon', 'φ': 'phi', 'χ': 'chi', 'ψ': 'psi', 'ω': 'omega',
}
# Explicit abbreviation cleanups applied before generic punctuation handling.
ABBREV = [
    ('(3R)', ' 3R '), ('(3S)', ' 3S '), ('(R)', ' R '), ('(S)', ' S '),
    ('NAD+', ' NAD '), ('NADP+', ' NADP '), ('Ca2+', ' Ca '), ('Mg2+', ' Mg '),
    ('Cu/Zn', ' Cu Zn '), ('alpha/beta', ' alpha beta '),
]
AA_LETTERS = set('ACDEFGHIKLMNPQRSTVWY')

# ── Quality-flag vocabulary (advisory only; never removes rows) ──────────────
VIRUS_TERMS = ['tobacco etch virus', 'mosaic virus', 'potyvirus',
               'virus', 'viral', 'viroid', 'phage']
BROAD_DESC_TERMS = {
    'receptor', 'receptors', 'transcription', 'factor', 'factors',
    'transporter', 'transporters', 'subunit', 'subunits', 'family',
    'domain', 'complex',
}
NON_LAND_PLANT = {'9DINO', '9EUKA'}

# Accession / cross-reference id patterns (uppercased input).
RE_UNIPROT_ACC = re.compile(
    r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$')
RE_REFSEQ      = re.compile(r'^(?:NP|XP|YP|WP|NM|XM|AP|ZP)_\d+(?:\.\d+)?$')
RE_TAIR        = re.compile(r'^AT[1-5MC]G\d{5}(?:\.\d+)?$')
RE_ENSEMBL     = re.compile(r'^[A-Z]{2,}\d{6,}(?:\.\d+)?$')
RE_PHYTO_LIKE  = re.compile(r'^[A-Z][a-z]+[A-Za-z]*\.?\d.*[A-Za-z].*$')


# ── Normalization helpers ────────────────────────────────────────────────────

def normalize_id(s):
    """UniProt/Phytozome lookup key convention (Phase 5): strip + uppercase."""
    return (s or '').strip().upper()


def normalize_text(s):
    """
    Lowercase, expand Greek/abbreviations, drop bracketed organism notes and
    punctuation while preserving alphanumeric IDs, collapse whitespace.
    Returns (normalized_string, informative_token_list).
    informative tokens = non-stopword, non-label, alphabetic-bearing tokens.
    """
    s = s or ''
    for g, name in GREEK.items():
        s = s.replace(g, ' ' + name + ' ')
    for a, b in ABBREV:
        s = s.replace(a, b)
    s = s.lower()
    s = re.sub(r'\[[^\]]*\]', ' ', s)       # bracketed clutter / organism notes
    s = re.sub(r'\([^)]*\)', ' ', s)        # parenthetical clutter
    s = s.replace('-', ' ').replace('_', ' ').replace('/', ' ')
    s = re.sub(r'[^a-z0-9 ]+', ' ', s)
    s = re.sub(r'\s+', ' ', s).strip()
    tokens = s.split()
    informative = [
        t for t in tokens
        if t not in STOPWORDS and t not in DESC_LABELS
        and not t.isdigit() and any(c.isalpha() for c in t)
    ]
    return s, informative


def informative_token_set(s):
    _, toks = normalize_text(s)
    return set(toks), toks


def description_tokens(desc):
    """Normalized informative token set for a UniProt description string."""
    toks, _ = informative_token_set(desc)
    return toks


def is_id_like(entity):
    n = normalize_id(entity)
    if not n or ' ' in n:
        return False
    return bool(
        RE_UNIPROT_ACC.match(n) or RE_REFSEQ.match(n) or RE_TAIR.match(n)
        or RE_ENSEMBL.match(n) or RE_PHYTO_LIKE.match(entity.strip())
    )


def looks_like_sequence_motif(entity):
    e = entity.strip()
    if ' ' in e or len(e) < 8 or any(c.isdigit() for c in e):
        return False
    return all(c.upper() in AA_LETTERS for c in e)


def quality_flags(entity, organism, candidate_summary, is_description_row):
    """
    Advisory quality assessment for a Phase 8 / merged row. Returns
    (phase8_quality_flag, phase8_quality_notes). Never alters matching.
      - broad_taxon_warning      : organism mnemonic begins with '9' (genus-level)
      - non_land_plant_warning   : organism is 9DINO / 9EUKA
      - virus_term_warning       : virus/phage term in entity or candidates
      - generic_description_warning : description-row entity dominated by broad terms
      - multiple_warnings        : >1 of the above (all listed in notes)
    """
    warns = []
    org = (organism or '').strip().upper()
    if org.startswith('9'):
        warns.append('broad_taxon_warning')
    if org in NON_LAND_PLANT:
        warns.append('non_land_plant_warning')
    text = f"{entity or ''} {candidate_summary or ''}".lower()
    if any(t in text for t in VIRUS_TERMS):
        warns.append('virus_term_warning')
    if is_description_row:
        norm, _ = normalize_text(entity)
        if (set(norm.split()) & BROAD_DESC_TERMS) or ('protein family' in (entity or '').lower()):
            warns.append('generic_description_warning')
    if not warns:
        return 'ok', ''
    if len(warns) > 1:
        return 'multiple_warnings', '; '.join(warns)
    return warns[0], warns[0]


# ── Pre-filter (Step A) ──────────────────────────────────────────────────────

def hard_skip_reason(row):
    """Global hard skips that apply regardless of matching path. Returns reason or None."""
    entity = row['entity']
    if row['bucket'] == 'E':
        return 'bucket_E_noise'
    if not any(c.isalpha() for c in entity):
        return 'no_alphabetic_characters'
    if entity.strip().isdigit():
        return 'purely_numeric'
    if looks_like_sequence_motif(entity):
        return 'sequence_motif'
    return None


def description_prefilter_reason(entity, min_tokens):
    """
    Description-specific reject (Step A nuance): reject only if no SPECIFIC
    informative token remains after normalization, or too few informative tokens.
    """
    tok_set, _ = informative_token_set(entity)
    if len(tok_set) < min_tokens:
        return 'too_few_informative_tokens'
    specific = tok_set - GENERIC_TERMS
    if not specific:
        return 'generic_only'
    return None


# ── Crosswalk ────────────────────────────────────────────────────────────────

def load_crosswalk(path):
    """Returns sets/maps for validating expected codes/mnemonics."""
    codes, mnemonics = set(), set()
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r['phytozome_code'].strip():
                codes.add(r['phytozome_code'].strip())
            if r['uniprot_mnemonic'].strip():
                mnemonics.add(r['uniprot_mnemonic'].strip())
    return codes, mnemonics


def build_plant_allowlist(crosswalk_mn):
    """
    Plant-organism allowlist for fuzzy matches = every UniProt organism mnemonic
    already accepted by the validated Phase 6 datasets, plus the crosswalk
    mnemonics. The UniProt index has no taxonomy column and the plant-database
    cross-ref fields (ensembl_plants/gramene/tair) proved unreliable (a fungal
    organism carried ~29k such xrefs while a real plant carried 1). Bounding
    Phase 8 to baseline-accepted taxa is the most defensible automatic filter:
    Phase 8 cannot introduce organisms the conservative baseline never vouched
    for. Residual non-land-plant codes (e.g. 9DINO) remain only because Phase 6
    already included them; a stricter Viridiplantae filter needs external
    taxonomy that is not available here.
    """
    allow = set(crosswalk_mn)
    for path in (P6_CONFIDENCE, P6_COVERAGE):
        with open(path) as fh:
            for r in csv.DictReader(fh):
                org = r.get('uniprot_organism', '').strip()
                if org:
                    allow.add(org)
    return allow


# ── SQLite query helpers (indexed columns only, parameterized) ───────────────

def _uniprot_cands(cur, where_sql, params, limit=10):
    """(sprot-preferred candidate rows, total count) for an indexed WHERE clause."""
    cur.execute(f"SELECT COUNT(*) FROM uniprot_entries WHERE {where_sql}", params)
    count = cur.fetchone()[0]
    if count == 0:
        return [], 0
    cur.execute(
        f"SELECT accession, entry_name, organism_mnemonic, source_tier, "
        f"description, gene_name FROM uniprot_entries "
        f"WHERE {where_sql} AND source_tier='sprot' LIMIT {limit}", params)
    cands = cur.fetchall()
    if not cands:
        cur.execute(
            f"SELECT accession, entry_name, organism_mnemonic, source_tier, "
            f"description, gene_name FROM uniprot_entries "
            f"WHERE {where_sql} LIMIT {limit}", params)
        cands = cur.fetchall()
    return cands, count


def query_uniprot_accession(cur, accession, limit=10):
    return _uniprot_cands(cur, "accession=?", (normalize_id(accession),), limit)


def query_uniprot_gene(cur, norm_gene, mnemonic=None, limit=10):
    if mnemonic:
        return _uniprot_cands(
            cur, "gene_name_normalized=? AND organism_mnemonic=?",
            (norm_gene, mnemonic), limit)
    return _uniprot_cands(cur, "gene_name_normalized=?", (norm_gene,), limit)


def query_phytozome_global(cur, norm_id, limit=5):
    """Match a gene id against any species (indexed gene_id / base_gene_id)."""
    cur.execute(
        "SELECT phytozome_code, gene_id, base_gene_id, sequence_length "
        "FROM phytozome_sequences "
        "WHERE gene_id_normalized=? OR base_gene_id_normalized=? LIMIT ?",
        (norm_id, norm_id, limit))
    return cur.fetchall()


def dedup_uniprot(cands):
    seen, out = set(), []
    for c in cands:
        key = (c[0], c[1], c[2], c[3])
        if key not in seen:
            seen.add(key)
            out.append(c)
    return out


def fmt_uniprot_summary(cands, n=5):
    return '; '.join(f"{c[0]}|{c[1]}|{c[2]}|{c[3]}" for c in dedup_uniprot(cands)[:n])


# ── FTS5 description index (Step C) ──────────────────────────────────────────

def check_fts5_available():
    try:
        con = sqlite3.connect(":memory:")
        con.execute("CREATE VIRTUAL TABLE _t USING fts5(x)")
        con.close()
        return True
    except sqlite3.OperationalError:
        return False


def build_fts_index(fts_path, uniprot_db, reuse):
    """
    Build a TEMPORARY on-disk FTS5 index over UniProt descriptions, once.
    Stores metadata as UNINDEXED columns so BM25 hits map back to accessions.
    Returns an open connection to the FTS database (read path).
    The source index is attached read-only and never modified.
    """
    if reuse and fts_path.exists():
        con = sqlite3.connect(str(fts_path), uri=True)
        row = con.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='_fts_done'"
        ).fetchone()
        if row:
            n = con.execute("SELECT COUNT(*) FROM uniprot_fts").fetchone()[0]
            print(f"  Reusing cached FTS5 index at {fts_path} ({n:,} rows).")
            return con
        con.close()
        fts_path.unlink()

    print(f"  Building FTS5 description index at {fts_path} ...")
    con = sqlite3.connect(str(fts_path), uri=True)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA cache_size=-200000")
    con.execute("CREATE VIRTUAL TABLE uniprot_fts USING fts5("
                "description, accession UNINDEXED, entry_name UNINDEXED, "
                "organism_mnemonic UNINDEXED, source_tier UNINDEXED, "
                "gene_name UNINDEXED)")
    con.execute("ATTACH DATABASE ? AS uni", (f"file:{uniprot_db}?mode=ro",))
    max_id = con.execute("SELECT MAX(id) FROM uni.uniprot_entries").fetchone()[0] or 0
    chunk = 1_000_000
    inserted = 0
    for lo in range(1, max_id + 1, chunk):
        hi = lo + chunk - 1
        con.execute(
            "INSERT INTO uniprot_fts(description, accession, entry_name, "
            "organism_mnemonic, source_tier, gene_name) "
            "SELECT COALESCE(description,'') || ' ' || COALESCE(gene_name,''), "
            "accession, entry_name, organism_mnemonic, source_tier, gene_name "
            "FROM uni.uniprot_entries WHERE id BETWEEN ? AND ?", (lo, hi))
        con.commit()
        inserted = min(hi, max_id)
        print(f"    FTS build: {inserted:,}/{max_id:,} rows")
    con.execute("CREATE TABLE _fts_done(ok INTEGER)")
    con.execute("INSERT INTO _fts_done VALUES (1)")
    con.commit()
    con.execute("DETACH DATABASE uni")
    print(f"  FTS5 index built: {inserted:,} rows.")
    return con


def fts_query(fts_cur, tokens, top_k):
    """BM25 top-k candidates for OR-joined literal tokens."""
    if not tokens:
        return []
    match = ' OR '.join('"' + t.replace('"', '') + '"' for t in tokens[:12])
    try:
        fts_cur.execute(
            "SELECT description, accession, entry_name, organism_mnemonic, "
            "source_tier, gene_name FROM uniprot_fts WHERE uniprot_fts MATCH ? "
            "ORDER BY bm25(uniprot_fts) LIMIT ?", (match, top_k))
        return fts_cur.fetchall()
    except sqlite3.OperationalError:
        return []


# ── Match record ─────────────────────────────────────────────────────────────

def new_record(urow, raw):
    return {
        'entity': urow['entity'],
        'bucket': urow['bucket'],
        'sub_type': urow['sub_type'],
        'unresolved_reason': urow['unresolved_reason'],
        'species_resolution_status': urow['species_resolution_status'],
        'species_match_status': urow['species_match_status'],
        'expected_phytozome_code': raw.get('expected_phytozome_code', ''),
        'expected_uniprot_mnemonic': raw.get('expected_uniprot_mnemonic', ''),
        'phase8_match_status': 'no_match',
        'phase8_match_type': '',
        'phase8_match_category': '',
        'tier': '',                      # 'confidence' | 'coverage'
        'phase8_confidence_label': '',
        'phase8_confidence_score': '',
        'description_overlap_score': '',
        'u_acc': '', 'u_entry': '', 'u_org': '', 'u_source': '',
        'p_code': '', 'p_gene': '', 'p_base': '', 'seq_len': '',
        'has_sequence': 'false',
        'candidate_count': 0,
        'candidate_summary': '',
        'notes': '',
        'best_method': 'none',
    }


# Category -> (tier, label, score)
CATEGORY_TIER = {
    'phase8_sequence_confirmed_match':        ('confidence', 'high', 95),
    'manual_override_sequence_match':         ('confidence', 'high', 95),
    'manual_override_species_confirmed_match':('confidence', 'high', 88),
    'description_species_confirmed_match':    ('confidence', 'high', 85),
    'phase8_crossref_species_confirmed_match':('confidence', 'high', 85),
    'description_single_species_match':       ('coverage',   'medium', 60),
    'phase8_crossref_uniprot_match':          ('coverage',   'medium', 65),
    'unknown_prefix_single_species_candidate':('coverage',   'medium', 50),
    'manual_override_global_match':           ('coverage',   'medium', 55),
    'description_multi_species_ambiguous':    ('coverage',   'low', 45),
    'description_partial_match':              ('coverage',   'low', 40),
    'unknown_prefix_multi_species_ambiguous': ('coverage',   'low', 40),
}


def set_category(rec, category, method, overlap=''):
    tier, label, score = CATEGORY_TIER[category]
    rec['phase8_match_status'] = 'matched'
    rec['phase8_match_type'] = method
    rec['phase8_match_category'] = category
    rec['tier'] = tier
    rec['phase8_confidence_label'] = label
    rec['phase8_confidence_score'] = score
    rec['best_method'] = method
    if overlap != '':
        rec['description_overlap_score'] = round(overlap, 4)


def fill_uniprot(rec, cand, count, summary):
    rec['u_acc'], rec['u_entry'], rec['u_org'], rec['u_source'] = \
        cand[0], cand[1], cand[2], cand[3]
    rec['candidate_count'] = count
    rec['candidate_summary'] = summary


def fill_phytozome(rec, prow):
    rec['p_code'], rec['p_gene'], rec['p_base'], rec['seq_len'] = \
        prow[0], prow[1], prow[2], prow[3]
    rec['has_sequence'] = 'true'
    rec['candidate_count'] = max(rec['candidate_count'], 1)
    if not rec['candidate_summary']:
        rec['candidate_summary'] = f"{prow[0]}|{prow[1]}|len={prow[3]}"


# ── Matching paths ───────────────────────────────────────────────────────────

def try_crossref(rec, urow, raw, ucur, pcur, xref_map):
    """Step D: id-like entity -> Phytozome (sequence) / UniProt accession / xref."""
    entity = urow['entity']
    sym = raw.get('extracted_gene_symbol', '').strip()
    expected_mn = raw.get('expected_uniprot_mnemonic', '').strip()
    candidates_ids = [entity]
    if sym and is_id_like(sym):
        candidates_ids.append(sym)
    if not any(is_id_like(c) for c in candidates_ids):
        return False
    rec['best_method'] = 'crossref'

    # Phytozome sequence (strongest evidence).
    for cid in candidates_ids:
        prows = query_phytozome_global(pcur, normalize_id(cid))
        if prows:
            fill_phytozome(rec, prows[0])
            rec['notes'] = 'Phytozome gene-id cross-reference'
            set_category(rec, 'phase8_sequence_confirmed_match', 'crossref')
            return True

    # UniProt accession (indexed).
    for cid in candidates_ids:
        cands, count = query_uniprot_accession(ucur, cid)
        if cands:
            return _finish_uniprot_crossref(rec, cands, count, expected_mn,
                                            'UniProt accession cross-reference')

    # UniProt xref columns (refseq/ensembl/tair/gramene) via prebuilt map.
    for cid in candidates_ids:
        hit = xref_map.get(normalize_id(cid))
        if hit:
            cands, count = _uniprot_cands(ucur, "accession=?", (hit,))
            if cands:
                return _finish_uniprot_crossref(
                    rec, cands, count, expected_mn,
                    'UniProt cross-reference id (refseq/ensembl/tair/gramene)')
    return False


def _finish_uniprot_crossref(rec, cands, count, expected_mn, note):
    cands = dedup_uniprot(cands)
    summary = fmt_uniprot_summary(cands)
    if expected_mn:
        for c in cands:
            if c[2] == expected_mn:
                fill_uniprot(rec, c, count, summary)
                rec['notes'] = note + ' (species-confirmed)'
                set_category(rec, 'phase8_crossref_species_confirmed_match', 'crossref')
                return True
    fill_uniprot(rec, cands[0], count, summary)
    rec['notes'] = note
    set_category(rec, 'phase8_crossref_uniprot_match', 'crossref')
    return True


def try_description(rec, urow, raw, fts_cur, args, plant_allow):
    """Step C: FTS5 BM25 retrieval + token-overlap scoring (Bucket D)."""
    entity = urow['entity']
    rec['best_method'] = 'description'
    skip = description_prefilter_reason(entity, args.min_informative_tokens)
    if skip:
        rec['notes'] = f'description prefilter: {skip}'
        return False

    q_set, q_tokens = informative_token_set(entity)
    cand_rows = fts_query(fts_cur, q_tokens, args.description_top_k)
    if not cand_rows:
        rec['notes'] = 'no FTS candidates'
        return False

    scored = []
    seen = set()
    for desc, acc, entry, org, source, gene in cand_rows:
        key = (acc, entry, org, source)
        if key in seen:
            continue
        seen.add(key)
        dtoks = description_tokens(desc)
        if not dtoks:
            continue
        inter = q_set & dtoks
        overlap = len(inter) / len(q_set) if q_set else 0.0
        cand_overlap = len(inter) / len(dtoks)
        scored.append((overlap, cand_overlap, acc, entry, org, source))

    kept = [s for s in scored if s[0] >= args.description_min_partial_score]
    # Plant-organism filter: resolve/ambiguity-break among plant candidates only.
    if plant_allow is not None:
        plant_kept = [s for s in kept if s[4] in plant_allow]
        if not plant_kept and kept:
            rec['notes'] = 'description match rejected: no plant-organism candidate'
            return False
        kept = plant_kept
    if not kept:
        rec['notes'] = 'description_no_good_match'
        return False
    kept.sort(key=lambda s: (-s[0], -s[1]))
    best = kept[0]
    best_overlap = best[0]
    expected_mn = raw.get('expected_uniprot_mnemonic', '').strip()
    summary = '; '.join(f"{s[2]}|{s[3]}|{s[4]}|{s[5]}(ov={s[0]:.2f})" for s in kept[:5])
    high = [s for s in kept if s[0] >= args.description_high_score]
    high_orgs = {s[4] for s in high}

    if high and expected_mn:
        for s in high:
            if s[4] == expected_mn:
                fill_uniprot(rec, (s[2], s[3], s[4], s[5]), len(kept), summary)
                rec['notes'] = 'description match, species-confirmed via expected mnemonic'
                set_category(rec, 'description_species_confirmed_match',
                             'description', best_overlap)
                return True

    if high and len(high_orgs) == 1:
        s = high[0]
        fill_uniprot(rec, (s[2], s[3], s[4], s[5]), len(kept), summary)
        rec['notes'] = 'strong description match, single organism'
        set_category(rec, 'description_single_species_match', 'description', best_overlap)
        return True
    if high and len(high_orgs) > 1:
        s = high[0]
        fill_uniprot(rec, (s[2], s[3], s[4], s[5]), len(kept), summary)
        rec['notes'] = f'strong description match, {len(high_orgs)} organisms (ambiguous)'
        set_category(rec, 'description_multi_species_ambiguous', 'description', best_overlap)
        return True

    # partial: min_partial <= best < high
    fill_uniprot(rec, (best[2], best[3], best[4], best[5]), len(kept), summary)
    rec['notes'] = 'partial description match (coverage)'
    set_category(rec, 'description_partial_match', 'description', best_overlap)
    return True


def try_unknown_prefix(rec, urow, raw, ucur, plant_allow):
    """Step E: global UniProt gene-symbol lookup for unknown species prefixes."""
    sym = raw.get('extracted_gene_symbol', '').strip()
    if not sym:
        return False
    rec['best_method'] = 'unknown_prefix'
    norm = normalize_id(sym)
    cands, count = query_uniprot_gene(ucur, norm)
    if not cands:
        rec['notes'] = 'no global symbol match'
        return False
    cands = dedup_uniprot(cands)
    if plant_allow is not None:
        plant_cands = [c for c in cands if c[2] in plant_allow]
        if not plant_cands:
            rec['notes'] = 'global symbol match rejected: no plant-organism candidate'
            return False
        cands = plant_cands
        count = len(cands)
    orgs = {c[2] for c in cands}
    summary = fmt_uniprot_summary(cands)
    fill_uniprot(rec, cands[0], count, summary)
    if len(orgs) == 1:
        rec['notes'] = f'global symbol {sym}: single organism {cands[0][2]}'
        set_category(rec, 'unknown_prefix_single_species_candidate', 'unknown_prefix')
    else:
        rec['notes'] = f'global symbol {sym}: {len(orgs)} organisms (ambiguous)'
        set_category(rec, 'unknown_prefix_multi_species_ambiguous', 'unknown_prefix')
    return True


def try_manual_override(rec, urow, ov, ucur, pcur):
    """Step F: apply a validated manual override row."""
    rec['best_method'] = 'manual_override'
    sym = ov.get('target_gene_symbol', '').strip()
    exp_mn = ov.get('expected_uniprot_mnemonic', '').strip()
    exp_code = ov.get('expected_phytozome_code', '').strip()
    note = ov.get('notes', '').strip()

    # Phytozome sequence (filtered by code if provided).
    for key in (sym, urow['entity']):
        if not key:
            continue
        prows = query_phytozome_global(pcur, normalize_id(key))
        if exp_code:
            prows = [p for p in prows if p[0] == exp_code]
        if prows:
            fill_phytozome(rec, prows[0])
            rec['notes'] = ('manual override: ' + note) if note else 'manual override'
            set_category(rec, 'manual_override_sequence_match', 'manual_override')
            return True

    if sym:
        cands, count = query_uniprot_gene(ucur, normalize_id(sym), exp_mn or None)
        if cands:
            cands = dedup_uniprot(cands)
            summary = fmt_uniprot_summary(cands)
            fill_uniprot(rec, cands[0], count, summary)
            rec['notes'] = ('manual override: ' + note) if note else 'manual override'
            if exp_mn:
                set_category(rec, 'manual_override_species_confirmed_match', 'manual_override')
            else:
                set_category(rec, 'manual_override_global_match', 'manual_override')
            return True
    return False


def load_manual_overrides(crosswalk_codes, crosswalk_mn, unresolved_entities):
    """Returns (overrides_by_entity, warnings list)."""
    overrides, warnings = {}, []
    if not MANUAL_OVERRIDE.exists():
        return overrides, warnings
    required = {'entity', 'target_gene_symbol', 'expected_uniprot_mnemonic',
                'expected_phytozome_code', 'notes'}
    with open(MANUAL_OVERRIDE) as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            warnings.append(f"manual override file missing required columns "
                            f"{sorted(required)}; file ignored")
            return overrides, warnings
        for i, r in enumerate(reader, start=2):
            entity = (r.get('entity') or '').strip()
            if not entity:
                warnings.append(f"row {i}: empty entity; skipped")
                continue
            if entity not in unresolved_entities:
                warnings.append(f"row {i}: entity {entity!r} not in unresolved_entries; skipped")
                continue
            mn = (r.get('expected_uniprot_mnemonic') or '').strip()
            if mn and mn not in crosswalk_mn:
                warnings.append(f"row {i}: expected_uniprot_mnemonic {mn!r} not in crosswalk; skipped")
                continue
            code = (r.get('expected_phytozome_code') or '').strip()
            if code and code not in crosswalk_codes:
                warnings.append(f"row {i}: expected_phytozome_code {code!r} not in crosswalk; skipped")
                continue
            overrides[entity] = r
    return overrides, warnings


# ── Cross-reference map build (single targeted scan) ─────────────────────────

def build_xref_map(ucur, candidate_ids):
    """
    One pass over UniProt building {xref_value_upper -> accession} restricted to
    values present in candidate_ids. Single scan (not per-entry); membership test
    against a set keeps memory bounded to actual matches.
    """
    if not candidate_ids:
        return {}
    print(f"  Scanning UniProt for {len(candidate_ids):,} candidate cross-ref ids ...")
    xref_map = {}
    ucur.execute("SELECT accession, refseq, ensembl_plants, tair, gramene "
                 "FROM uniprot_entries")
    n = 0
    for acc, refseq, ens, tair, gram in ucur:
        n += 1
        for field in (refseq, ens, tair, gram):
            if not field:
                continue
            for tok in re.split(r'[;, ]+', field):
                u = tok.strip().upper()
                if u and u in candidate_ids and u not in xref_map:
                    xref_map[u] = acc
        if n % 5_000_000 == 0:
            print(f"    xref scan: {n:,} rows, {len(xref_map):,} hits")
    print(f"  Cross-ref scan complete: {len(xref_map):,} id->accession entries.")
    return xref_map


# ── Output row builders ──────────────────────────────────────────────────────

DATASET_FIELDS = [
    'dataset_tier', 'entity', 'bucket', 'sub_type',
    'final_match_category', 'final_confidence_label', 'final_confidence_score',
    'uniprot_accession', 'uniprot_entry_name', 'uniprot_organism', 'uniprot_source_tier',
    'phytozome_code', 'phytozome_gene_id', 'phytozome_base_gene_id',
    'sequence_length', 'has_sequence',
    'expected_phytozome_code', 'expected_uniprot_mnemonic',
    'original_match_status', 'original_match_type', 'original_confidence_tier',
    'ambiguity_count', 'matched_species_category', 'species_match_status',
    'species_resolution_status', 'candidate_summary', 'notes',
]
PHASE8_EXTRA = ['source_phase', 'recovered_by_phase8', 'phase8_match_type',
                'phase8_match_category', 'description_overlap_score', 'phase8_notes',
                'phase8_quality_flag', 'phase8_quality_notes']
V2_FIELDS = DATASET_FIELDS + PHASE8_EXTRA

EXPANDED_FIELDS = [
    'entity', 'original_bucket', 'original_sub_type', 'original_unresolved_reason',
    'phase8_match_status', 'phase8_match_type', 'phase8_match_category',
    'phase8_confidence_label', 'phase8_confidence_score', 'description_overlap_score',
    'selected_uniprot_accession', 'selected_uniprot_entry_name',
    'selected_uniprot_organism', 'selected_uniprot_source_tier',
    'selected_phytozome_code', 'selected_phytozome_gene_id',
    'selected_phytozome_base_gene_id', 'selected_sequence_length',
    'has_sequence', 'candidate_count', 'candidate_summary', 'notes',
    'phase8_quality_flag', 'phase8_quality_notes',
]
STILL_FIELDS = None  # set from unresolved header + extras


def expanded_row(rec):
    qflag, qnotes = quality_flags(
        rec['entity'], rec['u_org'], rec['candidate_summary'],
        rec['phase8_match_type'] == 'description')
    return {
        'entity': rec['entity'],
        'original_bucket': rec['bucket'],
        'original_sub_type': rec['sub_type'],
        'original_unresolved_reason': rec['unresolved_reason'],
        'phase8_match_status': rec['phase8_match_status'],
        'phase8_match_type': rec['phase8_match_type'],
        'phase8_match_category': rec['phase8_match_category'],
        'phase8_confidence_label': rec['phase8_confidence_label'],
        'phase8_confidence_score': rec['phase8_confidence_score'],
        'description_overlap_score': rec['description_overlap_score'],
        'selected_uniprot_accession': rec['u_acc'],
        'selected_uniprot_entry_name': rec['u_entry'],
        'selected_uniprot_organism': rec['u_org'],
        'selected_uniprot_source_tier': rec['u_source'],
        'selected_phytozome_code': rec['p_code'],
        'selected_phytozome_gene_id': rec['p_gene'],
        'selected_phytozome_base_gene_id': rec['p_base'],
        'selected_sequence_length': rec['seq_len'],
        'has_sequence': rec['has_sequence'],
        'candidate_count': rec['candidate_count'],
        'candidate_summary': rec['candidate_summary'],
        'notes': rec['notes'],
        'phase8_quality_flag': qflag,
        'phase8_quality_notes': qnotes,
    }


def recovered_v2_row(rec):
    qflag, qnotes = quality_flags(
        rec['entity'], rec['u_org'], rec['candidate_summary'],
        rec['phase8_match_type'] == 'description')
    """Recovered Phase 8 row mapped into Phase 6 dataset schema + Phase 8 extras."""
    return {
        'dataset_tier': rec['tier'],
        'entity': rec['entity'],
        'bucket': rec['bucket'],
        'sub_type': rec['sub_type'],
        'final_match_category': rec['phase8_match_category'],
        'final_confidence_label': rec['phase8_confidence_label'],
        'final_confidence_score': rec['phase8_confidence_score'],
        'uniprot_accession': rec['u_acc'],
        'uniprot_entry_name': rec['u_entry'],
        'uniprot_organism': rec['u_org'],
        'uniprot_source_tier': rec['u_source'],
        'phytozome_code': rec['p_code'],
        'phytozome_gene_id': rec['p_gene'],
        'phytozome_base_gene_id': rec['p_base'],
        'sequence_length': rec['seq_len'],
        'has_sequence': rec['has_sequence'],
        'expected_phytozome_code': rec['expected_phytozome_code'],
        'expected_uniprot_mnemonic': rec['expected_uniprot_mnemonic'],
        'original_match_status': 'unresolved',
        'original_match_type': rec['unresolved_reason'],
        'original_confidence_tier': '',
        'ambiguity_count': rec['candidate_count'],
        'matched_species_category': rec['phase8_match_category'],
        'species_match_status': rec['species_match_status'],
        'species_resolution_status': rec['species_resolution_status'],
        'candidate_summary': rec['candidate_summary'],
        'notes': rec['notes'],
        'source_phase': '8',
        'recovered_by_phase8': 'true',
        'phase8_match_type': rec['phase8_match_type'],
        'phase8_match_category': rec['phase8_match_category'],
        'description_overlap_score': rec['description_overlap_score'],
        'phase8_notes': rec['notes'],
        'phase8_quality_flag': qflag,
        'phase8_quality_notes': qnotes,
    }


def phase6_v2_row(p6row):
    out = {k: p6row.get(k, '') for k in DATASET_FIELDS}
    qflag, qnotes = quality_flags(
        p6row.get('entity', ''), p6row.get('uniprot_organism', ''),
        p6row.get('candidate_summary', ''), is_description_row=False)
    out.update({'source_phase': '6', 'recovered_by_phase8': 'false',
                'phase8_match_type': '', 'phase8_match_category': '',
                'description_overlap_score': '', 'phase8_notes': '',
                'phase8_quality_flag': qflag, 'phase8_quality_notes': qnotes})
    return out


# ── Pre-flight ───────────────────────────────────────────────────────────────

def precheck(args):
    required = [UNRESOLVED_CSV, RAW_CSV, P6_CONFIDENCE, P6_COVERAGE,
                CROSSWALK_CSV, UNIPROT_DB, PHYTOZOME_DB]
    for p in required:
        if not p.exists():
            sys.exit(f"ERROR: required input missing: {p}")

    if MANUAL_OVERRIDE.exists():
        print(f"  Manual override file present: {MANUAL_OVERRIDE}")
    else:
        print("  No manual override file (optional) — skipping overrides.")

    # SQLite tables.
    for db, table in ((UNIPROT_DB, 'uniprot_entries'),
                      (PHYTOZOME_DB, 'phytozome_sequences')):
        con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        row = con.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
            (table,)).fetchone()
        con.close()
        if not row:
            sys.exit(f"ERROR: table {table} not found in {db}")

    # Required CSV columns.
    need_unres = {'entity', 'bucket', 'sub_type', 'unresolved_reason',
                  'species_resolution_status', 'species_match_status'}
    with open(UNRESOLVED_CSV) as fh:
        header = set(next(csv.reader(fh)))
    missing = need_unres - header
    if missing:
        sys.exit(f"ERROR: {UNRESOLVED_CSV} missing columns: {sorted(missing)}")

    need_raw = {'entity', 'extracted_gene_symbol', 'extracted_species_prefix',
                'expected_phytozome_code', 'expected_uniprot_mnemonic'}
    with open(RAW_CSV) as fh:
        header = set(next(csv.reader(fh)))
    missing = need_raw - header
    if missing:
        sys.exit(f"ERROR: {RAW_CSV} missing columns: {sorted(missing)}")

    # FTS5 capability.
    if not check_fts5_available():
        sys.exit("ERROR: SQLite FTS5 is not available in this Python build. "
                 "Phase 8 description matching requires FTS5. Aborting.")
    print("  SQLite FTS5: available.")

    # Output guard.
    existing = [p for p in ALL_OUTPUTS if p.exists()]
    if existing and not args.overwrite:
        for p in existing:
            print(f"ERROR: output exists: {p}")
        sys.exit("Refusing to overwrite. Re-run with --overwrite.")
    OUT_EXPANDED.parent.mkdir(parents=True, exist_ok=True)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Phase 8: unresolved-entry expansion")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--test-limit', type=int, metavar='N',
                      help='Process only the first N unresolved entries')
    mode.add_argument('--full', action='store_true', help='Process all unresolved entries')
    ap.add_argument('--overwrite', action='store_true',
                    help='Replace existing Phase 8 output files')
    ap.add_argument('--description-top-k', type=int, default=10)
    # Default tightened from the spec's 0.50 to 0.67 after test-driven calibration:
    # 0.50 (one shared token of a 2-token name) recovered ~63% of Bucket D, far
    # above the 15-30% target and dominated by weak partials. 0.67 (>=2 of 3
    # tokens) restores a realistic, higher-precision recovery rate.
    ap.add_argument('--description-min-partial-score', type=float, default=0.67)
    ap.add_argument('--description-high-score', type=float, default=0.80)
    ap.add_argument('--audit-sample-size', type=int, default=50)
    ap.add_argument('--random-seed', type=int, default=42)
    ap.add_argument('--min-informative-tokens', type=int, default=2,
                    help='Min informative tokens for description matching '
                         '(see module docstring re: deviation from "<3").')
    ap.add_argument('--fts-cache', type=str, default=None,
                    help='Persist/reuse the FTS5 index at this path (recommended '
                         'for the test->full workflow to avoid rebuilds).')
    ap.add_argument('--no-plant-filter', action='store_true',
                    help='Disable the plant-organism allowlist filter on fuzzy '
                         '(description / unknown-prefix) matches.')
    args = ap.parse_args()

    started = datetime.now()
    mode_label = 'FULL' if args.full else f'TEST (limit={args.test_limit})'
    print(f"Phase 8 — mode: {mode_label}")

    print("Pre-flight checks ...")
    precheck(args)

    crosswalk_codes, crosswalk_mn = load_crosswalk(CROSSWALK_CSV)
    if args.no_plant_filter:
        plant_allow = None
        print("  Plant-organism filter: DISABLED (--no-plant-filter).")
    else:
        plant_allow = build_plant_allowlist(crosswalk_mn)
        print(f"  Plant-organism allowlist: {len(plant_allow):,} mnemonics "
              "(Phase 6 datasets + crosswalk).")

    # Load unresolved (driver) rows.
    with open(UNRESOLVED_CSV) as fh:
        unresolved_header = next(csv.reader(fh))
    with open(UNRESOLVED_CSV) as fh:
        unresolved = list(csv.DictReader(fh))
    if not args.full:
        unresolved = unresolved[:args.test_limit]
    unresolved_entities = {r['entity'] for r in unresolved}
    print(f"  Unresolved rows to process: {len(unresolved):,}")

    # Load raw_matches lookup for extracted/expected fields.
    raw_lookup = {}
    with open(RAW_CSV) as fh:
        for r in csv.DictReader(fh):
            if r['entity'] in unresolved_entities:
                raw_lookup[r['entity']] = {
                    'extracted_gene_symbol': r['extracted_gene_symbol'],
                    'extracted_species_prefix': r['extracted_species_prefix'],
                    'expected_phytozome_code': r['expected_phytozome_code'],
                    'expected_uniprot_mnemonic': r['expected_uniprot_mnemonic'],
                }

    overrides, ov_warnings = load_manual_overrides(
        crosswalk_codes, crosswalk_mn, unresolved_entities)
    for w in ov_warnings:
        print(f"  WARNING (manual override): {w}")

    # Candidate cross-reference id set (for the single targeted xref scan).
    candidate_ids = set()
    for r in unresolved:
        raw = raw_lookup.get(r['entity'], {})
        for c in (r['entity'], raw.get('extracted_gene_symbol', '')):
            if c and is_id_like(c):
                candidate_ids.add(normalize_id(c))

    # Open databases (read-only).
    ucon = sqlite3.connect(f"file:{UNIPROT_DB}?mode=ro", uri=True)
    pcon = sqlite3.connect(f"file:{PHYTOZOME_DB}?mode=ro", uri=True)
    ucur = ucon.cursor()
    pcur = pcon.cursor()

    # Cross-reference map (single scan, only for non-accession xref ids).
    xref_candidates = {c for c in candidate_ids if not RE_UNIPROT_ACC.match(c)}
    xref_scan_cur = ucon.cursor()
    xref_map = build_xref_map(xref_scan_cur, xref_candidates)

    # FTS5 description index.
    if args.fts_cache:
        fts_path = Path(args.fts_cache)
        fts_path.parent.mkdir(parents=True, exist_ok=True)
        reuse = True
        cleanup_fts = False
    else:
        tmpdir = tempfile.mkdtemp(prefix='phase8_fts_')
        fts_path = Path(tmpdir) / 'uniprot_fts.db'
        reuse = False
        cleanup_fts = True
    fts_con = build_fts_index(fts_path, UNIPROT_DB.resolve(), reuse)
    fts_cur = fts_con.cursor()

    # ── Process ──────────────────────────────────────────────────────────────
    records = []           # matched records
    still = []             # still-unresolved (urow, attempted, reason, method, overlap, notes)
    attempted_by_method = Counter()
    n = 0
    print("Matching unresolved entries ...")
    for urow in unresolved:
        n += 1
        if n % 1000 == 0:
            print(f"  {n:,}/{len(unresolved):,} processed — "
                  f"recovered={len(records)} still={len(still)}")
        raw = raw_lookup.get(urow['entity'], {})
        rec = new_record(urow, raw)

        # Step A: hard skips.
        hskip = hard_skip_reason(urow)
        if hskip:
            still.append((urow, 'false', 'skipped_unrecoverable', 'none', '', hskip))
            continue

        matched = False

        # Step F: manual override (highest priority).
        if urow['entity'] in overrides:
            attempted_by_method['manual_override'] += 1
            matched = try_manual_override(rec, urow, overrides[urow['entity']], ucur, pcur)

        # Step D: cross-reference.
        if not matched and (is_id_like(urow['entity'])
                            or is_id_like(raw.get('extracted_gene_symbol', ''))):
            attempted_by_method['crossref'] += 1
            matched = try_crossref(rec, urow, raw, ucur, pcur, xref_map)

        # Step C: description matching (Bucket D).
        if not matched and urow['bucket'] == 'D':
            attempted_by_method['description'] += 1
            matched = try_description(rec, urow, raw, fts_cur, args, plant_allow)

        # Step E: unknown species prefix.
        if not matched and urow['species_resolution_status'] == 'unknown_species_prefix':
            attempted_by_method['unknown_prefix'] += 1
            matched = try_unknown_prefix(rec, urow, raw, ucur, plant_allow)

        if matched:
            records.append(rec)
        else:
            reason = ('skipped_unrecoverable'
                      if rec['best_method'] == 'description'
                      and rec['notes'].startswith('description prefilter')
                      else 'no_phase8_match')
            still.append((urow, 'true', reason, rec['best_method'],
                          rec['description_overlap_score'], rec['notes']))

    ucon.close()
    pcon.close()
    fts_con.close()
    if cleanup_fts:
        try:
            os.remove(fts_path)
            os.rmdir(fts_path.parent)
        except OSError:
            pass

    # ── Split recovered ────────────────────────────────────────────────────
    rec_conf = [r for r in records if r['tier'] == 'confidence']
    rec_cov = [r for r in records if r['tier'] == 'coverage']

    # ── Write expanded matches ───────────────────────────────────────────────
    def write_csv(path, fields, rows):
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    write_csv(OUT_EXPANDED, EXPANDED_FIELDS, [expanded_row(r) for r in records])
    write_csv(OUT_REC_CONF, V2_FIELDS, [recovered_v2_row(r) for r in rec_conf])
    write_csv(OUT_REC_COV, V2_FIELDS, [recovered_v2_row(r) for r in rec_cov])

    # ── Still unresolved ─────────────────────────────────────────────────────
    still_fields = unresolved_header + [
        'phase8_attempted', 'phase8_reason_still_unresolved',
        'phase8_best_attempted_method', 'description_overlap_score', 'phase8_notes']
    still_rows = []
    for urow, attempted, reason, method, overlap, notes in still:
        d = dict(urow)
        d['phase8_attempted'] = attempted
        d['phase8_reason_still_unresolved'] = reason
        d['phase8_best_attempted_method'] = method
        d['description_overlap_score'] = overlap
        d['phase8_notes'] = notes
        still_rows.append(d)
    write_csv(OUT_STILL, still_fields, still_rows)

    # ── Merged v2 outputs ────────────────────────────────────────────────────
    def load_p6(path):
        with open(path) as fh:
            return list(csv.DictReader(fh))
    p6_conf = load_p6(P6_CONFIDENCE)
    p6_cov = load_p6(P6_COVERAGE)

    def merge(p6_rows, recovered):
        seen = set()
        out = []
        for r in p6_rows:
            key = (r['entity'], r['uniprot_accession'], r['phytozome_gene_id'])
            seen.add(key)
            out.append(phase6_v2_row(r))
        dups = 0
        for r in recovered:
            key = (r['entity'], r['u_acc'], r['p_gene'])
            if key in seen:
                dups += 1
                continue
            seen.add(key)
            out.append(recovered_v2_row(r))
        return out, dups

    v2_conf, dup_conf = merge(p6_conf, rec_conf)
    v2_cov, dup_cov = merge(p6_cov, rec_cov)
    write_csv(OUT_V2_CONF, V2_FIELDS, v2_conf)
    write_csv(OUT_V2_COV, V2_FIELDS, v2_cov)

    # ── Audit sample (stratified, deterministic) ─────────────────────────────
    rng = random.Random(args.random_seed)
    strata = defaultdict(list)
    for r in records:
        strata[(r['tier'], r['phase8_match_category'])].append(r)
    for k in strata:
        rng.shuffle(strata[k])
    audit, keys = [], sorted(strata.keys())
    idx = {k: 0 for k in keys}
    while len(audit) < args.audit_sample_size and keys:
        progressed = False
        for k in list(keys):
            if len(audit) >= args.audit_sample_size:
                break
            if idx[k] < len(strata[k]):
                audit.append(strata[k][idx[k]])
                idx[k] += 1
                progressed = True
        if not progressed:
            break

    audit_fields = [
        'entity', 'original_bucket', 'original_unresolved_reason',
        'phase8_match_category', 'phase8_confidence_label',
        'selected_uniprot_accession', 'selected_uniprot_entry_name',
        'selected_uniprot_organism', 'selected_uniprot_source_tier',
        'selected_phytozome_code', 'selected_phytozome_gene_id',
        'selected_sequence_length', 'description_overlap_score',
        'candidate_summary', 'notes',
        'phase8_quality_flag', 'phase8_quality_notes',
        'manual_review_label', 'manual_review_notes']
    audit_rows = []
    for r in audit:
        qflag, qnotes = quality_flags(
            r['entity'], r['u_org'], r['candidate_summary'],
            r['phase8_match_category'].startswith('description'))
        audit_rows.append({
            'entity': r['entity'],
            'original_bucket': r['bucket'],
            'original_unresolved_reason': r['unresolved_reason'],
            'phase8_match_category': r['phase8_match_category'],
            'phase8_confidence_label': r['phase8_confidence_label'],
            'selected_uniprot_accession': r['u_acc'],
            'selected_uniprot_entry_name': r['u_entry'],
            'selected_uniprot_organism': r['u_org'],
            'selected_uniprot_source_tier': r['u_source'],
            'selected_phytozome_code': r['p_code'],
            'selected_phytozome_gene_id': r['p_gene'],
            'selected_sequence_length': r['seq_len'],
            'description_overlap_score': r['description_overlap_score'],
            'candidate_summary': r['candidate_summary'],
            'notes': r['notes'],
            'phase8_quality_flag': qflag,
            'phase8_quality_notes': qnotes,
            'manual_review_label': '',
            'manual_review_notes': '',
        })
    write_csv(OUT_AUDIT, audit_fields, audit_rows)

    # ── Summary report ───────────────────────────────────────────────────────
    cat_counts = Counter(r['phase8_match_category'] for r in records)
    by_bucket_total = Counter(u['bucket'] for u in unresolved)
    by_bucket_rec = Counter(r['bucket'] for r in records)
    by_reason_total = Counter(u['unresolved_reason'] for u in unresolved)
    by_reason_rec = Counter(r['unresolved_reason'] for r in records)
    still_reason = Counter(s[2] for s in still)
    n_recovered = len(records)
    total = len(unresolved)
    plant_n = 0 if plant_allow is None else len(plant_allow)
    write_summary(
        started, mode_label, args, total, still_reason.get('skipped_unrecoverable', 0),
        attempted_by_method, len(rec_conf), len(rec_cov), len(still), n_recovered,
        by_bucket_total, by_bucket_rec, by_reason_total, by_reason_rec,
        cat_counts, still_reason, ov_warnings, dup_conf, dup_cov,
        len(v2_conf), len(v2_cov), len(p6_conf), len(p6_cov),
        plant_filter_on=(plant_allow is not None), plant_n=plant_n)

    print("\n=== PHASE 8 DONE ===")
    print(f"  Attempted (post-skip): {total - still_reason.get('skipped_unrecoverable', 0):,}")
    print(f"  Recovered confidence:  {len(rec_conf):,}")
    print(f"  Recovered coverage:    {len(rec_cov):,}")
    print(f"  Still unresolved:      {len(still):,}")
    print(f"  Recovery rate:         {n_recovered/total*100:.1f}% of {total:,}")
    print(f"  v2 confidence rows:    {len(v2_conf):,} (Phase 6: {len(p6_conf):,})")
    print(f"  v2 coverage rows:      {len(v2_cov):,} (Phase 6: {len(p6_cov):,})")
    for p in ALL_OUTPUTS:
        print(f"    {p}")


def write_summary(started, mode_label, args, total, n_skipped, attempted,
                  n_conf, n_cov, n_still, n_recovered,
                  bucket_total, bucket_rec, reason_total, reason_rec,
                  cat_counts, still_reason, ov_warnings, dup_conf, dup_cov,
                  v2_conf_n, v2_cov_n, p6_conf_n, p6_cov_n,
                  plant_filter_on, plant_n):
    L = []
    L.append("# Phase 8 — Unresolved-Entry Expansion Summary\n")
    L.append(f"**Generated:** {started:%Y-%m-%d %H:%M:%S}  ")
    L.append(f"**Mode:** {mode_label}\n")
    L.append("## Thresholds")
    L.append("")
    L.append(f"- description_top_k = {args.description_top_k}")
    L.append(f"- description_min_partial_score = {args.description_min_partial_score} "
             "(tightened from spec default 0.50 after test calibration)")
    L.append(f"- description_high_score = {args.description_high_score}")
    L.append(f"- min_informative_tokens = {args.min_informative_tokens} "
             "(deviation from literal '<3'; see script docstring)")
    if plant_filter_on:
        L.append(f"- plant-organism filter = ON ({plant_n:,} allowed mnemonics "
                 "from Phase 6 datasets + crosswalk)")
    else:
        L.append("- plant-organism filter = OFF (--no-plant-filter)")
    L.append(f"- audit_sample_size = {args.audit_sample_size}")
    L.append(f"- random_seed = {args.random_seed}\n")

    L.append("## Headline counts\n")
    L.append("| Metric | Count |")
    L.append("|--------|-------|")
    L.append(f"| Original unresolved processed | {total:,} |")
    L.append(f"| Pre-filter skipped (unrecoverable) | {n_skipped:,} |")
    L.append(f"| Recovered confidence | {n_conf:,} |")
    L.append(f"| Recovered coverage | {n_cov:,} |")
    L.append(f"| Total recovered | {n_recovered:,} |")
    L.append(f"| Still unresolved | {n_still:,} |")
    rate = n_recovered / total * 100 if total else 0
    L.append(f"| Overall recovery rate | {rate:.1f}% |\n")

    L.append("## Rows attempted by method\n")
    L.append("| Method | Attempted |")
    L.append("|--------|-----------|")
    for m in ('manual_override', 'crossref', 'description', 'unknown_prefix'):
        L.append(f"| {m} | {attempted.get(m, 0):,} |")
    L.append("")

    L.append("## Recovery rate by original bucket\n")
    L.append("| Bucket | Unresolved | Recovered | Rate |")
    L.append("|--------|-----------|-----------|------|")
    for b in sorted(bucket_total):
        t, r = bucket_total[b], bucket_rec.get(b, 0)
        L.append(f"| {b} | {t:,} | {r:,} | {(r/t*100 if t else 0):.1f}% |")
    L.append("")

    L.append("## Recovery rate by original unresolved_reason\n")
    L.append("| Reason | Unresolved | Recovered | Rate |")
    L.append("|--------|-----------|-----------|------|")
    for reason in sorted(reason_total, key=lambda k: -reason_total[k]):
        t, r = reason_total[reason], reason_rec.get(reason, 0)
        L.append(f"| {reason} | {t:,} | {r:,} | {(r/t*100 if t else 0):.1f}% |")
    L.append("")

    L.append("## Recovered counts by phase8_match_category\n")
    L.append("| Category | Count |")
    L.append("|----------|-------|")
    for c, k in sorted(cat_counts.items(), key=lambda x: -x[1]):
        L.append(f"| {c} | {k:,} |")
    L.append("")

    L.append("## Still-unresolved counts by reason\n")
    L.append("| Reason | Count |")
    L.append("|--------|-------|")
    for c, k in sorted(still_reason.items(), key=lambda x: -x[1]):
        L.append(f"| {c} | {k:,} |")
    L.append("")

    L.append("## Manual overrides\n")
    L.append(f"- Warnings (malformed/invalid rows): {len(ov_warnings)}")
    for w in ov_warnings[:20]:
        L.append(f"  - {w}")
    L.append("")

    L.append("## Merged final datasets (v2)\n")
    L.append("| Dataset | Phase 6 rows | v2 rows | Duplicates skipped |")
    L.append("|---------|-------------|---------|--------------------|")
    L.append(f"| final_confidence_dataset_v2.csv | {p6_conf_n:,} | {v2_conf_n:,} | {dup_conf:,} |")
    L.append(f"| final_coverage_dataset_v2.csv | {p6_cov_n:,} | {v2_cov_n:,} | {dup_cov:,} |")
    L.append("")

    L.append("## Notes\n")
    L.append("- Fuzzy description matching increases coverage but **must be audited**; "
             "review `outputs/phase8_audit_sample.csv` before treating recovered "
             "coverage as final.")
    L.append("- Most recovered matches are coverage-tier by design; only "
             "sequence-confirmed, species-confirmed, and validated manual-override "
             "matches enter recovered confidence.")
    L.append("- **Phase 6 remains the conservative baseline.** Phase 8 outputs are "
             "an expansion layer, not a replacement. The `source_phase` column in "
             "the v2 files distinguishes original (6) from recovered (8) rows.")
    L.append("- Realistic recovery target was ~2,500-4,500 additional matches, "
             "mostly coverage-tier.")
    OUT_SUMMARY.write_text('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
