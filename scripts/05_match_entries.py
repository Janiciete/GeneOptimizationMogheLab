#!/usr/bin/env python3
"""Phase 5: Match classified plant protein entries against UniProt and Phytozome indexes."""

import argparse
import csv
import re
import sqlite3
import sys
from collections import Counter
from pathlib import Path

# ── Constants ──────────────────────────────────────────────────────────────────
CLASSIFIED_CSV = Path("outputs/classified_entries.csv")
CROSSWALK_CSV  = Path("data/species_crosswalk.csv")
UNIPROT_DB     = Path("indexes/uniprot_index.db")
PHYTOZOME_DB   = Path("indexes/phytozome_index.db")
OUT_CSV        = Path("outputs/raw_matches.csv")

# Bucket B: literature species prefix → (phytozome_code, uniprot_mnemonic)
# None means not mapped in our index; global search still attempted.
PREFIX_TO_SPECIES = {
    'At':  ('Atha', 'ARATH'),  # Arabidopsis thaliana
    'Ath': ('Atha', 'ARATH'),
    'Os':  ('Osat', 'ORYSJ'),  # Oryza sativa japonica
    'Sl':  ('Slyc', 'SOLLC'),  # Solanum lycopersicum
    'Le':  ('Slyc', 'SOLLC'),  # older tomato prefix
    'Zm':  ('Zmay', 'MAIZE'),  # Zea mays
    'Gm':  ('Gmax', 'SOYBN'),  # Glycine max
    'Mt':  ('Mtru', 'MEDTR'),  # Medicago truncatula
    'Sb':  ('Sbic', 'SORBI'),  # Sorghum bicolor
    'Bd':  ('Bdis', 'BRADI'),  # Brachypodium distachyon
    'Vv':  ('Vvin', 'VITVI'),  # Vitis vinifera
    'Vvi': ('Vvin', 'VITVI'),
    'Cs':  ('Csin', None),     # Citrus sinensis (no standard UniProt mnemonic)
    'Pt':  ('Ptri', None),     # Populus trichocarpa
    'Ptr': ('Ptri', None),
    'Rc':  ('Rcom', None),     # Ricinus communis
    'Pp':  ('Ppat', 'PHYPA'),  # Physcomitrella patens
    'Eg':  ('Egra', None),     # Eucalyptus grandis
    'Mp':  ('Mpol', None),     # Marchantia polymorpha
    'Cp':  ('Cpap', 'CARPA'),  # Carica papaya
    'Dc':  ('Dcar', None),     # Daucus carota
    'Br':  (None,  'BRAOL'),   # Brassica oleracea
    'Bra': ('Bole', 'BRAOL'),
    'Nt':  (None,  'TOBAC'),   # Nicotiana tabacum
    'St':  (None,  'SOLTU'),   # Solanum tuberosum
    'Ta':  (None,  'WHEAT'),   # Triticum aestivum
    'Hv':  (None,  'HORVU'),   # Hordeum vulgare
    'Sm':  (None,  'SELML'),   # Selaginella moellendorffii
}

# Bucket A: locus sub_type → (phytozome_code, uniprot_mnemonic)
LOCUS_TO_SPECIES = {
    'arabidopsis_locus':  ('Atha', 'ARATH'),
    'tomato_locus':       ('Slyc', 'SOLLC'),
    'maize_locus':        ('Zmay', 'MAIZE'),
    'soybean_locus':      ('Gmax', 'SOYBN'),
    'sorghum_locus':      ('Sbic', 'SORBI'),
    'brachypodium_locus': ('Bdis', 'BRADI'),
    'medicago_locus':     ('Mtru', 'MEDTR'),
    'rice_locus_rap':     ('Osat', 'ORYSJ'),
    'rice_locus_msu':     ('Osat', 'ORYSJ'),
}

# Output columns — new species columns appended after notes
OUT_FIELDS = [
    'entity', 'bucket', 'sub_type', 'extracted_species_prefix', 'extracted_gene_symbol',
    'match_status', 'match_type', 'confidence_tier', 'ambiguity_count',
    'selected_uniprot_accession', 'selected_uniprot_entry_name',
    'selected_uniprot_organism', 'selected_uniprot_source_tier',
    'selected_phytozome_code', 'selected_phytozome_gene_id',
    'selected_phytozome_base_gene_id', 'selected_sequence_length', 'has_sequence',
    'candidate_summary', 'notes',
    # Species-awareness columns (new)
    'expected_phytozome_code',
    'expected_uniprot_mnemonic',
    'species_resolution_status',
    'species_match_status',
    'matched_species_category',
    'all_species_hit_count',
    'species_specific_hit_count',
]

# ── Helpers ────────────────────────────────────────────────────────────────────
def normalize(s):
    return (s or '').strip().upper()

def compute_base_id(gene_id):
    """Same conservative isoform-stripping used in Phase 4."""
    s = gene_id
    if len(s) > 2 and s[-2:].lower() == '.p':
        s = s[:-2]
    s = re.sub(r'_P\d+$', '', s)
    while True:
        m = re.match(r'^(.+)\.\d{1,2}$', s)
        if m:
            s = m.group(1)
        else:
            break
    return s

def dedup_candidates(candidates):
    """Remove duplicate (accession, entry_name, organism, source_tier) rows, preserving order."""
    seen = set()
    out  = []
    for c in candidates:
        key = (c[0], c[1], c[2], c[3])
        if key not in seen:
            seen.add(key)
            out.append(c)
    return out

def fmt_uniprot(candidates):
    """Deduplicated compact representation of up to 5 UniProt candidates."""
    return '; '.join(
        f"{c[0]}|{c[1]}|{c[2]}|{c[3]}"
        for c in dedup_candidates(candidates)[:5]
    )

def fmt_phytozome(candidates, pcode):
    """Compact representation of up to 5 Phytozome candidates."""
    seen = set()
    parts = []
    for c in candidates:
        key = (pcode, c[0])
        if key not in seen:
            seen.add(key)
            parts.append(f"{pcode}|{c[0]}|len={c[2]}")
        if len(parts) >= 5:
            break
    return '; '.join(parts)

def make_result(input_row, **kwargs):
    """Build output dict with safe defaults; override with kwargs."""
    r = {
        'entity':                       input_row['entity'],
        'bucket':                       input_row['bucket'],
        'sub_type':                     input_row['sub_type'],
        'extracted_species_prefix':     input_row['extracted_species_prefix'],
        'extracted_gene_symbol':        input_row['extracted_gene_symbol'],
        'match_status':                 '',
        'match_type':                   '',
        'confidence_tier':              6,
        'ambiguity_count':              0,
        'selected_uniprot_accession':   '',
        'selected_uniprot_entry_name':  '',
        'selected_uniprot_organism':    '',
        'selected_uniprot_source_tier': '',
        'selected_phytozome_code':         '',
        'selected_phytozome_gene_id':      '',
        'selected_phytozome_base_gene_id': '',
        'selected_sequence_length':        '',
        'has_sequence':                    'false',
        'candidate_summary':               '',
        'notes':                           '',
        # species columns
        'expected_phytozome_code':     '',
        'expected_uniprot_mnemonic':   '',
        'species_resolution_status':   '',
        'species_match_status':        '',
        'matched_species_category':    '',
        'all_species_hit_count':       0,
        'species_specific_hit_count':  0,
    }
    r.update(kwargs)
    return r


# ── SQLite query helpers ────────────────────────────────────────────────────────
# All WHERE clauses use indexed columns only.

def _uniprot_cands(ucur, where_sql, params, limit=5):
    """Return (sprot-preferred candidates, total count) for a given WHERE clause."""
    ucur.execute(f"SELECT COUNT(*) FROM uniprot_entries WHERE {where_sql}", params)
    count = ucur.fetchone()[0]
    ucur.execute(
        f"SELECT accession, entry_name, organism_mnemonic, source_tier "
        f"FROM uniprot_entries WHERE {where_sql} AND source_tier='sprot' LIMIT {limit}",
        params,
    )
    cands = ucur.fetchall()
    if not cands and count > 0:
        ucur.execute(
            f"SELECT accession, entry_name, organism_mnemonic, source_tier "
            f"FROM uniprot_entries WHERE {where_sql} LIMIT {limit}",
            params,
        )
        cands = ucur.fetchall()
    return cands, count

def count_uniprot_gene(ucur, norm_gene, organism_mnemonic=None):
    """Fast COUNT-only query; used separately when only the total is needed."""
    if organism_mnemonic:
        ucur.execute(
            "SELECT COUNT(*) FROM uniprot_entries "
            "WHERE gene_name_normalized=? AND organism_mnemonic=?",
            (norm_gene, organism_mnemonic),
        )
    else:
        ucur.execute(
            "SELECT COUNT(*) FROM uniprot_entries WHERE gene_name_normalized=?",
            (norm_gene,),
        )
    return ucur.fetchone()[0]

def query_uniprot_accession(ucur, accession):
    return _uniprot_cands(ucur, "accession=?", (normalize(accession),))

def query_uniprot_gene(ucur, norm_gene, organism_mnemonic=None):
    if organism_mnemonic:
        return _uniprot_cands(
            ucur, "gene_name_normalized=? AND organism_mnemonic=?",
            (norm_gene, organism_mnemonic),
        )
    return _uniprot_cands(ucur, "gene_name_normalized=?", (norm_gene,))

def query_phytozome(pcur, pcode, norm_entity, norm_base, limit=5):
    """Return (candidates, count); prefer exact gene_id_normalized over base_gene_id."""
    pcur.execute(
        "SELECT COUNT(*) FROM phytozome_sequences "
        "WHERE phytozome_code=? AND (gene_id_normalized=? OR base_gene_id_normalized=?)",
        (pcode, norm_entity, norm_base),
    )
    count = pcur.fetchone()[0]
    pcur.execute(
        "SELECT gene_id, base_gene_id, sequence_length "
        "FROM phytozome_sequences "
        "WHERE phytozome_code=? AND gene_id_normalized=? LIMIT ?",
        (pcode, norm_entity, limit),
    )
    cands = pcur.fetchall()
    if not cands and count > 0:
        pcur.execute(
            "SELECT gene_id, base_gene_id, sequence_length "
            "FROM phytozome_sequences "
            "WHERE phytozome_code=? AND base_gene_id_normalized=? LIMIT ?",
            (pcode, norm_base, limit),
        )
        cands = pcur.fetchall()
    return cands, count


# ── Pre-flight checks ──────────────────────────────────────────────────────────
def precheck():
    ok = True
    for path in [CLASSIFIED_CSV, CROSSWALK_CSV, UNIPROT_DB, PHYTOZOME_DB]:
        if not path.exists():
            print(f"ERROR: missing file: {path}")
            ok = False
    if not ok:
        sys.exit(1)
    for db_path, table in [(UNIPROT_DB, 'uniprot_entries'),
                           (PHYTOZOME_DB, 'phytozome_sequences')]:
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table,))
        if not cur.fetchone():
            print(f"ERROR: table '{table}' not found in {db_path}")
            ok = False
        con.close()
    if not ok:
        sys.exit(1)

def load_crosswalk(csv_path):
    pcode_to_mn = {}
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            pc = row['phytozome_code'].strip()
            mn = row['uniprot_mnemonic'].strip()
            if pc and mn:
                pcode_to_mn[pc] = mn
    return pcode_to_mn


# ── Bucket matchers ────────────────────────────────────────────────────────────

def match_e(row):
    return make_result(row,
        match_status='skipped_noise',
        match_type='noise_skipped',
        confidence_tier=6,
        notes=row['sub_type'],
        species_resolution_status='not_applicable',
        species_match_status='not_applicable',
        matched_species_category='noise_skipped',
    )

def match_d(row):
    return make_result(row,
        match_status='deferred_description',
        match_type='descriptive_deferred',
        confidence_tier=5,
        notes='descriptive matching deferred to later version',
        species_resolution_status='not_applicable',
        species_match_status='not_applicable',
        matched_species_category='descriptive_deferred',
    )

# ── Bucket A ───────────────────────────────────────────────────────────────────

def match_a(row, ucur, pcur):
    sub    = row['sub_type']
    entity = row['entity'].strip()

    if sub == 'uniprot_accession':
        return _match_accession(row, ucur, entity)
    if sub in LOCUS_TO_SPECIES:
        return _match_locus(row, ucur, pcur, entity, sub)
    # refseq_id, ncbi_gi, and other unhandled sub_types
    return make_result(row,
        match_status='deferred_description',
        match_type=f'{sub}_deferred',
        confidence_tier=5,
        notes=f'{sub} matching not implemented in Phase 5 v1',
        species_resolution_status='not_applicable',
        species_match_status='not_applicable',
        matched_species_category='descriptive_deferred',
    )

def _match_accession(row, ucur, entity):
    cands, count = query_uniprot_accession(ucur, entity)
    if count == 0:
        return make_result(row,
            match_status='no_match',
            match_type='uniprot_accession_exact',
            confidence_tier=6,
            notes='accession not found in UniProt index',
            species_resolution_status='not_applicable',
            species_match_status='species_no_match',
            matched_species_category='no_species_no_match',
        )
    best = cands[0]
    return make_result(row,
        match_status='matched',
        match_type='uniprot_accession_exact',
        confidence_tier=1,
        ambiguity_count=count,
        selected_uniprot_accession=best[0],
        selected_uniprot_entry_name=best[1],
        selected_uniprot_organism=best[2],
        selected_uniprot_source_tier=best[3],
        candidate_summary=fmt_uniprot(cands),
        species_resolution_status='not_applicable',
        species_match_status='species_confirmed',
        matched_species_category='species_confirmed_uniprot_match',
        all_species_hit_count=count,
        species_specific_hit_count=count,
    )

def _match_locus(row, ucur, pcur, entity, sub_type):
    pcode, mnemonic = LOCUS_TO_SPECIES[sub_type]
    norm_entity = normalize(entity)
    norm_base   = normalize(compute_base_id(entity))

    # Phytozome: locus IDs like GRMZM6G175135 and Solyc* match well here
    phy_cands, phy_count = query_phytozome(pcur, pcode, norm_entity, norm_base)
    has_phy = bool(phy_cands)

    # UniProt: species-filtered (TAIR/RAP IDs are stored as gene-name synonyms)
    species_cands, species_count = query_uniprot_gene(ucur, norm_entity, mnemonic)
    uni_filtered = bool(species_cands)

    # UniProt: global count (for all_species_hit_count)
    all_uni_count = count_uniprot_gene(ucur, norm_entity)

    # Display candidates: species-filtered preferred, else global
    if species_cands:
        display_uni = species_cands
    elif all_uni_count > 0:
        display_uni, _ = query_uniprot_gene(ucur, norm_entity)
    else:
        display_uni = []
    has_uni = bool(display_uni)

    # species_specific_hit_count = species-filtered UniProt + Phytozome
    spec_specific = species_count + phy_count

    # Species match status
    if has_phy or uni_filtered:
        smatch = 'species_confirmed'
    elif has_uni:
        smatch = 'species_mismatch'
    else:
        smatch = 'species_no_match'

    # matched_species_category
    if has_phy:
        msc = 'sequence_confirmed_match'
    elif uni_filtered:
        msc = 'species_confirmed_uniprot_match'
    elif has_uni:
        msc = 'species_mismatched_match'
    else:
        msc = 'species_confirmed_no_match'

    if not has_phy and not has_uni:
        return make_result(row,
            match_status='no_match',
            match_type='phytozome_gene_exact',
            confidence_tier=6,
            notes=f'no match in Phytozome ({pcode}) or UniProt ({mnemonic}) for {entity}',
            expected_phytozome_code=pcode,
            expected_uniprot_mnemonic=mnemonic,
            species_resolution_status='known_locus_species',
            species_match_status=smatch,
            matched_species_category=msc,
            all_species_hit_count=all_uni_count,
            species_specific_hit_count=spec_specific,
        )

    # Confidence tier
    if has_phy and phy_count == 1:
        tier = 1
    elif has_phy:
        tier = 2
    elif uni_filtered and species_count == 1:
        tier = 2
    elif species_count == 1 or all_uni_count == 1:
        tier = 3
    else:
        tier = 4

    ambig    = max(phy_count, all_uni_count)
    best_uni = display_uni[0] if has_uni else None
    best_phy = phy_cands[0]   if has_phy else None

    summaries = []
    if has_phy:
        summaries.append(fmt_phytozome(phy_cands, pcode))
    if has_uni:
        summaries.append(fmt_uniprot(display_uni))

    return make_result(row,
        match_status='matched' if ambig == 1 else 'ambiguous',
        match_type='phytozome_gene_exact' if has_phy else 'species_prefixed_uniprot_symbol',
        confidence_tier=tier,
        ambiguity_count=ambig,
        selected_uniprot_accession=best_uni[0] if best_uni else '',
        selected_uniprot_entry_name=best_uni[1] if best_uni else '',
        selected_uniprot_organism=best_uni[2] if best_uni else '',
        selected_uniprot_source_tier=best_uni[3] if best_uni else '',
        selected_phytozome_code=pcode if has_phy else '',
        selected_phytozome_gene_id=best_phy[0] if has_phy else '',
        selected_phytozome_base_gene_id=best_phy[1] if has_phy else '',
        selected_sequence_length=best_phy[2] if has_phy else '',
        has_sequence='true' if has_phy and best_phy[2] and best_phy[2] > 0 else 'false',
        candidate_summary=' | '.join(summaries),
        notes=f'phy={phy_count} uni_species={species_count} uni_global={all_uni_count}',
        expected_phytozome_code=pcode,
        expected_uniprot_mnemonic=mnemonic,
        species_resolution_status='known_locus_species',
        species_match_status=smatch,
        matched_species_category=msc,
        all_species_hit_count=all_uni_count,
        species_specific_hit_count=spec_specific,
    )


# ── Bucket B ───────────────────────────────────────────────────────────────────

def match_b(row, ucur, pcur):
    prefix  = row['extracted_species_prefix']
    symbol  = row['extracted_gene_symbol']
    entity  = row['entity'].strip()

    mapping      = PREFIX_TO_SPECIES.get(prefix)
    prefix_known = mapping is not None
    pcode, mnemonic = mapping if mapping else (None, None)

    norm_symbol = normalize(symbol)
    norm_entity = normalize(entity)
    norm_base   = normalize(compute_base_id(entity))

    # Global UniProt count (always — needed for all_species_hit_count)
    all_uni_count = count_uniprot_gene(ucur, norm_symbol)

    # Species-filtered UniProt candidates (if mnemonic is known)
    species_cands, species_count = [], 0
    if mnemonic:
        species_cands, species_count = query_uniprot_gene(ucur, norm_symbol, mnemonic)

    # Global UniProt candidates (for display when species-filtered is empty)
    global_cands = []
    if not species_cands and all_uni_count > 0:
        global_cands, _ = query_uniprot_gene(ucur, norm_symbol)

    # Phytozome (gene symbols rarely match locus IDs, but try for completeness)
    phy_cands, phy_count = [], 0
    if pcode:
        phy_cands, phy_count = query_phytozome(pcur, pcode, norm_entity, norm_base)

    has_phy      = bool(phy_cands)
    display_uni  = species_cands if species_cands else global_cands
    has_uni      = bool(display_uni)
    spec_specific = species_count + phy_count

    # ── species_resolution_status ────────────────────────────────────────────
    sres = 'known_species_prefix' if prefix_known else 'unknown_species_prefix'

    # ── species_match_status ─────────────────────────────────────────────────
    if not prefix_known:
        smatch = 'unknown_species_prefix'
    elif mnemonic:
        if species_count > 0 or has_phy:
            smatch = 'species_confirmed'
        elif all_uni_count > 0:
            smatch = 'species_mismatch'
        else:
            smatch = 'species_no_match'
    else:
        # prefix known but no UniProt mnemonic
        if has_phy:
            smatch = 'species_confirmed'
        elif all_uni_count > 0:
            smatch = 'multi_species_ambiguous'
        else:
            smatch = 'species_no_match'

    # ── matched_species_category ─────────────────────────────────────────────
    if has_phy:
        msc = 'sequence_confirmed_match'
    elif species_cands:
        msc = 'species_confirmed_uniprot_match'
    elif global_cands:
        if not prefix_known:
            msc = 'unknown_species_prefix_global_match'
        elif mnemonic:
            msc = 'species_mismatched_match'
        else:
            msc = 'multi_species_ambiguous'
    else:
        msc = 'unknown_species_prefix_no_match' if not prefix_known else 'species_confirmed_no_match'

    if not has_uni and not has_phy:
        return make_result(row,
            match_status='no_match',
            match_type='species_prefixed_uniprot_symbol',
            confidence_tier=6,
            notes=f'no hit for symbol={norm_symbol} prefix={prefix} mnemonic={mnemonic}',
            expected_phytozome_code=pcode or '',
            expected_uniprot_mnemonic=mnemonic or '',
            species_resolution_status=sres,
            species_match_status=smatch,
            matched_species_category=msc,
            all_species_hit_count=all_uni_count,
            species_specific_hit_count=spec_specific,
        )

    # ── Confidence tier ───────────────────────────────────────────────────────
    if mnemonic and species_count == 1:
        tier = 2   # single species-specific UniProt hit
    elif mnemonic and species_count > 1:
        tier = 3   # multiple species-specific hits
    elif all_uni_count == 1:
        tier = 3   # single global hit (no species filter possible)
    else:
        tier = 4   # multiple global hits, ambiguous

    ambig   = max(phy_count, all_uni_count)
    best_u  = display_uni[0] if has_uni else None
    best_p  = phy_cands[0]   if has_phy else None

    summaries = []
    if has_phy:
        summaries.append(fmt_phytozome(phy_cands, pcode))
    if has_uni:
        summaries.append(fmt_uniprot(display_uni))

    return make_result(row,
        match_status='matched' if ambig == 1 else 'ambiguous',
        match_type='species_prefixed_uniprot_symbol',
        confidence_tier=tier,
        ambiguity_count=ambig,
        selected_uniprot_accession=best_u[0] if best_u else '',
        selected_uniprot_entry_name=best_u[1] if best_u else '',
        selected_uniprot_organism=best_u[2] if best_u else '',
        selected_uniprot_source_tier=best_u[3] if best_u else '',
        selected_phytozome_code=pcode if has_phy else '',
        selected_phytozome_gene_id=best_p[0] if has_phy else '',
        selected_phytozome_base_gene_id=best_p[1] if has_phy else '',
        selected_sequence_length=best_p[2] if has_phy else '',
        has_sequence='true' if has_phy and best_p[2] and best_p[2] > 0 else 'false',
        candidate_summary=' | '.join(summaries),
        notes=(f'prefix={prefix} pcode={pcode} mnemonic={mnemonic} '
               f'species_count={species_count} global_count={all_uni_count}'),
        expected_phytozome_code=pcode or '',
        expected_uniprot_mnemonic=mnemonic or '',
        species_resolution_status=sres,
        species_match_status=smatch,
        matched_species_category=msc,
        all_species_hit_count=all_uni_count,
        species_specific_hit_count=spec_specific,
    )


# ── Bucket C ───────────────────────────────────────────────────────────────────

def match_c(row, ucur):
    entity = row['entity'].strip()
    norm   = normalize(entity)

    display_uni, all_count = query_uniprot_gene(ucur, norm)
    has_uni = bool(display_uni)

    # No species is provided for Bucket C — never call it species-confirmed.
    sres   = 'no_species_provided'
    smatch = 'no_species_provided'
    msc    = 'multi_species_ambiguous' if has_uni else 'no_species_no_match'

    if all_count == 0:
        return make_result(row,
            match_status='no_match',
            match_type='bare_symbol_uniprot',
            confidence_tier=6,
            species_resolution_status=sres,
            species_match_status=smatch,
            matched_species_category=msc,
            all_species_hit_count=0,
            species_specific_hit_count=0,
        )

    best = display_uni[0]
    tier = 3 if all_count == 1 else 4

    return make_result(row,
        match_status='matched' if all_count == 1 else 'ambiguous',
        match_type='bare_symbol_uniprot',
        confidence_tier=tier,
        ambiguity_count=all_count,
        selected_uniprot_accession=best[0],
        selected_uniprot_entry_name=best[1],
        selected_uniprot_organism=best[2],
        selected_uniprot_source_tier=best[3],
        candidate_summary=fmt_uniprot(display_uni),
        species_resolution_status=sres,
        species_match_status=smatch,
        matched_species_category=msc,
        all_species_hit_count=all_count,
        species_specific_hit_count=0,
    )


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Phase 5: Match classified entries against UniProt and Phytozome indexes"
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--test-limit', type=int, metavar='N',
                      help='Process only the first N classified entries')
    mode.add_argument('--full', action='store_true',
                      help='Process all classified entries')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing output CSV')
    args = parser.parse_args()

    # ── Pre-flight ────────────────────────────────────────────────────────────
    precheck()

    if OUT_CSV.exists():
        if args.overwrite:
            print(f"WARNING: Overwriting {OUT_CSV}")
            OUT_CSV.unlink()
        else:
            print(f"ERROR: {OUT_CSV} already exists. Use --overwrite to replace it.")
            sys.exit(1)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    pcode_to_mn = load_crosswalk(CROSSWALK_CSV)
    print(f"Crosswalk: {len(pcode_to_mn)} species with both codes")

    ucon = sqlite3.connect(f"file:{UNIPROT_DB}?mode=ro", uri=True)
    pcon = sqlite3.connect(f"file:{PHYTOZOME_DB}?mode=ro", uri=True)
    ucur = ucon.cursor()
    pcur = pcon.cursor()

    limit = None if args.full else args.test_limit
    print(f"Mode: {'FULL' if args.full else f'TEST (limit={limit})'}")

    # ── Stream and match ──────────────────────────────────────────────────────
    stats_bucket  = Counter()
    stats_status  = Counter()
    stats_sres    = Counter()
    stats_smatch  = Counter()
    stats_msc     = Counter()
    has_seq_count = 0
    processed     = 0

    with open(CLASSIFIED_CSV) as inf, open(OUT_CSV, 'w', newline='') as outf:
        reader = csv.DictReader(inf)
        writer = csv.DictWriter(outf, fieldnames=OUT_FIELDS)
        writer.writeheader()

        for row in reader:
            if limit is not None and processed >= limit:
                break

            bucket = row['bucket']
            if   bucket == 'E': result = match_e(row)
            elif bucket == 'D': result = match_d(row)
            elif bucket == 'A': result = match_a(row, ucur, pcur)
            elif bucket == 'B': result = match_b(row, ucur, pcur)
            elif bucket == 'C': result = match_c(row, ucur)
            else:
                result = make_result(row, match_status='no_match',
                                     notes=f'unknown bucket {bucket}',
                                     species_resolution_status='not_applicable',
                                     species_match_status='not_applicable')

            writer.writerow(result)
            processed += 1
            stats_bucket[bucket]               += 1
            stats_status[result['match_status']] += 1
            stats_sres[result['species_resolution_status']] += 1
            stats_smatch[result['species_match_status']]    += 1
            stats_msc[result['matched_species_category']]   += 1
            if result['has_sequence'] == 'true':
                has_seq_count += 1

            if processed % 1000 == 0:
                print(f"  {processed:>6,}  matched={stats_status['matched']}  "
                      f"ambiguous={stats_status['ambiguous']}  "
                      f"no_match={stats_status['no_match']}  "
                      f"skipped/deferred="
                      f"{stats_status['skipped_noise'] + stats_status['deferred_description']}")

    ucon.close()
    pcon.close()

    # ── Final report ──────────────────────────────────────────────────────────
    print(f"\n=== FINAL REPORT ===")
    print(f"  Processed:                    {processed:,}")
    print(f"  By bucket:                    {dict(stats_bucket)}")
    print(f"  By match_status:              {dict(stats_status)}")
    print(f"  By species_resolution_status: {dict(stats_sres)}")
    print(f"  By species_match_status:      {dict(stats_smatch)}")
    print(f"  By matched_species_category:  {dict(stats_msc)}")
    print(f"  Has sequence:                 {has_seq_count:,}")
    print(f"  Output:                       {OUT_CSV}")
    print("\nDone.")


if __name__ == '__main__':
    main()
