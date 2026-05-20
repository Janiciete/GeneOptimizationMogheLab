#!/usr/bin/env python3
"""Phase 7: Export mapped plant proteins to nested JSON/FASTA dataset format."""

import argparse
import csv
import json
import re
import shutil
import sqlite3
import sys
from datetime import datetime
from pathlib import Path

# ── Path constants ─────────────────────────────────────────────────────────────
RAW_MATCHES_CSV = Path("outputs/raw_matches.csv")
CONFIDENCE_CSV  = Path("outputs/plant_protein_dataset_confidence.csv")
COVERAGE_CSV    = Path("outputs/plant_protein_dataset_coverage.csv")
UNRESOLVED_CSV  = Path("outputs/unresolved_entries.csv")
PHYTOZOME_DB    = Path("indexes/phytozome_index.db")
UNIPROT_DB      = Path("indexes/uniprot_index.db")
JSON_OUT_DIR    = Path("outputs/json_dataset")
FULL_DIR        = JSON_OUT_DIR / "full"
BEST_DIR        = JSON_OUT_DIR / "best_only"
FASTA_PATH      = JSON_OUT_DIR / "sequences.fasta"
MANIFEST_PATH   = JSON_OUT_DIR / "manifest.csv"
SUMMARY_PATH    = JSON_OUT_DIR / "json_dataset_summary.md"

MANIFEST_FIELDS = [
    "protein_json_id", "entity", "full_json_path", "best_only_json_path",
    "dataset_membership", "has_sequence", "sequence_ref",
    "match_count", "best_match_category",
]

TIER_MAP = {
    "1": ("high",      90),
    "2": ("high",      85),
    "3": ("medium",    70),
    "4": ("low",       50),
    "5": ("deferred",   0),
    "6": ("no_match",   0),
}


# ── Argument parsing ───────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Phase 7: Export plant protein JSON/FASTA dataset"
    )
    scope = p.add_mutually_exclusive_group(required=True)
    scope.add_argument("--test-limit", type=int, metavar="N",
                       help="Process first N rows of raw_matches.csv")
    scope.add_argument("--full", action="store_true",
                       help="Process all rows")
    p.add_argument("--candidate-mode", choices=["summary", "expanded"], default="summary",
                   help="summary: use candidate_summary field only; "
                        "expanded: re-query UniProt index for full candidate list")
    p.add_argument("--max-candidates-per-protein", type=int, default=1000,
                   metavar="N",
                   help="Cap on expanded UniProt candidates per protein (default 1000)")
    p.add_argument("--overwrite", action="store_true",
                   help="Delete and recreate outputs/json_dataset if it already exists")
    return p.parse_args()


# ── Pre-flight checks ──────────────────────────────────────────────────────────

def preflight(args):
    errors = []
    for f in [RAW_MATCHES_CSV, CONFIDENCE_CSV, COVERAGE_CSV, UNRESOLVED_CSV,
              PHYTOZOME_DB, UNIPROT_DB]:
        if not f.exists():
            errors.append(f"Missing required file: {f}")
    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    for db_path, table in [(PHYTOZOME_DB, "phytozome_sequences"),
                            (UNIPROT_DB,   "uniprot_entries")]:
        con = sqlite3.connect(db_path)
        cur = con.cursor()
        cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name=?", (table,)
        )
        if not cur.fetchone():
            print(f"ERROR: Table '{table}' missing in {db_path}", file=sys.stderr)
            sys.exit(1)
        con.close()

    if JSON_OUT_DIR.exists():
        if args.overwrite:
            shutil.rmtree(JSON_OUT_DIR)
            print(f"Removed existing {JSON_OUT_DIR}")
        else:
            print(
                f"ERROR: {JSON_OUT_DIR} already exists. Use --overwrite to replace.",
                file=sys.stderr,
            )
            sys.exit(1)

    FULL_DIR.mkdir(parents=True)
    BEST_DIR.mkdir(parents=True)
    print("Pre-flight checks passed.")


# ── Dataset membership lookup ──────────────────────────────────────────────────

def load_dataset_membership():
    """Return three dicts keyed by entity for confidence, coverage, and unresolved CSVs."""
    def load(path):
        d = {}
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                d[row["entity"]] = row
        return d

    return load(CONFIDENCE_CSV), load(COVERAGE_CSV), load(UNRESOLVED_CSV)


# ── candidate_summary parsing ──────────────────────────────────────────────────

def parse_candidate_summary(summary):
    """
    Parse Phase 5 candidate_summary into two lists.

    Phytozome segment format:  pcode|gene_id|len=N
    UniProt segment format:    acc|entry_name|org|tier
    Segments separated by ' | '; entries within a segment by '; '

    Returns (phy_cands, uni_cands) where each item is a dict.
    """
    phy_cands, uni_cands = [], []
    if not summary or not summary.strip():
        return phy_cands, uni_cands

    for section in summary.split(" | "):
        for entry in section.split("; "):
            entry = entry.strip()
            if not entry:
                continue
            parts = entry.split("|")
            if len(parts) >= 3 and parts[2].startswith("len="):
                # Phytozome: pcode | gene_id | len=N
                try:
                    seq_len = int(parts[2][4:])
                except ValueError:
                    seq_len = 0
                phy_cands.append({
                    "phytozome_code": parts[0].strip(),
                    "gene_id":        parts[1].strip(),
                    "sequence_length": seq_len,
                })
            elif len(parts) >= 4 and parts[3].strip() in ("sprot", "trembl"):
                # UniProt: acc | entry_name | org | tier
                uni_cands.append({
                    "accession":         parts[0].strip(),
                    "entry_name":        parts[1].strip(),
                    "organism_mnemonic": parts[2].strip(),
                    "source_tier":       parts[3].strip(),
                })
    return phy_cands, uni_cands


# ── Sequence helpers ───────────────────────────────────────────────────────────

def make_seq_ref(phytozome_code, gene_id):
    """Stable, filesystem-safe sequence reference key."""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", gene_id)
    return f"seq_{phytozome_code}_{safe}"


def query_sequence(pcur, phytozome_code, gene_id):
    """
    Query phytozome_sequences using indexed columns.
    Returns (sequence, sequence_length, actual_gene_id) or None.
    Falls back to base_gene_id_normalized if gene_id_normalized has no hit.
    """
    norm = gene_id.upper().strip()
    pcur.execute(
        "SELECT sequence, sequence_length, gene_id "
        "FROM phytozome_sequences "
        "WHERE phytozome_code=? AND gene_id_normalized=? LIMIT 1",
        (phytozome_code, norm),
    )
    row = pcur.fetchone()
    if row:
        return row
    # Strip isoform suffix (e.g. AT1G01120.1 → AT1G01120)
    base = re.sub(r"\.\d+$", "", norm)
    if base == norm:
        return None
    pcur.execute(
        "SELECT sequence, sequence_length, gene_id "
        "FROM phytozome_sequences "
        "WHERE phytozome_code=? AND base_gene_id_normalized=? LIMIT 1",
        (phytozome_code, base),
    )
    return pcur.fetchone()


def lookup_and_write_sequence(pcur, fasta_f, written_seqs, counters,
                              phytozome_code, gene_id, entity):
    """
    Look up sequence from DB; write to FASTA if not already written.
    Returns sequence_ref string, or '' if not found.
    """
    if not phytozome_code or not gene_id:
        return ""
    db_row = query_sequence(pcur, phytozome_code, gene_id)
    if not db_row:
        return ""
    sequence, seq_len, actual_gene_id = db_row
    seq_ref = make_seq_ref(phytozome_code, actual_gene_id)
    if seq_ref not in written_seqs:
        written_seqs.add(seq_ref)
        safe_entity = entity.replace("|", "_")
        fasta_f.write(
            f">{seq_ref}|entity={safe_entity}|phytozome={phytozome_code}"
            f"|gene_id={actual_gene_id}|len={seq_len}\n"
        )
        for i in range(0, len(sequence), 60):
            fasta_f.write(sequence[i : i + 60] + "\n")
        counters["seq_written"] += 1
    return seq_ref


# ── Expanded UniProt query ─────────────────────────────────────────────────────

def query_uniprot_expanded(ucur, gene_name_normalized, max_candidates):
    """
    Re-query UniProt index by normalized gene name.
    Fetches max_candidates+1 rows to detect truncation.
    Returns (rows, truncated, total_count).
    total_count is exact when truncated (additional COUNT query).
    """
    norm = gene_name_normalized.upper().strip()
    if not norm:
        return [], False, 0
    ucur.execute(
        "SELECT accession, entry_name, organism_mnemonic, source_tier "
        "FROM uniprot_entries "
        "WHERE gene_name_normalized=? LIMIT ?",
        (norm, max_candidates + 1),
    )
    rows = ucur.fetchall()
    truncated = len(rows) > max_candidates
    if truncated:
        rows = rows[:max_candidates]
        ucur.execute(
            "SELECT COUNT(*) FROM uniprot_entries WHERE gene_name_normalized=?",
            (norm,),
        )
        total = ucur.fetchone()[0]
    else:
        total = len(rows)
    return rows, truncated, total


# ── Match object builders ──────────────────────────────────────────────────────

def int_safe(v, default=0):
    try:
        return int(v)
    except (TypeError, ValueError):
        return default


def tier_label_score(tier_str):
    label, score = TIER_MAP.get(str(tier_str), ("unknown", 0))
    return label, score


def build_phy_block(code, gene_id, base_gene_id, seq_ref, seq_len):
    if not code:
        return None
    return {
        "phytozome_code": code,
        "gene_id":        gene_id,
        "base_gene_id":   base_gene_id,
        "sequence_ref":   seq_ref,
        "sequence_length": int_safe(seq_len),
    }


def build_uni_block(accession, entry_name, organism_mnemonic, source_tier):
    if not accession:
        return None
    return {
        "accession":         accession,
        "entry_name":        entry_name,
        "organism_mnemonic": organism_mnemonic,
        "source_tier":       source_tier,
    }


def build_match_entry(rank, is_best, match_source, match_category,
                      uni_block, phy_block, conf_label, conf_score, conf_reason):
    return {
        "rank":           rank,
        "is_best":        is_best,
        "match_source":   match_source,
        "match_category": match_category,
        "uniprot":        uni_block,
        "phytozome":      phy_block,
        "confidence": {
            "label":  conf_label,
            "score":  conf_score,
            "reason": conf_reason,
        },
    }


# ── Build matches list ─────────────────────────────────────────────────────────

def build_matches_summary(row, pcur, fasta_f, written_seqs, counters,
                          phy_cands, uni_cands,
                          sel_phy_code, sel_phy_gene, sel_phy_base, sel_phy_len,
                          sel_uni_acc, sel_uni_name, sel_uni_org, sel_uni_tier,
                          best_match_source, best_seq_ref,
                          conf_label, conf_score, ds_match_category):
    """
    Build the full matches list using parsed candidate_summary only.
    Rank-1 is always the selected best match.
    Remaining candidates are appended in order: phytozome then uniprot.
    """
    matches = []
    rank = 1
    entity = row["entity"]
    raw_category = row["matched_species_category"]

    if best_match_source:
        uni_blk = build_uni_block(sel_uni_acc, sel_uni_name, sel_uni_org, sel_uni_tier)
        phy_blk = build_phy_block(sel_phy_code, sel_phy_gene, sel_phy_base,
                                   best_seq_ref, sel_phy_len)
        matches.append(build_match_entry(
            rank, True, best_match_source, ds_match_category,
            uni_blk, phy_blk, conf_label, conf_score, row["notes"],
        ))
        rank += 1

    # Remaining phytozome candidates (skip the one already used as best)
    for phy_c in phy_cands:
        if phy_c["gene_id"] == sel_phy_gene and phy_c["phytozome_code"] == sel_phy_code:
            continue
        seq_ref_c = lookup_and_write_sequence(
            pcur, fasta_f, written_seqs, counters,
            phy_c["phytozome_code"], phy_c["gene_id"], entity,
        )
        phy_blk = build_phy_block(
            phy_c["phytozome_code"], phy_c["gene_id"], "",
            seq_ref_c, phy_c["sequence_length"],
        )
        matches.append(build_match_entry(
            rank, False, "phytozome", raw_category,
            None, phy_blk, conf_label, conf_score, "",
        ))
        rank += 1

    # Remaining uniprot candidates (skip the one already used as best)
    for uni_c in uni_cands:
        if uni_c["accession"] == sel_uni_acc:
            continue
        uni_blk = build_uni_block(
            uni_c["accession"], uni_c["entry_name"],
            uni_c["organism_mnemonic"], uni_c["source_tier"],
        )
        matches.append(build_match_entry(
            rank, False, "uniprot", raw_category,
            uni_blk, None, conf_label, conf_score, "",
        ))
        rank += 1

    return matches


def build_matches_expanded(row, ucur, pcur, fasta_f, written_seqs, counters,
                            phy_cands,
                            sel_phy_code, sel_phy_gene, sel_phy_base, sel_phy_len,
                            sel_uni_acc, sel_uni_name, sel_uni_org, sel_uni_tier,
                            best_match_source, best_seq_ref,
                            conf_label, conf_score, ds_match_category,
                            max_candidates):
    """
    Build the full matches list with expanded UniProt re-query.
    Returns (matches, truncated, total_candidate_count).
    """
    entity = row["entity"]
    raw_category = row["matched_species_category"]

    gene_sym = row["extracted_gene_symbol"] or row["entity"]
    exp_rows, truncated, total_count = query_uniprot_expanded(
        ucur, gene_sym, max_candidates
    )

    matches = []
    rank = 1

    if best_match_source:
        uni_blk = build_uni_block(sel_uni_acc, sel_uni_name, sel_uni_org, sel_uni_tier)
        phy_blk = build_phy_block(sel_phy_code, sel_phy_gene, sel_phy_base,
                                   best_seq_ref, sel_phy_len)
        matches.append(build_match_entry(
            rank, True, best_match_source, ds_match_category,
            uni_blk, phy_blk, conf_label, conf_score, row["notes"],
        ))
        rank += 1

    # Remaining phytozome candidates from summary
    for phy_c in phy_cands:
        if phy_c["gene_id"] == sel_phy_gene and phy_c["phytozome_code"] == sel_phy_code:
            continue
        seq_ref_c = lookup_and_write_sequence(
            pcur, fasta_f, written_seqs, counters,
            phy_c["phytozome_code"], phy_c["gene_id"], entity,
        )
        phy_blk = build_phy_block(
            phy_c["phytozome_code"], phy_c["gene_id"], "",
            seq_ref_c, phy_c["sequence_length"],
        )
        matches.append(build_match_entry(
            rank, False, "phytozome", raw_category,
            None, phy_blk, conf_label, conf_score, "",
        ))
        rank += 1

    # Expanded UniProt candidates (deduplicate against best)
    added_accs = {sel_uni_acc} if sel_uni_acc else set()
    for r in exp_rows:
        acc = r[0]
        if acc in added_accs:
            continue
        added_accs.add(acc)
        uni_blk = build_uni_block(acc, r[1], r[2], r[3])
        matches.append(build_match_entry(
            rank, False, "uniprot", raw_category,
            uni_blk, None, conf_label, conf_score, "",
        ))
        rank += 1

    return matches, truncated, total_count


# ── Per-protein export ─────────────────────────────────────────────────────────

def process_row(idx, row, args,
                pcur, ucur, fasta_f, written_seqs, counters,
                conf_rows, cov_rows, unres_rows):
    protein_id = f"protein_{idx + 1:06d}"
    entity = row["entity"]

    # ── Classification block ──
    classification = {
        "bucket":                    row["bucket"],
        "sub_type":                  row["sub_type"],
        "extracted_species_prefix":  row["extracted_species_prefix"],
        "extracted_gene_symbol":     row["extracted_gene_symbol"],
    }

    # ── Species context block ──
    species_context = {
        "expected_phytozome_code":   row["expected_phytozome_code"],
        "expected_uniprot_mnemonic": row["expected_uniprot_mnemonic"],
        "species_resolution_status": row["species_resolution_status"],
        "species_match_status":      row["species_match_status"],
    }

    # ── Phase 5 status block ──
    phase5_status = {
        "match_status":            row["match_status"],
        "match_type":              row["match_type"],
        "confidence_tier":         int_safe(row["confidence_tier"]),
        "matched_species_category": row["matched_species_category"],
        "ambiguity_count":         int_safe(row["ambiguity_count"]),
        "all_species_hit_count":   int_safe(row["all_species_hit_count"]),
        "species_specific_hit_count": int_safe(row["species_specific_hit_count"]),
        "notes":                   row["notes"],
    }

    # ── Dataset membership block ──
    in_conf  = entity in conf_rows
    in_cov   = entity in cov_rows
    in_unres = entity in unres_rows

    ds_row    = conf_rows.get(entity) or cov_rows.get(entity)
    unres_row = unres_rows.get(entity)

    ds_match_category = (
        ds_row["final_match_category"] if ds_row
        else row["matched_species_category"]
    )

    dataset_membership = {
        "in_confidence_dataset": in_conf,
        "in_coverage_dataset":   in_cov,
        "in_unresolved":         in_unres,
        "final_match_category":  ds_row["final_match_category"] if ds_row else "",
        "final_confidence_label": ds_row["final_confidence_label"] if ds_row else "",
        "final_confidence_score": int_safe(ds_row["final_confidence_score"]) if ds_row else 0,
        "unresolved_reason":     unres_row["unresolved_reason"] if unres_row else "",
    }

    # ── Selected (best) match fields ──
    sel_phy_code = row["selected_phytozome_code"]
    sel_phy_gene = row["selected_phytozome_gene_id"]
    sel_phy_base = row["selected_phytozome_base_gene_id"]
    sel_phy_len  = int_safe(row["selected_sequence_length"])
    sel_uni_acc  = row["selected_uniprot_accession"]
    sel_uni_name = row["selected_uniprot_entry_name"]
    sel_uni_org  = row["selected_uniprot_organism"]
    sel_uni_tier = row["selected_uniprot_source_tier"]

    if sel_phy_code and sel_uni_acc:
        best_match_source = "combined"
    elif sel_phy_code:
        best_match_source = "phytozome"
    elif sel_uni_acc:
        best_match_source = "uniprot"
    else:
        best_match_source = ""

    conf_label, conf_score = tier_label_score(row["confidence_tier"])

    # ── Sequence lookup for best match ──
    best_seq_ref = lookup_and_write_sequence(
        pcur, fasta_f, written_seqs, counters,
        sel_phy_code, sel_phy_gene, entity,
    )
    # Correct the sequence_length from DB if available (summary may have stale value)
    if best_seq_ref and sel_phy_gene:
        db_row = query_sequence(pcur, sel_phy_code, sel_phy_gene)
        if db_row:
            sel_phy_len = db_row[1]

    if best_seq_ref:
        counters["has_seq_ref"] += 1

    # ── Parse candidate_summary ──
    phy_cands, uni_cands = parse_candidate_summary(row["candidate_summary"])

    # ── Build full matches list ──
    truncated = False
    total_candidate_count = 0

    if args.candidate_mode == "summary":
        matches = build_matches_summary(
            row, pcur, fasta_f, written_seqs, counters,
            phy_cands, uni_cands,
            sel_phy_code, sel_phy_gene, sel_phy_base, sel_phy_len,
            sel_uni_acc, sel_uni_name, sel_uni_org, sel_uni_tier,
            best_match_source, best_seq_ref,
            conf_label, conf_score, ds_match_category,
        )
    else:  # expanded
        matches, truncated, total_candidate_count = build_matches_expanded(
            row, ucur, pcur, fasta_f, written_seqs, counters,
            phy_cands,
            sel_phy_code, sel_phy_gene, sel_phy_base, sel_phy_len,
            sel_uni_acc, sel_uni_name, sel_uni_org, sel_uni_tier,
            best_match_source, best_seq_ref,
            conf_label, conf_score, ds_match_category,
            args.max_candidates_per_protein,
        )
        if truncated:
            counters["truncated"] += 1

    match_count = len(matches)
    # For best_only, match_count uses ambiguity_count (Phase 5's authoritative total)
    best_match_count = int_safe(row["ambiguity_count"]) or match_count
    has_alternates = (match_count > 1) or (best_match_count > 1)
    if has_alternates:
        counters["has_alternates"] += 1

    # ── Build full JSON object ──
    full_obj = {
        "protein_json_id": protein_id,
        "entity":          entity,
        "classification":  classification,
        "species_context": species_context,
        "phase5_status":   phase5_status,
        "dataset_membership": dataset_membership,
        "matches":         matches,
        "match_count":     match_count,
    }
    if truncated:
        full_obj["truncated"]               = True
        full_obj["candidate_limit"]         = args.max_candidates_per_protein
        full_obj["total_candidate_count"]   = total_candidate_count

    # ── Build best_only JSON object ──
    best_match_obj = None
    if best_match_source:
        uni_blk = build_uni_block(sel_uni_acc, sel_uni_name, sel_uni_org, sel_uni_tier)
        phy_blk = build_phy_block(sel_phy_code, sel_phy_gene, sel_phy_base,
                                   best_seq_ref, sel_phy_len)
        best_match_obj = build_match_entry(
            1, True, best_match_source, ds_match_category,
            uni_blk, phy_blk, conf_label, conf_score, row["notes"],
        )

    best_only_obj = {
        "protein_json_id":    protein_id,
        "entity":             entity,
        "classification":     classification,
        "species_context":    species_context,
        "dataset_membership": dataset_membership,
        "best_match":         best_match_obj,
        "match_count":        best_match_count,
        "has_alternates":     has_alternates,
    }

    # ── Write JSON files ──
    full_path = FULL_DIR / f"{protein_id}.json"
    best_path = BEST_DIR / f"{protein_id}.json"

    with open(full_path, "w", encoding="utf-8") as jf:
        json.dump(full_obj, jf, indent=2, ensure_ascii=False)
    with open(best_path, "w", encoding="utf-8") as jf:
        json.dump(best_only_obj, jf, indent=2, ensure_ascii=False)

    counters["full_written"] += 1
    counters["best_written"] += 1

    # ── Manifest row ──
    if in_conf:
        ds_tier = "confidence"
    elif in_cov:
        ds_tier = "coverage"
    elif in_unres:
        ds_tier = "unresolved"
    else:
        ds_tier = "none"

    return {
        "protein_json_id":     protein_id,
        "entity":              entity,
        "full_json_path":      str(full_path),
        "best_only_json_path": str(best_path),
        "dataset_membership":  ds_tier,
        "has_sequence":        "true" if best_seq_ref else "false",
        "sequence_ref":        best_seq_ref,
        "match_count":         match_count,
        "best_match_category": ds_match_category,
    }


# ── Summary markdown ───────────────────────────────────────────────────────────

def write_summary(args, counters, elapsed_sec):
    mode = "full" if args.full else f"test-limit={args.test_limit}"
    now  = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    trunc_note = (
        f"\n- **Truncated candidate lists:** {counters['truncated']} proteins "
        f"exceeded the {args.max_candidates_per_protein}-candidate limit in expanded mode."
        if counters["truncated"] else ""
    )
    content = f"""# JSON Dataset Export Summary

## Purpose

This export converts the plant protein mapping dataset (produced by Phases 1–6)
into a nested JSON-per-protein directory format plus a shared `sequences.fasta`
file. The goal is to provide a downstream-ready format for MTGNN-style ML
pipelines without replacing the canonical CSV outputs.

## Export Details

| Field | Value |
|---|---|
| Date / time | {now} |
| Mode | {mode} |
| Candidate mode | {args.candidate_mode} |
| Max candidates per protein | {args.max_candidates_per_protein} |
| Elapsed | {elapsed_sec:.1f}s |

## Output Counts

| Metric | Count |
|---|---|
| Proteins exported | {counters['processed']} |
| Full JSON files | {counters['full_written']} |
| Best-only JSON files | {counters['best_written']} |
| Unique sequences in FASTA | {counters['seq_written']} |
| Proteins with sequence_ref | {counters['has_seq_ref']} |
| Proteins with >1 match | {counters['has_alternates']} |
| Proteins with truncated lists | {counters['truncated']} |

## Directory Structure

```
outputs/json_dataset/
├── full/           — one JSON per protein with ALL candidate matches
├── best_only/      — one JSON per protein with ONLY the selected best match
├── sequences.fasta — all unique Phytozome amino acid sequences
├── manifest.csv    — index mapping protein_json_id → entity → files
└── json_dataset_summary.md
```

## full/ vs best_only/

**`full/`** contains every candidate match that the system can retrieve for each
entry. In `--candidate-mode summary` (default), this is parsed from the
`candidate_summary` field of `raw_matches.csv` (up to 5 candidates per type).
In `--candidate-mode expanded`, UniProt is re-queried by `gene_name_normalized`
to recover the complete candidate set.

**`best_only/`** contains only the highest-confidence selected match (the match
already chosen in Phases 5–6), plus a `match_count` field indicating how many
alternative candidates existed.

## sequences.fasta

Amino acid sequences are stored **only** in `sequences.fasta`, never inside JSON
files. Each JSON references its sequence via a stable `sequence_ref` key of the
form `seq_<phytozome_code>_<gene_id>`. This keeps individual JSON files small
and avoids duplicating sequences across proteins that share the same gene.
{trunc_note}

## Notes

- The CSV outputs from Phases 1–6 are preserved as validation/reference outputs.
- This export adds the JSON/FASTA format as a parallel downstream artifact.
"""
    SUMMARY_PATH.write_text(content, encoding="utf-8")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    preflight(args)

    print("Loading dataset membership tables...")
    conf_rows, cov_rows, unres_rows = load_dataset_membership()
    print(f"  confidence: {len(conf_rows):,}  coverage: {len(cov_rows):,}  "
          f"unresolved: {len(unres_rows):,}")

    pcon = sqlite3.connect(PHYTOZOME_DB)
    ucon = sqlite3.connect(UNIPROT_DB)
    pcur = pcon.cursor()
    ucur = ucon.cursor()

    written_seqs = set()
    counters = {
        "processed":    0,
        "full_written": 0,
        "best_written": 0,
        "seq_written":  0,
        "has_seq_ref":  0,
        "has_alternates": 0,
        "truncated":    0,
    }
    manifest_rows = []
    start_time = datetime.now()

    mode_label = "full" if args.full else f"test-limit={args.test_limit}"
    print(f"Starting export — mode: {mode_label}, "
          f"candidate-mode: {args.candidate_mode}")

    with open(RAW_MATCHES_CSV, newline="", encoding="utf-8") as raw_f, \
         open(FASTA_PATH, "w", encoding="utf-8") as fasta_f:

        reader = csv.DictReader(raw_f)
        for idx, row in enumerate(reader):
            if not args.full and idx >= args.test_limit:
                break

            manifest_row = process_row(
                idx, row, args,
                pcur, ucur, fasta_f, written_seqs, counters,
                conf_rows, cov_rows, unres_rows,
            )
            manifest_rows.append(manifest_row)
            counters["processed"] += 1

            if counters["processed"] % 1000 == 0:
                print(f"  ... {counters['processed']:,} proteins processed")

    pcon.close()
    ucon.close()

    # ── Write manifest ──
    print(f"Writing manifest ({len(manifest_rows):,} rows)...")
    with open(MANIFEST_PATH, "w", newline="", encoding="utf-8") as mf:
        writer = csv.DictWriter(mf, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(manifest_rows)

    # ── Write summary markdown ──
    elapsed = (datetime.now() - start_time).total_seconds()
    write_summary(args, counters, elapsed)

    # ── Final report ──
    print()
    print("── Phase 7 Export Complete " + "─" * 46)
    print(f"  proteins processed:            {counters['processed']:>8,}")
    print(f"  full JSON files written:       {counters['full_written']:>8,}")
    print(f"  best-only JSON files written:  {counters['best_written']:>8,}")
    print(f"  unique sequences written:      {counters['seq_written']:>8,}")
    print(f"  records with sequence_ref:     {counters['has_seq_ref']:>8,}")
    print(f"  records with match_count > 1:  {counters['has_alternates']:>8,}")
    print(f"  records with truncated lists:  {counters['truncated']:>8,}")
    print(f"  output directory:              {JSON_OUT_DIR}")
    print("─" * 73)


if __name__ == "__main__":
    main()
