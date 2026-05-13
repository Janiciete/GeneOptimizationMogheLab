#!/usr/bin/env python3
"""Phase 6: Generate final filtered datasets and a summary report from raw Phase 5 matches."""

import argparse
import csv
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

# ── Constants ──────────────────────────────────────────────────────────────────
RAW_CSV        = Path("outputs/raw_matches.csv")
OUT_CONFIDENCE = Path("outputs/plant_protein_dataset_confidence.csv")
OUT_COVERAGE   = Path("outputs/plant_protein_dataset_coverage.csv")
OUT_UNRESOLVED = Path("outputs/unresolved_entries.csv")
OUT_SUMMARY    = Path("outputs/dataset_summary.md")

REQUIRED_COLUMNS = {
    'entity', 'bucket', 'sub_type', 'match_status', 'match_type', 'confidence_tier',
    'ambiguity_count', 'selected_uniprot_accession', 'selected_uniprot_entry_name',
    'selected_uniprot_organism', 'selected_uniprot_source_tier',
    'selected_phytozome_code', 'selected_phytozome_gene_id', 'selected_phytozome_base_gene_id',
    'selected_sequence_length', 'has_sequence', 'candidate_summary', 'notes',
    'expected_phytozome_code', 'expected_uniprot_mnemonic', 'species_resolution_status',
    'species_match_status', 'matched_species_category', 'all_species_hit_count',
    'species_specific_hit_count',
}

DATASET_FIELDS = [
    'dataset_tier',
    'entity', 'bucket', 'sub_type',
    'final_match_category', 'final_confidence_label', 'final_confidence_score',
    'uniprot_accession', 'uniprot_entry_name', 'uniprot_organism', 'uniprot_source_tier',
    'phytozome_code', 'phytozome_gene_id', 'phytozome_base_gene_id',
    'sequence_length', 'has_sequence',
    'expected_phytozome_code', 'expected_uniprot_mnemonic',
    'original_match_status', 'original_match_type', 'original_confidence_tier',
    'ambiguity_count', 'matched_species_category', 'species_match_status',
    'species_resolution_status',
    'candidate_summary', 'notes',
]

UNRESOLVED_FIELDS = [
    'entity', 'bucket', 'sub_type',
    'match_status', 'match_type', 'confidence_tier', 'ambiguity_count',
    'matched_species_category', 'species_match_status', 'species_resolution_status',
    'selected_uniprot_accession', 'selected_phytozome_gene_id',
    'candidate_summary', 'notes',
    'unresolved_reason',
]

# Never included in either dataset.
_HARD_EXCLUDE_MSC = {
    'noise_skipped', 'descriptive_deferred', 'unknown_species_prefix_no_match',
    'no_species_no_match', 'species_confirmed_no_match', 'species_mismatched_match',
}

# Excluded from coverage rule 4 (matched, tier ≤ 3).
_RULE4_EXCLUDE_MSC = {
    'species_mismatched_match', 'unknown_species_prefix_global_match',
    'unknown_species_prefix_no_match', 'descriptive_deferred', 'noise_skipped',
}


# ── Classification helpers ──────────────────────────────────────────────────────

def _fmc(msc, row):
    """Map matched_species_category + row evidence → final_match_category."""
    if (msc == 'sequence_confirmed_match'
            and row['has_sequence'].lower() == 'true'
            and row['selected_phytozome_gene_id'].strip()):
        return 'sequence_confirmed'
    if msc in ('sequence_confirmed_match', 'species_confirmed_uniprot_match'):
        return 'species_confirmed_uniprot'
    if msc == 'multi_species_ambiguous':
        return 'ambiguous_but_useful'
    return 'lower_confidence_symbol'


def _unresolved_reason(msc, match_status):
    _MAP = {
        'descriptive_deferred':                'needs_description_matching',
        'noise_skipped':                       'noise_or_not_gene',
        'species_mismatched_match':            'species_mismatch',
        'unknown_species_prefix_no_match':     'unknown_prefix_no_match',
        'unknown_species_prefix_global_match': 'unknown_prefix_ambiguous',
        'multi_species_ambiguous':             'multi_species_ambiguous',
        'species_confirmed_no_match':          'species_confirmed_no_match',
        'no_species_no_match':                 'no_species_no_match',
    }
    if msc in _MAP:
        return _MAP[msc]
    if match_status == 'no_match':
        return 'no_match'
    return 'excluded_by_filter'


def classify_row(row):
    """
    Returns one of:
      ('confidence', fmc, label, score)
      ('coverage',   fmc, label, score)
      ('unresolved', reason)
    """
    msc     = row['matched_species_category']
    tier    = int(row['confidence_tier'])
    ambig   = int(row['ambiguity_count'])
    mstatus = row['match_status']
    sms     = row['species_match_status']
    source  = row['selected_uniprot_source_tier']
    acc     = row['selected_uniprot_accession'].strip()
    phy     = row['selected_phytozome_gene_id'].strip()
    cand    = row['candidate_summary'].strip()

    # ── Confidence rule 1: sequence_confirmed_match ──────────────────────────
    if msc == 'sequence_confirmed_match':
        fmc   = _fmc(msc, row)
        score = 100 if ambig <= 1 else 95
        return 'confidence', fmc, 'high', score

    # ── Hard excludes — never useful ─────────────────────────────────────────
    if msc in _HARD_EXCLUDE_MSC or mstatus == 'no_match':
        return 'unresolved', _unresolved_reason(msc, mstatus)

    # ── Safety: no identifiers → unresolved ──────────────────────────────────
    if not acc and not phy:
        return 'unresolved', _unresolved_reason(msc, mstatus)

    # ── Confidence rule 2: species-confirmed UniProt, tier 1–2 ───────────────
    if (msc == 'species_confirmed_uniprot_match'
            and mstatus == 'matched'
            and tier in (1, 2)
            and sms == 'species_confirmed'):
        fmc   = _fmc(msc, row)
        score = 90 if source == 'sprot' else 85
        return 'confidence', fmc, 'high', score

    # ── Coverage rule 2: species-confirmed UniProt, tier ≤ 3 ─────────────────
    if msc == 'species_confirmed_uniprot_match' and tier <= 3:
        fmc   = _fmc(msc, row)
        score = 75 if source == 'sprot' else 70
        return 'coverage', fmc, 'medium', score

    # ── Coverage rule 4: any matched row, tier ≤ 3, not excluded msc ─────────
    if mstatus == 'matched' and tier <= 3 and msc not in _RULE4_EXCLUDE_MSC:
        fmc   = _fmc(msc, row)
        score = 65 if fmc == 'species_confirmed_uniprot' else 55
        return 'coverage', fmc, 'medium', score

    # ── Coverage rule 5: multi-species ambiguous with evidence ────────────────
    if (msc == 'multi_species_ambiguous'
            and mstatus == 'ambiguous'
            and tier == 4
            and cand):
        return 'coverage', 'ambiguous_but_useful', 'low', 40

    return 'unresolved', _unresolved_reason(msc, mstatus)


# ── Row builders ───────────────────────────────────────────────────────────────

def make_dataset_row(row, tier_label, fmc, clabel, cscore):
    return {
        'dataset_tier':              tier_label,
        'entity':                    row['entity'],
        'bucket':                    row['bucket'],
        'sub_type':                  row['sub_type'],
        'final_match_category':      fmc,
        'final_confidence_label':    clabel,
        'final_confidence_score':    cscore,
        'uniprot_accession':         row['selected_uniprot_accession'],
        'uniprot_entry_name':        row['selected_uniprot_entry_name'],
        'uniprot_organism':          row['selected_uniprot_organism'],
        'uniprot_source_tier':       row['selected_uniprot_source_tier'],
        'phytozome_code':            row['selected_phytozome_code'],
        'phytozome_gene_id':         row['selected_phytozome_gene_id'],
        'phytozome_base_gene_id':    row['selected_phytozome_base_gene_id'],
        'sequence_length':           row['selected_sequence_length'],
        'has_sequence':              row['has_sequence'],
        'expected_phytozome_code':   row['expected_phytozome_code'],
        'expected_uniprot_mnemonic': row['expected_uniprot_mnemonic'],
        'original_match_status':     row['match_status'],
        'original_match_type':       row['match_type'],
        'original_confidence_tier':  row['confidence_tier'],
        'ambiguity_count':           row['ambiguity_count'],
        'matched_species_category':  row['matched_species_category'],
        'species_match_status':      row['species_match_status'],
        'species_resolution_status': row['species_resolution_status'],
        'candidate_summary':         row['candidate_summary'],
        'notes':                     row['notes'],
    }


def make_unresolved_row(row, reason):
    return {
        'entity':                    row['entity'],
        'bucket':                    row['bucket'],
        'sub_type':                  row['sub_type'],
        'match_status':              row['match_status'],
        'match_type':                row['match_type'],
        'confidence_tier':           row['confidence_tier'],
        'ambiguity_count':           row['ambiguity_count'],
        'matched_species_category':  row['matched_species_category'],
        'species_match_status':      row['species_match_status'],
        'species_resolution_status': row['species_resolution_status'],
        'selected_uniprot_accession': row['selected_uniprot_accession'],
        'selected_phytozome_gene_id': row['selected_phytozome_gene_id'],
        'candidate_summary':         row['candidate_summary'],
        'notes':                     row['notes'],
        'unresolved_reason':         reason,
    }


# ── Pre-flight ─────────────────────────────────────────────────────────────────

def precheck(overwrite):
    if not RAW_CSV.exists():
        print(f"ERROR: missing input file: {RAW_CSV}")
        sys.exit(1)

    with open(RAW_CSV) as fh:
        header = set(next(csv.reader(fh)))
    missing = REQUIRED_COLUMNS - header
    if missing:
        print(f"ERROR: missing columns in {RAW_CSV}: {sorted(missing)}")
        sys.exit(1)

    ok = True
    for path in [OUT_CONFIDENCE, OUT_COVERAGE, OUT_UNRESOLVED, OUT_SUMMARY]:
        if path.exists():
            if overwrite:
                print(f"WARNING: Overwriting {path}")
                path.unlink()
            else:
                print(f"ERROR: {path} already exists. Use --overwrite to replace.")
                ok = False
    if not ok:
        sys.exit(1)

    OUT_CONFIDENCE.parent.mkdir(parents=True, exist_ok=True)


# ── Summary report ─────────────────────────────────────────────────────────────

def write_summary(generated_at, input_file, n_raw, mode_label,
                  n_conf, n_cov_total, n_unresolved,
                  fmc_conf, fmc_cov, ureason_counts, raw_msc_counts,
                  seq_conf, seq_cov):
    lines = [
        "# Plant Protein Dataset — Phase 6 Summary",
        "",
        f"**Generated:** {generated_at}",
        f"**Input:** {input_file}",
        f"**Mode:** {mode_label}",
        f"**Raw rows processed:** {n_raw:,}",
        "",
        "## Output counts",
        "",
        "| Dataset | Rows |",
        "|---------|------|",
        f"| plant_protein_dataset_confidence.csv | {n_conf:,} |",
        f"| plant_protein_dataset_coverage.csv   | {n_cov_total:,} |",
        f"| unresolved_entries.csv               | {n_unresolved:,} |",
        "",
        "## Filtering rules",
        "",
        "### Maximum-confidence dataset",
        "",
        "Include only rows meeting one of:",
        "1. `matched_species_category = sequence_confirmed_match`",
        "   → score = 100 (ambiguity ≤ 1) or 95 (ambiguity > 1)",
        "2. `matched_species_category = species_confirmed_uniprot_match`",
        "   AND `match_status = matched` AND `confidence_tier ∈ {1, 2}`",
        "   AND `species_match_status = species_confirmed`",
        "   → score = 90 (Swiss-Prot) or 85 (TrEMBL)",
        "",
        "### Maximum-coverage dataset",
        "",
        "All confidence rows plus:",
        "3. `species_confirmed_uniprot_match` with `confidence_tier ≤ 3`",
        "   → score = 75 (Swiss-Prot) or 70 (TrEMBL)",
        "4. Any `match_status = matched` row with `confidence_tier ≤ 3`",
        "   where `matched_species_category` is not in:",
        "   species_mismatched_match, unknown_species_prefix_global_match,",
        "   unknown_species_prefix_no_match, descriptive_deferred, noise_skipped",
        "   → score = 65 (species-confirmed) or 55 (other)",
        "5. `multi_species_ambiguous` with `match_status = ambiguous`,",
        "   `confidence_tier = 4`, non-empty `candidate_summary`",
        "   → score = 40",
        "",
        "Always excluded from both datasets:",
        "- noise_skipped, descriptive_deferred, unknown_species_prefix_no_match,",
        "  no_species_no_match, species_confirmed_no_match, species_mismatched_match",
        "- Any row with `match_status = no_match`",
        "- Rows with no UniProt accession AND no Phytozome gene ID",
        "",
        "## Raw input — matched_species_category distribution",
        "",
        "| Category | Count |",
        "|----------|-------|",
    ]
    for cat, cnt in sorted(raw_msc_counts.items(), key=lambda x: -x[1]):
        lines.append(f"| {cat} | {cnt:,} |")

    lines += [
        "",
        "## Confidence dataset — final_match_category distribution",
        "",
        "| Category | Count |",
        "|----------|-------|",
    ]
    for cat, cnt in sorted(fmc_conf.items(), key=lambda x: -x[1]):
        lines.append(f"| {cat} | {cnt:,} |")
    lines.append(f"| *(has_sequence = true)* | {seq_conf:,} |")

    lines += [
        "",
        "## Coverage dataset — final_match_category distribution",
        "",
        "| Category | Count |",
        "|----------|-------|",
    ]
    for cat, cnt in sorted(fmc_cov.items(), key=lambda x: -x[1]):
        lines.append(f"| {cat} | {cnt:,} |")
    lines.append(f"| *(has_sequence = true)* | {seq_cov:,} |")

    lines += [
        "",
        "## Unresolved entries — unresolved_reason distribution",
        "",
        "| Reason | Count |",
        "|--------|-------|",
    ]
    for reason, cnt in sorted(ureason_counts.items(), key=lambda x: -x[1]):
        lines.append(f"| {reason} | {cnt:,} |")

    lines += [
        "",
        "## Notes",
        "",
        "- **Confidence dataset** is the safest for automatic downstream use (GNN/ML).",
        "- **Coverage dataset** is broader and should be reviewed before high-stakes use.",
        "  It is a strict superset of the confidence dataset (`dataset_tier` column distinguishes rows).",
        "- **Descriptive names** (~14k entries) were deferred, not lost — they await",
        "  description-based matching in a future phase.",
        "- **Ambiguous gene symbols** are common in plants because many genes share names",
        "  across species; these land in `ambiguous_but_useful` or unresolved.",
        "- **Sequence-confirmed rows** (Phytozome hits) are the strongest for future ML work.",
        "- UniProt source tier: `sprot` = manually curated Swiss-Prot;",
        "  `trembl` = computationally annotated TrEMBL.",
    ]

    OUT_SUMMARY.write_text('\n'.join(lines) + '\n')


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Phase 6: Generate final filtered datasets from raw Phase 5 matches"
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--test-limit', type=int, metavar='N',
                      help='Process only the first N rows of raw_matches.csv')
    mode.add_argument('--full', action='store_true',
                      help='Process all rows')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing output files')
    args = parser.parse_args()

    precheck(args.overwrite)

    limit      = None if args.full else args.test_limit
    mode_label = 'FULL' if args.full else f'TEST (limit={limit})'
    print(f"Mode: {mode_label}")

    # ── Read and classify ──────────────────────────────────────────────────────
    confidence_rows = []   # → confidence.csv
    coverage_rows   = []   # → coverage.csv only (not confidence)
    unresolved_rows = []
    raw_msc_counts  = Counter()
    n_raw = 0

    with open(RAW_CSV) as fh:
        reader = csv.DictReader(fh)
        for i, row in enumerate(reader):
            if limit is not None and i >= limit:
                break

            raw_msc_counts[row['matched_species_category']] += 1
            n_raw += 1

            result = classify_row(row)

            if result[0] == 'confidence':
                _, fmc, clabel, cscore = result
                confidence_rows.append(make_dataset_row(row, 'confidence', fmc, clabel, cscore))
            elif result[0] == 'coverage':
                _, fmc, clabel, cscore = result
                coverage_rows.append(make_dataset_row(row, 'coverage', fmc, clabel, cscore))
            else:
                _, reason = result
                unresolved_rows.append(make_unresolved_row(row, reason))

            if n_raw % 5000 == 0:
                print(f"  {n_raw:>6,} processed — "
                      f"conf={len(confidence_rows)}  "
                      f"cov_extra={len(coverage_rows)}  "
                      f"unresolved={len(unresolved_rows)}")

    # coverage.csv = confidence rows + coverage-only rows
    all_coverage = confidence_rows + coverage_rows

    # ── Write CSV outputs ──────────────────────────────────────────────────────
    def write_csv(path, fields, rows):
        with open(path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    write_csv(OUT_CONFIDENCE, DATASET_FIELDS, confidence_rows)
    write_csv(OUT_COVERAGE,   DATASET_FIELDS, all_coverage)
    write_csv(OUT_UNRESOLVED, UNRESOLVED_FIELDS, unresolved_rows)

    # ── Stats ──────────────────────────────────────────────────────────────────
    fmc_conf    = Counter(r['final_match_category'] for r in confidence_rows)
    fmc_cov     = Counter(r['final_match_category'] for r in all_coverage)
    ureason     = Counter(r['unresolved_reason']    for r in unresolved_rows)
    seq_conf    = sum(r['has_sequence'].lower() == 'true' for r in confidence_rows)
    seq_cov     = sum(r['has_sequence'].lower() == 'true' for r in all_coverage)

    # ── Write summary.md ───────────────────────────────────────────────────────
    write_summary(
        generated_at=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        input_file=str(RAW_CSV),
        n_raw=n_raw,
        mode_label=mode_label,
        n_conf=len(confidence_rows),
        n_cov_total=len(all_coverage),
        n_unresolved=len(unresolved_rows),
        fmc_conf=fmc_conf,
        fmc_cov=fmc_cov,
        ureason_counts=ureason,
        raw_msc_counts=raw_msc_counts,
        seq_conf=seq_conf,
        seq_cov=seq_cov,
    )

    # ── Final report ───────────────────────────────────────────────────────────
    total_check = len(confidence_rows) + len(coverage_rows) + len(unresolved_rows)
    print(f"\n=== FINAL REPORT ===")
    print(f"  Raw rows processed:              {n_raw:,}")
    print(f"  Confidence rows:                 {len(confidence_rows):,}")
    print(f"  Coverage rows (total, superset): {len(all_coverage):,}")
    print(f"  Unresolved rows:                 {len(unresolved_rows):,}")
    print(f"  Sum check (conf+cov_only+unres): {total_check:,}  "
          f"({'OK' if total_check == n_raw else 'MISMATCH — check classify_row'})")
    print(f"\n  Confidence by final_match_category:")
    for cat, cnt in sorted(fmc_conf.items(), key=lambda x: -x[1]):
        print(f"    {cat}: {cnt:,}")
    print(f"\n  Coverage by final_match_category:")
    for cat, cnt in sorted(fmc_cov.items(), key=lambda x: -x[1]):
        print(f"    {cat}: {cnt:,}")
    print(f"\n  Unresolved by reason:")
    for reason, cnt in sorted(ureason.items(), key=lambda x: -x[1]):
        print(f"    {reason}: {cnt:,}")
    print(f"\n  has_sequence=true in confidence: {seq_conf:,}")
    print(f"  has_sequence=true in coverage:   {seq_cov:,}")
    print(f"\n  Outputs:")
    for p in [OUT_CONFIDENCE, OUT_COVERAGE, OUT_UNRESOLVED, OUT_SUMMARY]:
        print(f"    {p}")
    print("\nDone.")


if __name__ == '__main__':
    main()
