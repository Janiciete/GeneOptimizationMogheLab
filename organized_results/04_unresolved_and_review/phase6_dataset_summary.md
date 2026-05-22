# Plant Protein Dataset — Phase 6 Summary

**Generated:** 2026-05-13 17:24:02
**Input:** outputs/raw_matches.csv
**Mode:** FULL
**Raw rows processed:** 23,976

## Output counts

| Dataset | Rows |
|---------|------|
| plant_protein_dataset_confidence.csv | 376 |
| plant_protein_dataset_coverage.csv   | 4,868 |
| unresolved_entries.csv               | 19,108 |

## Filtering rules

### Maximum-confidence dataset

Include only rows meeting one of:
1. `matched_species_category = sequence_confirmed_match`
   → score = 100 (ambiguity ≤ 1) or 95 (ambiguity > 1)
2. `matched_species_category = species_confirmed_uniprot_match`
   AND `match_status = matched` AND `confidence_tier ∈ {1, 2}`
   AND `species_match_status = species_confirmed`
   → score = 90 (Swiss-Prot) or 85 (TrEMBL)

### Maximum-coverage dataset

All confidence rows plus:
3. `species_confirmed_uniprot_match` with `confidence_tier ≤ 3`
   → score = 75 (Swiss-Prot) or 70 (TrEMBL)
4. Any `match_status = matched` row with `confidence_tier ≤ 3`
   where `matched_species_category` is not in:
   species_mismatched_match, unknown_species_prefix_global_match,
   unknown_species_prefix_no_match, descriptive_deferred, noise_skipped
   → score = 65 (species-confirmed) or 55 (other)
5. `multi_species_ambiguous` with `match_status = ambiguous`,
   `confidence_tier = 4`, non-empty `candidate_summary`
   → score = 40

Always excluded from both datasets:
- noise_skipped, descriptive_deferred, unknown_species_prefix_no_match,
  no_species_no_match, species_confirmed_no_match, species_mismatched_match
- Any row with `match_status = no_match`
- Rows with no UniProt accession AND no Phytozome gene ID

## Raw input — matched_species_category distribution

| Category | Count |
|----------|-------|
| descriptive_deferred | 14,535 |
| multi_species_ambiguous | 3,581 |
| unknown_species_prefix_global_match | 1,352 |
| species_confirmed_uniprot_match | 1,058 |
| species_mismatched_match | 947 |
| no_species_no_match | 919 |
| species_confirmed_no_match | 598 |
| unknown_species_prefix_no_match | 425 |
| noise_skipped | 295 |
| sequence_confirmed_match | 266 |

## Confidence dataset — final_match_category distribution

| Category | Count |
|----------|-------|
| sequence_confirmed | 266 |
| species_confirmed_uniprot | 110 |
| *(has_sequence = true)* | 266 |

## Coverage dataset — final_match_category distribution

| Category | Count |
|----------|-------|
| ambiguous_but_useful | 3,581 |
| species_confirmed_uniprot | 1,021 |
| sequence_confirmed | 266 |
| *(has_sequence = true)* | 266 |

## Unresolved entries — unresolved_reason distribution

| Reason | Count |
|--------|-------|
| needs_description_matching | 14,535 |
| unknown_prefix_ambiguous | 1,352 |
| species_mismatch | 947 |
| no_species_no_match | 919 |
| species_confirmed_no_match | 598 |
| unknown_prefix_no_match | 425 |
| noise_or_not_gene | 295 |
| excluded_by_filter | 37 |

## Notes

- **Confidence dataset** is the safest for automatic downstream use (GNN/ML).
- **Coverage dataset** is broader and should be reviewed before high-stakes use.
  It is a strict superset of the confidence dataset (`dataset_tier` column distinguishes rows).
- **Descriptive names** (~14k entries) were deferred, not lost — they await
  description-based matching in a future phase.
- **Ambiguous gene symbols** are common in plants because many genes share names
  across species; these land in `ambiguous_but_useful` or unresolved.
- **Sequence-confirmed rows** (Phytozome hits) are the strongest for future ML work.
- UniProt source tier: `sprot` = manually curated Swiss-Prot;
  `trembl` = computationally annotated TrEMBL.
