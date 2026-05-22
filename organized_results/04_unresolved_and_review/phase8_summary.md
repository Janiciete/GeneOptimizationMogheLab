# Phase 8 — Unresolved-Entry Expansion Summary

**Generated:** 2026-05-21 18:49:14  
**Mode:** FULL

## Thresholds

- description_top_k = 10
- description_min_partial_score = 0.67 (tightened from spec default 0.50 after test calibration)
- description_high_score = 0.8
- min_informative_tokens = 2 (deviation from literal '<3'; see script docstring)
- plant-organism filter = ON (330 allowed mnemonics from Phase 6 datasets + crosswalk)
- audit_sample_size = 50
- random_seed = 42

## Headline counts

| Metric | Count |
|--------|-------|
| Original unresolved processed | 19,108 |
| Pre-filter skipped (unrecoverable) | 4,354 |
| Recovered confidence | 49 |
| Recovered coverage | 5,696 |
| Total recovered | 5,745 |
| Still unresolved | 13,363 |
| Overall recovery rate | 30.1% |

## Rows attempted by method

| Method | Attempted |
|--------|-----------|
| manual_override | 0 |
| crossref | 1,068 |
| description | 14,313 |
| unknown_prefix | 1,772 |

## Recovery rate by original bucket

| Bucket | Unresolved | Recovered | Rate |
|--------|-----------|-----------|------|
| A | 96 | 38 | 39.6% |
| B | 3,276 | 1,324 | 40.4% |
| C | 913 | 13 | 1.4% |
| D | 14,528 | 4,370 | 30.1% |
| E | 295 | 0 | 0.0% |

## Recovery rate by original unresolved_reason

| Reason | Unresolved | Recovered | Rate |
|--------|-----------|-----------|------|
| needs_description_matching | 14,535 | 4,371 | 30.1% |
| unknown_prefix_ambiguous | 1,352 | 1,323 | 97.9% |
| species_mismatch | 947 | 0 | 0.0% |
| no_species_no_match | 919 | 13 | 1.4% |
| species_confirmed_no_match | 598 | 0 | 0.0% |
| unknown_prefix_no_match | 425 | 1 | 0.2% |
| noise_or_not_gene | 295 | 0 | 0.0% |
| excluded_by_filter | 37 | 37 | 100.0% |

## Recovered counts by phase8_match_category

| Category | Count |
|----------|-------|
| description_multi_species_ambiguous | 3,001 |
| description_single_species_match | 894 |
| unknown_prefix_multi_species_ambiguous | 749 |
| unknown_prefix_single_species_candidate | 570 |
| description_partial_match | 387 |
| phase8_crossref_uniprot_match | 95 |
| phase8_crossref_species_confirmed_match | 37 |
| phase8_sequence_confirmed_match | 12 |

## Still-unresolved counts by reason

| Reason | Count |
|--------|-------|
| no_phase8_match | 9,009 |
| skipped_unrecoverable | 4,354 |

## Manual overrides

- Warnings (malformed/invalid rows): 0

## Merged final datasets (v2)

| Dataset | Phase 6 rows | v2 rows | Duplicates skipped |
|---------|-------------|---------|--------------------|
| final_confidence_dataset_v2.csv | 376 | 425 | 0 |
| final_coverage_dataset_v2.csv | 4,868 | 10,564 | 0 |

## Notes

- Fuzzy description matching increases coverage but **must be audited**; review `outputs/phase8_audit_sample.csv` before treating recovered coverage as final.
- Most recovered matches are coverage-tier by design; only sequence-confirmed, species-confirmed, and validated manual-override matches enter recovered confidence.
- **Phase 6 remains the conservative baseline.** Phase 8 outputs are an expansion layer, not a replacement. The `source_phase` column in the v2 files distinguishes original (6) from recovered (8) rows.
- Realistic recovery target was ~2,500-4,500 additional matches, mostly coverage-tier.
