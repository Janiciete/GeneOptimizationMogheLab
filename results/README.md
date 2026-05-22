# Results Directory

This folder contains selected linkage output files from the plant protein mapping pipeline.

## Folder Structure

- confidence_linkages/
  - final_confidence_dataset_v2.csv: safest combined high-confidence linkage dataset, 425 rows.
  - phase8_recovered_confidence.csv: 49 new high-confidence links recovered by Phase 8.
- candidate_linkages/
  - phase8_recovered_coverage.csv: 5,696 Phase 8 recovered candidate links, not all 100% confident.
  - phase8_audit_sample.csv: 50-row audit sample for manual review.
- all_linkages/
  - final_coverage_dataset_v2.csv: broader combined linkage dataset, 10,564 rows.

## Notes

- Confidence files are safest for downstream use.
- Candidate and all-linkage files are broader and should be reviewed before final biological interpretation.
- Use source_phase to tell whether a row came from Phase 6 or Phase 8.
- Use phase8_quality_flag and phase8_quality_notes to identify rows needing review.
- Large generated files such as SQLite indexes, JSON exports, and FTS caches are not stored here.
