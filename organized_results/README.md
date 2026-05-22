# Organized Results Directory

This folder is a clean, human-readable copy of the important output files from the plant protein mapping pipeline.

The original `outputs/` folder is still kept for pipeline execution. This folder is only for easier navigation, sharing, and review.

## Folder Guide

### 01_confidence_linkages/

Safest linkage files.

Use these when you want the most reliable results.

Includes:
- Phase 6 conservative high-confidence matches
- Phase 8 recovered high-confidence matches
- Final confidence dataset combining Phase 6 + Phase 8

Recommended main file:
- `phase6_plus_phase8_final_confidence_dataset_v2.csv`

### 02_candidate_linkages/

Useful candidate linkages that are not all 100% confident.

These include broader coverage matches and Phase 8 fuzzy/expanded matches.

Use these for:
- manual review
- exploratory analysis
- finding additional possible linkages

Important:
Rows with `phase8_quality_flag` other than `ok` should be reviewed before downstream use.

### 03_all_linkages/

Broadest combined linkage datasets.

Recommended main file:
- `phase6_plus_phase8_final_coverage_dataset_v2.csv`

This combines Phase 6 coverage rows with Phase 8 recovered candidate rows.

Use:
- `source_phase` to tell whether a row came from Phase 6 or Phase 8
- `phase8_quality_flag` to identify rows needing review

### 04_unresolved_and_review/

Files for review and improvement.

Includes:
- unresolved entries
- still unresolved after Phase 8
- audit sample
- summary reports

Use this folder to decide what could be improved in future phases.

### 05_phase_outputs/

Intermediate phase outputs.

These are useful for debugging and provenance, but are not the main final datasets.

### 06_json_dataset/

Metadata and manifest files for the JSON export.

The many individual JSON files may not be copied here because there are many of them. The manifest explains the JSON dataset structure.

### 07_sequences/

FASTA sequence outputs from the JSON/sequence export.

### 08_archives/

Compressed archives, if any exist.

## Main Files to Share With PI

Most reliable:
- `01_confidence_linkages/phase6_plus_phase8_final_confidence_dataset_v2.csv`

Expanded candidates:
- `02_candidate_linkages/phase8_recovered_candidate_coverage.csv`

Broadest all-linkage dataset:
- `03_all_linkages/phase6_plus_phase8_final_coverage_dataset_v2.csv`

Manual review:
- `04_unresolved_and_review/phase8_audit_sample_for_manual_review.csv`
- `04_unresolved_and_review/phase8_summary.md`

## Interpretation

- Confidence files = safest confirmed linkages.
- Candidate files = additional possible linkages, not all fully confirmed.
- All-linkage files = broader merged datasets.
- Unresolved/review files = entries needing more work or human review.
- Phase 8 is an expansion layer, not a replacement for the conservative Phase 6 baseline.
