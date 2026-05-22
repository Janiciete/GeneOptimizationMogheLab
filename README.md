# GeneOptimizationMogheLab

This project builds a plant protein mapping pipeline that converts messy gene/protein entries into structured, confidence-scored datasets linked to UniProt metadata and Phytozome protein sequences.

The goal is to create a reliable plant protein mapping dataset that can later support downstream computational biology or machine learning work, such as adapting MTGNN-style approaches to plant protein research.

This project does not train MTGNN yet. It focuses on building the plant protein mapping dataset foundation.

---

## Project Overview

The pipeline takes a raw list of gene/protein names and:

1. Standardizes species information across Phytozome and UniProt.
2. Classifies each input entry by identifier type.
3. Builds searchable UniProt and Phytozome SQLite indexes.
4. Matches entries to UniProt metadata and/or Phytozome protein sequences.
5. Separates high-confidence linkages from broader candidate linkages.
6. Exports CSV, JSON, FASTA, audit, unresolved, and organized result files.
7. Adds Phase 8 expansion logic to recover additional linkages from previously unresolved entries.

The current project state includes Phases 1–8.

---

## Current Validated Results

Original input size:

- 23,976 original gene/protein entries

Phase 6 conservative baseline:

- 376 high-confidence matches
- 4,868 broader coverage matches
- 19,108 unresolved entries
- 266 sequence-linked records
- 265 unique FASTA sequences after sequence deduplication
- 23,976 full JSON files generated locally
- 23,976 best-only JSON files generated locally

Phase 8 expanded recovery:

- 49 additional high-confidence recovered rows
- 5,696 additional candidate/coverage recovered rows
- 5,745 total additional recovered rows
- 13,363 entries still unresolved after Phase 8
- 425 final confidence rows in final_confidence_dataset_v2.csv
- 10,564 final coverage rows in final_coverage_dataset_v2.csv

The difference between 266 sequence-linked records and 265 unique FASTA sequences means that two records reference the same deduplicated Phytozome sequence.

---

## Repository Result Organization

The repository now includes selected organized result files in:

- results/

This folder is intended for easy GitHub viewing and sharing.

### GitHub Results Folder

The GitHub results are organized as:

- results/confidence_linkages/
- results/candidate_linkages/
- results/all_linkages/
- results/README.md

### confidence_linkages/

Use this folder when you want the safest confirmed linkages.

Files:

- results/confidence_linkages/final_confidence_dataset_v2.csv
  - Phase 6 high-confidence matches plus Phase 8 recovered high-confidence matches.
  - Contains 425 final high-confidence rows.
  - This is the safest linkage file to use for downstream computational work.

- results/confidence_linkages/phase8_recovered_confidence.csv
  - Only the new high-confidence matches recovered by Phase 8.
  - Contains 49 Phase 8 recovered confidence rows.
  - Useful for seeing exactly what Phase 8 added to the confirmed dataset.

### candidate_linkages/

Use this folder when you want additional possible linkages that are not all 100% confident.

Files:

- results/candidate_linkages/phase8_recovered_coverage.csv
  - Phase 8 recovered coverage-tier/candidate matches.
  - Contains 5,696 additional recovered candidate rows.
  - Includes phase8_quality_flag and phase8_quality_notes columns.
  - Rows with warnings should be reviewed before downstream use.

- results/candidate_linkages/phase8_audit_sample.csv
  - A 50-row stratified audit sample for manual review.
  - Used to assess whether Phase 8 fuzzy matches are biologically reasonable.
  - Includes blank manual review columns.

### all_linkages/

Use this folder when you want the broadest combined linkage dataset.

Files:

- results/all_linkages/final_coverage_dataset_v2.csv
  - Phase 6 coverage dataset plus Phase 8 recovered coverage dataset.
  - Contains 10,564 broader linkage rows.
  - Includes high-confidence and candidate-level linkages.
  - Use source_phase to tell whether a row came from Phase 6 or Phase 8.
  - Use phase8_quality_flag to identify rows that need additional review.

### GitHub CSV Viewing Note

Large CSVs may not preview fully in the GitHub web interface. If a file preview looks blank or slow, use one of these options:

- Click Raw
- Click Download raw file
- Run git pull locally and open the file from the cloned repository

The CSVs are still present even if GitHub preview is limited.

---

## BioHPC Result Organization

On BioHPC, the full working project is located at:

- /local/storage/jedric/1_GeneOptimization/

The original generated outputs are in:

- /local/storage/jedric/1_GeneOptimization/outputs/

A cleaner human-readable copy of important results is in:

- /local/storage/jedric/1_GeneOptimization/organized_results/

### Best BioHPC Files to Use

Most reliable confidence dataset:

- organized_results/01_confidence_linkages/phase6_plus_phase8_final_confidence_dataset_v2.csv

Phase 8-only recovered confidence rows:

- organized_results/01_confidence_linkages/phase8_recovered_confidence_only.csv

Phase 8 recovered candidate rows:

- organized_results/02_candidate_linkages/phase8_recovered_candidate_coverage.csv

Broadest combined coverage dataset:

- organized_results/03_all_linkages/phase6_plus_phase8_final_coverage_dataset_v2.csv

Audit sample for manual review:

- organized_results/04_unresolved_and_review/phase8_audit_sample_for_manual_review.csv

Phase 8 summary:

- organized_results/04_unresolved_and_review/phase8_summary.md

JSON manifest:

- organized_results/06_json_dataset/json_manifest.csv

FASTA sequences:

- organized_results/07_sequences/phase7_sequences.fasta

### Original BioHPC Output Files

The original pipeline outputs remain in outputs/.

Important files include:

- outputs/plant_protein_dataset_confidence.csv
- outputs/plant_protein_dataset_coverage.csv
- outputs/unresolved_entries.csv
- outputs/raw_matches.csv
- outputs/phase8_recovered_confidence.csv
- outputs/phase8_recovered_coverage.csv
- outputs/final_confidence_dataset_v2.csv
- outputs/final_coverage_dataset_v2.csv
- outputs/phase8_audit_sample.csv
- outputs/phase8_summary.md

JSON and FASTA outputs are stored in:

- outputs/json_dataset/full/
- outputs/json_dataset/best_only/
- outputs/json_dataset/manifest.csv
- outputs/json_dataset/sequences.fasta

---

## Pipeline Phases

### Phase 1: Species Crosswalk

Script:

- scripts/01_build_species_crosswalk.py

Output:

- data/species_crosswalk.csv

This phase links species identifiers across Phytozome and UniProt.

Examples:

- Atha ↔ ARATH ↔ Arabidopsis thaliana
- Slyc ↔ SOLLC ↔ Solanum lycopersicum
- Gmax ↔ SOYBN ↔ Glycine max

This matters because different biological databases use different species codes for the same organism.

---

### Phase 2: Entry Classification

Script:

- scripts/02_classify_entries.py

Output:

- outputs/classified_entries.csv

This phase classifies each raw gene/protein entry into one of five buckets:

- A = clean database IDs / locus IDs / UniProt accessions
- B = species-prefixed gene symbols
- C = bare gene symbols
- D = descriptive protein names
- E = noise / unmappable entries

Examples:

- A0A178WN35 → Bucket A, UniProt accession
- AT1G01120 → Bucket A, Arabidopsis locus
- AtAPX7 → Bucket B, species-prefixed symbol
- CHS → Bucket C, bare gene symbol
- (+)-larreatricin hydroxylase 1 → Bucket D, descriptive name
- (C/T)ACGTGTC → Bucket E, noise

Phase 2 is a triage step. It does not decide the final biological match; it decides what matching strategy should be used later.

---

### Phase 3: UniProt SQLite Index

Script:

- scripts/03_build_uniprot_index.py

Local output:

- indexes/uniprot_index.db

This phase builds a searchable SQLite index from UniProt plant protein files.

It stores:

- Gene name / synonym
- Normalized gene name
- UniProt accession
- UniProt entry name
- Organism mnemonic
- Source tier: Swiss-Prot or TrEMBL
- Cross-reference metadata

The full Phase 3 build inserted over 21 million searchable UniProt rows.

This database is used to answer questions such as:

- Does this gene symbol exist in UniProt?
- Which species does the UniProt match belong to?
- Is the match reviewed Swiss-Prot or unreviewed TrEMBL?
- Does an input ID appear as a UniProt cross-reference?

---

### Phase 4: Phytozome Sequence Index

Script:

- scripts/04_build_phytozome_index.py

Local output:

- indexes/phytozome_index.db

This phase builds a searchable SQLite index from Phytozome plant protein FASTA files.

It stores:

- Phytozome species code
- Scientific name
- Raw FASTA header
- Gene ID
- Base gene ID
- Amino acid sequence
- Sequence length
- Source FASTA file

The full Phase 4 build indexed 1,085,219 protein sequences across 37 plant species.

This database is used to answer questions such as:

- Does this plant locus ID have a protein sequence?
- Which Phytozome species does it belong to?
- What is the amino acid sequence length?
- Can this input be sequence-confirmed?

---

### Phase 5: Species-Aware Raw Matching

Script:

- scripts/05_match_entries.py

Output:

- outputs/raw_matches.csv

This phase matches the classified entries against the UniProt and Phytozome indexes.

It uses species-aware logic so that a gene found in the wrong plant species is not treated as a high-confidence match.

Important match categories include:

- sequence_confirmed_match
- species_confirmed_uniprot_match
- multi_species_ambiguous
- unknown_species_prefix_global_match
- unknown_species_prefix_no_match
- species_mismatched_match
- species_confirmed_no_match
- descriptive_deferred
- noise_skipped

This matters because a gene symbol like CHS, PAL, or HSP70 may exist in many plants. The system avoids falsely assigning a gene to the wrong species by checking whether the expected species matches the database result.

---

### Phase 6: Final CSV Dataset Generation

Script:

- scripts/06_generate_datasets.py

Outputs:

- outputs/plant_protein_dataset_confidence.csv
- outputs/plant_protein_dataset_coverage.csv
- outputs/unresolved_entries.csv
- outputs/dataset_summary.md

This phase filters raw matches into usable datasets.

Confidence dataset:

- Strictest dataset
- Best for reliable downstream computational work
- Contains sequence-confirmed and species-confirmed UniProt matches

Coverage dataset:

- Broader dataset
- Includes useful but less certain entries
- Useful for exploration and manual review

Unresolved dataset:

- Contains entries that were not safely resolved
- Includes descriptive protein names, unknown prefixes, no-match rows, species mismatches, and noise

Phase 6 remains the conservative baseline of the project.

---

### Phase 7: JSON and FASTA Export

Script:

- scripts/07_export_json_dataset.py

Local output directory:

- outputs/json_dataset/

This phase exports the mapped protein dataset into a downstream-friendly JSON/FASTA format.

It creates:

- outputs/json_dataset/full/
- outputs/json_dataset/best_only/
- outputs/json_dataset/sequences.fasta
- outputs/json_dataset/manifest.csv
- outputs/json_dataset/json_dataset_summary.md

The full/ folder contains one JSON file per original entry. Each JSON stores all candidate matches that can be reconstructed from the matching output.

Use full/ for:

- Manual review
- Comparative analysis
- Debugging ambiguous matches
- Maximum-completeness workflows

The best_only/ folder contains one JSON file per original entry, but only the selected best match is stored.

Each best-only file also includes:

- match_count
- has_alternates

Use best_only/ for:

- Cleaner downstream workflows
- Faster per-protein lookup
- Machine learning preparation
- MTGNN-style input preparation

The sequences.fasta file stores amino acid sequences separately from the JSON files.

The JSON files do not contain full amino acid sequences. Instead, they reference sequences by ID, such as:

- sequence_ref: seq_Atha_AT1G01120.1

The actual sequence is stored in FASTA format in outputs/json_dataset/sequences.fasta.

This keeps JSON files small and avoids duplicating long biological sequences.

---

### Phase 8: Unresolved-Entry Expansion

Script:

- scripts/08_expand_unresolved_matches.py

Outputs:

- outputs/phase8_expanded_matches.csv
- outputs/phase8_recovered_confidence.csv
- outputs/phase8_recovered_coverage.csv
- outputs/phase8_still_unresolved.csv
- outputs/phase8_audit_sample.csv
- outputs/final_confidence_dataset_v2.csv
- outputs/final_coverage_dataset_v2.csv
- outputs/phase8_summary.md

Phase 8 attempts to recover additional linkages from entries that were unresolved after Phase 6.

It uses four recovery paths:

1. Manual override support, if data/manual_name_overrides.csv exists
2. Cross-reference matching against UniProt and Phytozome IDs
3. FTS5 description matching over UniProt descriptions
4. Unknown-prefix global gene-symbol lookup

Phase 8 includes quality controls:

- Pre-filtering of unrecoverable/noisy entries
- Plant-organism filtering
- Description-overlap score thresholds
- Confidence vs coverage separation
- phase8_quality_flag and phase8_quality_notes warning columns
- 50-row audit sample for manual review

Quality flags include:

- ok
- broad_taxon_warning
- generic_description_warning
- virus_term_warning
- non_land_plant_warning
- multiple_warnings

Phase 8 is an expansion layer, not a replacement for the conservative Phase 6 baseline.

---

## How to Generate Outputs on BioHPC

Run all commands from the repository directory:

- cd /local/storage/jedric/1_GeneOptimization

Run phases in order:

- python scripts/01_build_species_crosswalk.py
- python scripts/02_classify_entries.py
- python scripts/03_build_uniprot_index.py --full --overwrite
- python scripts/04_build_phytozome_index.py --full --overwrite
- python scripts/05_match_entries.py --full --overwrite
- python scripts/06_generate_datasets.py --full --overwrite
- python scripts/07_export_json_dataset.py --full --candidate-mode summary --overwrite
- python scripts/08_expand_unresolved_matches.py --full --overwrite --fts-cache logs/phase8_fts_cache.db

---

## How to Access Particular Results

### Safest Confirmed Linkages

Use this when you want only the most reliable results.

GitHub:

- results/confidence_linkages/final_confidence_dataset_v2.csv

BioHPC organized copy:

- organized_results/01_confidence_linkages/phase6_plus_phase8_final_confidence_dataset_v2.csv

BioHPC original output:

- outputs/final_confidence_dataset_v2.csv

Meaning:

- Phase 6 high-confidence matches
- Plus Phase 8 recovered high-confidence matches
- 425 rows total

---

### Phase 8-Only Newly Confirmed Linkages

Use this when you want to see only what Phase 8 added to the confirmed set.

GitHub:

- results/confidence_linkages/phase8_recovered_confidence.csv

BioHPC organized copy:

- organized_results/01_confidence_linkages/phase8_recovered_confidence_only.csv

BioHPC original output:

- outputs/phase8_recovered_confidence.csv

Meaning:

- 49 new high-confidence recovered rows from Phase 8

---

### Candidate Linkages

Use this when you want additional possible linkages that are not all 100% confident.

GitHub:

- results/candidate_linkages/phase8_recovered_coverage.csv

BioHPC organized copy:

- organized_results/02_candidate_linkages/phase8_recovered_candidate_coverage.csv

BioHPC original output:

- outputs/phase8_recovered_coverage.csv

Meaning:

- 5,696 additional candidate rows recovered by Phase 8
- Includes description matches, unknown-prefix symbol matches, and cross-reference candidates
- Use phase8_quality_flag and phase8_quality_notes to identify rows needing review

---

### Broadest All-Linkage Dataset

Use this when you want the broadest available dataset.

GitHub:

- results/all_linkages/final_coverage_dataset_v2.csv

BioHPC organized copy:

- organized_results/03_all_linkages/phase6_plus_phase8_final_coverage_dataset_v2.csv

BioHPC original output:

- outputs/final_coverage_dataset_v2.csv

Meaning:

- Phase 6 coverage rows plus Phase 8 recovered coverage rows
- 10,564 rows total
- Includes both high-confidence and candidate-level linkages
- Use source_phase to tell whether a row came from Phase 6 or Phase 8
- Use phase8_quality_flag to identify rows needing manual review

---

### Manual Review and Audit Files

Use these when you want to inspect uncertain matches or tune future recovery logic.

GitHub:

- results/candidate_linkages/phase8_audit_sample.csv

BioHPC organized copy:

- organized_results/04_unresolved_and_review/phase8_audit_sample_for_manual_review.csv
- organized_results/04_unresolved_and_review/phase8_summary.md

BioHPC original output:

- outputs/phase8_audit_sample.csv
- outputs/phase8_summary.md

Meaning:

- phase8_audit_sample.csv contains a 50-row manual review sample
- phase8_summary.md contains Phase 8 thresholds, counts, and recovery summaries

---

### JSON and FASTA Outputs

JSON and FASTA outputs are generated locally on BioHPC and are not fully committed to GitHub because they include tens of thousands of generated files.

BioHPC original output:

- outputs/json_dataset/full/
- outputs/json_dataset/best_only/
- outputs/json_dataset/manifest.csv
- outputs/json_dataset/sequences.fasta

BioHPC organized copy:

- organized_results/06_json_dataset/json_manifest.csv
- organized_results/07_sequences/phase7_sequences.fasta

Use full/ when you want all candidate matches per protein.

Use best_only/ when you want only the selected best match per protein.

Use sequences.fasta when you need amino acid sequences.

---

## GitHub Viewing Notes

Large CSV files may not preview fully in the GitHub web interface.

If a CSV preview looks blank, use one of these options:

- Click Raw
- Click Download raw file
- Run git pull locally and open the file from the cloned repository

The files can still be present even if the GitHub preview is limited.

---

## Files Generated Locally but Not Fully Tracked in GitHub

Large generated files are intentionally ignored by Git and are created locally on BioHPC.

These include:

- indexes/uniprot_index.db
- indexes/phytozome_index.db
- outputs/raw_matches.csv
- outputs/json_dataset/
- logs/phase8_fts_cache.db

Selected CSV results are tracked in results/, while the full local working outputs remain in outputs/ on BioHPC.

---

## Validation Commands

### Phase 6 Validation

Run:

- wc -l outputs/plant_protein_dataset_confidence.csv
- wc -l outputs/plant_protein_dataset_coverage.csv
- wc -l outputs/unresolved_entries.csv
- wc -l outputs/raw_matches.csv

Expected full-run line counts, including header row:

- 377 outputs/plant_protein_dataset_confidence.csv
- 4869 outputs/plant_protein_dataset_coverage.csv
- 19109 outputs/unresolved_entries.csv
- 23977 outputs/raw_matches.csv

---

### Phase 7 JSON/FASTA Validation

Run:

- find outputs/json_dataset/full -name "*.json" | wc -l
- find outputs/json_dataset/best_only -name "*.json" | wc -l
- wc -l outputs/json_dataset/manifest.csv
- grep -c "^>" outputs/json_dataset/sequences.fasta

Expected full-run counts:

- 23976 full JSON files
- 23976 best-only JSON files
- 23977 manifest lines including header
- 265 FASTA sequence records

---

### Phase 8 Validation

Run:

- wc -l outputs/phase8_expanded_matches.csv
- wc -l outputs/phase8_recovered_confidence.csv
- wc -l outputs/phase8_recovered_coverage.csv
- wc -l outputs/phase8_still_unresolved.csv
- wc -l outputs/final_confidence_dataset_v2.csv
- wc -l outputs/final_coverage_dataset_v2.csv
- wc -l outputs/phase8_audit_sample.csv

Expected full-run line counts, including header row:

- 5746 outputs/phase8_expanded_matches.csv
- 50 outputs/phase8_recovered_confidence.csv
- 5697 outputs/phase8_recovered_coverage.csv
- 13364 outputs/phase8_still_unresolved.csv
- 426 outputs/final_confidence_dataset_v2.csv
- 10565 outputs/final_coverage_dataset_v2.csv
- 51 outputs/phase8_audit_sample.csv

---

## Output Interpretation

### Confidence Linkages

Use confidence files when you want the safest matches.

They include:

- Sequence-confirmed matches
- Species-confirmed UniProt matches
- Phase 8 cross-reference/sequence-confirmed recovered rows

### Candidate Linkages

Use candidate files when you want more possible linkages and are willing to review uncertainty.

They may include:

- Multi-species ambiguous description matches
- Single-species description candidates
- Unknown-prefix global symbol matches
- Cross-reference matches without full species confirmation

Rows with warning flags should be manually reviewed.

### All Linkages

Use all-linkage files when you want the broadest available dataset.

The all-linkage dataset includes Phase 6 coverage rows plus Phase 8 recovered coverage rows.

### Unresolved Entries

Use unresolved files to improve the pipeline later.

Common unresolved reasons include:

- needs_description_matching
- unknown_prefix_ambiguous
- species_mismatch
- no_species_no_match
- species_confirmed_no_match
- unknown_prefix_no_match
- noise_or_not_gene

---

## Notes

- This project currently builds the plant protein mapping dataset only.
- It does not train MTGNN yet.
- It does not create embeddings yet.
- It does not build graph edges yet.
- Phase 6 remains the conservative baseline.
- Phase 8 adds broader recovered linkages but keeps warnings and confidence labels separate.
- The most reliable GitHub result file is results/confidence_linkages/final_confidence_dataset_v2.csv.
- The broadest GitHub result file is results/all_linkages/final_coverage_dataset_v2.csv.
