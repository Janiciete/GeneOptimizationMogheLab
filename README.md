# GeneOptimizationMogheLab

This project builds a plant protein mapping pipeline that converts messy gene/protein entries into structured, confidence-scored datasets linked to UniProt metadata and Phytozome protein sequences.

## Project Goal

The goal is to create a reliable plant protein mapping dataset that can later support downstream computational biology or machine learning work, such as adapting MTGNN-style approaches to plant protein research.

This project does not train MTGNN yet. It focuses only on building the plant protein mapping dataset foundation.

## What the Pipeline Does

The pipeline takes a raw list of gene/protein names and:

1. Standardizes species information across biological databases.
2. Classifies each entry by identifier type.
3. Builds searchable UniProt and Phytozome SQLite indexes.
4. Matches entries to UniProt metadata and/or Phytozome protein sequences.
5. Assigns confidence, ambiguity, and species-match categories.
6. Produces final confidence, coverage, unresolved, JSON, and FASTA outputs.

## Final Validated Results

From 23,976 original entries, the pipeline produced:

- 376 high-confidence matches
- 4,868 broader coverage matches
- 19,108 unresolved entries
- 266 sequence-linked records
- 265 unique FASTA sequences after sequence deduplication
- 23,976 full JSON files
- 23,976 best-only JSON files

The difference between 266 sequence-linked records and 265 unique FASTA sequences means that two records reference the same deduplicated Phytozome sequence.

## Pipeline Phases

### Phase 1: Species Crosswalk

Script: `scripts/01_build_species_crosswalk.py`

Output: `data/species_crosswalk.csv`

This phase links species identifiers across Phytozome and UniProt.

Examples:

- Atha ↔ ARATH ↔ Arabidopsis thaliana
- Slyc ↔ SOLLC ↔ Solanum lycopersicum
- Gmax ↔ SOYBN ↔ Glycine max

This matters because different biological databases use different species codes for the same organism.

### Phase 2: Entry Classification

Script: `scripts/02_classify_entries.py`

Output: `outputs/classified_entries.csv`

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

Phase 2 is a first-pass triage step. It does not decide the final biological match; it decides what matching strategy should be used later.

### Phase 3: UniProt SQLite Index

Script: `scripts/03_build_uniprot_index.py`

Local output: `indexes/uniprot_index.db`

This phase builds a searchable SQLite index from UniProt plant protein files.

It stores:

- gene name / synonym
- normalized gene name
- UniProt accession
- UniProt entry name
- organism mnemonic
- source tier: Swiss-Prot or TrEMBL
- cross-reference metadata

The full Phase 3 build inserted over 21 million searchable UniProt rows.

This database is used to answer questions such as:

- Does this gene symbol exist in UniProt?
- Which species does the UniProt match belong to?
- Is the match reviewed Swiss-Prot or unreviewed TrEMBL?

### Phase 4: Phytozome Sequence Index

Script: `scripts/04_build_phytozome_index.py`

Local output: `indexes/phytozome_index.db`

This phase builds a searchable SQLite index from Phytozome plant protein FASTA files.

It stores:

- Phytozome species code
- scientific name
- raw FASTA header
- gene ID
- base gene ID
- amino acid sequence
- sequence length
- source FASTA file

The full Phase 4 build indexed 1,085,219 protein sequences across 37 plant species.

This database is used to answer:

- Does this plant locus ID have a protein sequence?
- Which Phytozome species does it belong to?
- What is the amino acid sequence length?

### Phase 5: Species-Aware Raw Matching

Script: `scripts/05_match_entries.py`

Output: `outputs/raw_matches.csv`

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

### Phase 6: Final CSV Dataset Generation

Script: `scripts/06_generate_datasets.py`

Outputs:

- `outputs/plant_protein_dataset_confidence.csv`
- `outputs/plant_protein_dataset_coverage.csv`
- `outputs/unresolved_entries.csv`
- `outputs/dataset_summary.md`

This phase filters raw matches into usable datasets.

Confidence dataset:

- `outputs/plant_protein_dataset_confidence.csv`
- Strictest dataset
- Best for reliable downstream computational work
- Contains sequence-confirmed and species-confirmed UniProt matches

Coverage dataset:

- `outputs/plant_protein_dataset_coverage.csv`
- Broader dataset
- Includes useful but less certain entries
- Useful for exploration and manual review

Unresolved entries:

- `outputs/unresolved_entries.csv`
- Contains entries that were not safely resolved
- Includes descriptive protein names, unknown prefixes, no-match rows, species mismatches, and noise

Summary report:

- `outputs/dataset_summary.md`
- Human-readable summary of dataset counts, filtering rules, and unresolved categories

### Phase 7: JSON and FASTA Export

Script: `scripts/07_export_json_dataset.py`

Local output directory: `outputs/json_dataset/`

This phase exports the mapped protein dataset into a downstream-friendly JSON/FASTA format.

It creates:

- `outputs/json_dataset/full/`
- `outputs/json_dataset/best_only/`
- `outputs/json_dataset/sequences.fasta`
- `outputs/json_dataset/manifest.csv`
- `outputs/json_dataset/json_dataset_summary.md`

The `full/` folder contains one JSON file per original entry. Each JSON stores all candidate matches that can be reconstructed from the matching output.

Use `full/` for:

- manual review
- comparative analysis
- debugging ambiguous matches
- maximum-completeness workflows

The `best_only/` folder contains one JSON file per original entry, but only the selected best match is stored.

Each best-only file also includes:

- match_count
- has_alternates

Use `best_only/` for:

- cleaner downstream workflows
- faster per-protein lookup
- machine learning preparation
- MTGNN-style input preparation

The `sequences.fasta` file stores amino acid sequences separately from the JSON files.

The JSON files do not contain full amino acid sequences. Instead, they reference sequences by ID, such as:

`"sequence_ref": "seq_Atha_AT1G01120.1"`

The actual sequence is stored in FASTA format in `outputs/json_dataset/sequences.fasta`.

This keeps JSON files small and avoids duplicating long biological sequences.

The `manifest.csv` file maps each protein JSON ID to its output files. It includes:

- protein_json_id
- entity
- full_json_path
- best_only_json_path
- dataset_membership
- has_sequence
- sequence_ref
- match_count
- best_match_category

## How to Generate the Outputs

Run phases in order.

Build species crosswalk:

`python scripts/01_build_species_crosswalk.py`

Classify input entries:

`python scripts/02_classify_entries.py`

Build UniProt index in test mode:

`python scripts/03_build_uniprot_index.py --test-limit 1000 --overwrite`

Build UniProt index in full mode:

`python scripts/03_build_uniprot_index.py --full --overwrite`

Build Phytozome sequence index in test mode:

`python scripts/04_build_phytozome_index.py --test-files 3 --test-records 1000 --overwrite`

Build Phytozome sequence index in full mode:

`python scripts/04_build_phytozome_index.py --full --overwrite`

Run raw matching in test mode:

`python scripts/05_match_entries.py --test-limit 1000 --overwrite`

Run raw matching in full mode:

`python scripts/05_match_entries.py --full --overwrite`

Generate final CSV datasets in test mode:

`python scripts/06_generate_datasets.py --test-limit 1000 --overwrite`

Generate final CSV datasets in full mode:

`python scripts/06_generate_datasets.py --full --overwrite`

Export JSON and FASTA outputs in test mode:

`python scripts/07_export_json_dataset.py --test-limit 1000 --candidate-mode summary --overwrite`

Export JSON and FASTA outputs in full mode:

`python scripts/07_export_json_dataset.py --full --candidate-mode summary --overwrite`

## Where to Find the Data

### Files visible in GitHub

The repository tracks the scripts and selected small outputs.

Tracked code files include:

- `scripts/01_build_species_crosswalk.py`
- `scripts/02_classify_entries.py`
- `scripts/03_build_uniprot_index.py`
- `scripts/04_build_phytozome_index.py`
- `scripts/05_match_entries.py`
- `scripts/06_generate_datasets.py`
- `scripts/07_export_json_dataset.py`

`outputs/classified_entries.csv` may also be visible in GitHub if it was previously committed.

### Files generated locally but not shown in GitHub

Large generated files are intentionally ignored by Git and are created locally on BioHPC.

These include:

- `indexes/uniprot_index.db`
- `indexes/phytozome_index.db`
- `outputs/raw_matches.csv`
- `taset_confidence.csv`
- `outputs/plant_protein_dataset_coverage.csv`
- `outputs/unresolved_entries.csv`
- `outputs/dataset_summary.md`
- `outputs/json_dataset/`

They are ignored because they are large generated artifacts and can be rebuilt from the scripts.

### Where to access JSON files locally

After running Phase 7, access JSON outputs here:

- `outputs/json_dataset/full/`
- `outputs/json_dataset/best_only/`

Examples:

`ls outputs/json_dataset/best_only | head`

Pretty-print one JSON file:

`python -m json.tool outputs/json_dataset/best_only/protein_002243.json | head -80`

Access the shared FASTA file:

`outputs/json_dataset/sequences.fasta`

Count FASTA records:

`grep -c "^>" outputs/json_dataset/sequences.fasta`

Access the manifest:

`outputs/json_dataset/manifest.csv`

Preview the manifest:

`head outputs/json_dataset/manifest.csv`

## Why the JSON Files Are Not Shown on GitHub

The JSON export creates tens of thousands of generated files:

- 23,976 files in `outputs/json_dataset/full/`
- 23,976 files in `outputs/json_dataset/best_only/`

These are intentionally not committed to GitHub because they are generated artifacts. Instead, they can be regenerated locally with:

`python scripts/07_export_json_dataset.py --full --candidate-mode summary --overwrite`

This keeps the repository smaller, cleaner, and easier to maintain.

## Validation Commands

Check Phase 6 final CSV outputs:

`wc -l outputs/plant_protein_dataset_confidence.csv`

`wc -l outputs/plant_protein_dataset_coverage.csv`

`wc -l outputs/unresolved_entries.csv`

`wc -l outputs/raw_matches.csv`

Expected full-run counts:

- 377 `outputs/plant_protein_dataset_confidence.csv`
- 4869 `outputs/plant_protein_dataset_coverage.csv`
- 19109 `outputs/unresolved_entries.csv`
- 23977 `outputs/raw_matches.csv`

CSV files include one header row.

Check Phase 7 JSON/FASTA outputs:

`find outputs/json_dataset/full -name "*.json" | wc -l`

`find outputs/json_dataset/best_only -name "*.json" | wc -l`

`wc -l outputs/json_dataset/manifest.csv`

`grep -c "^>" outputs/json_dataset/sequences.fasta`

Expected full-run counts:

- 23976
- 23976
- 23977
- 265

Validate JSON files parse correctly:

`python - <<'PY'
import json, glob

bad = []
for folder in ["full", "best_only"]:
    for p in glob.glob(f"outputs/json_dataset/{folder}/*.json"):
        try:
            with open(p) as f:
                json.load(f)
        except Exception as e:
            bad.append((p, str(e)))

print("bad JSON files:", len(bad))
PY`

Expected:

`bad JSON files: 0`

Validate JSON sequence references:

`python - <<'PY'
import json, glob

fasta_ids = set()
with open("outputs/json_dataset/sequences.fasta") as f:
    for line in f:
        if line.startswith(">"):
            fasta_ids.add(line[1:].split("|")[0])

missing = []
refs_checked = 0

for p in glob.glob("outputs/json_dataset/**/*.json", recursive=True):
    with open(p) as f:
        obj = json.load(f)

    bm = obj.get("best_match")
    if isinstance(bm, dict):
        phy = bm.get("phytozome")
        if isinstance(phy, dict) and phy.get("sequence_ref"):
            refs_checked += 1
            if phy["sequence_ref"] not in fasta_ids:
                missing.append((p, phy["sequence_ref"]))

    matches = obj.get("matches", [])
    if isinstance(matches, list):
        for m in matches:
            if not isinstance(m, dict):
                continue
            phy = m.get("phytozome")
            if isinstance(phy, dict) and phy.get("sequence_ref"):
                refs_checked += 1
                if phy["sequence_ref"] not in fasta_ids:
                    missing.append((p, phy["sequence_ref"]))

print("FASTA sequence IDs:", len(fasta_ids))
print("JSON sequence refs checked:", refs_checked)
print("missing sequence refs:", len(missing))
PY`

Expected:

`missing sequence refs: 0`

## Output Interpretation

### Confidence Dataset

Use this when you want the safest matches.

Contains:

- sequence_confirmed
- species_confirmed_uniprot

### Coverage Dataset

Use this when you want more rows and are willing to review ambiguity.

Contains:

- sequence_confirmed
- species_confirmed_uniprot
- ambiguous_but_useful

### Unresolved Dataset

Use this to improve the pipeline later.

Common unresolved reasons include:

- needs_description_matching
- unknown_prefix_ambiguous
- species_mismatch
- no_species_no_match
- species_confirmed_no_match
- unknown_prefix_no_match
- noise_or_not_gene

## Notes

- This project currently builds the plant protein mapping dataset only.
- It does not train MTGNN yet.
- It does not create embeddings yet.
- It does not build graph edges yet.
- Large generated files are intentionally ignored by GitHub.
- The JSON and FASTA outputs are local generated artifacts and can be recreated by running Phase 7.
