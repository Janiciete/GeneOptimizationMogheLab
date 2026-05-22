# JSON Dataset Export Summary

## Purpose

This export converts the plant protein mapping dataset (produced by Phases 1–6)
into a nested JSON-per-protein directory format plus a shared `sequences.fasta`
file. The goal is to provide a downstream-ready format for MTGNN-style ML
pipelines without replacing the canonical CSV outputs.

## Export Details

| Field | Value |
|---|---|
| Date / time | 2026-05-20 18:32:01 |
| Mode | full |
| Candidate mode | summary |
| Max candidates per protein | 1000 |
| Elapsed | 5.5s |

## Output Counts

| Metric | Count |
|---|---|
| Proteins exported | 23976 |
| Full JSON files | 23976 |
| Best-only JSON files | 23976 |
| Unique sequences in FASTA | 265 |
| Proteins with sequence_ref | 266 |
| Proteins with >1 match | 6278 |
| Proteins with truncated lists | 0 |

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


## Notes

- The CSV outputs from Phases 1–6 are preserved as validation/reference outputs.
- This export adds the JSON/FASTA format as a parallel downstream artifact.
