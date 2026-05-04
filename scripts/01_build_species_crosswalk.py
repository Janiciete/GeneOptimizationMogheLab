#!/usr/bin/env python3
"""Phase 1: Build species crosswalk linking Phytozome codes <-> UniProt mnemonics."""

import csv
import sys
from pathlib import Path

PHYTOZOME_TAB  = Path("/local/storage/0_databases/2_protein/2_phytozome/proteome_files.tab")
UNIPROT_ORGMAP = Path("/local/storage/0_databases/2_protein/3_uniprot_2024/uniprot_sprot_plants_OrganismMap.txt")
OUT_CSV        = Path("/local/storage/jedric/1_GeneOptimization/data/species_crosswalk.csv")

# Manual corrections: phytozome scientific name (after underscore→space) -> canonical name for UniProt lookup
MANUAL_FIXES = {
    "Mallus domesticus": "Malus domestica",
}

# ── 1. Read Phytozome table ───────────────────────────────────────────────────
phytozome_species = []   # list of (phytozome_code, scientific_name_for_lookup, original_name)
with open(PHYTOZOME_TAB) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        filename   = parts[0]                          # e.g. "Atha_prot.fa"
        raw_name   = parts[1].replace("_", " ")        # e.g. "Arabidopsis thaliana"
        pcode      = filename[:4]                      # e.g. "Atha"
        lookup_name = MANUAL_FIXES.get(raw_name, raw_name)
        phytozome_species.append((pcode, lookup_name, raw_name))

print(f"Phytozome species loaded: {len(phytozome_species)}")

# ── 2. Read UniProt OrganismMap — build {scientific_name -> mnemonic} ─────────
uniprot_name_to_mnemonic = {}   # {canonical scientific name -> 5-letter mnemonic}
with open(UNIPROT_ORGMAP) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        protein  = row["Protein"].strip()          # e.g. "10HGO_CATRO"
        organism = row["Organism"].strip()         # e.g. "Catharanthus roseus"
        mnemonic = protein.rsplit("_", 1)[-1]     # e.g. "CATRO"
        if organism and mnemonic and organism not in uniprot_name_to_mnemonic:
            uniprot_name_to_mnemonic[organism] = mnemonic

print(f"UniProt unique organisms loaded: {len(uniprot_name_to_mnemonic)}")

# ── 3. Join & build crosswalk ─────────────────────────────────────────────────
rows = []
matched   = 0
unmatched = []

for pcode, lookup_name, original_name in phytozome_species:
    mnemonic = uniprot_name_to_mnemonic.get(lookup_name, "")
    notes    = ""

    if not mnemonic:
        # Try prefix match (e.g. "Oryza sativa" matches "Oryza sativa subsp. japonica")
        prefix_hits = [
            (name, mn) for name, mn in uniprot_name_to_mnemonic.items()
            if name.startswith(lookup_name)
        ]
        if len(prefix_hits) == 1:
            mnemonic = prefix_hits[0][1]
            notes = f"prefix match on '{prefix_hits[0][0]}'"
        elif len(prefix_hits) > 1:
            # Pick shortest (most specific base name)
            prefix_hits.sort(key=lambda x: len(x[0]))
            mnemonic = prefix_hits[0][1]
            notes = f"prefix match (ambiguous, picked shortest): '{prefix_hits[0][0]}'"

    if original_name != lookup_name:
        note_fix = f"typo corrected: '{original_name}' -> '{lookup_name}'"
        notes = (notes + "; " + note_fix).lstrip("; ") if notes else note_fix

    if mnemonic:
        matched += 1
    else:
        unmatched.append((pcode, lookup_name))

    rows.append({
        "phytozome_code":   pcode,
        "uniprot_mnemonic": mnemonic,
        "scientific_name":  lookup_name,
        "notes":            notes,
    })

# ── 4. Save output ────────────────────────────────────────────────────────────
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_CSV, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["phytozome_code", "uniprot_mnemonic", "scientific_name", "notes"])
    writer.writeheader()
    writer.writerows(rows)

print(f"\nSaved: {OUT_CSV}  ({len(rows)} rows)")

# ── 5. Report ─────────────────────────────────────────────────────────────────
print(f"\n=== REPORT ===")
print(f"Total Phytozome species:  {len(rows)}")
print(f"Matched to UniProt:       {matched}")
print(f"Unmatched:                {len(unmatched)}")

if unmatched:
    print("\nUnmatched species:")
    for pcode, name in unmatched:
        print(f"  {pcode}  {name}")

print("\nSample matched rows (first 5 with a mnemonic):")
shown = 0
for r in rows:
    if r["uniprot_mnemonic"]:
        print(f"  {r['phytozome_code']:6}  {r['uniprot_mnemonic']:8}  {r['scientific_name']}")
        shown += 1
        if shown >= 5:
            break

notes_rows = [r for r in rows if r["notes"]]
if notes_rows:
    print(f"\nRows with notes ({len(notes_rows)}):")
    for r in notes_rows:
        print(f"  {r['phytozome_code']:6}  {r['uniprot_mnemonic']:8}  {r['scientific_name']:30}  [{r['notes']}]")
