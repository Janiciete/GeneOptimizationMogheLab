#!/usr/bin/env python3
"""Phase 3: Build SQLite index of UniProt plant gene/protein names for fast lookup."""

import argparse
import csv
import sqlite3
import sys
from pathlib import Path

# ── Constants ──────────────────────────────────────────────────────────────────
SPROT_FILE  = Path("/local/storage/0_databases/2_protein/3_uniprot_2024/uniprot_sprot_plants_output.txt.geneNames")
TREMBL_FILE = Path("/local/storage/0_databases/2_protein/3_uniprot_2024/uniprot_trembl_plants_output.txt.geneNames")
OUT_DB      = Path("/local/storage/jedric/1_GeneOptimization/indexes/uniprot_index.db")

BATCH_SIZE     = 10_000   # rows per executemany call
COMMIT_EVERY   = 100_000  # rows read between commits / progress prints

# ── Helpers ────────────────────────────────────────────────────────────────────
def parse_accession(raw):
    """Return first accession before the first semicolon."""
    return raw.split(";")[0].strip()

def parse_accession_description(raw_accession_field, description_col):
    """
    Use the Description column if non-empty; otherwise fall back to the text
    after the first semicolon in the Accession field (legacy layout).
    """
    if description_col:
        return description_col
    parts = raw_accession_field.split(";", 1)
    return parts[1].strip() if len(parts) > 1 else ""

def parse_organism_mnemonic(entry_name):
    """Extract organism mnemonic from suffix after last underscore. e.g. 10HGO_CATRO -> CATRO."""
    if "_" in entry_name:
        return entry_name.rsplit("_", 1)[-1]
    return ""

def parse_gene_names(raw):
    """
    Return list of gene name strings from the Gene_Name field.
    Splits on ';', strips whitespace, filters ECO evidence tokens ({ECO:...}).
    """
    if not raw or not raw.strip():
        return []
    tokens = [t.strip() for t in raw.split(";")]
    return [t for t in tokens if t and not t.startswith("{")]

def normalize_gene_name(name):
    """Strip whitespace, remove surrounding quotes, uppercase for lookup."""
    return name.strip().strip('"').strip("'").upper()

def get_col(row, col_index, key):
    """Safely get a named column from a raw csv.reader row list; return '' if absent."""
    idx = col_index.get(key)
    if idx is None or idx >= len(row):
        return ""
    return (row[idx] or "").strip()


# ── Database setup ─────────────────────────────────────────────────────────────
def setup_db(db_path):
    con = sqlite3.connect(db_path)
    cur = con.cursor()

    # Speed pragmas — safe for a one-time generated index; do not use on transactional DBs
    cur.executescript("""
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous  = OFF;
        PRAGMA cache_size   = -1000000;
        PRAGMA temp_store   = MEMORY;
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS uniprot_entries (
            id                   INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_name            TEXT,
            gene_name_normalized TEXT,
            accession            TEXT,
            entry_name           TEXT,
            organism_mnemonic    TEXT,
            source_tier          TEXT,
            description          TEXT,
            refseq               TEXT,
            ensembl_plants       TEXT,
            tair                 TEXT,
            gramene              TEXT,
            aa_count             TEXT
        )
    """)

    cur.execute("CREATE INDEX IF NOT EXISTS idx_gene_name_norm ON uniprot_entries(gene_name_normalized)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_accession      ON uniprot_entries(accession)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_organism       ON uniprot_entries(organism_mnemonic)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_source_tier    ON uniprot_entries(source_tier)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_gene_organism  ON uniprot_entries(gene_name_normalized, organism_mnemonic)")

    con.commit()
    return con


# ── File streaming ─────────────────────────────────────────────────────────────
def open_genenames_file(filepath):
    """
    Open a .geneNames file, skip the leading '#' comment line, read the header,
    and return (file_handle, csv.reader, col_index).

    Gene names are physically appended as the last field of each row regardless of
    how many DB cross-reference columns precede them.  Rows range from ~43 to 106
    fields; only those with exactly 106 fields would align with the 'Gene_Name'
    header column.  Use row[-1] for the gene name and col_index for everything else.
    """
    fh = open(filepath, newline="", encoding="utf-8")
    reader = csv.reader(fh, delimiter="\t")
    first = next(reader)
    if first and first[0].startswith("#"):
        header = next(reader)
    else:
        header = first
    col_index = {name: i for i, name in enumerate(header)}
    return fh, reader, col_index


# ── Core processing ────────────────────────────────────────────────────────────
def process_file(filepath, source_tier, con, limit=None):
    """Stream a .geneNames file, parse rows, insert into SQLite.
    Returns (rows_read, rows_with_genes, rows_inserted, rows_skipped).
    rows_read     = total source rows processed
    rows_with_genes = source rows that had at least one gene name (rows_read - rows_skipped)
    rows_inserted = total DB rows written (>= rows_with_genes due to synonym expansion)
    rows_skipped  = source rows with no usable gene name
    """
    cur = con.cursor()
    batch = []
    rows_read = rows_with_genes = rows_inserted = rows_skipped = 0

    fh, reader, col_index = open_genenames_file(filepath)
    try:
        for row in reader:
            if not row:
                continue
            if limit is not None and rows_read >= limit:
                break

            rows_read += 1

            entry_name      = get_col(row, col_index, "Protein_Name")
            raw_accession   = get_col(row, col_index, "Accession")
            accession       = parse_accession(raw_accession)
            description     = parse_accession_description(raw_accession, get_col(row, col_index, "Description"))
            organism_mnem   = parse_organism_mnemonic(entry_name)
            refseq          = get_col(row, col_index, "RefSeq")
            ensembl_plants  = get_col(row, col_index, "EnsemblPlants")
            tair            = get_col(row, col_index, "TAIR")
            gramene         = get_col(row, col_index, "Gramene")
            aa_count        = get_col(row, col_index, "Sequence_AA_Count")

            # Gene name is always the physical last field, appended after however
            # many DB cross-reference columns the row happens to have.
            gene_names = parse_gene_names(row[-1])

            if not gene_names:
                rows_skipped += 1
                continue

            rows_with_genes += 1
            for gname in gene_names:
                gname_norm = normalize_gene_name(gname)
                if not gname_norm:
                    continue
                batch.append((
                    gname, gname_norm, accession, entry_name,
                    organism_mnem, source_tier, description,
                    refseq, ensembl_plants, tair, gramene, aa_count
                ))
                rows_inserted += 1

            if len(batch) >= BATCH_SIZE:
                cur.executemany("""
                    INSERT INTO uniprot_entries
                        (gene_name, gene_name_normalized, accession, entry_name,
                         organism_mnemonic, source_tier, description,
                         refseq, ensembl_plants, tair, gramene, aa_count)
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
                """, batch)
                batch.clear()

            if rows_read % COMMIT_EVERY == 0:
                con.commit()
                print(f"  [{source_tier}] {rows_read:>10,} read | {rows_with_genes:>10,} with gene | {rows_inserted:>10,} DB rows | {rows_skipped:>7,} skipped")

    finally:
        fh.close()

    # Flush remaining batch
    if batch:
        cur.executemany("""
            INSERT INTO uniprot_entries
                (gene_name, gene_name_normalized, accession, entry_name,
                 organism_mnemonic, source_tier, description,
                 refseq, ensembl_plants, tair, gramene, aa_count)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
        """, batch)
    con.commit()

    return rows_read, rows_with_genes, rows_inserted, rows_skipped


# ── Sample queries ─────────────────────────────────────────────────────────────
def run_sample_queries(con):
    cur = con.cursor()
    test_genes = ["CYP76AH1", "RBCS1", "CHS", "PAL", "ADH1", "10HGO", "GRF1"]
    print("\n=== SAMPLE QUERIES ===")
    for gene in test_genes:
        norm = normalize_gene_name(gene)
        cur.execute(
            "SELECT gene_name, accession, organism_mnemonic, source_tier "
            "FROM uniprot_entries WHERE gene_name_normalized = ? LIMIT 3",
            (norm,)
        )
        hits = cur.fetchall()
        if hits:
            first = hits[0]
            print(f"  {gene:<12}  {len(hits)} hit(s)  →  acc={first[1]}  org={first[2]}  [{first[3]}]")
        else:
            print(f"  {gene:<12}  no hits")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Phase 3: Build UniProt SQLite index for plant gene/protein lookup"
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--test-limit", type=int, metavar="N",
        help="Process only the first N data rows from each file (safe test mode)"
    )
    mode.add_argument(
        "--full", action="store_true",
        help="Process all rows from both Swiss-Prot and TrEMBL files"
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Overwrite the existing database if it exists"
    )
    args = parser.parse_args()

    limit = None if args.full else args.test_limit

    # ── Safety check ────────────────────────────────────────────────────────────
    if OUT_DB.exists():
        if args.overwrite:
            print(f"WARNING: Overwriting existing database: {OUT_DB}")
            OUT_DB.unlink()
        else:
            print(f"ERROR: {OUT_DB} already exists.")
            print("Use --overwrite to replace it, or delete it manually.")
            sys.exit(1)

    OUT_DB.parent.mkdir(parents=True, exist_ok=True)

    print("Phase 3: Building UniProt SQLite index")
    print(f"  Swiss-Prot:  {SPROT_FILE}")
    print(f"  TrEMBL:      {TREMBL_FILE}")
    print(f"  Output DB:   {OUT_DB}")
    print(f"  Mode:        {'TEST (limit {:,} rows per file)'.format(limit) if limit else 'FULL'}")

    con = setup_db(OUT_DB)

    # ── Swiss-Prot ──────────────────────────────────────────────────────────────
    print(f"\n[Swiss-Prot] Streaming {SPROT_FILE.name} ...")
    sp_read, sp_with_genes, sp_inserted, sp_skipped = process_file(SPROT_FILE, "sprot", con, limit=limit)
    print(f"[Swiss-Prot] Done: {sp_read:,} read | {sp_with_genes:,} with gene name | {sp_inserted:,} DB rows | {sp_skipped:,} skipped")

    # ── TrEMBL ─────────────────────────────────────────────────────────────────
    print(f"\n[TrEMBL] Streaming {TREMBL_FILE.name} ...")
    tr_read, tr_with_genes, tr_inserted, tr_skipped = process_file(TREMBL_FILE, "trembl", con, limit=limit)
    print(f"[TrEMBL] Done: {tr_read:,} read | {tr_with_genes:,} with gene name | {tr_inserted:,} DB rows | {tr_skipped:,} skipped")

    # ── Final report ────────────────────────────────────────────────────────────
    total_read       = sp_read + tr_read
    total_with_genes = sp_with_genes + tr_with_genes
    total_inserted   = sp_inserted + tr_inserted
    total_skipped    = sp_skipped + tr_skipped

    print("\n=== FINAL REPORT ===")
    print(f"                          {'Swiss-Prot':>12}  {'TrEMBL':>12}  {'Total':>12}")
    print(f"  Source rows read:       {sp_read:>12,}  {tr_read:>12,}  {total_read:>12,}")
    print(f"  Rows with gene name:    {sp_with_genes:>12,}  {tr_with_genes:>12,}  {total_with_genes:>12,}")
    print(f"  Skipped (no gene name): {sp_skipped:>12,}  {tr_skipped:>12,}  {total_skipped:>12,}")
    print(f"  DB rows inserted:       {sp_inserted:>12,}  {tr_inserted:>12,}  {total_inserted:>12,}")
    print(f"  (DB rows > source rows with genes when synonyms exist)")
    print(f"\nOutput database: {OUT_DB}")

    run_sample_queries(con)
    con.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
