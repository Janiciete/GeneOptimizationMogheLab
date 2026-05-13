#!/usr/bin/env python3
"""Phase 4: Build SQLite index of Phytozome protein sequences for fast lookup."""

import argparse
import re
import sqlite3
import sys
from pathlib import Path

# ── Constants ──────────────────────────────────────────────────────────────────
PHYTOZOME_DIR = Path("/local/storage/0_databases/2_protein/2_phytozome")
PROTEOME_TAB  = PHYTOZOME_DIR / "proteome_files.tab"
OUT_DB        = Path("/local/storage/jedric/1_GeneOptimization/indexes/phytozome_index.db")

BATCH_SIZE = 5_000   # rows per executemany call

# ── Helpers ────────────────────────────────────────────────────────────────────
def load_species_map(tab_path):
    """Return {filename -> scientific_name} from proteome_files.tab."""
    species_map = {}
    with open(tab_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                fname    = parts[0].strip()
                sci_name = parts[1].strip().replace("_", " ")
                species_map[fname] = sci_name
    return species_map

def compute_base_gene_id(gene_id):
    """
    Compute a stable base gene ID by removing isoform/protein suffixes.
    Conservative rules (applied in order):
      1. Strip trailing .p or .P  (Phytozome protein marker)
      2. Strip trailing _P<digits> (maize-style protein isoform: _P01, _P02)
      3. Repeatedly strip trailing .<1-2 digits>  (isoform numbers: .1, .1.1, .3)
         Only 1-2 digit numbers to avoid touching IDs like supercontig_0.100
    gene_id is never modified; base_gene_id is a separate derived value.
    """
    s = gene_id
    if len(s) > 2 and s[-2:].lower() == ".p":
        s = s[:-2]
    s = re.sub(r"_P\d+$", "", s)
    while True:
        m = re.match(r"^(.+)\.\d{1,2}$", s)
        if m:
            s = m.group(1)
        else:
            break
    return s

def normalize(name):
    """Uppercase and strip for index lookups."""
    return name.strip().upper()


# ── Database setup ─────────────────────────────────────────────────────────────
def setup_db(db_path):
    con = sqlite3.connect(db_path)
    cur = con.cursor()

    # Speed pragmas for one-time index build
    cur.executescript("""
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous  = OFF;
        PRAGMA cache_size   = -512000;
        PRAGMA temp_store   = MEMORY;
    """)

    cur.execute("""
        CREATE TABLE IF NOT EXISTS phytozome_sequences (
            id                      INTEGER PRIMARY KEY AUTOINCREMENT,
            phytozome_code          TEXT,
            scientific_name         TEXT,
            raw_header              TEXT,
            gene_id                 TEXT,
            gene_id_normalized      TEXT,
            base_gene_id            TEXT,
            base_gene_id_normalized TEXT,
            sequence                TEXT,
            sequence_length         INTEGER,
            source_file             TEXT
        )
    """)

    cur.execute("CREATE INDEX IF NOT EXISTS idx_pcode          ON phytozome_sequences(phytozome_code)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_gene_id_norm   ON phytozome_sequences(gene_id_normalized)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_base_gene_norm ON phytozome_sequences(base_gene_id_normalized)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_pcode_gene     ON phytozome_sequences(phytozome_code, gene_id_normalized)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_pcode_base     ON phytozome_sequences(phytozome_code, base_gene_id_normalized)")

    con.commit()
    return con


# ── FASTA streaming ────────────────────────────────────────────────────────────
def stream_fasta(filepath, phytozome_code, scientific_name, record_limit=None):
    """
    Yield one tuple per FASTA record, streaming line by line.
    Tuple column order matches the INSERT below.

    Header handling:
      - If the header starts with "{phytozome_code}_", strip that prefix for gene_id.
      - Otherwise keep the header as-is (e.g. Slyc, Pabi, Ptae, Afil, Scuc).
    Sequence:
      - Concatenates continuation lines; strips trailing '*'.
    """
    source_file = filepath.name
    prefix      = phytozome_code + "_"
    records_yielded = 0

    current_header = None
    seq_parts      = []

    def make_tuple():
        raw = current_header
        gene_id = raw[len(prefix):] if raw.startswith(prefix) else raw
        seq     = "".join(seq_parts).rstrip("*")
        base_id = compute_base_gene_id(gene_id)
        return (
            phytozome_code,
            scientific_name,
            raw,
            gene_id,
            normalize(gene_id),
            base_id,
            normalize(base_id),
            seq,
            len(seq),
            source_file,
        )

    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    yield make_tuple()
                    records_yielded += 1
                    if record_limit is not None and records_yielded >= record_limit:
                        return
                current_header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)

        # Emit the final record
        if current_header is not None:
            if record_limit is None or records_yielded < record_limit:
                yield make_tuple()


# ── Core processing ────────────────────────────────────────────────────────────
INSERT_SQL = """
    INSERT INTO phytozome_sequences
        (phytozome_code, scientific_name, raw_header, gene_id,
         gene_id_normalized, base_gene_id, base_gene_id_normalized,
         sequence, sequence_length, source_file)
    VALUES (?,?,?,?,?,?,?,?,?,?)
"""

def process_file(fasta_path, phytozome_code, scientific_name, con, record_limit=None):
    """Parse one FASTA file and batch-insert into SQLite. Returns count inserted."""
    cur     = con.cursor()
    batch   = []
    inserted = 0

    for row_tuple in stream_fasta(fasta_path, phytozome_code, scientific_name, record_limit):
        batch.append(row_tuple)
        inserted += 1

        if len(batch) >= BATCH_SIZE:
            cur.executemany(INSERT_SQL, batch)
            batch.clear()
            con.commit()

    if batch:
        cur.executemany(INSERT_SQL, batch)
        con.commit()

    return inserted


# ── Sample queries ─────────────────────────────────────────────────────────────
def run_sample_queries(con):
    cur = con.cursor()
    print("\n=== SAMPLE QUERIES ===")
    for code in ["Atha", "Slyc", "Osat", "Zmay", "Gmax"]:
        cur.execute(
            "SELECT gene_id, base_gene_id, sequence_length "
            "FROM phytozome_sequences WHERE phytozome_code = ? LIMIT 2",
            (code,)
        )
        rows = cur.fetchall()
        if rows:
            for r in rows:
                print(f"  {code}  gene_id={r[0]}  base={r[1]}  len={r[2]}")
        else:
            print(f"  {code}: not in index")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Phase 4: Build Phytozome SQLite sequence index"
    )
    parser.add_argument("--test-files", type=int, metavar="N",
                        help="Process only the first N FASTA files (alphabetical order)")
    parser.add_argument("--test-records", type=int, metavar="N",
                        help="Process only the first N records per file")
    parser.add_argument("--full", action="store_true",
                        help="Process all *_prot.fa files completely")
    parser.add_argument("--overwrite", action="store_true",
                        help="Delete and rebuild an existing database")
    args = parser.parse_args()

    if not args.full and args.test_files is None:
        parser.error("Specify --full, or --test-files N (optionally with --test-records N)")

    # ── Safety check ────────────────────────────────────────────────────────────
    if OUT_DB.exists():
        if args.overwrite:
            print(f"WARNING: Overwriting existing database: {OUT_DB}")
            OUT_DB.unlink()
        else:
            print(f"ERROR: {OUT_DB} already exists. Use --overwrite to replace it.")
            sys.exit(1)

    OUT_DB.parent.mkdir(parents=True, exist_ok=True)

    # ── Load species map ────────────────────────────────────────────────────────
    species_map = load_species_map(PROTEOME_TAB)
    print(f"Loaded {len(species_map)} species from {PROTEOME_TAB.name}")

    # ── Discover FASTA files ────────────────────────────────────────────────────
    fasta_files = sorted(PHYTOZOME_DIR.glob("*_prot.fa"))
    if not fasta_files:
        print(f"ERROR: No *_prot.fa files found in {PHYTOZOME_DIR}")
        sys.exit(1)

    if args.test_files is not None:
        fasta_files = fasta_files[: args.test_files]

    record_limit = None if args.full else args.test_records

    print(f"\nPhase 4: Building Phytozome SQLite sequence index")
    print(f"  FASTA directory: {PHYTOZOME_DIR}")
    print(f"  Output DB:       {OUT_DB}")
    print(f"  Files to process: {len(fasta_files)}")
    if args.full:
        print(f"  Mode: FULL")
    else:
        print(f"  Mode: TEST  (files={args.test_files}, records/file={record_limit or 'all'})")

    con = setup_db(OUT_DB)

    # ── Process each FASTA file ─────────────────────────────────────────────────
    total_inserted = 0
    for fasta_path in fasta_files:
        pcode    = fasta_path.name.replace("_prot.fa", "")
        sci_name = species_map.get(fasta_path.name, "")
        inserted = process_file(fasta_path, pcode, sci_name, con, record_limit=record_limit)
        total_inserted += inserted
        print(f"  [{pcode:6}]  {fasta_path.name:25}  {sci_name:35}  {inserted:>7,} seqs")

    # ── Final report ────────────────────────────────────────────────────────────
    print(f"\n=== FINAL REPORT ===")
    print(f"  FASTA files processed: {len(fasta_files)}")
    print(f"  Total sequences:       {total_inserted:,}")
    print(f"  Output database:       {OUT_DB}")

    run_sample_queries(con)
    con.close()
    print("\nDone.")


if __name__ == "__main__":
    main()
