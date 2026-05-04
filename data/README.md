## genes_dedup.csv location

Put `genes_dedup.csv` in this folder so other scripts can always reference:

- `1_GeneOptimization/data/genes_dedup.csv`

### Option A: copy the file here

```bash
cp "/Users/user/Downloads/genes_dedup.csv" "1_GeneOptimization/data/genes_dedup.csv"
```

### Option B: symlink (keeps a single copy)

```bash
ln -sf "/Users/user/Downloads/genes_dedup.csv" "1_GeneOptimization/data/genes_dedup.csv"
```

