# Cell Filtering — Low-RAM Execution Guide

Use this two-step workflow when the full QC checkpoint (`checkpoint_02_seurat_qc.rds`)
cannot be loaded into memory all at once (typically on machines with ≤ 16 GB RAM).

## Why this is needed

`checkpoint_02_seurat_qc.rds` is saved as a single list of all Seurat objects.
Loading it with `readRDS()` deserialises the entire list at once — a 2.6 GB file
can expand to ~12–15 GB in RAM. The standard `03_filter_cells.sh` runs under a
13 GB cgroup limit and is killed before any filtering begins.

The solution splits the list into individual per-library files so that step 03
loads one library at a time, keeping peak RAM to roughly the size of a single
Seurat object (~1–2 GB).

---

## Step 1 — Split checkpoint_02 into per-library files (run once)

```bash
bash code/03a_split_qc_checkpoint.sh
```

- Runs under a **22 GB** `systemd` memory limit (enough for the full list).
- Reads `checkpoint_02_seurat_qc.rds` once, then saves each library separately:
  - `RData/qc_libs/qc_<lib_name>.rds` — one file per library
  - `RData/lib_names.rds` — ordered list of library names
- Nulls each slot after saving, so RAM is freed incrementally.
- Only needs to run **once**. If `qc_libs/qc_<lib>.rds` already exists it is skipped.

---

## Step 2 — Filter cells one library at a time

```bash
bash code/03_filter_cells.sh
```

The updated `03_filter_cells.R` now:

1. Reads `lib_names.rds` to get the list of libraries (negligible RAM).
2. Iterates sequentially: for each library —
   - loads `qc_libs/qc_<lib>.rds`
   - applies QC thresholds (see `global_variables.R`)
   - saves `filtered_libs/filtered_<lib>.rds`
   - calls `rm()` + `gc()` before loading the next library
3. Assembles all filtered libraries into a named list and saves
   `checkpoint_03_seurat_filtered.rds`.

Peak RAM during filtering ≈ size of one Seurat object. The 13 GB cgroup limit
is sufficient.

---

## QC thresholds applied (set in `code/global_variables.R`)

| Parameter | Variable | Default |
|---|---|---|
| Min genes per cell | `QC_MIN_FEATURES` | 200 |
| Max genes per cell | `QC_MAX_FEATURES` | 5 000 |
| Min UMI counts | `QC_MIN_COUNTS` | 500 |
| Max UMI counts | `QC_MAX_COUNTS` | 25 000 |
| Max % mitochondrial | `QC_MAX_MT` | 2 % |
| Min log10 genes/UMI | `QC_MIN_COMPLEXITY` | see `global_variables.R` |

Demuxlet singlet filtering (Layer 1) is available but commented out in
`03_filter_cells.R`. Uncomment the `DemuxletType` block if your experiment
includes multiplexed samples.

---

## Output files

```
RData/
├── lib_names.rds                        # library name vector (step 1)
├── qc_libs/
│   ├── qc_<lib1>.rds                    # pre-filter Seurat objects (step 1)
│   └── qc_<lib2>.rds
├── filtered_libs/
│   ├── filtered_<lib1>.rds              # post-filter Seurat objects (step 2)
│   └── filtered_<lib2>.rds
└── checkpoint_03_seurat_filtered.rds    # assembled list (step 2)
```

## Resuming after interruption

Both scripts are idempotent. If the process is killed mid-run:

- **Step 1**: re-run `03a_split_qc_checkpoint.sh` — already-saved `qc_libs/` files are skipped.
- **Step 2**: re-run `03_filter_cells.sh` — already-saved `filtered_libs/` files are skipped, and the final checkpoint is only written once all libraries are present.
