# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## CRITICAL: File Modification Safety Rules

**NEVER modify a file that was not explicitly named in the user's request.**

This is non-negotiable. Violations have already caused loss of work (deleted lines from `.qmd` notebooks that had to be manually reconstructed).

Specific rules:
- **Only touch the file(s) the user explicitly asks to edit.** If a task involves `00_packages.R`, do not open or write to any `.Rmd`, `.qmd`, or other script.
- **Never reformat, clean up, reorganize, or "optimize" code in a file unless the user asks for it explicitly.**
- **For `.Rmd` and `.qmd` notebooks**: only modify the specific chunk or section requested. Never restructure YAML headers, chunk options, prose, or surrounding code.
- **When in doubt, ask** which file to edit rather than making assumptions and touching multiple files.
- **Read before writing**: always read the full file content before making any edit, and use targeted `Edit` operations (never `Write` to overwrite an existing notebook).

---

## What this repo is

End-to-end snRNA-seq analysis pipeline built on **Seurat v5** (R). Primary use case: multi-sample WT vs KO mouse brain experiments. The `code_claude/` subdirectory contains a parallel adaptation for a human PBMC CITE-seq public dataset (GSE194315) that serves as a reference/example.

## Running the pipeline

There is no build or test runner. Work is done through RStudio / R interactively:

```r
# Restore the full R environment first (R 4.5.2, renv-managed)
renv::restore()

# Render the main analysis notebook
rmarkdown::render("rmds/Single_Cell_10X_Integrated_functions_SCT - PV_Cre-chacon22.Rmd")

# Run individual enrichment modules
source("code/03_GO.R")       # GO over-representation
source("code/05_gse.R")      # GSEA (GO + KEGG ranked lists)
source("code/09_gseKEGG.R")  # KEGG pathway analysis + pathview
source("code/04_strings.R")  # STRING PPI networks
```

Logs from script runs land in `logs/` (gitignored pattern: `*_<timestamp>.log`).

## Architecture

The pipeline is **modular and script-sourced**: RMarkdown notebooks orchestrate the analysis by `source()`-ing numbered R scripts from `code/`. Scripts are independent and can be run standalone after the upstream Seurat object exists.

**Execution order:**
1. `global_variables.R` — must be sourced first; sets all thresholds and organism parameters
2. `00_packages.R` — installs/loads dependencies via `pak`
3. `01_sc_functions.R` — QC utilities (`library_summary`, `generate_qc_plots`, `save_plot`)
4. Notebook: QC → SCTransform → RPCA integration → PCA → UMAP → Louvain clustering
5. `02_vulcano_plots.R` — downstream DE visualisation
6. `03_GO.R` / `05_gse.R` / `09_gseKEGG.R` / `08_EnrichR.R` — enrichment modules (run after DE)

**Key entry points:**

| Path | Purpose |
|---|---|
| `rmds/snRNAseq_pipeline.qmd` | **Active primary notebook** — new Quarto pipeline being built from scratch |
| `rmds/Single_Cell_10X_Integrated_functions_SCT - PV_Cre-chacon22.Rmd` | Legacy Rmd template for multi-sample experiments (SCT + RPCA) |
| `rmds/old/` | Archive of pre-SCTransform and experimental notebooks — do not modify |
| `code_claude/GSE194315_PBMC_SCT_Analysis.Rmd` | Human PBMC reference analysis (full annotated example) |

**For a new dataset:** use `rmds/snRNAseq_pipeline.qmd` as the main entry point. Set `output_path` and `data_path` in the setup chunk, tune thresholds in `global_variables.R`.

## Configuration

All analysis parameters live in `code/global_variables.R` (mouse) or `code_claude/global_variables_GSE194315.R` (human PBMC). **Never hardcode thresholds in notebooks.**

Key parameters to adjust per experiment:
- `kegg_organism` / `species` / `organism` / `keyType` — change together when switching species
- `QC_MAX_MT`: brain snRNA-seq → 1–5%; PBMC → 20%
- `PARALLEL_WORKERS` + `plan()` — set to 1/sequential for ≤16 GB RAM

## Checkpoint system (code_claude only)

`01_sc_functions_GSE194315.R` implements crash recovery via `save_checkpoint()` / `load_checkpoint()` / `load_or_run()`. Checkpoints are saved as `.rds` files under `output/RData/checkpoint_<step>.rds`. Use `load_or_run("03_filtered", compute_fn)` pattern in notebooks to skip already-completed steps on restart.

## Output structure

`create_directories(output_path)` in `global_variables.R` creates the full output tree on load:
```
output/<experiment>/
├── figures/        # enrichGO/, gseGO/, KEGG/, string/, panther/ subdirs
├── tables/
└── RData/          # Seurat objects, checkpoints
```

Figures are saved in dual format: TIFF (300 dpi) + PDF via `save_plot()`.

## Data

- `rawdata/` — gitignored; expects 10X CellRanger `filtered_feature_bc_matrix/` per sample
- `output/` — gitignored; all generated figures, tables, and Seurat objects
- `renv.lock` — fully specifies the R 4.5.2 environment; always commit changes to this file

## Available skills

### `single-cell-rna-qc-v1.0.0` — scverse-based QC pipeline (Python/AnnData)

Invoke via `/single-cell-rna-qc-v1.0.0` when the user wants QC on `.h5ad` or 10X `.h5` files following MAD-based scverse best practices. This skill is **complementary** to the Seurat R pipeline in this repo: use it when working with AnnData-format data or when the user explicitly wants a Python/scanpy QC workflow.

**Quick usage:**
```bash
python3 scripts/qc_analysis.py input.h5ad          # AnnData input
python3 scripts/qc_analysis.py raw_feature_bc_matrix.h5  # 10X .h5 input
```

Key differences from the Seurat R pipeline in this repo:
- Uses **MAD-based outlier detection** (adaptive per-dataset) vs. fixed thresholds in `global_variables.R`
- Outputs filtered `.h5ad` files; to continue in R use `SeuratDisk::Convert()` or `anndata::read_h5ad()`
- MT gene prefix: `mt-` for mouse, `MT-` for human (vs. regex in `PercentageFeatureSet`)

**Species-specific parameters to set:**
- Mouse: `--mt-pattern "^mt-"` (lowercase)
- Human: `--mt-pattern "^MT-"` (uppercase, default)

## Common issues

- **Seurat v5 assay errors**: `options(Seurat.object.assay.version = "v3")` or `as(obj[["RNA"]], "Assay")`
- **OOM during integration**: reduce `N_INTEGRATION_FEATURES` (3000 → 1500) or switch to `plan("sequential")`
- **clusterProfiler key mismatch**: `keyType` must match gene ID format in data (mouse = `"UNIPROT"`, human PBMC = `"SYMBOL"`)
- **renv restore fails for Seurat v5**: install from GitHub first: `remotes::install_github("satijalab/seurat", "seurat5")`
