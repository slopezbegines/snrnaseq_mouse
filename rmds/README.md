# RMarkdown Notebooks — Usage Guide

This directory contains all analysis notebooks. They share the modular scripts in `../code/` and `../code_claude/`.

---

## Which notebook to use

| Notebook | Type | When to use |
|---|---|---|
| `Single_Cell_10X_Integrated_functions_SCT - PV_Cre-chacon22.Rmd` | ⭐ **Template** | Full snRNA-seq pipeline (SCTransform + RPCA integration). **Start here for any new multi-sample experiment.** Covers QC → normalisation → integration → UMAP → clustering → DE. |
| `Single_Cell_10X_Integrated_functions_SCT -UBC_Cre.Rmd` | ⭐ **Template** | Same pipeline adapted for UBC-Cre CSP-Flox experimental model. Use as secondary template reference. |
| `Clustering Association_FindAllMarkers.Rmd` | Downstream | Cluster annotation (SingleR + Azimuth) and marker gene discovery with `FindAllMarkers`. Run after the main pipeline has produced a clustered Seurat object. |
| `Clustering Association.Rmd` | Downstream | Clustering association analysis without `FindAllMarkers`. Lighter alternative for quick annotation checks. |
| `Single_Cell_10X_Integrated_functions.Rmd` | Legacy | Pre-SCTransform version using log-normalisation. Superseded by the SCT notebooks. Kept for reference only. |
| `Single_Cell_10X_Integrated_functions_UBC.Rmd` | Legacy | UBC-Cre variant of the log-normalisation pipeline. Superseded. |
| `Single_Cell_10X_Integrated.Rmd` | Legacy | Older integrated pipeline without modular functions. Superseded. |
| `Single_Cell_10X_merged.Rmd` | Legacy | Simple sample merge without batch correction or integration. Use only when samples have no batch effects. |
| `Single_Cell_10X_SCTransform.Rmd` | Legacy | SCTransform applied to a single sample. No multi-sample integration. |

---

## Recommended workflow for a new dataset

```
1. Copy the SCT template:
   cp "Single_Cell_10X_Integrated_functions_SCT - PV_Cre-chacon22.Rmd" \
      "MyExperiment_SCT.Rmd"

2. Edit the setup chunk (~lines 20–40):
   - Set output_path to your output directory
   - Point data loading to your CellRanger filtered_feature_bc_matrix/ folders
   - Adjust QC thresholds in ../code/global_variables.R (or ../code_claude/global_variables_GSE194315.R as reference)

3. Run Clustering Association_FindAllMarkers.Rmd after clustering is complete.
```

---

## GSE194315 adaptation (`code_claude/`)

The `code_claude/` directory contains scripts adapted for the **GSE194315 PBMC CITE-seq** reference dataset (human, 3 conditions: Healthy / PSA / PSO). These scripts demonstrate how to adapt the pipeline to a new organism and dataset type, and serve as an annotated example for clients.

| Script | Purpose |
|---|---|
| `global_variables_GSE194315.R` | QC thresholds, organism parameters, hardware config, PBMC marker genes |
| `00_packages_GSE194315.R` | Dependency loading with auto-install helpers |
| `01_sc_functions_GSE194315.R` | Checkpoint system, crash recovery, extended QC utilities |
| `GSE194315_PBMC_SCT_Analysis.Rmd` | Full end-to-end analysis notebook for the GSE194315 dataset |

---

## Notes

- All notebooks expect `../code/global_variables.R` to be sourced first.
- Outputs are written to `../output/<experiment_name>/` (gitignored).
- Raw 10X data (`rawdata/`) is gitignored — see the main README for download instructions.
