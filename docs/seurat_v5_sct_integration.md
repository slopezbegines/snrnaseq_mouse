# Seurat v5 SCTransform + RPCA Integration — Workflow Reference

## Overview

This document explains the Seurat v5 native integration workflow implemented in
`code/06_integration.R`, why it replaces the v4 approach, and the key conceptual
differences a user needs to understand before running downstream analysis.

---

## 1. v4 vs v5: What Changed and Why

### Seurat v4 workflow (replaced)

```
SCTransform (per object)
  → PrepSCTIntegration
    → FindIntegrationAnchors
      → IntegrateData
          └─ creates a new "integrated" Assay (stores corrected expression values)
```

The corrected expression matrix was stored in a dedicated assay slot. All
downstream steps (UMAP, clustering, DE) ran on the `integrated` assay.

**Problem:** `IntegrateData` builds the full anchor weight matrix for all
samples simultaneously. On datasets with many cells this exhausts RAM before
completion.

### Seurat v5 workflow (current)

```
merge()
  → SCTransform (per-layer / per-sample)
    → RunPCA
      → IntegrateLayers
          └─ creates a new "integrated.rpca" REDUCTION (not an assay)
```

`IntegrateLayers` processes one pair of datasets at a time, releasing memory
between iterations. The correction is stored as a low-dimensional embedding
(reduction), not as a full expression matrix.

**Key consequence:** there is no `integrated` assay in the v5 workflow. The
`SCT` assay is used for all gene-expression tasks; `integrated.rpca` is used
only for graph construction and embedding.

---

## 2. The Seurat v5 Layer System

In Seurat v5, `Assay5` replaces the v4 `Assay` class. The core change:

| | Seurat v4 `Assay` | Seurat v5 `Assay5` |
|---|---|---|
| Data storage | Fixed slots: `@counts`, `@data`, `@scale.data` | Named layers: any number, any name |
| Multi-sample | One matrix (cells from all samples concatenated) | Separate layers per sample: `counts.sample1`, `counts.sample2`, … |
| Normalization | Must be applied to entire concatenated matrix | Applied per-layer (per-sample) independently |

When `merge()` is called on a list of Seurat v5 objects, the resulting `RNA`
assay automatically gets one layer per source object (e.g.
`counts.Ctrl_forebrain_1`, `counts.Ctrl_forebrain_2`). This split is what
allows `SCTransform` to fit a separate model per sample.

**Note on `SCTAssay`:** `SCTransform` always produces an `SCTAssay` object,
which inherits from the v4 `Assay` class (not `Assay5`). This is intentional:
`SCTAssay` needs to carry per-cell regression model parameters in
`@SCTModel.list`, a slot that has no equivalent in `Assay5`. The
`Layer counts isn't present` warning seen in v4-style `IntegrateData` calls is
a consequence of mixing APIs; it does not appear with `IntegrateLayers`.

---

## 3. Step-by-Step Rationale

### Step 1 — `merge()`

```r
merged <- merge(x = seurat_singlets[[1]], y = seurat_singlets[-1],
                add.cell.ids = lib_names)
```

- Concatenates all filtered/singlet objects into a single Seurat object.
- `add.cell.ids` prefixes each barcode with the library name, preventing
  collisions across samples.
- The resulting `RNA` `Assay5` has layers `counts.<lib>` and `data.<lib>` for
  every source object.

### Step 2 — `SCTransform()` on the merged object

```r
merged <- SCTransform(merged, vst.flavor = "v2", method = "glmGamPoi",
                      vars.to.regress = "percent.MT",
                      variable.features.n = N_INTEGRATION_FEATURES)
```

- Seurat v5 `SCTransform` detects split layers and fits one NB regression model
  per sample independently.
- `vst.flavor = "v2"` uses regularised regression (recommended since Seurat
  4.2 / SCTransform v2).
- `method = "glmGamPoi"` is substantially faster than the default Poisson
  method on large datasets.
- `vars.to.regress = "percent.MT"` removes mitochondrial content as a
  confounding factor from the residuals.
- `variable.features.n` controls how many variable features the SCT model
  retains. Reducing from the 3000 default limits peak RAM during `scale.data`
  construction.

### Step 3 — `RunPCA()`

```r
merged <- RunPCA(merged, npcs = N_PCS_MAX)
```

- PCA is run on the `SCT` assay `scale.data` (the per-sample-normalised,
  per-feature-scaled residuals).
- This produces the `pca` reduction, which is the starting point for RPCA
  integration in the next step.

### Step 4 — `IntegrateLayers()` with RPCA

```r
merged <- IntegrateLayers(
  object               = merged,
  method               = RPCAIntegration,
  orig.reduction       = "pca",
  normalization.method = "SCT",
  dims                 = 1:N_PCS_INTEGRATION,
  verbose              = FALSE
)
```

- `RPCAIntegration`: applies Reciprocal PCA, a reciprocal projection approach
  that is faster and more conservative than CCA. Recommended for datasets that
  share broadly similar cell type composition.
- `orig.reduction = "pca"`: the PCA computed in step 3 is projected onto each
  reference dataset to find mutual nearest neighbours.
- `normalization.method = "SCT"`: instructs `IntegrateLayers` to use the SCT
  model for re-scaling during anchor scoring (equivalent to the v4
  `PrepSCTIntegration` step, now handled internally).
- Output: a new reduction `integrated.rpca` is added to `merged@reductions`.
  **No new assay is created.**

Alternative methods available in `IntegrateLayers`:

| Method | Use case |
|---|---|
| `RPCAIntegration` | Large datasets, similar cell-type composition (default choice) |
| `CCAIntegration` | Heterogeneous datasets; finds anchors across disparate cell types |
| `HarmonyIntegration` | Very large datasets; requires the `harmony` package |
| `FastMNNIntegration` | Requires `batchelor`; fastest option |
| `scVIIntegration` | Deep-learning approach; requires `reticulate` + `scvi-tools` |

### Step 5 — `JoinLayers()`

```r
merged[["SCT"]] <- JoinLayers(merged[["SCT"]])
```

- Collapses per-sample `counts.<lib>` / `data.<lib>` layers into single
  `counts` and `data` layers.
- Required for `FindMarkers`, `DotPlot`, `FeaturePlot`, and most other
  downstream functions that expect a single expression matrix.
- **Does not modify `integrated.rpca`.**

---

## 4. Downstream Usage

The integrated reduction must be specified explicitly in all graph-based steps:

```r
# Correct v5 downstream workflow
merged <- RunUMAP(merged,
                  reduction = "integrated.rpca",
                  dims      = 1:N_PCS_INTEGRATION)

merged <- FindNeighbors(merged,
                        reduction = "integrated.rpca",
                        dims      = 1:N_PCS_INTEGRATION)

merged <- FindClusters(merged, resolution = 0.5)
```

For differential expression, use the `SCT` assay (the DefaultAssay after
`JoinLayers`):

```r
DefaultAssay(merged) <- "SCT"
markers <- FindAllMarkers(merged, only.pos = TRUE)
```

---

## 5. RAM Considerations (16 GB machine)

| Step | v4 peak RAM | v5 peak RAM |
|---|---|---|
| Normalization | N × SCTransform per object | 1 × SCTransform on merged (sequential layers) |
| Integration | Full anchor weight matrix (all samples) | Pairwise anchor matrices (freed between pairs) |
| Output assay | Full corrected expression matrix | Only low-dim reduction (kb, not GB) |

Knobs to reduce RAM in this script:
- `variable.features.n`: lower → less `scale.data` → less RAM in SCT and PCA.
- `dims` in `IntegrateLayers`: fewer → smaller pairwise matrices.
- `plan("sequential")` in `global_variables.R`: avoids forking overhead and
  RAM duplication from `{future}` workers.

---

## 6. References

| Resource | URL |
|---|---|
| Seurat v5 Integration Introduction | https://satijalab.org/seurat/articles/seurat5_integration |
| SCTransform v2 vignette | https://satijalab.org/seurat/articles/sctransform_vignette |
| `IntegrateLayers` function reference | https://satijalab.org/seurat/reference/integratelayers |
| SCTransform original paper (Hafemeister & Satija 2019) | https://doi.org/10.1186/s13059-019-1874-1 |
| SCTransform v2 paper (Choudhary & Satija 2022) | https://doi.org/10.1016/j.cels.2021.12.010 |
| RPCA integration paper (Stuart et al. 2019) | https://doi.org/10.1016/j.cell.2019.05.031 |

The authoritative API reference for all Seurat v5 functions is at
`https://satijalab.org/seurat/reference/` — each function page includes the
full parameter list, default values, and worked examples.
