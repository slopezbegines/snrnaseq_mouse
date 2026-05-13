# Seurat v5 + BPCells — Full Pipeline Workflow Reference

## Purpose

This document describes the recommended end-to-end snRNA-seq pipeline using
Seurat v5 and BPCells, from raw data import to UMAP visualisation. For each
step it explains which function to use, why it is correct for a v5/BPCells
pipeline, and why the Assay v3 counterpart fails or is insufficient.

For integration-specific rationale (v4 vs v5 API differences, RAM table,
`IntegrateLayers` parameter guide) see `seurat_v5_sct_integration.md`.

---

## The Central Architectural Distinction

Before the workflow: understanding why Assay5 matters.

| | Assay v3 (`Assay` / `SCTAssay`) | Assay v5 (`Assay5`) |
|---|---|---|
| Storage | Fixed S4 slots: `@counts`, `@data`, `@scale.data` | Named layers: any number, any name |
| Slot type validation | S4 `checkSlotAssignment` enforces `dgCMatrix` | `LayerData<-` accepts `IterableMatrix` (BPCells) |
| BPCells compatible | No — `slot<-` rejects `MatrixDir` at the type-check level | Yes — designed for on-disk `IterableMatrix` |
| Multi-sample | One concatenated matrix | One layer per sample: `counts.sample1`, `counts.sample2` |
| Memory model | Always fully in RAM | On-disk; R holds a lazy pointer, not the data |

**Why BPCells only works with Assay5:** BPCells matrices are instances of
`IterableMatrix`. Assay v3's S4 slot validator (`checkSlotAssignment`) requires
`dgCMatrix` or `dgTMatrix` — it rejects `IterableMatrix` with a type error.
Assay5's `LayerData<-` setter explicitly accepts `IterableMatrix`, so no
validation is triggered and the matrix stays on disk.

Any Seurat operation that internally calls `slot<-` on an Assay v3 object
(such as `RenameCells()`, which is called by `merge()`) will fail with a
`checkSlotAssignment` error if BPCells data is present — even if the initial
write succeeded via an `attr<-` bypass.

---

## Workflow

### Step 1 — Load Raw Data as BPCells (on-disk)

**Functions:** `BPCells::open_matrix_10x_hdf5()` → `BPCells::write_matrix_dir()`
→ `BPCells::open_matrix_dir()` → `CreateSeuratObject()`

`open_matrix_10x_hdf5()` reads the count matrix from a CellBender/CellRanger
HDF5 file and returns a lazy `IterableMatrix` — no data is loaded into RAM.
`write_matrix_dir()` serialises that matrix to a compressed bit-packed directory
on disk (BPCells native format). `open_matrix_dir()` reopens it as an
`IterableMatrix` backed by those files.

`CreateSeuratObject()` accepts an `IterableMatrix` as the `counts` argument and
automatically creates an **Assay5** RNA assay with a `counts` layer pointing to
the on-disk files. The cell metadata (nCount_RNA, nFeature_RNA) is computed
lazily from the on-disk matrix without loading it.

**Why not the v3 path:** The v3 path (`Read10X()` → `CreateSeuratObject()`)
loads the full count matrix as a `dgCMatrix` into RAM immediately. For a
16 GB machine with large datasets this is the bottleneck. `Read10X()` also
returns an in-memory sparse matrix that cannot be passed to a BPCells workflow.

---

### Step 2 — Quality Control

**Functions:** `PercentageFeatureSet()`, `VlnPlot()` / `ggplot2`, `subset()`

`PercentageFeatureSet()` works on Assay5 objects with on-disk BPCells layers.
It computes the percentage of counts matching a gene pattern (e.g. `"^mt-"`)
by streaming the BPCells matrix one chunk at a time, so only a small chunk of
the count matrix is in RAM at any moment.

`subset()` filters cells by the computed metadata columns (nFeature_RNA,
percent.MT). It reindexes the BPCells layer to the retained cell barcodes —
the on-disk files are not rewritten; only the in-memory row/column index
changes.

**Doublet detection** (`scDblFinder::scDblFinder()`): BPCells matrices must be
materialised to `dgCMatrix` before passing to `scDblFinder`, because
`SingleCellExperiment` does not accept `IterableMatrix`. This is done
per-sample (one sample in RAM at a time), which is safe at 16 GB.

**Why not the v3 path:** `subset()` on a v3 `Assay` with BPCells layers fails
for the same reason `merge()` does — internally it calls `slot<-` to assign
the subsetted matrix, which triggers `checkSlotAssignment` and rejects
`IterableMatrix`.

---

### Step 3 — SCTransform Normalisation (per sample)

**Function:** `SCTransform()` with `vst.flavor = "v2"`, `method = "glmGamPoi"`

**Why per sample, not on the merged object:** `SCTransform()` cannot run on a
BPCells-backed Assay5 object. Internally, the function calls `dim(x) <-
length(x)` to reshape matrix data — this operation is invalid for a lazy
`IterableMatrix`, which is read-only and cannot be reshaped. The only
workaround is to run `SCTransform()` on one sample at a time, where the
in-memory cost is bounded to a single library.

**What `SCTransform()` produces:** An `SCTAssay` object, not `Assay5`. `SCTAssay`
is a subclass of Assay v3 that carries `@SCTModel.list` — a per-cell record of
the fitted negative-binomial regression parameters (sequencing depth offset,
dispersion estimate, gene-specific coefficients). This slot has no equivalent
in Assay5.  The three data representations stored are:

- `scale.data` — Pearson residuals (non-sparse; stored for variable genes only). This is the PCA input.
- `counts` — depth-corrected UMI counts (what we would observe at uniform sequencing depth).
- `data` — log-normalised corrected counts (for visualisation).

**Why `SCTAssay` cannot become `Assay5` without consequence:** Converting
`SCTAssay` to `Assay5` discards `@SCTModel.list`. This breaks
`PrepSCTFindMarkers()`, which recorrects counts before DE testing, and it
breaks the internal anchor-rescoring step inside
`IntegrateLayers(normalization.method = "SCT")`. In practice, the loss of
`@SCTModel.list` means DE results after integration are less accurate. This is
a known open issue in the Seurat GitHub (issues #7542, #9707) and as of 2025
there is no official resolution.

**`vst.flavor = "v2"`** uses regularised regression across genes (avoids
overfitting to low-count genes). **`method = "glmGamPoi"`** fits the
negative-binomial GLM with the Gamma-Poisson approximation, ~10× faster than
the default.

---

### Step 4 — Offload SCT Object to BPCells (per sample, after SCTransform)

**Function:** `sct_to_assay5_bpcells()` (project-internal, `code/02_sc_functions.R`)

After `SCTransform()` produces an `SCTAssay` (in-memory), the `counts` and
`data` layers are large sparse matrices. This function:

1. Coerces `SCTAssay` → `Assay5` with `as(sct_assay, "Assay5")`.
2. Writes `counts` and `data` to BPCells on-disk directories via
   `BPCells::write_matrix_dir()`.
3. Reassigns those layers via `LayerData<-` (valid because Assay5 accepts
   `IterableMatrix`).

The resulting object is saved as an RDS. The RDS file is lightweight (contains
only pointers to the BPCells directories and the in-memory `scale.data`
residuals for variable genes).

**Trade-off:** Coercing to Assay5 discards `@SCTModel.list`. See Step 3 for
consequences. This is accepted here in favour of RAM efficiency during the
merge step. `PrepSCTFindMarkers()` cannot be used after this conversion;
`FindMarkers()` on the SCT `counts` layer directly is the alternative.

**Why not `sct_to_bpcells()`:** That function keeps the assay as `SCTAssay`
and uses `attr<-` to bypass S4 validation. It succeeds at write time, but any
subsequent Seurat operation that reassigns slots through the proper S4 path
(including `merge()` → `RenameCells()` → `slot<-`) will fail with a
`checkSlotAssignment` error.

---

### Step 5 — Merge

**Function:** `merge()` — single call with `y = list(obj2, ..., objN)`

A single `merge()` call is mandatory (not iterative). Iterative merging
(`merge(a, b)` then `merge(result, c)`) appends `.1`, `.2`, … suffixes to
layer names at every iteration, so the first sample accumulates N-1 nested
suffixes after N-1 merges. A single call produces clean layer names:
`counts.sample1`, `counts.sample2`, etc.

`merge()` internally calls `RenameCells()` to disambiguate barcodes across
samples. This calls `slot<-` on each assay — which is why the SCT assay must
be **Assay5** (not `SCTAssay`) by this point. With Assay5, `slot<-` targets
the layer index (a `list`), not a typed S4 slot, so `IterableMatrix` values
are accepted without type validation.

After merge, the object has one layer per sample under the SCT assay:
`counts.sample1`, `data.sample1`, `counts.sample2`, `data.sample2`, etc.
`scale.data` is not split per sample — it is the concatenated residual matrix
for all variable genes.

---

### Step 6 — PCA

**Function:** `RunPCA()` with `assay = "SCT"`

PCA is computed on the `scale.data` slot of the SCT assay — the Pearson
residual matrix. This matrix is already variance-stabilised and mean-centred by
SCTransform, so no additional `ScaleData()` call is needed (unlike the
LogNormalize path, where `ScaleData()` is a required prior step).

`RunPCA()` works on Assay5 because it reads `scale.data` via `GetAssayData()`,
which works equally for Assay5 and Assay v3. The PCA reduction is stored in
`object@reductions[["pca"]]` — not in the assay — so assay version is
irrelevant for the output.

**Why not `ScaleData()` + `RunPCA()` as in the v3/LogNormalize path:**
`ScaleData()` would overwrite the Pearson residuals in `scale.data` with
z-scores of log-normalised data. That is correct for the LogNormalize pathway,
but incorrect after SCTransform, where `scale.data` already contains the
proper input for PCA.

---

### Step 7 — Integration

#### Option A: RPCA — `IntegrateLayers()` with `method = RPCAIntegration`

**Function:** `IntegrateLayers()` with `normalization.method = "SCT"`,
`orig.reduction = "pca"`, `method = RPCAIntegration`

RPCA finds mutual nearest neighbours by projecting each dataset onto the PCA
space of a reference dataset (reciprocal projection). It is faster than CCA
and more conservative — it corrects only genuine batch effects, not
biological variance. Best suited when datasets share broad cell-type
composition.

`normalization.method = "SCT"` tells `IntegrateLayers` to use the SCT model
parameters during anchor rescoring (equivalent to the v4 `PrepSCTIntegration`
step, now internalised). This requires that the SCT assay still carries
enough model information — which is the reason the model-parameter loss from
Step 4 is a documented compromise.

Output: a new reduction `integrated.rpca` is added to `@reductions`. No new
assay is created. All downstream graph-building and embedding uses this
reduction.

**Why not the v4 `FindIntegrationAnchors()` + `IntegrateData()` path:**
`IntegrateData()` constructs the full anchor weight matrix for all sample pairs
simultaneously and stores a corrected expression matrix as a new assay. Both
operations exhaust RAM on large datasets. `IntegrateLayers()` processes one
pair at a time and stores only a low-dimensional reduction (kilobytes, not
gigabytes).

#### Option B: Harmony — `IntegrateLayers()` with `method = HarmonyIntegration`

**Function:** `IntegrateLayers()` with `method = HarmonyIntegration`,
`normalization.method = "SCT"`, `group.by.vars = "orig.ident"`

Harmony iteratively adjusts cell embeddings in PCA space to remove
batch-specific variance while preserving biological variance. It is faster
than RPCA for very large datasets (> 100k cells) and handles more complex
batch structures (e.g. multiple confounders). It requires the `harmony` R
package.

Output: a new reduction `harmony` is added to `@reductions`. Use
`reduction = "harmony"` in all subsequent steps.

**When to choose RPCA vs Harmony:**

| | RPCA | Harmony |
|---|---|---|
| Speed | Moderate | Fast |
| Datasets with very different compositions | May over-correct | More robust |
| Requires anchor finding | Yes | No |
| Multiple batch variables | No (one batch variable) | Yes |
| Best for | WT/KO, small N samples | Many donors, complex designs |

---

### Step 8 — Join Layers

**Function:** `JoinLayers()` on the SCT assay

After integration, the SCT assay still has per-sample layers
(`counts.sample1`, `data.sample1`, …). `JoinLayers()` collapses them into
single `counts` and `data` layers. This is required before `FindMarkers()`,
`DotPlot()`, `FeaturePlot()`, and most downstream functions that expect a
single expression matrix.

`JoinLayers()` is an Assay5-only function. There is no v3 equivalent because
v3 assays never had split layers — they always stored one concatenated matrix.

**Important:** `JoinLayers()` on a BPCells-backed Assay5 will materialise the
on-disk layers into RAM (because concatenating on-disk matrices requires
reading all data). Plan for the resulting RAM cost.

---

### Step 9 — Graph Construction and Clustering

**Functions:** `FindNeighbors()`, `FindClusters()`

`FindNeighbors()` builds a k-nearest-neighbour graph in the integrated
reduction space (`reduction = "integrated.rpca"` or `reduction = "harmony"`).
It operates on the low-dimensional embedding, not on the assay — so assay
version is irrelevant here.

`FindClusters()` applies the Louvain (default) or Leiden algorithm to the KNN
graph. The `resolution` parameter controls cluster granularity (higher =
more clusters). Explore a range (0.2–1.5) and compare to known biology.

---

### Step 10 — UMAP

**Function:** `RunUMAP()` with `reduction = "integrated.rpca"` (or `"harmony"`)

UMAP embeds the integrated low-dimensional space into 2D for visualisation.
Using the integrated reduction (not `"pca"`) is mandatory — using `"pca"` would
visualise unintegrated variance, making batch effects visible instead of
biological structure.

`return.model = TRUE` is recommended if you intend to project new data onto the
same embedding later (e.g. reference mapping or sketch-to-full projection).

---

### Step 11 — Visualisation

**Functions:** `DimPlot()`, `FeaturePlot()`, `VlnPlot()`, `DotPlot()`

All these functions work on Assay5 objects. They access expression values via
`GetAssayData()`, which is version-agnostic. `DimPlot()` reads from
`@reductions`, not from any assay, so it is fully independent of assay version.

`FeaturePlot()` and `VlnPlot()` read from the `data` layer of the default
assay (SCT after `JoinLayers()`). If layers are still split (Step 8 not yet
done), these functions will only see `data.sample1` and fail to find a unified
`data` layer.

---

## Summary Table

| Step | Key Function(s) | Assay produced | BPCells on disk? |
|---|---|---|---|
| 1 — Import | `open_matrix_10x_hdf5`, `write_matrix_dir`, `CreateSeuratObject` | RNA (Assay5) | Yes |
| 2 — QC | `PercentageFeatureSet`, `subset` | RNA (Assay5) | Yes |
| 3 — SCTransform | `SCTransform` | SCT (SCTAssay v3) | No (in-memory) |
| 4 — Offload SCT | `sct_to_assay5_bpcells` | SCT (Assay5) | Yes |
| 5 — Merge | `merge` | SCT (Assay5, split layers) | Yes |
| 6 — PCA | `RunPCA` | reduction: `pca` | — |
| 7 — Integration | `IntegrateLayers` | reduction: `integrated.rpca` / `harmony` | — |
| 8 — Join layers | `JoinLayers` | SCT (Assay5, unified) | No (materialised) |
| 9 — Clustering | `FindNeighbors`, `FindClusters` | graph + cluster labels | — |
| 10 — UMAP | `RunUMAP` | reduction: `umap` | — |
| 11 — Visualise | `DimPlot`, `FeaturePlot`, `DotPlot` | reads SCT data layer | — |

---

## Known Limitations of This Workflow

**`@SCTModel.list` loss at Step 4.** Converting `SCTAssay` → `Assay5` discards
the per-cell regression model. Downstream consequences:

- `PrepSCTFindMarkers()` cannot be used. Use `FindMarkers()` directly on the
  SCT `counts` layer instead.
- Anchor rescoring in `IntegrateLayers(normalization.method = "SCT")` uses an
  approximation (no per-cell model available). Results are still good but not
  identical to the full-model path.

This is an active limitation of the BPCells + SCT path and is tracked in Seurat
GitHub issues [#7814](https://github.com/satijalab/seurat/issues/7814) and
[#7542](https://github.com/satijalab/seurat/issues/7542). If the Seurat team
officially supports a BPCells-compatible `SCTransform` in a future release,
Step 4 can be replaced by the official function.

**SCTransform cannot run on BPCells-backed objects.** Running `SCTransform()`
directly on an Assay5 object backed by BPCells fails because the function
internally calls `dim(x) <- length(x)` to reshape the matrix — an operation
that is invalid for a lazy, read-only `IterableMatrix`. Per-sample execution
(Step 3) is the only viable workaround.

---

## References

| Resource | URL |
|---|---|
| Seurat v5 BPCells vignette | https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette |
| Seurat v5 integration vignette | https://satijalab.org/seurat/articles/seurat5_integration |
| SCTransform vignette | https://satijalab.org/seurat/articles/sctransform_vignette |
| BPCells package | https://bnprks.github.io/BPCells |
| SCT + BPCells GitHub issue | https://github.com/satijalab/seurat/issues/7814 |
| SCT + IntegrateLayers GitHub issue | https://github.com/satijalab/seurat/issues/7542 |
| Companion doc: integration rationale | `seurat_v5_sct_integration.md` |
