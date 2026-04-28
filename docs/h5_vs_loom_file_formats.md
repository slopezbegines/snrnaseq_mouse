# H5 vs Loom: File Formats and Origins in the 10X snRNA-seq Workflow

**Dataset:** GSE262881 — Wang et al., BMC Biology (2025) 23:124  
**Paper:** A single-cell transcriptomic atlas of all cell types in the brain of 5xFAD Alzheimer mice in response to dietary inulin supplementation

---

## Full Processing Workflow (per Methods, p. 19)

```
FASTQ  →  CellRanger 7.1.0  →  raw_feature_bc_matrix/
                                        ↓
                              CellBender 0.2.2
                                        ↓
                    *_CellBender_feature_bc_matrix_filtered.h5   ← .h5 files on GEO

FASTQ  →  CellRanger 7.1.0 (BAM)  →  velocyto 0.17.17
                                        ↓
                                    *.loom.gz   ← .loom files on GEO
```

---

## `.h5` Files — Gene Expression Matrix (main analysis)

**Produced by: CellRanger → CellBender** (two sequential steps)

- **CellRanger 7.1.0** aligns reads to the mouse reference genome mm10 (introns included) and generates the initial barcode × gene count matrix.
- **CellBender 0.2.2** performs a second, critical step: it models and removes **ambient RNA contamination** — transcripts from lysed cells that leak into other droplets. This is particularly important in snRNA-seq from frozen brain tissue.
- The resulting `.h5` file is what is deposited on GEO.

**Contents:** sparse gene × nucleus matrix with corrected RNA counts.  
**Pipeline usage:** loaded directly into Seurat with `Read10X_h5()`.

---

## The `filtered` Suffix

CellRanger produces two types of output:

| Suffix | Contents | Typical size |
|--------|----------|--------------|
| `raw_feature_bc_matrix` | **All** detected barcodes (~700k), including empty droplets | Very large |
| `filtered_feature_bc_matrix` | Only barcodes classified as real cells (EmptyDrops / knee-point detection) | ~10–50k nuclei |

CellBender takes the CellRanger `filtered` output (already pre-filtered) and applies its probabilistic ambient RNA model to it. The `filtered` suffix in the GEO filenames indicates the object **excludes empty droplets** and is ready for QC in Seurat.

**QC thresholds applied after loading (per paper):**
- 500–6,000 identified genes per nucleus
- 1,000–50,000 unique molecular identifiers (UMIs)
- < 0.2% mitochondrial genes
- Top 8% highest doublet-score nuclei removed (DoubletFinder v2.0.3)

---

## `.loom.gz` Files — Spliced/Unspliced Counts (RNA velocity only)

**Produced by: velocyto 0.17.17** (run on CellRanger BAM files)

- `velocyto` scans the CellRanger BAM files and counts, per barcode, how many reads map to **exons** (spliced RNA) vs **introns** (unspliced pre-mRNA).
- These two layers are stored in **loom format** (an HDF5 specialization for single-cell data).
- In the paper, loom files are used exclusively for **RNA velocity analysis** with scVelo 0.3.2: imported into Python, combined with Seurat object metadata, and used to infer the direction and rate of cellular differentiation (Figs. 3Q, 4K, 5I–J).

**Contents:** two matrices (spliced + unspliced) × barcodes × genes.  
**Not used** for the main gene expression analysis — RNA velocity only.

---

## Why Forebrain Samples Lack `.loom.gz` Files

Files `GSM8181019–024` (forebrain) have only `.h5` files, while interbrain, brainstem, and cerebellum samples have both. This is consistent with the paper's selective RNA velocity analyses:

| Brain region | RNA velocity applied to |
|---|---|
| Cerebellum | Granule cell subpopulations (Fig. 5I–J) |
| Interbrain | Astrocyte subpopulations (Fig. 4K) |
| Forebrain | Microglia subpopulations (Fig. 3Q) |

The forebrain loom files were likely generated locally by the authors but not deposited on GEO, or were omitted from the submission.

---

## Summary for Pipeline Reproduction

| Goal | Files to use | Tool |
|------|-------------|------|
| Gene expression analysis (clustering, DE, enrichment) | `*_CellBender_feature_bc_matrix_filtered.h5` | Seurat (`Read10X_h5()`) |
| RNA velocity (differentiation trajectories) | `*.loom.gz` | velocyto + scVelo (Python) |
