# 06_integration.R — Merge + SCTransform + RPCA integration (Seurat v5 native)
#
# WORKFLOW CHANGE FROM v4:
#   v4: SCTransform per object → FindIntegrationAnchors → IntegrateData → new "integrated" assay
#   v5: merge all → SCTransform on merged (per-layer) → PCA → IntegrateLayers → new reduction
#
# CRITICAL OUTPUT DISTINCTION:
#   IntegrateLayers outputs a new REDUCTION ("integrated.rpca"), NOT a new assay.
#   Use integrated.rpca for FindNeighbors / UMAP / FindClusters.
#   The SCT assay is retained for gene expression and DE analysis.
#
# Input:  checkpoint_06_seurat_singlets.rds  (named list of QC-filtered Seurat objects)
# Output: checkpoint_07_data_integrated.rds  (merged object with integrated.rpca reduction)
#
# See: docs/seurat_v5_sct_integration.md for workflow rationale and Seurat references.
# Run: bash code/06_integration.sh

rm(list = ls())
gc()

image_number <- 1
output_path <- "./rmds/output/GSE262881/SCT/"
set.seed(123)
options(scipen = 3, digits = 3)

library(future)
plan("sequential")

source("code/00_packages.R")
source("code/global_variables.R")
source("code/01_aux_functions.R")
source("code/02_sc_functions.R")

# --- Logging ------------------------------------------------------------------
log_dir <- "logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(
  log_dir,
  paste0("06_integration_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)
log_con <- file(log_file, open = "wt")

log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(...))
  message(line)
  writeLines(line, log_con)
}

on.exit(
  {
    log_msg("[EXIT] Script finished.")
    close(log_con)
  },
  add = TRUE
)

log_msg("[START] 06_integration.R (Seurat v5 native) — log: ", log_file)
log_msg(
  "[INFO]  Seurat ", packageVersion("Seurat"),
  " | SeuratObject ", packageVersion("SeuratObject")
)
log_msg("[INFO]  output_path = ", output_path)

# --- Load singlets checkpoint -------------------------------------------------
if (!check_checkpoint("06_seurat_singlets")) {
  stop(
    "[ABORT] checkpoint_06_seurat_singlets.rds not found in: ",
    output_path, "RData/  — run steps 01–05 first."
  )
}
log_msg("[LOAD] Loading checkpoint_06_seurat_singlets...")
seurat_singlets <- load_checkpoint("06_seurat_singlets")

lib_names <- names(seurat_singlets)
total_cells <- sum(vapply(seurat_singlets, ncol, integer(1)))
log_msg("[OK] ", length(lib_names), " libraries | ", total_cells, " total cells")
log_msg("[LIBS] ", paste(lib_names, collapse = ", "))
log_msg("[RAM] After load: ", ram_mb(), " MB")

# --- STEP 1: Merge all objects into one Seurat object -------------------------
# merge() in Seurat v5 creates an Assay5 with per-sample layers (counts.<lib>, data.<lib>),
# replacing the need to keep objects separate for integration.
if (check_checkpoint("07_merged")) {
  log_msg("[SKIP] 07_merged checkpoint found.")
  merged <- load_checkpoint("07_merged")
} else {
  log_msg("[MERGE] Merging ", length(lib_names), " objects...")
  tic()

  merged <- merge(
    x            = seurat_singlets[[1]],
    y            = seurat_singlets[-1],
    add.cell.ids = lib_names
  )

  rm(seurat_singlets)
  gc()
  elapsed <- toc(quiet = TRUE)

  log_msg("[MERGE] ", ncol(merged), " cells x ", nrow(merged), " features")
  log_msg("[MERGE] RNA class: ", class(merged[["RNA"]])[1])
  log_msg("[MERGE] RNA layers: ", paste(Layers(merged[["RNA"]]), collapse = ", "))
  log_msg("[TIME] Merge: ", round(elapsed$toc - elapsed$tic, 1), " s")
  log_msg("[RAM] After merge + gc: ", ram_mb(), " MB")

  save_checkpoint(merged, "07_merged")
}

# --- STEP 2: SCTransform (per-layer / per-sample) ----------------------------
# Seurat v5 SCTransform detects split layers and fits a separate NB regression
# model per sample, producing a single SCTAssay with per-sample scale.data.
# variable.features.n = N_INTEGRATION_FEATURES reduces peak RAM vs default 3000.
if (check_checkpoint("07_sct")) {
  log_msg("[SKIP] 07_sct checkpoint found.")
  merged <- load_checkpoint("07_sct")
} else {
  log_msg(
    "[SCT] SCTransform | vst.flavor=v2 | method=glmGamPoi | features=",
    N_INTEGRATION_FEATURES
  )
  tic()

  merged <- SCTransform(
    merged,
    vst.flavor          = "v2",
    method              = "glmGamPoi",
    vars.to.regress     = "percent.MT",
    variable.features.n = N_INTEGRATION_FEATURES,
    verbose             = TRUE
  )

  gc()
  elapsed <- toc(quiet = TRUE)
  log_msg("[SCT] assay class: ", class(merged[["SCT"]])[1])
  log_msg("[SCT] variable features: ", length(VariableFeatures(merged)))
  log_msg("[TIME] SCTransform: ", round(elapsed$toc - elapsed$tic, 1), " s")
  log_msg("[RAM] After SCT + gc: ", ram_mb(), " MB")

  save_checkpoint(merged, "07_sct")
}

# --- STEP 3: PCA on SCT-normalised merged data --------------------------------
if (check_checkpoint("07_pca")) {
  log_msg("[SKIP] 07_pca checkpoint found.")
  merged <- load_checkpoint("07_pca")
} else {
  log_msg("[PCA] RunPCA | npcs = ", N_PCS_MAX, " | DefaultAssay: ", DefaultAssay(merged))
  tic()

  merged <- RunPCA(merged, npcs = N_PCS_MAX, verbose = TRUE)

  gc()
  elapsed <- toc(quiet = TRUE)
  log_msg("[PCA] Reductions: ", paste(names(merged@reductions), collapse = ", "))
  log_msg("[TIME] PCA: ", round(elapsed$toc - elapsed$tic, 1), " s")
  log_msg("[RAM] After PCA + gc: ", ram_mb(), " MB")

  save_checkpoint(merged, "07_pca")
}

# --- STEP 4: IntegrateLayers — RPCA + SCT ------------------------------------
# Seurat v5 native integration API. Instead of creating a new assay (v4 IntegrateData),
# IntegrateLayers adds a new reduction: "integrated.rpca".
# All downstream steps (FindNeighbors, UMAP, FindClusters) must use this reduction.
if (check_checkpoint("07_integrated")) {
  log_msg("[SKIP] 07_integrated checkpoint found.")
  merged <- load_checkpoint("07_integrated")
} else {
  log_msg("[INTEGRATE] IntegrateLayers | method=RPCA | SCT | dims=1:", N_PCS_INTEGRATION)
  tic()

  merged <- IntegrateLayers(
    object               = merged,
    method               = RPCAIntegration,
    orig.reduction       = "pca",
    normalization.method = "SCT",
    dims                 = 1:N_PCS_INTEGRATION,
    verbose              = TRUE
  )

  gc()
  elapsed <- toc(quiet = TRUE)
  log_msg("[INTEGRATE] Reductions: ", paste(names(merged@reductions), collapse = ", "))
  log_msg("[TIME] IntegrateLayers: ", round(elapsed$toc - elapsed$tic, 1), " s")
  log_msg("[RAM] After integration + gc: ", ram_mb(), " MB")

  save_checkpoint(merged, "07_integrated")
}

# --- STEP 5: Join SCT layers --------------------------------------------------
# Collapses per-sample SCT layers into a single layer for convenient access
# in downstream functions (DimPlot, FeaturePlot, FindMarkers, etc.).
# Does NOT modify integrated.rpca.
log_msg("[JOIN] Joining SCT layers...")
merged[["SCT"]] <- JoinLayers(merged[["SCT"]])
DefaultAssay(merged) <- "SCT"

log_msg("[JOIN] SCT layers after join: ", paste(Layers(merged[["SCT"]]), collapse = ", "))
log_msg("[JOIN] DefaultAssay: ", DefaultAssay(merged))
log_msg("[JOIN] All assays: ", paste(Assays(merged), collapse = ", "))

# --- Save final integrated object ---------------------------------------------
save_checkpoint(merged, "07_data_integrated")
log_msg("[CHECKPOINT] Saved: checkpoint_07_data_integrated.rds")

summary_line <- paste0(
  "Integrated: ", ncol(merged), " cells | ",
  nrow(merged), " features | ",
  "Reductions: ", paste(names(merged@reductions), collapse = ", ")
)
cat(summary_line, "\n")
log_msg(summary_line)
log_msg("[NEXT] Use reduction = 'integrated.rpca' for RunUMAP, FindNeighbors, FindClusters.")
message("\n[DONE] Step 06 complete.")
