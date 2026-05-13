# 08_lognormalize.R — LogNormalize pipeline: merge → normalise → integrate → cluster → UMAP
#
# Alternative to the SCTransform path (05_sct.R + integration).
# Works entirely within Assay5 + BPCells — no SCTAssay conversion needed.
# Officially supported BPCells pathway per Seurat v5 documentation.
#
# Strategy:
#   1. Merge all per-sample singlet objects in a single merge() call.
#   2. NormalizeData + FindVariableFeatures on split layers (per-sample, streaming).
#   3. ScaleData on variable genes only (loads HVGs into RAM — bounded).
#   4. RunPCA.
#   5. IntegrateLayers: RPCA and Harmony (separate checkpoints).
#   6. JoinLayers → FindNeighbors → FindClusters → RunUMAP.
#
# Input:  output_path/RData/checkpoint_06_seurat_singlets.rds
# Output: output_path/RData/checkpoint_lognorm_01_merged.rds
#         output_path/RData/checkpoint_lognorm_02_rpca.rds
#         output_path/RData/checkpoint_lognorm_03_harmony.rds
#
# Run from project root:
#   bash code/08_lognormalize.sh

rm(list = ls())
gc()

output_path <- "./rmds/output/GSE262881/LogNorm/"
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
  paste0("08_lognormalize_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
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

image_number <- 1

log_msg("[START] 08_lognormalize.R — log: ", log_file)
log_msg("[INFO]  output_path = ", output_path)
log_msg("[INFO]  N_INTEGRATION_FEATURES = ", N_INTEGRATION_FEATURES)
log_msg("[INFO]  N_PCS_INTEGRATION = ", N_PCS_INTEGRATION)

create_directories(output_path)

# --- Guard: skip if final checkpoint already exists ---------------------------
if (check_checkpoint("lognorm_03_harmony")) {
  log_msg("[SKIP] checkpoint_lognorm_03_harmony found — nothing to do.")
  quit(save = "no", status = 0)
}

# --- Load checkpoint_06_seurat_singlets ---------------------------------------
sct_output_path <- "./rmds/output/GSE262881/SCT/"

singlets_path <- file.path(
  sct_output_path, "RData",
  paste0(CHECKPOINT_PREFIX, "06_seurat_singlets.rds")
)
if (!file.exists(singlets_path)) {
  stop(
    "[ABORT] checkpoint_06_seurat_singlets.rds not found at: ", singlets_path,
    "\n  Run QC + doublet removal steps first."
  )
}

log_msg("[LOAD] Loading checkpoint_06_seurat_singlets from: ", singlets_path)
seurat_singlets <- readRDS(singlets_path)
sample_names <- names(seurat_singlets)
log_msg("[LOAD] ", length(sample_names), " samples: ", paste(sample_names, collapse = ", "))
log_msg("[RAM]  After load: ", ram_mb(), " MB")


# BLOCK 1: MERGE + NORMALISE + VARIABLE FEATURES + SCALE + PCA


if (check_checkpoint("lognorm_01_merged")) {
  log_msg("[SKIP] lognorm_01_merged found — loading...")
  merged <- load_checkpoint("lognorm_01_merged")
} else {
  # --- Merge (single call) ----------------------------------------------------
  # Single merge() avoids ".1"/".2" suffix accumulation from iterative merging.
  # RNA assay stays Assay5 + BPCells on disk throughout.
  log_msg("[MERGE] Single merge call — ", length(sample_names), " objects...")
  merged <- merge(
    seurat_singlets[[1]],
    y = seurat_singlets[-1]
  )
  rm(seurat_singlets)
  gc()

  log_msg(sprintf(
    "[MERGE] Done: %d cells | %d genes | RAM: %d MB",
    ncol(merged), nrow(merged), ram_mb()
  ))

  # --- Normalise (per-layer, streaming over BPCells) --------------------------
  # Seurat v5 detects split layers and applies log-normalisation independently
  # per sample — equivalent to normalising each sample separately before merging.
  log_msg("[NORM] NormalizeData (LogNormalize, scale.factor = 10000)...")
  merged <- NormalizeData(merged,
    normalization.method = "LogNormalize",
    scale.factor         = 10000,
    verbose              = FALSE
  )
  log_msg("[RAM]  After NormalizeData: ", ram_mb(), " MB")

  # --- Variable features (consensus across samples) ---------------------------
  # FindVariableFeatures on split-layer objects computes HVGs per sample and
  # returns a consensus ranked set. Avoids one sample dominating HVG selection.
  log_msg("[HVG] FindVariableFeatures (nfeatures = ", N_INTEGRATION_FEATURES, ")...")
  merged <- FindVariableFeatures(merged,
    selection.method = "vst",
    nfeatures        = N_INTEGRATION_FEATURES,
    verbose          = FALSE
  )
  log_msg("[HVG] ", length(VariableFeatures(merged)), " variable features selected.")
  log_msg("[RAM]  After FindVariableFeatures: ", ram_mb(), " MB")

  # --- Scale (variable genes only, loads into RAM) ----------------------------
  # Only ~N_INTEGRATION_FEATURES genes are scaled — manageable RAM cost.
  # Required before PCA because NormalizeData does not centre or scale variance.
  log_msg("[SCALE] ScaleData on variable features...")
  merged <- ScaleData(merged,
    features = VariableFeatures(merged),
    verbose  = FALSE
  )
  log_msg("[RAM]  After ScaleData: ", ram_mb(), " MB")

  # --- PCA --------------------------------------------------------------------
  log_msg("[PCA] RunPCA (npcs = ", N_PCS_MAX, ")...")
  merged <- RunPCA(merged,
    assay   = "RNA",
    npcs    = N_PCS_MAX,
    verbose = FALSE
  )
  log_msg("[RAM]  After RunPCA: ", ram_mb(), " MB")

  # Elbow plot to guide choice of N_PCS_INTEGRATION
  p_elbow <- ElbowPlot(merged, ndims = N_PCS_MAX)
  save_plot("lognorm_elbow", p_elbow, width = 8, height = 5, subdir = "integration/")

  save_checkpoint(merged, "lognorm_01_merged")
  log_msg("[CHECKPOINT] lognorm_01_merged saved.")
}


# BLOCK 2: RPCA INTEGRATION


if (check_checkpoint("lognorm_02_rpca")) {
  log_msg("[SKIP] lognorm_02_rpca found — loading...")
  obj_rpca <- load_checkpoint("lognorm_02_rpca")
} else {
  log_msg("[RPCA] IntegrateLayers — method = RPCAIntegration...")
  log_msg("[INFO]  group.by.vars = sample_id | dims = 1:", N_PCS_INTEGRATION)

  obj_rpca <- IntegrateLayers(
    object               = merged,
    method               = RPCAIntegration,
    orig.reduction       = "pca",
    normalization.method = "LogNormalize",
    dims                 = 1:N_PCS_INTEGRATION,
    verbose              = TRUE
  )
  gc()
  log_msg("[RAM]  After RPCA integration: ", ram_mb(), " MB")

  # --- Join layers (required before FindMarkers and most visualisation) -------
  log_msg("[LAYERS] JoinLayers on RNA assay...")
  obj_rpca[["RNA"]] <- JoinLayers(obj_rpca[["RNA"]])
  gc()
  log_msg("[RAM]  After JoinLayers: ", ram_mb(), " MB")

  # --- Clustering + UMAP (RPCA) -----------------------------------------------
  log_msg("[CLUSTER] FindNeighbors + FindClusters (RPCA)...")
  obj_rpca <- FindNeighbors(obj_rpca,
    reduction = "integrated.rpca",
    dims      = 1:N_PCS_INTEGRATION,
    verbose   = FALSE
  )

  for (res in CLUSTERING_RESOLUTIONS) {
    obj_rpca <- FindClusters(obj_rpca,
      resolution = res,
      algorithm  = 1,
      verbose    = FALSE
    )
    log_msg(sprintf(
      "[CLUSTER] resolution %.2f → %d clusters",
      res, length(unique(obj_rpca@meta.data[[paste0("RNA_snn_res.", res)]]))
    ))
  }

  log_msg("[UMAP] RunUMAP on integrated.rpca (dims = 1:", N_PCS_INTEGRATION, ")...")
  obj_rpca <- RunUMAP(obj_rpca,
    reduction      = "integrated.rpca",
    dims           = 1:N_PCS_INTEGRATION,
    umap.method    = "uwot",
    reduction.name = "umap.rpca",
    verbose        = FALSE
  )

  # Plots per resolution
  for (res in CLUSTERING_RESOLUTIONS) {
    col_nm <- paste0("RNA_snn_res.", res)
    Idents(obj_rpca) <- col_nm
    p <- DimPlot(obj_rpca,
      reduction = "umap.rpca",
      label     = TRUE,
      repel     = TRUE
    ) +
      ggtitle(paste0("RPCA | resolution ", res))
    save_plot(
      paste0("lognorm_rpca_umap_res", gsub("\\.", "", as.character(res))),
      p,
      width = 10, height = 8, subdir = "clustering/"
    )
  }

  p_sample <- DimPlot(obj_rpca,
    reduction  = "umap.rpca",
    group.by   = "sample_id",
    label      = FALSE
  ) +
    ggtitle("RPCA | by sample")
  save_plot("lognorm_rpca_umap_sample", p_sample,
    width = 10, height = 8,
    subdir = "clustering/"
  )

  save_checkpoint(obj_rpca, "lognorm_02_rpca")
  log_msg("[CHECKPOINT] lognorm_02_rpca saved.")
  log_msg("[RAM]  Final RPCA block: ", ram_mb(), " MB")
}


# BLOCK 3: HARMONY INTEGRATION


if (check_checkpoint("lognorm_03_harmony")) {
  log_msg("[SKIP] lognorm_03_harmony found — nothing to do.")
} else {
  log_msg("[HARMONY] IntegrateLayers — method = HarmonyIntegration...")
  log_msg("[INFO]    group.by.vars = sample_id | dims = 1:", N_PCS_INTEGRATION)

  # Work from the pre-joined merged object so Harmony gets the split layers
  # (required for IntegrateLayers — JoinLayers must not be done beforehand).
  obj_harmony <- IntegrateLayers(
    object               = merged,
    method               = HarmonyIntegration,
    orig.reduction       = "pca",
    group.by.vars        = "sample_id",
    normalization.method = "LogNormalize",
    dims                 = 1:N_PCS_INTEGRATION,
    verbose              = TRUE
  )
  gc()
  log_msg("[RAM]  After Harmony integration: ", ram_mb(), " MB")

  # --- Join layers ------------------------------------------------------------
  log_msg("[LAYERS] JoinLayers on RNA assay...")
  obj_harmony[["RNA"]] <- JoinLayers(obj_harmony[["RNA"]])
  gc()

  # --- Clustering + UMAP (Harmony) --------------------------------------------
  log_msg("[CLUSTER] FindNeighbors + FindClusters (Harmony)...")
  obj_harmony <- FindNeighbors(obj_harmony,
    reduction = "harmony",
    dims      = 1:N_PCS_INTEGRATION,
    verbose   = FALSE
  )

  for (res in CLUSTERING_RESOLUTIONS) {
    obj_harmony <- FindClusters(obj_harmony,
      resolution = res,
      algorithm  = 1,
      verbose    = FALSE
    )
    log_msg(sprintf(
      "[CLUSTER] resolution %.2f → %d clusters",
      res, length(unique(obj_harmony@meta.data[[paste0("RNA_snn_res.", res)]]))
    ))
  }

  log_msg("[UMAP] RunUMAP on harmony (dims = 1:", N_PCS_INTEGRATION, ")...")
  obj_harmony <- RunUMAP(obj_harmony,
    reduction      = "harmony",
    dims           = 1:N_PCS_INTEGRATION,
    umap.method    = "uwot",
    reduction.name = "umap.harmony",
    verbose        = FALSE
  )

  # Plots per resolution
  for (res in CLUSTERING_RESOLUTIONS) {
    col_nm <- paste0("RNA_snn_res.", res)
    Idents(obj_harmony) <- col_nm
    p <- DimPlot(obj_harmony,
      reduction = "umap.harmony",
      label     = TRUE,
      repel     = TRUE
    ) +
      ggtitle(paste0("Harmony | resolution ", res))
    save_plot(
      paste0("lognorm_harmony_umap_res", gsub("\\.", "", as.character(res))),
      p,
      width = 10, height = 8, subdir = "clustering/"
    )
  }

  p_sample <- DimPlot(obj_harmony,
    reduction  = "umap.harmony",
    group.by   = "sample_id",
    label      = FALSE
  ) +
    ggtitle("Harmony | by sample")
  save_plot("lognorm_harmony_umap_sample", p_sample,
    width = 10, height = 8,
    subdir = "clustering/"
  )

  save_checkpoint(obj_harmony, "lognorm_03_harmony")
  log_msg("[CHECKPOINT] lognorm_03_harmony saved.")
  log_msg("[RAM]  Final Harmony block: ", ram_mb(), " MB")
}

message("\n[DONE] Step 08 complete.")
message("  RPCA result  : checkpoint_lognorm_02_rpca.rds")
message("  Harmony result: checkpoint_lognorm_03_harmony.rds")
