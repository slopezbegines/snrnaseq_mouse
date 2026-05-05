# 03_filter_cells.R — Cell filtering (step 3 of GSE262881 pipeline)
# Run from project root via: bash code/03_filter_cells.sh
#
# Prerequisite: run bash code/03a_split_qc_checkpoint.sh once first.
# That script splits checkpoint_02_seurat_qc.rds into per-library files under
# RData/qc_libs/ so this script never loads the full list into RAM.

rm(list = ls())
gc()

image_number <- 1
output_path <- "rmds/output/GSE262881/SCT/"
set.seed(123)
options(scipen = 3, digits = 3)

library(future)
plan("sequential")

source("code/00_packages.R")
source("code/global_variables.R")
source("code/01_aux_functions.R")
source("code/02_sc_functions.R")

create_directories(output_path)

message("\n=== STEP 03: Cell Filtering ===\n")

if (check_checkpoint("03_seurat_filtered")) {
  seurat_filtered <- load_checkpoint("03_seurat_filtered")
  message("[SKIP] Checkpoint 03_seurat_filtered already exists.")
} else {
  libs_dir    <- paste0(output_path, "RData/filtered_libs/")
  qc_libs_dir <- paste0(output_path, "RData/qc_libs/")
  dir.create(libs_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Resolve library names from the per-library split produced by 03a --------
  lib_names_file <- paste0(output_path, "RData/lib_names.rds")
  if (!file.exists(lib_names_file)) {
    stop(
      "[ERROR] lib_names.rds not found.\n",
      "        Run bash code/03a_split_qc_checkpoint.sh first to split\n",
      "        checkpoint_02 into individual per-library files."
    )
  }
  lib_names <- readRDS(lib_names_file)
  message(sprintf("[INFO] Libraries to filter: %s", paste(lib_names, collapse = ", ")))

  # --- Sequential loop: load → filter → save → free (one library at a time) ---
  for (lib_name in lib_names) {
    lib_rds <- paste0(libs_dir, "filtered_", lib_name, ".rds")

    if (file.exists(lib_rds)) {
      message(sprintf("[SKIP] Already filtered: %s", lib_name))
      next
    }

    qc_rds <- paste0(qc_libs_dir, "qc_", lib_name, ".rds")
    if (!file.exists(qc_rds)) {
      stop(sprintf(
        "[ERROR] QC file not found: %s\n        Run 03a_split_qc_checkpoint.sh first.",
        qc_rds
      ))
    }

    message(sprintf("[LOAD]   %s", qc_rds))
    obj <- readRDS(qc_rds)
    n0  <- ncol(obj)

    # Layer 1: Demuxlet singlets only (remove multiplets and unassigned)
    # NA DemuxletType means cell was not in the metadata — exclude
    # singlet_pass <- !is.na(obj$DemuxletType) & obj$DemuxletType == "SNG"
    # obj <- subset(obj, cells = colnames(obj)[singlet_pass])

    # Layer 2: Standard QC thresholds
    obj <- subset(obj,
      subset = nFeature_RNA >= QC_MIN_FEATURES &
        nFeature_RNA <= QC_MAX_FEATURES &
        nCount_RNA   >= QC_MIN_COUNTS   &
        nCount_RNA   <= QC_MAX_COUNTS   &
        percent.MT   <= QC_MAX_MT       &
        log10GenesPerUMI >= QC_MIN_COMPLEXITY
    )
    n2 <- ncol(obj)

    message(sprintf(
      "[FILTER] %s: %d cells → QC %d (kept %.1f%%)",
      lib_name, n0, n2, n2 / n0 * 100
    ))

    saveRDS(obj, file = lib_rds)
    message(sprintf("[SAVED]  %s", lib_rds))

    rm(obj)
    gc()
  }

  # --- Assemble final list from individual filtered checkpoints ----------------
  message("\n[ASSEMBLE] Loading individual filtered libraries...")
  seurat_filtered <- lapply(lib_names, function(lib_name) {
    readRDS(paste0(libs_dir, "filtered_", lib_name, ".rds"))
  })
  names(seurat_filtered) <- lib_names

  save_checkpoint(seurat_filtered, "03_seurat_filtered")

  # To remove individual per-library files after the combined checkpoint is
  # confirmed good, uncomment and run manually:
  # unlink(paste0(libs_dir, "filtered_", lib_names, ".rds"))
}

message("\n[DONE] Step 03 complete.")
message(sprintf("[MEMORY] Final RSS: %.1f MB", sum(gc()[, 2])))
