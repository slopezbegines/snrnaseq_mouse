# 05_sct.R — SCT normalization per library (step 5 of GSE262881 pipeline)
# Run from project root via: bash code/05_sct.sh
#
# Strategy: split the merged object by sample_id, SCTransform each library in
# isolation (peak ~4-5 GB/library instead of ~14 GB for the full merged object),
# cache each result to disk, then re-merge with merge.data = TRUE.
#
# Per-library caching means a crash on library N lets you resume from library N
# without reprocessing 1..N-1.
#
# Input:  output_path/RData/checkpoint_07_seurat_merged.rds
# Output: output_path/RData/sct_objects/sct_<sample_id>.rds  (per-library)
#         output_path/RData/checkpoint_08_seurat_merged_sct.rds

rm(list = ls())
gc()

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
  paste0("05_sct_normalization_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
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

ram_mb <- function() round(sum(gc(verbose = FALSE)[, 2]), 0)

log_msg("[START] 05_sct.R — log: ", log_file)
log_msg("[INFO]  output_path = ", output_path)

create_directories(output_path)
sct_dir <- file.path(output_path, "RData", "sct_objects")
bpcells_sct_root <- file.path(output_path, "RData", "BPCells_SCT")
dir.create(sct_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bpcells_sct_root, recursive = TRUE, showWarnings = FALSE)

# --- Skip if already done -----------------------------------------------------
if (check_checkpoint("08_seurat_merged_sct")) {
  log_msg("[SKIP] checkpoint_08 found — nothing to do.")
  quit(save = "no", status = 0)
}

# --- Load checkpoint_07 and split by sample -----------------------------------
if (!check_checkpoint("06_seurat_singlets")) {
  log_msg(
    "[ABORT] checkpoint_06_seurat_singlets.rds not found in: ",
    output_path, "RData/"
  )
  stop("Cannot proceed: run step 06 first, or update output_path.")
}

log_msg("[LOAD] Loading checkpoint_06_seurat_singlets...")
seurat_singlets <- load_checkpoint("06_seurat_singlets")
log_msg("[RAM]  After load: ", ram_mb(), " MB")

gc()

sample_names <- names(seurat_singlets)
log_msg("[SPLIT] ", length(sample_names), " samples: ", paste(sample_names, collapse = ", "))
log_msg("[RAM]  After split: ", ram_mb(), " MB")

# Convert SCTAssay (v3) to Assay5 and offload counts + data to BPCells on disk.
# SCTransform always creates SCTAssay (v3), which rejects BPCells MatrixDir layers.
# Coercing to Assay5 first makes LayerData<- accept BPCells matrices.
# Idempotent: if BPCells dirs already exist, they are reused (safe to re-run).
sct_to_assay5_bpcells <- function(obj, nm) {
  obj[["SCT"]] <- as(obj[["SCT"]], "Assay5")
  gc()

  for (layer_nm in c("counts", "data")) {
    bp_dir <- file.path(bpcells_sct_root, nm, layer_nm)
    if (!dir.exists(bp_dir)) {
      mat <- LayerData(obj, assay = "SCT", layer = layer_nm)
      bp_mat <- BPCells::write_matrix_dir(mat, dir = bp_dir)
      rm(mat)
      gc()
    } else {
      bp_mat <- BPCells::open_matrix_dir(bp_dir)
    }
    LayerData(obj, assay = "SCT", layer = layer_nm) <- bp_mat
    rm(bp_mat)
  }
  gc()
  obj
}

# --- SCT per library (with per-library disk cache) ----------------------------
for (nm in sample_names) {
  sct_file <- file.path(sct_dir, paste0("sct_", nm, ".rds"))

  if (file.exists(sct_file)) {
    log_msg("[SCT] [SKIP] ", nm, " — per-library cache found")
    seurat_singlets[[nm]] <- NULL
    gc()
    next
  }

  obj <- seurat_singlets[[nm]]
  seurat_singlets[[nm]] <- NULL # free list slot before SCT to minimise peak RAM
  gc()

  log_msg("[SCT] [START] ", nm, " | ", ncol(obj), " cells")
  log_msg("[RAM]  Before SCT: ", ram_mb(), " MB")

  obj <- SCTransform(
    obj,
    vst.flavor          = "v2",
    variable.features.n = N_INTEGRATION_FEATURES,
    seed.use            = 123,
    verbose             = TRUE
  )

  # SCTransform creates SCTAssay (v3). Convert to Assay5 + BPCells so the
  # saved .rds is lightweight and the merge loop requires no retroactive conversion.
  obj <- sct_to_assay5_bpcells(obj, nm)
  log_msg("[SCT] SCTAssay→Assay5 + BPCells done | RAM: ", ram_mb(), " MB")

  saveRDS(obj, sct_file)
  rm(obj)
  gc()

  log_msg("[SCT] [DONE] ", nm, " — saved: ", sct_file)
  log_msg("[RAM]  After SCT + gc: ", ram_mb(), " MB")
}

rm(seurat_singlets)
gc()
log_msg("[RAM]  After all libraries: ", ram_mb(), " MB")

# --- Merge per-library SCT objects --------------------------------------------
# Load ALL objects first, then call merge() once with y = list(obj2, ..., objN).
# Iterative merge(x, y=obj_i) adds ".1"/".2" suffixes on every iteration, so
# after N-1 merges the first sample accumulates N-1 nested suffixes.
# A single merge call produces clean layer names: counts.1, counts.2, ..., counts.N.
log_msg("[MERGE] Loading all ", length(sample_names), " SCT objects (Assay5 + BPCells)...")

all_sct <- lapply(seq_along(sample_names), function(i) {
  nm <- sample_names[i]
  obj <- readRDS(file.path(sct_dir, paste0("sct_", nm, ".rds")))
  obj <- sct_to_assay5_bpcells(obj, nm)
  log_msg(sprintf("[LOAD] [%d/%d] %s | RAM: %d MB", i, length(sample_names), nm, ram_mb()))
  obj
})
names(all_sct) <- sample_names

log_msg("[SCT] Reconciling variable features (intersection across all libraries)...")
vf_per_sample <- lapply(all_sct, VariableFeatures)
vf_shared     <- Reduce("intersect", vf_per_sample)
rm(vf_per_sample)
log_msg(sprintf("[SCT] Variable features in intersection: %d", length(vf_shared)))

log_msg("[MERGE] Single merge call...")
merged_sct <- merge(
  all_sct[[1]],
  y          = all_sct[-1],
  merge.data = TRUE
)
rm(all_sct)
gc()

VariableFeatures(merged_sct) <- vf_shared
rm(vf_shared)

log_msg(sprintf(
  "[MERGE] Done: %d cells | %d genes | RAM: %d MB",
  ncol(merged_sct), nrow(merged_sct), ram_mb()
))
# --- Save checkpoint ----------------------------------------------------------
save_checkpoint(merged_sct, "08_seurat_merged_sct")
log_msg("[RAM]  Final: ", ram_mb(), " MB")

message("\n[DONE] Step 05 complete.")
message("  Next: bash code/06_integration.sh")
