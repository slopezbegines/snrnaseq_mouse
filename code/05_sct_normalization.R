# 05_sct_normalization.R — SCT normalization per library (step 5 of GSE262881 pipeline)
# Run from project root via: bash code/05_sct_normalization.sh
#
# Input:  output_path/RData/checkpoint_04_seurat_singlets.rds
# Output: output_path/RData/sct_objects/sct_<lib_name>.rds  (one per library)
#
# Memory strategy: process ONE library at a time; null each slot after loading it
# to keep peak RAM bounded. No combined checkpoint_05 is written here — the
# integration script (06) assembles from individual sct_*.rds files.

rm(list = ls())
gc()

image_number <- 1
output_path <- "./rmds/output/GSE262881/SCT/" # edit this to point at your data
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

log_msg("[START] 05_sct_normalization.R — log: ", log_file)
log_msg("[INFO]  output_path = ", output_path)

create_directories(output_path)
sct_dir <- paste0(output_path, "RData/sct_objects/")
dir.create(sct_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load checkpoint_04 -------------------------------------------------------
if (!check_checkpoint("06_seurat_singlets")) {
  log_msg(
    "[ABORT] checkpoint_06_seurat_singlets.rds not found in: ",
    output_path, "RData/"
  )
  stop("Cannot proceed: run step 06 first, or update output_path to match existing files.")
}
log_msg("[LOAD] Loading checkpoint_06_seurat_singlets...")
seurat_singlets <- load_checkpoint("06_seurat_singlets")
lib_names <- names(seurat_singlets)
log_msg("[OK] ", length(lib_names), " libraries: ", paste(lib_names, collapse = ", "))
log_msg("[RAM] After loading checkpoint_05: ", ram_mb(), " MB")

# --- SCT per library ----------------------------------------------------------
message("\n=== STEP 05: SCTransform Normalization ===\n")

for (lib_name in lib_names) {
  sct_file <- paste0(sct_dir, "sct_", lib_name, ".rds")

  if (file.exists(sct_file)) {
    log_msg("[SKIP] Already exists: sct_", lib_name, ".rds")
    seurat_singlets[[lib_name]] <- NULL
    gc()
    next
  }

  n_cells <- ncol(seurat_singlets[[lib_name]])
  log_msg("[SCT] Processing: ", lib_name, "  (", n_cells, " cells)")

  obj <- seurat_singlets[[lib_name]]
  seurat_singlets[[lib_name]] <- NULL # free slot before SCT to cap peak RAM
  gc()
  log_msg("[RAM] Freed slot, before SCT: ", ram_mb(), " MB")

  obj <- SCTransform(
    obj,
    vst.flavor      = "v2",
    method          = "glmGamPoi",
    vars.to.regress = "percent.MT",
    verbose         = FALSE
  )

  # PCA run immediately after SCT — required for RPCA integration downstream
  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs     = N_PCS_MAX,
    verbose  = FALSE
  )

  # Atomic write (temp → rename) prevents corrupt files on mid-save crash
  tmp_file <- paste0(sct_file, ".tmp")
  saveRDS(obj, file = tmp_file)
  file.rename(tmp_file, sct_file)

  rm(obj)
  gc()
  log_msg("[SAVED] sct_", lib_name, ".rds  RAM after gc: ", ram_mb(), " MB")
}

rm(seurat_singlets)
gc()
log_msg("[RAM] After all libraries: ", ram_mb(), " MB")

# --- Summary ------------------------------------------------------------------
finished <- list.files(sct_dir, pattern = "^sct_.*\\.rds$", full.names = FALSE)
log_msg("[DONE] ", length(finished), " SCT objects in ", sct_dir)
for (f in finished) log_msg("  ", f)
message("\n[DONE] Step 05 complete.")
message("  Next: bash code/06_integration.sh  (>= 32 GB RAM recommended)")
