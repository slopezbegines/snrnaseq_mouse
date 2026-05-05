# 03a_split_qc_checkpoint.R — One-time conversion: split checkpoint_02 into
# per-library RDS files so that step 03 can load one library at a time.
#
# WHY: checkpoint_02_seurat_qc.rds expands to ~12-15 GB in RAM (2.6 GB on disk),
# exceeding the 13 GB systemd limit used in the main pipeline. Loading the full
# list is unavoidable here, but once split we never need to load it all at once
# again. Run this script ONCE with the higher-limit shell wrapper, then proceed
# normally with bash code/03_filter_cells.sh.

rm(list = ls())
gc()

output_path <- "rmds/output/GSE262881/SCT/"
set.seed(123)
options(scipen = 3, digits = 3)

source("code/00_packages.R")
source("code/global_variables.R")
source("code/01_aux_functions.R")

qc_libs_dir <- paste0(output_path, "RData/qc_libs/")
dir.create(qc_libs_dir, recursive = TRUE, showWarnings = FALSE)

message("\n=== STEP 03a: Splitting QC checkpoint into per-library files ===\n")

checkpoint_path <- paste0(output_path, "RData/", CHECKPOINT_PREFIX, "02_seurat_qc.rds")
if (!file.exists(checkpoint_path)) {
  stop("[ERROR] checkpoint_02_seurat_qc.rds not found at: ", checkpoint_path)
}

message("[LOAD] Reading checkpoint_02_seurat_qc.rds (this will use ~12-15 GB RAM)...")
seurat_qc <- readRDS(checkpoint_path)
lib_names  <- names(seurat_qc)
message(sprintf("[INFO] Found %d libraries: %s", length(lib_names), paste(lib_names, collapse = ", ")))

# Save the library name list so step 03 can iterate without loading anything
saveRDS(lib_names, paste0(output_path, "RData/lib_names.rds"))
message("[SAVED] lib_names.rds")

for (lib_name in lib_names) {
  out_file <- paste0(qc_libs_dir, "qc_", lib_name, ".rds")
  if (file.exists(out_file)) {
    message(sprintf("[SKIP]  Already exists: %s", out_file))
  } else {
    saveRDS(seurat_qc[[lib_name]], out_file)
    message(sprintf("[SAVED] %s  (%.0f MB on disk)", out_file, file.size(out_file) / 1e6))
  }
  seurat_qc[[lib_name]] <- NULL
  gc()
}

message("\n[DONE] All QC libraries saved to: ", qc_libs_dir)
message("       You can now run: bash code/03_filter_cells.sh")
message(sprintf("[MEMORY] Final RSS: %.1f MB", sum(gc()[, 2])))
