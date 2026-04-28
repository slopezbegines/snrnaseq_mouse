# 1. CHECKPOINT SYSTEM ####
# Enables crash recovery: save progress -> reload on restart
# Pattern: check_checkpoint() -> if not found, compute -> save_checkpoint()


#' Save a named R object as a checkpoint .rds file (atomic write)
#' Uses temp-file + rename to prevent corrupt files on mid-save crash
save_checkpoint <- function(obj, name, base = output_path) {
  dir <- paste0(base, "RData/")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dir, CHECKPOINT_PREFIX, name, "_TMP.rds")
  final <- paste0(dir, CHECKPOINT_PREFIX, name, ".rds")
  saveRDS(obj, file = tmp)
  file.rename(tmp, final)
  message(sprintf(
    "[CHECKPOINT] Saved: %s  (%s)",
    final, format(Sys.time(), "%H:%M:%S")
  ))
}

#' Check whether a checkpoint file exists
check_checkpoint <- function(name, base = output_path) {
  file.exists(paste0(base, "RData/", CHECKPOINT_PREFIX, name, ".rds"))
}

#' Load a checkpoint file and return the object
load_checkpoint <- function(name, base = output_path) {
  path <- paste0(base, "RData/", CHECKPOINT_PREFIX, name, ".rds")
  if (!file.exists(path)) stop("[CHECKPOINT] Not found: ", path)
  message(sprintf(
    "[CHECKPOINT] Loaded: %s  (%s)",
    path, format(Sys.time(), "%H:%M:%S")
  ))
  readRDS(path)
}

# 2. Plot Helpers ####

#' Save a ggplot to TIFF and PDF; increment global image_number counter
#' @param plotname  Base name of the file (no extension)
#' @param plot      A ggplot object
#' @param width     Width in inches
#' @param height    Height in inches
#' @param subdir    Subdirectory under figures/ (e.g. "QC/", "DE/")
save_plot <- function(plotname, plot, width = 8, height = 6, subdir = "") {
  if (!inherits(plot, c("ggplot", "patchwork", "Heatmap"))) {
    warning(sprintf("[SAVE_PLOT] '%s' is not a ggplot/patchwork object — skipping.", plotname))
    return(invisible(NULL))
  }
  dir <- paste0(output_path, "figures/", subdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0(dir, sprintf("%03d", image_number), "_", plotname)
  tryCatch(
    {
      ggsave(paste0(filename, tiff_extension), plot, width = width, height = height, units = "in", dpi = 300)
      ggsave(paste0(filename, pdf_extension), plot, width = width, height = height, units = "in")
      image_number <<- image_number + 1
    },
    error = function(e) {
      warning(sprintf("[SAVE_PLOT] Failed to save '%s': %s", filename, e$message))
    }
  )
}
