# Load CellBender-filtered h5 files as a named list of Seurat objects.
# File naming convention: GSM{ID}_{Ctrl|Inulin}_{region}_{replicate}_CellBender_feature_bc_matrix_filtered.h5
# Requires: output_path, CHECKPOINT_PREFIX (set in global_variables.R)


# 2. DATA LOADING ####


#' Load CellBender h5 files as a named list of Seurat objects
#'
#' @param data_path    Directory containing the h5 files
#' @param regions      Character vector of regions to load; NULL loads all four
#' @param min.cells    Passed to CreateSeuratObject()
#' @param min.features Passed to CreateSeuratObject()
#' @param checkpoint   Checkpoint name (e.g. "01_raw_data"); NULL skips checkpointing
#'
#' @return Named list of Seurat objects (names: Ctrl_forebrain_1, ...)
#'
#' @examples
#' # Load all regions (with crash recovery)
#' seurat_list <- load_h5_samples(data_path, checkpoint = "01_raw_data")
#'
#' # Load only forebrain (no checkpoint)
#' seurat_list <- load_h5_samples(data_path, regions = "forebrain")
load_h5_samples <- function(
  data_path = DATA_DIR,
  regions = c("forebrain", "interbrain", "brainstem", "cerebellum"),
  min.cells = 3,
  min.features = 200,
  checkpoint = NULL
) {
  if (!is.null(checkpoint) && check_checkpoint(checkpoint)) {
    return(load_checkpoint(checkpoint))
  }

  h5_files <- list.files(
    path = data_path,
    pattern = "_CellBender_feature_bc_matrix_filtered\\.h5$",
    full.names = TRUE
  )

  if (length(h5_files) == 0) {
    stop("No CellBender h5 files found in: ", data_path)
  }

  fnames <- basename(h5_files)
  pattern <- paste0(
    "^(GSM\\d+)",
    "_(Ctrl|Inulin)",
    "_(forebrain|interbrain|brainstem|cerebellum)",
    "_(\\d+)",
    "_CellBender.*\\.h5$"
  )
  matches <- regmatches(fnames, regexec(pattern, fnames))

  sample_table <- do.call(rbind, lapply(seq_along(matches), function(i) {
    m <- matches[[i]]
    if (length(m) == 0) {
      warning("Skipping unparseable filename: ", fnames[i])
      return(NULL)
    }
    data.frame(
      filepath = h5_files[i],
      filename = fnames[i],
      gsm_id = m[2],
      condition = m[3],
      region = m[4],
      replicate = as.integer(m[5]),
      sample_id = paste(m[3], m[4], m[5], sep = "_"),
      stringsAsFactors = FALSE
    )
  }))

  sample_table <- sample_table[sample_table$region %in% regions, ]

  if (nrow(sample_table) == 0) {
    stop("No samples found for region(s): ", paste(regions, collapse = ", "))
  }

  message(sprintf(
    "Loading %d sample(s) | regions: %s",
    nrow(sample_table), paste(regions, collapse = ", ")
  ))

  seurat_list <- lapply(seq_len(nrow(sample_table)), function(i) {
    s <- sample_table[i, ]
    message(sprintf("  [%d/%d] %s (%s)", i, nrow(sample_table), s$sample_id, s$gsm_id))

    counts <- Read10X_h5(s$filepath)

    obj <- CreateSeuratObject(
      counts       = counts,
      project      = s$sample_id,
      min.cells    = min.cells,
      min.features = min.features
    )

    obj$condition <- s$condition
    obj$region <- s$region
    obj$replicate <- s$replicate
    obj$sample_id <- s$sample_id
    obj$orig_file <- s$filename # GSM ID traceability

    obj
  })

  names(seurat_list) <- sample_table$sample_id

  if (!is.null(checkpoint)) save_checkpoint(seurat_list, checkpoint)

  seurat_list
}
