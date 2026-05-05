# 1. DATA LOADING ####

# Load CellBender-filtered h5 files as a named list of Seurat objects.
# File naming convention: GSM{ID}_{Ctrl|Inulin}_{region}_{replicate}_CellBender_feature_bc_matrix_filtered.h5
# Requires: output_path, CHECKPOINT_PREFIX (set in global_variables.R)

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
    data_path,
    pattern    = "_CellBender_feature_bc_matrix_filtered\\.h5$",
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


# 2. QUALITY CONTROL ####
#' Add standard QC metrics to a Seurat object
#' MT- pattern uppercase for HUMAN genes (differs from mouse ^mt-)
add_qc_metrics <- function(obj, pattern_mito = c("^MT-", "^mt-"), pattern_ribo = c("^RP[LS]", "^Rp[ls]")) {
  obj$percent.MT <- PercentageFeatureSet(obj, pattern = pattern_mito)
  obj$percent.ribosomal <- PercentageFeatureSet(obj, pattern = pattern_ribo)
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)
  obj
}


#' Build combined metadata data.frame from a named list of Seurat objects
build_combined_meta <- function(seurat_list) {
  do.call(dplyr::bind_rows, lapply(names(seurat_list), function(nm) {
    seurat_list[[nm]]@meta.data %>%
      dplyr::mutate(library = nm)
  }))
}


# 3. DOUBLETFINDER (per-library, with checkpoint recovery) ####


#' Detect doublets on one Seurat object using scDblFinder; returns a data.frame
#' with barcode, doublet classification, and score.
#' Results saved as: output_path/RData/doublets/df_<lib_name>.rds
#'
#' Replaces DoubletFinder (abandoned, recurring xtfrm.data.frame bug with
#' Seurat v5). scDblFinder is the Bioconductor standard and is already used
#' in code/Doublets_Finders.R.
#'
#' @param seurat_obj  Filtered (but not yet SCT-normalized) Seurat object
#' @param lib_name    Library identifier string
#' @param out_path    Base output path (output_path variable)
#' @return data.frame with columns: barcode, df_classification, pANN
run_doubletfinder <- function(seurat_obj, lib_name, out_path = output_path) {
  df_dir <- paste0(out_path, "RData/doublets/")
  df_file <- paste0(df_dir, "df_", lib_name, ".rds")
  dir.create(df_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(df_file)) {
    message(sprintf("[DF] Checkpoint found — loading: %s", lib_name))
    return(readRDS(df_file))
  }

  message(sprintf("[DF] Running scDblFinder: %s  (%d cells)", lib_name, ncol(seurat_obj)))
  tic()

  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder::scDblFinder(sce)

  result <- data.frame(
    barcode           = colnames(sce),
    df_classification = ifelse(sce$scDblFinder.class == "doublet", "Doublet", "Singlet"),
    pANN              = sce$scDblFinder.score,
    row.names         = colnames(sce),
    stringsAsFactors  = FALSE
  )

  saveRDS(result, file = df_file)
  toc()
  gc()

  message(sprintf(
    "  [DF] Singlets: %d | Doublets: %d",
    sum(result$df_classification == "Singlet"),
    sum(result$df_classification == "Doublet")
  ))
  result
}
