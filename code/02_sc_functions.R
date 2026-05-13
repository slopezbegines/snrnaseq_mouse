# 0. UTILITIES ####

#' Convert Ensembl IDs to gene symbols in a BPCells (or standard) matrix.
#' Replaces Azimuth:::ConvertEnsembleToSymbol without the Signac dependency.
#' Unmapped IDs are kept as-is (no data lost).
convert_ensembl_to_symbol <- function(mat, species = "mouse") {
  pkg <- switch(species,
    mouse = "org.Mm.eg.db",
    human = "org.Hs.eg.db",
    stop("species must be 'mouse' or 'human'")
  )
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install with: BiocManager::install('", pkg, "')")
  }
  db <- getExportedValue(pkg, pkg)
  ids <- rownames(mat)
  symbols <- suppressMessages(
    AnnotationDbi::mapIds(db,
      keys = ids, column = "SYMBOL",
      keytype = "ENSEMBL", multiVals = "first"
    )
  )
  symbols[is.na(symbols)] <- ids[is.na(symbols)]
  rownames(mat) <- unname(symbols)
  # For duplicate symbols keep the Ensembl entry with the highest total counts
  # (mirrors Azimuth:::ConvertEnsembleToSymbol; pseudogenes tend to rank lower)
  if (anyDuplicated(rownames(mat))) {
    rs <- rowSums(mat)
    ord <- order(rownames(mat), -rs)
    mat <- mat[ord, ]
    mat <- mat[!duplicated(rownames(mat)), ]
  }
  mat
}


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
#' @param species      Species name for gene symbol conversion (e.g. "human" or "mouse")
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
  species = species,
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

# BPCells-backed variant of load_h5_samples (code/02_sc_functions.R).
# Count matrices are written to disk under output_path/BPCells/<sample_id>/
# so that Seurat operates on-disk, keeping peak RAM low.
#
# Requires:
#   - BPCells package  (install: BiocManager::install("BPCells"))
#   - output_path, DATA_DIR  (set in global_variables.R)
#   - check_checkpoint(), load_checkpoint(), save_checkpoint()  (01_aux_functions.R)
#
# NOTE: the .rds checkpoint only stores the Seurat object skeleton and
# metadata; the actual count data lives in the BPCells directories on disk.
# Both must be accessible when loading a saved checkpoint.


#' Load CellBender h5 files as BPCells-backed Seurat objects
#'
#' Each sample's count matrix is stored on disk under
#' bpcells_root/<sample_id>/. If that directory already exists it is
#' reused without re-writing (safe to call across sessions).
#'
#' File naming convention (same as load_h5_samples):
#'   GSM{ID}_{Ctrl|Inulin}_{region}_{replicate}_CellBender_feature_bc_matrix_filtered.h5
#'
#' @param data_path     Directory containing the CellBender h5 files
#' @param regions       Regions to load; NULL loads all four
#' @param min.cells     Passed to CreateSeuratObject()
#' @param min.features  Passed to CreateSeuratObject()
#' @param checkpoint    Checkpoint name (e.g. "01_raw_data_bp"); NULL skips
#' @param bpcells_root  Root for on-disk BPCells matrices
#'
#' @return Named list of Seurat objects backed by on-disk BPCells matrices
#'
#' @examples
#' seurat_list <- load_h5_samples_bpcells(data_path, checkpoint = "01_raw_data_bp")
load_h5_samples_bpcells <- function(
  data_path = DATA_DIR,
  regions = c("forebrain", "interbrain", "brainstem", "cerebellum"),
  min.cells = 3,
  min.features = 200,
  checkpoint = NULL,
  bpcells_root = file.path(output_path, "RData/BPCells"),
  specie = c("mouse", "human")
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

  dir.create(bpcells_root, recursive = TRUE, showWarnings = FALSE)

  message(sprintf(
    "Loading %d sample(s) [BPCells] | regions: %s",
    nrow(sample_table), paste(regions, collapse = ", ")
  ))

  seurat_list <- lapply(seq_len(nrow(sample_table)), function(i) {
    s <- sample_table[i, ]
    bp_dir <- file.path(bpcells_root, s$sample_id)

    message(sprintf("  [%d/%d] %s (%s)", i, nrow(sample_table), s$sample_id, s$gsm_id))

    if (dir.exists(bp_dir)) {
      message(sprintf("    BPCells dir found — reusing: %s", bp_dir))
      mat <- BPCells::open_matrix_dir(bp_dir)
    } else {
      message(sprintf("    Writing BPCells matrix to: %s", bp_dir))
      mat_h5 <- BPCells::open_matrix_10x_hdf5(s$filepath)
      mat <- BPCells::write_matrix_dir(mat_h5, dir = bp_dir)
    }

    mat <- convert_ensembl_to_symbol(mat, species = specie)

    obj <- CreateSeuratObject(
      counts       = mat,
      project      = s$sample_id,
      min.cells    = min.cells,
      min.features = min.features
    )

    obj$condition <- s$condition
    obj$region <- s$region
    obj$replicate <- s$replicate
    obj$sample_id <- s$sample_id
    obj$orig_file <- s$filename

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


#' Run scDblFinder on a single BPCells-backed Seurat object.
#' BPCells counts are materialised to dgCMatrix per-sample — safe at 16 GB
#' because each sample is processed independently before merging.
#' Results cached under out_path/RData/doublets/df_<lib_name>.rds.
#'
#' @return data.frame with columns: barcode, df_classification, pANN
run_doubletfinder_bp <- function(seurat_obj, lib_name, out_path = output_path) {
  df_dir <- file.path(out_path, "RData", "doublets")
  df_file <- file.path(df_dir, paste0("df_", lib_name, ".rds"))
  dir.create(df_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(df_file)) {
    message(sprintf("[DF] Checkpoint found — loading: %s", lib_name))
    return(readRDS(df_file))
  }

  message(sprintf("[DF] scDblFinder: %s (%d cells)", lib_name, ncol(seurat_obj)))
  tic()

  # Materialise BPCells → dgCMatrix for scDblFinder (SCE does not support BPCells).
  counts_mat <- as(LayerData(seurat_obj, layer = "counts"), "dgCMatrix")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_mat)
  )
  rm(counts_mat)
  gc()

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


#' Add doublet classification and score to Seurat object metadata
apply_doublet_labels <- function(seurat_obj, df_result) {
  shared <- intersect(colnames(seurat_obj), rownames(df_result))
  seurat_obj$df_classification <- NA_character_
  seurat_obj$doublet_score <- NA_real_
  seurat_obj$df_classification[shared] <- df_result[shared, "df_classification"]
  seurat_obj$doublet_score[shared] <- df_result[shared, "pANN"]
  seurat_obj
}


#' Remove doublets from a labelled Seurat object
filter_doublets <- function(seurat_obj) {
  n_before <- ncol(seurat_obj)
  seurat_obj <- seurat_obj[, seurat_obj$df_classification == "Singlet"]
  message(sprintf(
    "  [DOUBLETS] %s: %d → %d cells (-%d doublets)",
    unique(seurat_obj$sample_id), n_before, ncol(seurat_obj), n_before - ncol(seurat_obj)
  ))
  seurat_obj
}



#' Full doublet-removal step: checkpoint → label → filter → save → summary.
#' Caller should run rm(seurat_filtered, df_results); gc() after this returns
#' to release the original references from the calling environment.
build_seurat_singlets <- function(seurat_filtered, df_results) {
  if (check_checkpoint("06_seurat_singlets")) {
    return(load_checkpoint("06_seurat_singlets"))
  }

  seurat_singlets <- mapply(
    function(obj, df) filter_doublets(apply_doublet_labels(obj, df)),
    seurat_filtered, df_results,
    SIMPLIFY = FALSE
  )

  rm(seurat_filtered, df_results)
  gc()
  save_checkpoint(seurat_singlets, "06_seurat_singlets")

  cat("Cells after doublet removal:", sum(sapply(seurat_singlets, ncol)), "\n")
  seurat_singlets
}


# 4. MERGE ----------------------------------------------------------------

#' Merge a named list of Seurat objects into one object.
#' In Seurat v5, each input becomes a separate layer in the merged object.
merge_samples <- function(seurat_list) {
  message(sprintf("[MERGE] Merging %d samples...", length(seurat_list)))

  merged <- merge(
    seurat_list[[1]],
    y            = seurat_list[-1],
    add.cell.ids = names(seurat_list),
    merge.data   = FALSE
  )

  gc()
  message(sprintf("[MERGE] %d cells | %d genes", ncol(merged), nrow(merged)))
  merged
}

# Convert SCTAssay from SCTransform to BPCells on disk.
# SCTransform always creates SCTAssay (v3), which rejects BPCells MatrixDir layers.
# Idempotent: if BPCells dirs already exist, they are reused (safe to re-run).
sct_to_bpcells <- function(obj, nm) {
  sct_assay <- obj[["SCT"]]
  for (layer_nm in c("counts", "data")) {
    bp_dir <- file.path(bpcells_sct_root, nm, layer_nm)
    if (!dir.exists(bp_dir)) {
      mat <- slot(sct_assay, layer_nm)
      bp_mat <- BPCells::write_matrix_dir(mat, dir = bp_dir)
      rm(mat)
      gc()
    } else {
      bp_mat <- BPCells::open_matrix_dir(bp_dir)
    }
    # slot<- calls checkSlotAssignment: SCTAssay rejects MatrixDir (not AnyMatrix).
    # attr<- writes directly to the attribute list, bypassing S4 type validation.
    attr(sct_assay, layer_nm) <- bp_mat
    rm(bp_mat)
    gc()
  }
  slot(obj, "assays")[["SCT"]] <- sct_assay # bypass [[<- coercion to Assay5
  rm(sct_assay)
  gc()
  obj
}

# Offload SCTAssay (v3) counts + data to BPCells on disk, keeping it as SCTAssay.
# Uses slot<- to bypass the LayerData<- validator that rejects BPCells in v3 assays.
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
# 4. HARMONY INTEGRATION --------------------------------------------------

#' PCA → Harmony integration → JoinLayers.
#' Requires VariableFeatures to already be set on the merged object
#' (done in the SCT merge step via lapply(all_sct, VariableFeatures) + table).
#'
#' @param group_vars Metadata column(s) used as Harmony batch variable(s)
#' @param npcs       Number of PCs to compute
run_harmony_integration <- function(
  obj,
  group_vars = "sample_id",
  npcs = N_PCS_INTEGRATION
) {
  message(sprintf("[PCA] Running PCA (%d PCs)...", npcs))
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  gc()

  message(sprintf("[HARMONY] Integrating by: %s", paste(group_vars, collapse = ", ")))
  obj <- IntegrateLayers(
    obj,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    group.by.vars = group_vars,
    normalization.method = "SCT",
    verbose = TRUE
  )
  gc()

  # Required in Seurat v5: collapse split layers after integration
  message("[LAYERS] Joining layers...")
  obj <- JoinLayers(obj)
  gc()

  obj
}
