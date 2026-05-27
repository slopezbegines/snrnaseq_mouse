# Heatmap utilities for snRNAseq cluster-level visualisation
#
# Two public entry points:
#
#   heatmap_de(seurat_obj, cluster_markers, ...)
#       One heatmap per cluster — rows = DE genes for that cluster.
#
#   heatmap_genes(seurat_obj, genes, ...)
#       Rows = user-supplied gene vector.
#       per_cluster = FALSE  → single combined heatmap across selected clusters
#       per_cluster = TRUE   → one heatmap per selected cluster
#
# Both functions:
#   • accept a `clusters` argument to restrict which clusters are plotted
#   • share identical palette / scaling / layout options
#   • call save_plot() to write TIFF + PDF to figures/<subdir>/
#
# Usage:
#   source("code/heatmap.R")
#   heatmap_de(data_select_bp, cluster_markers_by_condition)
#   heatmap_genes(data_select_bp, c("Gfap","Cx3cr1","Mbp"), clusters = c(0, 3))

# ── Private helpers ────────────────────────────────────────────────────────────

# Build a circlize colour function from an RColorBrewer or viridis palette.
.hm_col_fun <- function(mat, palette, rev_palette, scale_rows) {
  viridis_opts <- c("viridis", "magma", "plasma", "inferno", "cividis")
  if (tolower(palette) %in% viridis_opts) {
    cols <- viridisLite::viridis(100, option = tolower(palette))
  } else {
    n    <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    cols <- RColorBrewer::brewer.pal(n, palette)
  }
  if (rev_palette) cols <- rev(cols)
  key <- cols[c(1L, ceiling(length(cols) / 2L), length(cols))]
  if (scale_rows) {
    circlize::colorRamp2(c(-2.5, 0, 2.5), key)
  } else {
    q <- quantile(mat, c(0.02, 0.5, 0.98), na.rm = TRUE)
    circlize::colorRamp2(as.numeric(q), key)
  }
}

# Row-wise z-score capped at ±2.5.
.hm_scale <- function(mat) {
  mat <- t(scale(t(mat)))
  mat[is.nan(mat)] <- 0
  pmin(pmax(mat, -2.5), 2.5)
}

# Generate exactly n discrete colours by interpolating a RColorBrewer palette.
# Avoids NA entries when n > brewer.pal maxcolors (e.g. Set2 caps at 8).
.hm_discrete_colors <- function(n, palette) {
  max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  base  <- RColorBrewer::brewer.pal(max_n, palette)
  colorRampPalette(base)(n)
}

# Top annotation bar for a single condition vector.
.hm_cond_anno <- function(cond_vals, split_by) {
  levs   <- unique(cond_vals)
  colors <- setNames(.hm_discrete_colors(length(levs), "Set1"), levs)
  ComplexHeatmap::HeatmapAnnotation(
    condition               = cond_vals,
    col                     = list(condition = colors),
    annotation_name_side    = "left",
    annotation_legend_param = list(title = split_by)
  )
}

# Two-row top annotation: cluster + condition (used by the combined heatmap).
.hm_combined_anno <- function(cluster_vals, cond_vals, all_clusters, split_by) {
  cl_levs     <- all_clusters
  cond_levs   <- unique(cond_vals)
  cl_colors   <- setNames(.hm_discrete_colors(length(cl_levs),   "Set2"), cl_levs)
  cond_colors <- setNames(.hm_discrete_colors(length(cond_levs), "Set1"), cond_levs)
  ComplexHeatmap::HeatmapAnnotation(
    cluster               = cluster_vals,
    condition             = cond_vals,
    col                   = list(cluster = cl_colors, condition = cond_colors),
    annotation_name_side  = "left"
  )
}

# Assemble a ComplexHeatmap::Heatmap object from a prepped matrix.
.hm_build <- function(expr, title, col_fun, top_annotation,
                       cluster_rows, cluster_cols, show_colnames, scale_rows) {
  ComplexHeatmap::Heatmap(
    expr,
    name              = if (scale_rows) "z-score" else "Expression",
    col               = col_fun,
    top_annotation    = top_annotation,
    cluster_rows      = cluster_rows,
    cluster_columns   = cluster_cols,
    show_column_names = show_colnames,
    show_row_names    = TRUE,
    row_names_gp      = grid::gpar(fontsize = 8),
    column_title      = title,
    column_title_gp   = grid::gpar(fontsize = 12, fontface = "bold"),
    use_raster        = ncol(expr) > 500
  )
}

# Return barcodes whose active ident equals `cl`.
.hm_cells <- function(seurat_obj, cl) {
  names(Idents(seurat_obj))[as.character(Idents(seurat_obj)) == as.character(cl)]
}

# Save a ComplexHeatmap via the function branch of save_plot()
# (ggsave cannot dispatch grid.draw on Heatmap objects).
.hm_save <- function(name, ht, width, height, subdir) {
  save_plot(name, function() ComplexHeatmap::draw(ht),
            width = width, height = height, subdir = subdir)
}

# Extract and optionally scale an expression matrix.
.hm_expr <- function(seurat_obj, genes, cells, assay, layer, scale_rows) {
  mat <- as.matrix(
    Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)[
      genes, cells, drop = FALSE
    ]
  )
  if (scale_rows) .hm_scale(mat) else mat
}

# Resolve and validate a `clusters` argument against a reference set.
.hm_resolve_clusters <- function(clusters, reference) {
  reference <- as.character(reference)
  if (is.null(clusters)) return(reference)
  sel <- as.character(clusters)
  found <- sel[sel %in% reference]
  if (length(found) == 0)
    stop("None of the requested clusters (", paste(sel, collapse = ", "),
         ") found in the Seurat object / DE table.")
  found
}

# ── heatmap_de ─────────────────────────────────────────────────────────────────

#' Heatmap of DE genes for each cluster.
#'
#' @param seurat_obj   Seurat object (active ident must match cluster column).
#' @param cluster_markers  Data frame with columns gene, avg_log2FC, p_val_adj, cluster.
#' @param clusters     NULL = all clusters; c(0,2,5) or c("Oligo","Micro") to subset.
#' @param fc           FC threshold. NULL → global FC.
#' @param pval         p_val_adj threshold. NULL → global p_val.
#' @param n_top        Top N genes per cluster by |avg_log2FC|. NULL = all passing.
#' @param assay        Seurat assay (default "RNA").
#' @param layer        "data" (normalised) or "scale.data" (pre-scaled).
#' @param split_by     Metadata column used to colour and order cells (default "condition").
#' @param palette      RColorBrewer palette name or viridis option ("viridis","magma",...).
#' @param rev_palette  Reverse the colour scale (default TRUE).
#' @param scale_rows   z-score genes across cells, capped ±2.5 (default TRUE).
#' @param cluster_rows Hierarchically cluster genes (default TRUE).
#' @param cluster_cols Hierarchically cluster cells (default FALSE).
#' @param show_colnames Show cell barcodes on x-axis (default FALSE).
#' @param width,height Figure dimensions in inches.
#' @param subdir       Subdirectory under figures/ for save_plot().
heatmap_de <- function(
  seurat_obj,
  cluster_markers,
  clusters      = NULL,
  fc            = NULL,
  pval          = NULL,
  n_top         = 20,
  assay         = "RNA",
  layer         = "data",
  split_by      = "condition",
  palette       = "RdBu",
  rev_palette   = TRUE,
  scale_rows    = TRUE,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  show_colnames = FALSE,
  width         = 10,
  height        = 8,
  subdir        = "heatmaps"
) {
  fc   <- if (is.null(fc))   get("FC",    envir = .GlobalEnv) else fc
  pval <- if (is.null(pval)) get("p_val", envir = .GlobalEnv) else pval

  de <- cluster_markers %>%
    dplyr::filter(abs(avg_log2FC) > fc, p_val_adj < pval)

  if (nrow(de) == 0) {
    message("No DE genes pass FC > ", fc, " / p_val_adj < ", pval, ". No heatmaps generated.")
    return(invisible(NULL))
  }

  cl_ids <- .hm_resolve_clusters(clusters, sort(unique(de$cluster)))

  for (cl in cl_ids) {
    cl_de <- de %>% dplyr::filter(as.character(cluster) == as.character(cl))
    if (!is.null(n_top))
      cl_de <- cl_de %>% dplyr::slice_max(abs(avg_log2FC), n = n_top, with_ties = FALSE)

    genes <- cl_de %>%
      dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
      dplyr::pull(gene) %>%
      intersect(rownames(seurat_obj))

    if (length(genes) == 0) { message("Cluster ", cl, ": no genes — skipping."); next }

    cells <- .hm_cells(seurat_obj, cl)
    if (length(cells) == 0) { message("Cluster ", cl, ": no cells — skipping."); next }

    expr      <- .hm_expr(seurat_obj, genes, cells, assay, layer, scale_rows)
    meta      <- seurat_obj@meta.data[cells, , drop = FALSE]
    cond_vals <- if (split_by %in% colnames(meta)) meta[[split_by]] else NULL

    if (!is.null(cond_vals)) {
      ord       <- order(cond_vals)
      expr      <- expr[, ord, drop = FALSE]
      cond_vals <- cond_vals[ord]
    }

    ha  <- if (!is.null(cond_vals)) .hm_cond_anno(cond_vals, split_by) else NULL
    ht  <- .hm_build(
      expr, title = paste0("Cluster ", cl, "  |  ", length(genes), " DE genes"),
      col_fun       = .hm_col_fun(expr, palette, rev_palette, scale_rows),
      top_annotation = ha,
      cluster_rows  = cluster_rows, cluster_cols  = cluster_cols,
      show_colnames = show_colnames, scale_rows    = scale_rows
    )
    .hm_save(paste0("heatmap_de_cluster_", cl), ht, width, height, subdir)
  }
  invisible(NULL)
}

# ── heatmap_genes ──────────────────────────────────────────────────────────────

#' Heatmap for a user-defined gene list across selected clusters.
#'
#' @param seurat_obj  Seurat object.
#' @param genes       Character vector of gene names.
#' @param clusters    NULL = all; c(0,2,5) or c("Oligo","Micro") to subset.
#' @param per_cluster FALSE → single combined heatmap annotated by cluster + condition.
#'                    TRUE  → one heatmap per selected cluster.
#' @inheritParams heatmap_de
heatmap_genes <- function(
  seurat_obj,
  genes,
  clusters      = NULL,
  per_cluster   = FALSE,
  assay         = "RNA",
  layer         = "data",
  split_by      = "condition",
  palette       = "RdBu",
  rev_palette   = TRUE,
  scale_rows    = TRUE,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  show_colnames = FALSE,
  width         = 10,
  height        = 8,
  subdir        = "heatmaps"
) {
  genes <- intersect(genes, rownames(seurat_obj))
  if (length(genes) == 0) {
    message("None of the supplied genes found in the Seurat object.")
    return(invisible(NULL))
  }

  cl_ids <- .hm_resolve_clusters(
    clusters,
    sort(unique(as.character(Idents(seurat_obj))))
  )

  if (per_cluster) {
    # ── One heatmap per selected cluster ──────────────────────────────────────
    for (cl in cl_ids) {
      cells <- .hm_cells(seurat_obj, cl)
      if (length(cells) == 0) { message("Cluster ", cl, ": no cells — skipping."); next }

      expr      <- .hm_expr(seurat_obj, genes, cells, assay, layer, scale_rows)
      meta      <- seurat_obj@meta.data[cells, , drop = FALSE]
      cond_vals <- if (split_by %in% colnames(meta)) meta[[split_by]] else NULL

      if (!is.null(cond_vals)) {
        ord       <- order(cond_vals)
        expr      <- expr[, ord, drop = FALSE]
        cond_vals <- cond_vals[ord]
      }

      ha <- if (!is.null(cond_vals)) .hm_cond_anno(cond_vals, split_by) else NULL
      ht <- .hm_build(
        expr, title = paste0("Cluster ", cl, "  |  ", length(genes), " genes"),
        col_fun       = .hm_col_fun(expr, palette, rev_palette, scale_rows),
        top_annotation = ha,
        cluster_rows  = cluster_rows, cluster_cols  = cluster_cols,
        show_colnames = show_colnames, scale_rows    = scale_rows
      )
      .hm_save(paste0("heatmap_genes_cluster_", cl), ht, width, height, subdir)
    }

  } else {
    # ── Single combined heatmap across all selected clusters ──────────────────
    cell_info <- do.call(rbind, lapply(cl_ids, function(cl) {
      cells <- .hm_cells(seurat_obj, cl)
      if (length(cells) == 0) return(NULL)
      meta  <- seurat_obj@meta.data[cells, , drop = FALSE]
      cond  <- if (split_by %in% colnames(meta)) meta[[split_by]] else NA_character_
      data.frame(cell = cells, cluster = cl, condition = as.character(cond),
                 stringsAsFactors = FALSE)
    }))

    if (is.null(cell_info) || nrow(cell_info) == 0) {
      message("No cells found for the selected clusters.")
      return(invisible(NULL))
    }

    # Order cells: cluster first, then condition within cluster
    cell_info <- cell_info[order(cell_info$cluster, cell_info$condition), ]

    expr <- .hm_expr(seurat_obj, genes, cell_info$cell, assay, layer, scale_rows)

    ha <- .hm_combined_anno(cell_info$cluster, cell_info$condition, cl_ids, split_by)
    ht <- .hm_build(
      expr,
      title = paste0(length(genes), " genes  |  clusters: ", paste(cl_ids, collapse = ", ")),
      col_fun       = .hm_col_fun(expr, palette, rev_palette, scale_rows),
      top_annotation = ha,
      cluster_rows  = cluster_rows, cluster_cols  = cluster_cols,
      show_colnames = show_colnames, scale_rows    = scale_rows
    )
    .hm_save("heatmap_genes_combined", ht, width, height, subdir)
  }

  invisible(NULL)
}

# ── heatmap_mean ───────────────────────────────────────────────────────────────

#' Heatmap of mean expression per cluster × condition.
#'
#' Each column is one (cluster, condition) group — the mean normalised expression
#' across all cells in that group. Much more compact than per-cell heatmaps and
#' directly comparable across conditions.
#'
#' @param seurat_obj Seurat object.
#' @param genes      Character vector of gene names.
#' @param clusters   NULL = all; c(0,2,5) or c("Oligo","Micro") to subset.
#' @param split_by   Metadata column that defines conditions (default "condition").
#' @param assay      Seurat assay (default "RNA").
#' @param layer      Layer to average: "data" (log-normalised, recommended).
#' @param scale_rows z-score each gene across all cluster×condition columns (default TRUE).
#' @param cluster_rows Hierarchically cluster genes (default TRUE).
#' @param cluster_cols Hierarchically cluster groups (default FALSE — keeps order).
#' @param palette    RColorBrewer palette name or viridis option.
#' @param rev_palette Reverse the colour scale (default TRUE).
#' @param width,height Figure dimensions in inches.
#' @param subdir     Subdirectory under figures/ for save_plot().
heatmap_mean <- function(
  seurat_obj,
  genes,
  clusters      = NULL,
  split_by      = "condition",
  assay         = "RNA",
  layer         = "data",
  scale_rows    = TRUE,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  palette       = "RdBu",
  rev_palette   = TRUE,
  width         = 8,
  height        = 6,
  subdir        = "heatmaps"
) {
  genes <- intersect(genes, rownames(seurat_obj))
  if (length(genes) == 0) {
    message("None of the supplied genes found in the Seurat object.")
    return(invisible(NULL))
  }

  cl_ids <- .hm_resolve_clusters(
    clusters,
    sort(unique(as.character(Idents(seurat_obj))))
  )

  # Build per-(cluster × condition) mean expression matrix ─────────────────────
  meta      <- seurat_obj@meta.data
  expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)

  if (!split_by %in% colnames(meta))
    stop("'split_by' column '", split_by, "' not found in Seurat metadata.")

  meta$..cl   <- as.character(Idents(seurat_obj))
  meta$..cond <- as.character(meta[[split_by]])
  meta_sel    <- meta[meta$..cl %in% cl_ids, , drop = FALSE]

  groups <- unique(meta_sel[, c("..cl", "..cond")])
  groups <- groups[order(groups$..cl, groups$..cond), ]

  if (nrow(groups) == 0) {
    message("No cells found for the selected clusters.")
    return(invisible(NULL))
  }

  # One column per group: mean expression across its cells
  avg_mat <- vapply(seq_len(nrow(groups)), function(i) {
    cells <- rownames(meta_sel)[
      meta_sel$..cl == groups$..cl[i] & meta_sel$..cond == groups$..cond[i]
    ]
    rowMeans(as.matrix(expr_data[genes, cells, drop = FALSE]))
  }, numeric(length(genes)))

  rownames(avg_mat) <- genes
  colnames(avg_mat) <- paste0(groups$..cl, " | ", groups$..cond)

  # Scale rows (z-score across all cluster×condition columns) ──────────────────
  if (scale_rows) avg_mat <- .hm_scale(avg_mat)

  # Column annotation: cluster + condition bars ─────────────────────────────────
  cl_levs     <- unique(groups$..cl)
  cond_levs   <- unique(groups$..cond)
  cl_colors   <- setNames(.hm_discrete_colors(length(cl_levs),   "Set2"), cl_levs)
  cond_colors <- setNames(.hm_discrete_colors(length(cond_levs), "Set1"), cond_levs)
  ha <- ComplexHeatmap::HeatmapAnnotation(
    cluster              = groups$..cl,
    condition            = groups$..cond,
    col                  = list(cluster = cl_colors, condition = cond_colors),
    annotation_name_side = "left"
  )

  ht <- .hm_build(
    avg_mat,
    title         = paste0("Mean expression  |  ", length(genes), " genes"),
    col_fun       = .hm_col_fun(avg_mat, palette, rev_palette, scale_rows),
    top_annotation = ha,
    cluster_rows  = cluster_rows,
    cluster_cols  = cluster_cols,
    show_colnames = TRUE,      # groups are few enough to always label
    scale_rows    = scale_rows
  )

  .hm_save("heatmap_mean", ht, width, height, subdir)
  invisible(NULL)
}
