# GSE (gseGO) analysis for cluster-level DE results
#
# Expected input columns in cluster_markers:
#   gene, avg_log2FC, p_val_adj, cluster

run_gse_analysis <- function(
  cluster_markers,
  fc = NULL, # FC threshold; NULL → global FC
  pval = NULL, # p-value threshold; NULL → global p_val
  org_db = org.Mm.eg.db,
  key_type = "SYMBOL",
  show_n = 10, # categories shown in dot plots
  width = 8,
  height = 6,
  subdir = "gseGO"
) {
  fc <- if (is.null(fc)) get("FC", envir = .GlobalEnv) else fc
  pval <- if (is.null(pval)) get("p_val", envir = .GlobalEnv) else pval

  # --- Prepare ranked gene lists (cluster × direction) -------------------------
  filtered <- cluster_markers %>%
    dplyr::filter(abs(avg_log2FC) > fc, p_val_adj < pval) %>%
    dplyr::mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"))

  cluster_split <- split(filtered, paste0("Cluster_", filtered$cluster))

  split_list <- unlist(lapply(cluster_split, function(df) {
    list(
      UP   = dplyr::filter(df, avg_log2FC > 0) %>% dplyr::mutate(direction = "UP"),
      DOWN = dplyr::filter(df, avg_log2FC < 0) %>% dplyr::mutate(direction = "DOWN")
    )
  }), recursive = FALSE)
  split_list <- Filter(function(x) nrow(x) > 0, split_list)

  gene_lists <- lapply(split_list, function(df) {
    gl <- df$avg_log2FC
    names(gl) <- df$gene
    gl <- na.omit(gl)
    sort(gl, decreasing = TRUE)
  })

  # --- Run gseGO ---------------------------------------------------------------
  results <- lapply(gene_lists, function(gene_set) {
    tryCatch(
      clusterProfiler::gseGO(
        geneList      = gene_set,
        ont           = "ALL",
        keyType       = key_type,
        minGSSize     = 10,
        maxGSSize     = 1000,
        pvalueCutoff  = 0.05,
        verbose       = FALSE,
        OrgDb         = org_db,
        pAdjustMethod = "none"
      ),
      error = function(e) {
        message("gseGO failed: ", e$message)
        NULL
      }
    )
  })
  results <- Filter(Negate(is.null), results)

  save(results, file = file.path(output_path, "RData", "results_list_gseGO_ALL.RData"))

  # Save result tables to xlsx
  results_df <- lapply(results, function(x) x@result)
  results_df <- Filter(function(x) nrow(x) > 0, results_df)
  writexl::write_xlsx(
    results_df,
    file.path(output_path, "tables", "results_df_gseGO_ALL.xlsx")
  )

  # --- Helper: clean plot label ------------------------------------------------
  .clean_name <- function(nm) gsub("_", " ", gsub("\\.[A-Z]+$", "", nm))

  # --- Dotplots ----------------------------------------------------------------
  for (nm in names(results)) {
    p <- tryCatch(
      enrichplot::dotplot(results[[nm]], showCategory = show_n, split = ".sign") +
        ggplot2::facet_grid(. ~ .sign) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::ggtitle(paste("GSE Dot plot —", .clean_name(nm))),
      error = function(e) NULL
    )
    if (!is.null(p)) {
      print(p)
      save_plot(paste0("gseGO_dotplot_", nm), p,
        width = width, height = height, subdir = subdir
      )
    }
  }

  # --- CNET plots --------------------------------------------------------------
  for (nm in names(results)) {
    p <- tryCatch(
      {
        x2 <- clusterProfiler::simplify(results[[nm]])
        enrichplot::cnetplot(x2, foldChange = results[[nm]]@geneList) +
          ggplot2::ggtitle(paste("GSE CNET plot —", .clean_name(nm)))
      },
      error = function(e) NULL
    )
    if (!is.null(p)) {
      print(p)
      save_plot(paste0("gseGO_cnet_", nm), p,
        width = width, height = height, subdir = subdir
      )
    }
  }

  # --- Heatmap plots -----------------------------------------------------------
  for (nm in names(results)) {
    p <- tryCatch(
      enrichplot::heatplot(results[[nm]], foldChange = results[[nm]]@geneList) +
        ggplot2::ggtitle(paste("GSE Heatmap plot —", .clean_name(nm))),
      error = function(e) NULL
    )
    if (!is.null(p)) {
      print(p)
      save_plot(paste0("gseGO_heatmap_", nm), p,
        width = width, height = height * 1.5, subdir = subdir
      )
    }
  }

  # --- Ridge plots -------------------------------------------------------------
  for (nm in names(results)) {
    p <- tryCatch(
      enrichplot::ridgeplot(results[[nm]]) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::ggtitle(paste("GSE Ridge plot —", .clean_name(nm))),
      error = function(e) NULL
    )
    if (!is.null(p)) {
      print(p)
      save_plot(paste0("gseGO_ridgeplot_", nm), p,
        width = width, height = height, subdir = subdir
      )
    }
  }
  # --- Lolliplots --------------------------------------------------------------
  results_df_all <- lapply(results, function(x) {
    x@result
  })

  for (nm in names(results_df_all)) {
    df <- results_df_all[[nm]]
    if (nrow(df) == 0) next
    p <- df %>%
      dplyr::mutate(Description = reorder(Description, enrichmentScore)) %>%
      dplyr::slice_max(enrichmentScore, n = 15) %>%
      ggplot2::ggplot(ggplot2::aes(x = Description, y = enrichmentScore, colour = p.adjust)) +
      ggplot2::geom_segment(ggplot2::aes(xend = Description, y = 0, yend = enrichmentScore)) +
      ggplot2::geom_point(ggplot2::aes(size = setSize)) +
      ggplot2::scale_color_viridis_c(option = "viridis") +
      ggplot2::facet_wrap(ONTOLOGY ~ ., scales = "free", ncol = 1) +
      ggplot2::coord_flip() +
      ggplot2::theme_gray() +
      ggplot2::labs(size = "N. of genes", x = "GO term", y = "Gene Ratio") +
      ggplot2::ggtitle(paste("Lolliplot —", .clean_name(nm)))
    print(p)
    save_plot(paste0("enrichGO_lolliplot_", nm), p,
      width = width, height = height * 1.5, subdir = subdir
    )
  }
  invisible(results)
}
