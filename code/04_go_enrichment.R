# GO over-representation analysis for cluster-level DE results
#
# Expected input columns in cluster_markers:
#   gene, avg_log2FC, p_val_adj, cluster

run_go_enrichment <- function(
  cluster_markers,
  fc = NULL, # FC threshold; NULL → global FC
  pval = NULL, # p-value threshold; NULL → global p_val
  org_db = "org.Mm.eg.db",
  key_type = "SYMBOL"
) {
  fc <- if (is.null(fc)) get("FC", envir = .GlobalEnv) else fc
  pval <- if (is.null(pval)) get("p_val", envir = .GlobalEnv) else pval

  # --- Prepare split list (cluster × direction) --------------------------------
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

  # --- Run enrichGO ------------------------------------------------------------
  go_results <- lapply(split_list, function(gene_set) {
    tryCatch(
      clusterProfiler::enrichGO(
        gene          = gene_set$gene,
        OrgDb         = org_db,
        keyType       = key_type,
        ont           = "ALL",
        pAdjustMethod = "none"
      ),
      error = function(e) {
        message("enrichGO failed: ", e$message)
        NULL
      }
    )
  })
  go_results <- Filter(Negate(is.null), go_results)

  save(go_results, file = file.path(output_path, "RData", "results_list_enrichGO_ALL.RData"))

  # Save result tables to xlsx
  results_df <- lapply(go_results, function(x) x@result)
  results_df <- Filter(function(x) nrow(x) > 0, results_df)
  writexl::write_xlsx(
    results_df,
    file.path(output_path, "tables", "results_df_enrichGO_ALL.xlsx")
  )
  return(go_results)
}

plot_go_enrichment <- function(go_results,
                               show_n = 10,
                               width = 8,
                               height = 6,
                               subdir = "enrichGO") {
  # --- Helper: clean plot label ------------------------------------------------
  .clean_name <- function(nm) gsub("_", " ", gsub("\\.[A-Z]+$", "", nm))

  rdata_path <- file.path(output_path, "RData", "results_list_enrichGO_ALL.RData")
  if (file.exists(normalizePath(rdata_path, mustWork = FALSE))) {
    load(normalizePath(rdata_path, mustWork = FALSE), envir = .GlobalEnv)
    go_results <- get("go_results", envir = .GlobalEnv)

    # --- Barplots ----------------------------------------------------------------
    for (nm in names(go_results)) {
      p <- tryCatch(
        barplot(go_results[[nm]], showCategory = show_n) +
          ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free") +
          ggplot2::ggtitle(paste("Bar plot —", .clean_name(nm))),
        error = function(e) NULL
      )
      if (!is.null(p)) {
        print(p)
        save_plot(paste0("enrichGO_barplot_", nm), p,
          width = width, height = height, subdir = subdir
        )
      }
    }

    # --- Dotplots ----------------------------------------------------------------
    for (nm in names(go_results)) {
      p <- tryCatch(
        enrichplot::dotplot(go_results[[nm]], showCategory = show_n) +
          ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free") +
          ggplot2::ggtitle(paste("Dot plot —", .clean_name(nm))),
        error = function(e) NULL
      )
      if (!is.null(p)) {
        print(p)
        save_plot(paste0("enrichGO_dotplot_", nm), p,
          width = width, height = height, subdir = subdir
        )
      }
    }

    # --- CNET plots --------------------------------------------------------------
    for (nm in names(go_results)) {
      p <- tryCatch(
        {
          x2 <- clusterProfiler::simplify(go_results[[nm]])
          enrichplot::cnetplot(x2) +
            ggplot2::ggtitle(paste("CNET plot —", .clean_name(nm)))
        },
        error = function(e) NULL
      )
      if (!is.null(p)) {
        print(p)
        save_plot(paste0("enrichGO_cnet_", nm), p,
          width = width, height = height, subdir = subdir
        )
      }
    }

    # --- UpSet plots -------------------------------------------------------------
    for (nm in names(go_results)) {
      p <- tryCatch(
        enrichplot::upsetplot(go_results[[nm]]) +
          ggplot2::ggtitle(paste("UpSet plot —", .clean_name(nm))),
        error = function(e) NULL
      )
      if (!is.null(p)) {
        print(p)
        save_plot(paste0("enrichGO_upset_", nm), p,
          width = width, height = height, subdir = subdir
        )
      }
    }

    # --- Heatmap plots -----------------------------------------------------------
    for (nm in names(go_results)) {
      p <- tryCatch(
        enrichplot::heatplot(go_results[[nm]]) +
          ggplot2::ggtitle(paste("Heatmap plot —", .clean_name(nm))),
        error = function(e) NULL
      )
      if (!is.null(p)) {
        print(p)
        save_plot(paste0("enrichGO_heatmap_", nm), p,
          width = width, height = height * 1.5, subdir = subdir
        )
      }
    }

    # --- Lolliplots --------------------------------------------------------------
    go_results_df_all <- lapply(go_results, function(x) {
      x@result %>%
        dplyr::mutate(GeneRatio = purrr::map_dbl(
          strsplit(GeneRatio, "/"),
          ~ as.numeric(.x[1]) / as.numeric(.x[2])
        ))
    })

    for (nm in names(go_results_df_all)) {
      df <- go_results_df_all[[nm]]
      if (nrow(df) == 0) next
      p <- df %>%
        dplyr::mutate(Description = reorder(Description, GeneRatio)) %>%
        dplyr::slice_max(GeneRatio, n = 15) %>%
        ggplot2::ggplot(ggplot2::aes(x = Description, y = GeneRatio, colour = p.adjust)) +
        ggplot2::geom_segment(ggplot2::aes(xend = Description, y = 0, yend = GeneRatio)) +
        ggplot2::geom_point(ggplot2::aes(size = Count)) +
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
  } else {
    warning("RData file not found: ", rdata_path, ". Please run 'run_go_enrichment()' first.")
  }

  invisible(go_results)
}
