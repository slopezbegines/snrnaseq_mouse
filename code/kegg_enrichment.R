# KEGG gene set enrichment analysis for cluster-level DE results
#
# Expected input columns in cluster_markers:
#   gene, avg_log2FC, p_val_adj, cluster

run_kegg_enrichment <- function(
  cluster_markers,
  fc = NULL, # FC threshold; NULL → global FC
  pval = NULL, # p-value threshold; NULL → global p_val
  org_db = "org.Mm.eg.db",
  kegg_org = "mmu", # KEGG organism code (mmu = mouse, hsa = human)
  n_perm = 10000,
  width = 8,
  height = 6,
  subdir = "KEGG"
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

  # --- Convert SYMBOL → ENTREZ -------------------------------------------------
  convert_to_entrez <- function(gene_set) {
    if (length(gene_set) == 0) return(NULL)   # bitr misbehaves on character(0)
    mapped <- clusterProfiler::bitr(names(gene_set),
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org_db
    )
    if (nrow(mapped) == 0) {
      warning("No ENTREZ IDs mapped. Returning NULL.")
      return(NULL)
    }
    mapped    <- mapped[!duplicated(mapped$SYMBOL), ]   # one ENTREZ per SYMBOL
    gl_entrez <- gene_set[mapped$SYMBOL]
    names(gl_entrez) <- mapped$ENTREZID
    sort(na.omit(gl_entrez), decreasing = TRUE)
  }

  gene_lists_entrez <- lapply(gene_lists, function(gl) {
    tryCatch(convert_to_entrez(gl),
      error = function(e) {
        message("ENTREZ conversion failed: ", e$message)
        NULL
      }
    )
  })
  gene_lists_entrez <- Filter(Negate(is.null), gene_lists_entrez)

  # --- Run gseKEGG -------------------------------------------------------------
  results <- lapply(gene_lists_entrez, function(gl) {
    if (length(gl) == 0) {
      return(NULL)
    }
    tryCatch(
      clusterProfiler::gseKEGG(
        geneList          = gl,
        organism          = kegg_org,
        keyType           = "ncbi-geneid",
        # nPerm             = n_perm,
        minGSSize         = 3,
        maxGSSize         = 800,
        pvalueCutoff      = 0.05,
        pAdjustMethod     = "none",
        use_internal_data = FALSE
      ),
      error = function(e) {
        message("gseKEGG failed: ", e$message)
        NULL
      }
    )
  })
  results <- Filter(Negate(is.null), results)

  save(results, file = file.path(output_path, "RData", "results_list_gseKEGG_ALL.RData"))

  # Save result tables to xlsx
  results_df <- lapply(results, function(x) x@result)
  results_df <- Filter(function(x) nrow(x) > 0, results_df)
  writexl::write_xlsx(
    results_df,
    file.path(output_path, "tables", "results_df_gseKEGG_ALL.xlsx")
  )

  # --- Helper: clean plot label ------------------------------------------------
  .clean_name <- function(nm) gsub("_", " ", gsub("\\.[A-Z]+$", "", nm))

  # --- Lolliplots --------------------------------------------------------------
  for (nm in names(results_df)) {
    df <- results_df[[nm]]
    if (nrow(df) == 0) next
    p <- df %>%
      dplyr::mutate(Description = reorder(Description, enrichmentScore)) %>%
      dplyr::slice_max(enrichmentScore, n = 15) %>%
      ggplot2::ggplot(ggplot2::aes(x = Description, y = enrichmentScore, colour = p.adjust)) +
      ggplot2::geom_segment(ggplot2::aes(xend = Description, y = 0, yend = enrichmentScore)) +
      ggplot2::geom_point(ggplot2::aes(size = setSize)) +
      ggplot2::scale_color_viridis_c(option = "viridis") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(size = "N. of genes", x = "KEGG pathway", y = "Enrichment Score") +
      ggplot2::ggtitle(paste("KEGG Lolliplot —", .clean_name(nm)))
    print(p)
    save_plot(paste0("gseKEGG_lolliplot_", nm), p,
      width = width, height = height, subdir = subdir
    )
  }

  invisible(results)
}
