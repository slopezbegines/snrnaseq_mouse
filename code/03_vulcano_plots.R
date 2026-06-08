# Volcano plot function for cluster-level DE results (Ctrl vs condition)
#
# Expected input columns in cluster_markers:
#   gene, avg_log2FC, p_val_adj, cluster

plot_volcano_clusters <- function(
  cluster_markers,
  fc = NULL, # FC threshold; NULL → uses global FC
  pval = NULL, # p-value threshold; NULL → uses global p_val
  top_n = 10, # genes to label per direction (UP / DOWN), ranked by p_val_adj
  colors = c("DOWN" = "#4DBBD5", "UP" = "#F39B7F", "NO" = "#A6A6A6"),
  width = 8,
  height = 6,
  subdir = "volcano"
) {
  fc <- if (is.null(fc)) get("FC", envir = .GlobalEnv) else fc
  pval <- if (is.null(pval)) get("p_val", envir = .GlobalEnv) else pval

  clusters <- sort(unique(cluster_markers$cluster))
  plots <- vector("list", length(clusters))
  names(plots) <- as.character(clusters)

  for (cl in clusters) {
    df <- cluster_markers[cluster_markers$cluster == cl, ]

    # -log10(0) → Inf; replace zeros with minimum non-zero value
    min_p <- min(df$p_val_adj[df$p_val_adj > 0], na.rm = TRUE)
    df$p_adj_plot <- ifelse(df$p_val_adj == 0 | is.na(df$p_val_adj), min_p, df$p_val_adj)

    df$direction <- dplyr::case_when(
      df$avg_log2FC > fc & df$p_val_adj < pval ~ "UP",
      df$avg_log2FC < -fc & df$p_val_adj < pval ~ "DOWN",
      TRUE ~ "NO"
    )
    # Ensure all three levels are present for consistent colour mapping
    df$direction <- factor(df$direction, levels = c("UP", "DOWN", "NO"))

    up_genes <- df[df$direction == "UP", ]
    down_genes <- df[df$direction == "DOWN", ]
    label_genes <- c(
      up_genes$gene[order(up_genes$p_adj_plot)[seq_len(min(top_n, nrow(up_genes)))]],
      down_genes$gene[order(down_genes$p_adj_plot)[seq_len(min(top_n, nrow(down_genes)))]]
    )
    df$label <- ifelse(df$gene %in% label_genes, df$gene, NA_character_)

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x     = avg_log2FC,
        y     = -log10(p_adj_plot),
        col   = direction,
        label = label
      )
    ) +
      ggplot2::geom_point(size = 1.5, alpha = 0.8) +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(paste("Volcano plot — Cluster", cl)) +
      ggplot2::labs(
        x = expression("avg log"[2] * "FC  (Ctrl vs condition)"),
        y = expression(-log[10] * "(p-value adj)")
      ) +
      ggplot2::geom_vline(xintercept = c(-fc, fc), col = "red", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(pval), col = "red", linetype = "dashed") +
      ggplot2::scale_colour_manual(values = colors, drop = FALSE) +
      ggrepel::geom_text_repel(na.rm = TRUE, max.overlaps = 15, size = 3) +
      ggplot2::guides(colour = ggplot2::guide_legend(title = NULL))

    save_plot(
      paste0("volcano_cluster_", cl), p,
      width = width,
      height = height,
      subdir = subdir
    )
    print(p) # Display the plot in the R session
    plots[[as.character(cl)]] <- p
  }

  invisible(plots)
}
