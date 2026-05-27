# STRING PPI network analysis for cluster-level DE results
#
# Expected input columns in cluster_markers:
#   gene, avg_log2FC, p_val_adj, cluster

run_string_network <- function(
  cluster_markers,
  fc = NULL, # FC threshold; NULL → global FC
  pval = NULL, # p-value threshold; NULL → global p_val
  species = 10090, # 10090 = mouse, 9606 = human
  score_threshold = 200,
  version = "12.0",
  subdir = "string"
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

  id_list <- lapply(split_list, function(df) data.frame(gene = df$gene))

  # --- Connect to STRING -------------------------------------------------------
  string_db <- STRINGdb$new(
    version         = version,
    species         = species,
    score_threshold = score_threshold,
    input_directory = ""
  )

  get_link <- function(x) {
    mapped <- string_db$map(x, "gene", removeUnmappedRows = TRUE)
    if (nrow(mapped) == 0) {
      return(NA_character_)
    }
    hits <- mapped$STRING_id[seq_len(nrow(mapped))]
    as.character(string_db$get_link(hits))
  }

  # --- Generate links and save -------------------------------------------------
  STRING_links <- lapply(names(id_list), function(nm) {
    link <- tryCatch(
      get_link(id_list[[nm]]),
      error = function(e) {
        message("STRING failed for ", nm, ": ", e$message)
        NA_character_
      }
    )
    data.frame(Name = nm, Link = link, stringsAsFactors = FALSE)
  })

  STRING_df <- do.call(rbind, STRING_links)
  writexl::write_xlsx(
    STRING_df,
    file.path(output_path, "tables", "STRING_plotnames_table.xlsx")
  )

  message("[STRING] Links saved to tables/STRING_plotnames_table.xlsx")
  print(knitr::kable(STRING_df, format = "simple"))

  invisible(STRING_df)
}
