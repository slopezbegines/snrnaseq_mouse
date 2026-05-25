

# libraries ####

#source("./code/00_packages.R")


#Filtering results
# Function to filter a single data frame
filter_dataframe <- function(df) {
  df %>%
    filter(abs(avg_log2FC) > FC, p_val_adj < 0.05)  %>%
    mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"))
}
filtered_list <- map(cluster_list, filter_dataframe)

# Function to split the filtered data frame into two based on avg_log2FC values
split_dataframe <- function(filtered_df) {
  up_df <- filtered_df %>%
    filter(avg_log2FC > 0) %>%
    mutate(direction = "UP")
  
  down_df <- filtered_df %>%
    filter(avg_log2FC < 0) %>%
    mutate(direction = "DOWN")
  
  return(list(UP = up_df, DOWN = down_df))
}
# Apply the splitting function to each filtered data frame
split_filtered_list <- map(filtered_list, split_dataframe)
split_filtered_list <- unlist(split_filtered_list, recursive = FALSE)# Remove one level in the list



# Gene Ontology ####
# GO functions ####
### EnrichGO functions ####
# GO terms list
# Function to enrichGO for the dataframes of a list
perform_enrichGO <- function(gene_set) {
  results <- clusterProfiler::enrichGO(gene = gene_set$gene, 
                                       OrgDb = "org.Mm.eg.db", 
                                       keyType = "SYMBOL", 
                                       ont = "ALL",
                                       pAdjustMethod = "none")
  return(results)
}

# Apply the enrichGO function to each gene set in differential_genes
#results_list_enrichGO_ALL <- lapply(filtered_list, perform_enrichGO)
results_list_enrichGO_ALL <- lapply(split_filtered_list, perform_enrichGO)
#results_list_enrichGO_ALL_names <- names(results_list_enrichGO_ALL)
# Eliminar elementos NULL de la lista
results_list_enrichGO_ALL <- Filter(Negate(is.null), results_list_enrichGO_ALL)

save(results_list_enrichGO_ALL, file = paste0(output_path, "RData/","results_list_enrichGO_ALL.RData"))

# load results_list_go_ALL RData file
load(paste0(output_path, "RData/","results_list_enrichGO_ALL.RData"))





# GO Barplots function
perform_barplots <- function(x) {
  if (length(x[1][10]) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    results <- barplot(x, showCategory = 30) + facet_grid(ONTOLOGY ~ ., scale = "free")
    #results <- enrichplot::dotplot(x, showCategory = 30)
    return(results)
  }
}

# GO Dotplots function
perform_dotplots <- function(x) {
  if (length(x[1][10]) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    results <- enrichplot::dotplot(x, showCategory = 30) + facet_grid(ONTOLOGY ~., scale = "free")
    #results <- enrichplot::dotplot(x, showCategory = 30)
    return(results)
  }
}

# GO cnetplots function
perform_cnetplots <- function(x){
  if (length(x@result$Count) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(0)
  } else {
    ## remove redundent GO terms
    x2 <- simplify(x)
    #results <- enrichplot::cnetplot(x2)
    #For circular Gene Concept map
    results <- enrichplot::cnetplot(x2, circular = TRUE, colorEdge = TRUE)
    return(results)
  }
}


#UpSet plot function
perform_upsetplot <- function(x){
  if (length(x@result$Count) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(0)
  } else {
    results <- enrichplot::upsetplot(x)
    return(results)
  }
}
#image_number <-129

#HeatMap plot function
perform_heatmapplot <- function(x){
  if (length(x@result$Count) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(0)
  } else {
    results <- enrichplot::heatplot(x)
    return(results)
  }
}

#image_number <- 000

### Barplot ####
barplot_results <- lapply(results_list_enrichGO_ALL, perform_barplots)
# Iterate over dotplot_results and print each plot individually
for (i in seq_along(barplot_results)) {
  if (nrow(barplot_results[[i]][[1]][10]) == 0) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_enrichGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.ID", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- barplot_results[[i]] +
      ggplot2::ggtitle(paste0("Bar plot for ", plot_name))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/","0",image_number +i,"_barplot_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Tiff format
  }
}


image_number <- image_number +i

### Dotplot ####

dotplot_results <- lapply(results_list_enrichGO_ALL, perform_dotplots)
# Iterate over dotplot_results and print each plot individually
for (i in seq_along(dotplot_results)) {
  if (nrow(dotplot_results[[i]][[1]][10]) == 0) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_enrichGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.ID", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- dotplot_results[[i]] +
      ggplot2::ggtitle(paste0("Dot plot for ", plot_name))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/","0",image_number +i,"_Dotplot_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Tiff format
  }
}
image_number <- image_number +i

### CNET plot ####

cnetplot_results <- lapply(results_list_enrichGO_ALL, perform_cnetplots)

for (i in seq_along(cnetplot_results)) {
  if (is.list(cnetplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_enrichGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.ID", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- cnetplot_results[[i]]+
      ggplot2::ggtitle(paste0("CNET plot for ", plot_name))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/","0",image_number+i,"_cnet_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Tiff format
  }
}
image_number <- image_number+i



### UPSET plot ####
upsetplot_results <- lapply(results_list_enrichGO_ALL, perform_upsetplot)

for (i in seq_along(upsetplot_results)) {
  if (is.list(upsetplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_enrichGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.ID", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- upsetplot_results[[i]]+
      ggplot2::ggtitle(paste0("UPSET plot for ", plot_name))+
      geom_bar(aes(stat = "identity", width=0.2)) 
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/","0",image_number+i,"_upset_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Tiff format
  }
}
image_number <- image_number+i


### Heatmap plot ####
heatmapplot_results <- lapply(results_list_enrichGO_ALL, perform_heatmapplot)

for (i in seq_along(heatmapplot_results)) {
  if (is.list(heatmapplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_enrichGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- heatmapplot_results[[i]]+
      ggplot2::ggtitle(paste0("HeatMap plot for ", plot_name))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/","0",image_number+i,"_heatmap_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in", dpi = 300)  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in", dpi = 300)  # Tiff format
  }
}
image_number <- image_number+i

#Look at this page follow it.
#https://bioc.ism.ac.jp/packages/3.7/bioc/vignettes/enrichplot/inst/doc/enrichplot.html#bar-plot

## Lolliplots ####

# Make dataframe list from GO results #
results_df_enrichGO_ALL <- map(results_list_enrichGO_ALL, ~ .x@result)

results_df_enrichGO_ALL <- lapply(results_df_enrichGO_ALL, function(df) {
  df <- df %>%
    mutate(GeneRatio = strsplit(GeneRatio, "/") %>%
             map_dbl(~ as.numeric(.x[1]) / as.numeric(.x[2])))
  return(df)
})

#Save xlsx
results_df_enrichGO_ALL %>% 
  write_xlsx(paste0(output_path,"tables/","results_df_enrichGO_ALL.xlsx"))

lolliplot <- function(data_name, df_list, file_prefix = NULL) {
  non_empty_indices <- which(sapply(df_list, function(x) nrow(x) > 0))
  
  if (length(non_empty_indices) == 0) {
    message("All data frames are empty. No plots generated.")
    return(NULL)
  }
  
  if (length(non_empty_indices) < data_name) {
    message("Selected index is out of range. No plot generated.")
    return(NULL)
  }
  
  df <- df_list[[non_empty_indices[data_name]]]
  
  if (nrow(df) == 0) {
    message("Selected data frame is empty. No plot generated.")
    return(NULL)
  }
  
  plot_name <- gsub("^names_", "", names(df_list)[non_empty_indices[data_name]])  # Remove "names_"
  plot_name <- gsub("_vs_", " vs ", plot_name)  # Replace "_vs_" with " vs "
  plot_name <- gsub("_", " ", plot_name)  # Replace "_" with " "
  plot_name <- gsub("UP$", "UP", plot_name)  # Remove "UP" from the end
  
  plot <- df %>%
    dplyr::mutate(Description = reorder(Description, GeneRatio)) %>%
    top_n(15, GeneRatio) %>%
    ggplot(aes(x = Description,
               y = GeneRatio,
               colour = p.adjust)) +
    geom_segment(aes(x = Description,
                     xend = Description,
                     y = 0,
                     yend = GeneRatio)) +
    geom_point(aes(size = Count), show.legend = TRUE) +
    scale_color_viridis_c(option = "viridis", direction = 1) +
    facet_wrap(ONTOLOGY ~ ., scale = "free", ncol=1)+
    coord_flip() +
    theme_gray() +
    labs(size = "N. of genes",
         x = "GO term",
         y = "Gene Ratio") + 
    ggtitle(paste("Plot for", plot_name))
  
  if (!is.null(file_prefix)) {
    tiff_filename <- paste0(file_prefix, "_", gsub(" ", "_", plot_name), ".tiff")
    pdf_filename <- paste0(file_prefix, "_", gsub(" ", "_", plot_name), ".pdf")
    ggsave(filename = tiff_filename, plot = plot, device = "tiff")
    ggsave(filename = pdf_filename, plot = plot, device = "pdf")
  }
  print(plot)
  return(plot) 
}

# Apply lolliplot function
for (i in seq_along(results_df_enrichGO_ALL)) {
  lolliplot_result <- lolliplot(i, results_df_enrichGO_ALL, file_prefix = paste0(output_path, "figures/enrichGO/","Lolliplot_0",i))
}

