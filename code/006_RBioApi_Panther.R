



# libraries ####

source("./code/00_packages.R")
# Setting thresholds ####
# Setting threshold for p-value and Fold-Change
p_val <- 0.05
FC <- 0.5

# Global variables
tiff_extension <- ".tiff"
pdf_extension <- ".pdf" #Vectorial format


# Load data ####

#Load data results tables
excel_sheet_reader <- function(filename) {
  # Get the names of all the sheets in the Excel file
  sheets <- readxl::excel_sheets(filename)
  # Read in each sheet, rename the first column to "ID", and store it in a list
  sheet_data <- lapply(sheets, function(sheet_name) {
    df <- readxl::read_excel(filename, sheet = sheet_name) %>%
      dplyr::rename(ID = ...1) %>%  # Rename the first column to "ID"
      dplyr::select(ID, everything())  # Reorder the columns with "ID" as the first column
    return(df)
  })
  # Combine the list of data frames into a single named list
  result <- setNames(sheet_data, sheets)
  # Return the list
  return(result)
}
# Call the function with the filename as the argument
input_data <- excel_sheet_reader("./data/input/res0.5_KO_1_vs_WT_1_response_Minpct0.1_LFC0.1.xlsx")

#Filtering results
# Function to filter a single data frame
filter_dataframe <- function(df) {
  df %>%
    filter(abs(avg_log2FC) > FC, p_val_adj < 0.05)  %>%
    mutate(direction = ifelse(avg_log2FC > 0, "UP", "DOWN"))
}
filtered_list <- map(input_data, filter_dataframe)

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



# Enrichment ####
#Statistical overrepresentation test

#species: 9606 for Human, 10090 for mouse, 7955 for zebrafish
perform_enrich <- function(gene_set) {
  enriched <- rba_panther_enrich(genes = gene_set$ID,
                                 organism = 10090,
                                 annot_dataset = "GO:0008150", #from annots result
                                 cutoff = 0.05)
  
  return(enriched)
}
# Apply the enrichGO function to each gene set in differential_genes
results_panther <- lapply(split_filtered_list, perform_enrich)

#Save results
save(results_panther, file = "data/output/RData/SC/results_panther.RData")

# load results_list_go_ALL RData file
load("data/output/RData/SC/results_panther.RData")

# Make dataframe list from Panther results #
GO_panther <- map(results_panther, ~ .x$result)

# Plots ####

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
    dplyr::mutate(term.label = reorder(term.label, fold_enrichment)) %>%
    top_n(15, fold_enrichment) %>%
    ggplot(aes(x = term.label,
               y = fold_enrichment,
               colour = -log10(fdr))) +
    geom_segment(aes(x = term.label,
                     xend = term.label,
                     y = 0,
                     yend = fold_enrichment)) +
    geom_point(aes(size = number_in_list), show.legend = TRUE) +
    scale_color_viridis_c(option = "viridis", direction = 1) +
    facet_wrap(ONTOLOGY ~ ., scale = "free", ncol=1)+
    coord_flip() +
    theme_gray() +
    labs(size = "N. of genes",
         x = "GO term",
         y = "Fold enrichment") + 
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
for (i in seq_along(GO_panther)) {
  lolliplot_result <- lolliplot(i, GO_panther, file_prefix = paste0("data/output/SC/panther/split/Lolliplot_0",i))
}
