
source("./code/00_packages.R")
library(tidyr)
library(dplyr)
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
# load results_list_go_ALL RData file
load("./data/output/RData/SC/results_list_enrichGO_ALL.RData")


# Make dataframe list from GO results #
results_df_enrichGO_ALL <- map(results_list_enrichGO_ALL, ~ .x@result)
#Calculate GeneRatio from character fraction
results_df_enrichGO_ALL <- lapply(results_df_enrichGO_ALL, function(df) {
  df <- df %>%
    mutate(GeneRatio = strsplit(GeneRatio, "/") %>%
             map_dbl(~ as.numeric(.x[1]) / as.numeric(.x[2])))
  return(df)
})
# Extract the names of the dataframes without the ".up" or ".down" suffix
df_names <- sub("\\.UP$|\\.DOWN$", "", names(results_df_enrichGO_ALL))

# Function to add a column with the dataframe name to each dataframe
add_df_name <- function(df, name) {
  if (nrow(df) > 0) {
    df$cluster <- name
  }
  return(df)
}

# Use lapply to apply the function to each dataframe in the list
list_of_dfs <- mapply(add_df_name, results_df_enrichGO_ALL, df_names, SIMPLIFY = FALSE)
# Combine all dataframes into a single dataframe
combined_df <- do.call(rbind, list_of_dfs)

selected_clusters <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster11", "cluster14", "cluster18")
heat_df <- combined_df %>% 
  dplyr::filter(ONTOLOGY == "BP", cluster %in% selected_clusters ) %>% 
  dplyr::select(ID, Description,GeneRatio, p.adjust, Count, cluster)


heatmap_GO_function <- function(x){
# Aggregate the data to find the lowest p-value for each combination of GO term and cluster
agg_df <- aggregate(p.adjust ~ Description + cluster, data = heat_df, FUN = min)


# Reshape the data to create a matrix suitable for heatmap
heatmap_data <- reshape(agg_df, idvar = "Description", timevar = "cluster", direction = "wide")
heatmap_data <- heatmap_data[order(heatmap_data$Description),]
# Remove "p.adjust." prefix from column names
colnames(heatmap_data)[-1] <- gsub("^p.adjust\\.", "", colnames(heatmap_data)[-1])
row.names(heatmap_data) <- heatmap_data$Description



# Extract the data to normalize
heatmap_data_matrix <- as.matrix(heatmap_data[, -1])  # Exclude the GO_term column

# Z-score normalization
heatmap_data_scaled <- scale(heatmap_data_matrix)

# Close any open graphics devices
#dev.off()

# Create heatmap using heatmap.2 with Z-score normalized data
p <- heatmap.2(heatmap_data_scaled, 
          Rowv = FALSE, 
          Colv = FALSE, 
          dendrogram = "none", 
          trace = "none", 
          col = cm.colors(256),
          key = TRUE, 
          keysize = 1.0, 
          density.info = "none", 
          margins = c(5, 10))

print(p)
#tiff_filename <- paste0("data/output/SC/heatmap/", plot_name, ".tiff")
#pdf_filename <- paste0("data/output/SC/heatmap/", plot_name, ".pdf")
#ggsave(filename = tiff_filename, plot = p, device = "tiff")
#ggsave(filename = pdf_filename, plot = p, device = "pdf")


}
#All clusters
heat_df <- combined_df %>% 
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::select(ID, Description,GeneRatio, p.adjust, Count, cluster)

heatmap_GO_function(heat_df)


# Cluster only for excitatory neurons
selected_clusters <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster11", "cluster14", "cluster18")
heat_df <- combined_df %>% 
  dplyr::filter(ONTOLOGY == "BP", cluster %in% selected_clusters ) %>% 
  dplyr::select(ID, Description,GeneRatio, p.adjust, Count, cluster)

heatmap_GO_function(heat_df)


# Cluster only for excitatory neurons
selected_clusters <- c()
heat_df <- combined_df %>% 
  dplyr::filter(ONTOLOGY == "BP") %>% 
  dplyr::select(ID, Description,GeneRatio, p.adjust, Count, cluster) %>% 
  group_by(Description, cluster) %>%
  slice_min(p.adjust) %>% 
  pivot_wider(names_from = cluster, values_from = p.adjust)


# Cluster only for non-neuronal cells
# Filter and select rows with the lowest p-value for each cluster



# Reshape the filtered data to create a matrix suitable for heatmap
heatmap_data <- filtered_df %>%
  pivot_wider(names_from = cluster, values_from = p.adjust)



