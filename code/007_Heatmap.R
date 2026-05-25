

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
      #dplyr::select(ID, everything())  # Reorder the columns with "ID" as the first column
      dplyr::select(ID, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)  # Reorder the columns with "ID" as the first column
        return(df)
  })
  # Combine the list of data frames into a single named list
  result <- setNames(sheet_data, sheets)
  # Return the list
  return(result)
}
# Call the function with the filename as the argument
#input_data <- excel_sheet_reader("./data/input/res0.5_KO_1_vs_WT_1_response_Minpct0.1_LFC0.1.xlsx")
input_data <- excel_sheet_reader("./data/input/final_raw_DEGs_PostAnnotation_KOs_vs_WTs.xlsx")


# Gene list
gene_list <- c("Hmgcr", "Srebf2", "Apoe","Lrp1", "Lrp1b", "Abca1")


# Initialize an empty dataframe to store the final result
final_df <- data.frame(ID = factor(), avg_log2FC = numeric(), p_val_adj = numeric())


process_df <- function(df, name) {
  df <- df %>%
    dplyr::filter(ID %in% gene_list) %>%
    dplyr::select(ID, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
    dplyr::mutate(cluster_index = name)  # Add cluster index as a column
  return(df)
}

sheet_names <- names(input_data)

# Check if input_data is a named list
if (!is.null(names(input_data))) {
  # Apply the function to each dataframe in the list and bind them together
  final_df <- bind_rows(lapply(names(input_data), function(name) process_df(input_data[[name]], name)))
  
  # Print the final dataframe
  print(final_df)
} else {
  print("The input_data list does not seem to be properly named.")
}

final_df$ID <- as.factor(final_df$ID)
final_df$cluster_index <- as.factor(final_df$cluster_index)

'
# Check if input_data is a named list
if(!is.null(names(input_data))) {
  # Apply the function to each dataframe in the list and bind them together
  #final_df <- bind_rows(lapply(0:20, function(i) process_df(input_data[[paste0("cluster", i)]], i)))
 
  # Remove rows where avg_log2FC is 0
 # final_df <- final_df %>%
  #  filter(p != 0)
  
  # Reorder columns as required
  final_df <- final_df[, c("ID", "avg", "p", "cluster_index")] %>% 
    #filter_("p" >0)
  #final_df <- final_df[, c("ID", "p", "cluster_index")]
  # Print the final dataframe
  print(final_df)
} else {
  print("The input_data list does not seem to be properly named.")
}

'
# Convert final_df to a dataframe
final_df <- as.data.frame(final_df)

# Obtener los nombres de los genes presentes en final_df
genes_presentes <- unique(final_df$ID)

subset_df <- final_df %>% 
  dplyr::select(ID, avg_log2FC, cluster_index) %>% 
  #dplyr::mutate(logp = -log10(p)) %>% 
  na.omit()
  

# make a df with the gene names and avg_log2FC values for clusters

heatmap_data <- data.frame(Gene = genes_presentes)  # Crear columna de genes
for (name in unique(subset_df$cluster_index)) {
  cluster_values <- subset_df %>% 
    dplyr::filter(cluster_index == name) %>% 
    dplyr::select(ID, avg_log2FC) %>% 
    setNames(c("Gene", name))  # Renombrar la columna avg_log2FC
    heatmap_data <- merge(heatmap_data, cluster_values, by = "Gene", all.x = TRUE)  # Combinar datos
}
# Substitute NA for 0
heatmap_data <- replace(heatmap_data, is.na(heatmap_data), 0)

# Remove the "Gene" column as it's not needed for the heatmap plot
heatmap_data_no_gene <- heatmap_data[, -1]

# Convert the data to a matrix for plotting
heatmap_matrix <- as.matrix(heatmap_data_no_gene)

rownames(heatmap_matrix) <- heatmap_data$Gene

specific_order <- c("Cluster_0","Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4", "Cluster_5", 
                    "Cluster_6", "Cluster_7", "Cluster_8", "Cluster_9", "Cluster_10", 
                    "Cluster_11", "Cluster_12", "Cluster_13", "Cluster_14", "Cluster_15", 
                    "Cluster_16", "Cluster_17", "Cluster_18", "Cluster_19", "Cluster_20")  # Add the names in the desired order

specific_order <- c("ExN1",  "ExN2",  "ExN3",  "ExN4",  "ExN5",
                    "InN1",  "InN2",  "InN3",  "InN4",
                    "Oligo", "OPCs",  "Micro", "Astro", "Vasc")


heatmap <- pheatmap(heatmap_matrix, cutree_rows = 6, row_names = TRUE,
         cluster_cols = FALSE, scale = "none", 
         labels_col = specific_order,  # Add cluster names to x-axis
         labels_row = heatmap_data$Gene,  # Add gene names to y-axis
         cluster_rows = TRUE,
         legend_title = "Log2 FC")  # Add legend title


plot_name <- "1_heatmap_6genes"
tiff_filename <- paste0("data/output/SC/heatmap/0", plot_name, ".tiff")
pdf_filename <- paste0("data/output/SC/heatmap/0", plot_name, ".pdf")
ggsave(filename = tiff_filename, plot = heatmap, device = "tiff")
ggsave(filename = pdf_filename, plot = heatmap, device = "pdf")

# Define custom color palette
my_color_palette <- colorRampPalette(viridis(12))(100)#define your color scale

# Plot heatmap with specific order of cluster names on x-axis
pheatmap(heatmap_matrix, cutree_rows = 6, row_names = TRUE, 
         cluster_cols = FALSE, scale = "none",
         #annotation_col = heatmap_data$Gene_Present,
         #annotation_colors = list(Gene_Present = c("No" = "white", "Yes" = "black")),
         #annotation_legend = TRUE,
         labels_col = specific_order,  # Add cluster names to x-axis
         labels_row = heatmap_data$Gene,  # Add gene names to y-axis
         #annotation_names_col = FALSE,
         cluster_rows = TRUE,
         #legend_title = "Log2 FC",  # Add legend title
         color = my_color_palette)


# Convert cluster_index column to a factor with specific order
final_df$cluster_index <- factor(final_df$cluster_index, levels = rev(specific_order))

bar_plot_genes <- final_df %>% 
  filter(p_val_adj < 0.05) %>% 
  ggplot() +
  aes(x = cluster_index, y = avg_log2FC, fill = cluster_index) +
  geom_col() +
  #scale_fill_viridis_d(option = "viridis", direction = -1) +
  coord_flip() +
  theme_minimal() +
  facet_wrap(vars(ID), scales = "free_y", nrow = 2L)
  #facet_wrap(vars(ID), nrow = 2L)
bar_plot_genes
plot_name <- "bar_plot_genes_pval_filter_default"
tiff_filename <- paste0("data/output/SC/heatmap/0", plot_name, ".tiff")
pdf_filename <- paste0("data/output/SC/heatmap/0", plot_name, ".pdf")
ggsave(filename = tiff_filename, plot = bar_plot_genes, device = "tiff")
ggsave(filename = pdf_filename, plot = bar_plot_genes, device = "pdf")




final_df %>%
  filter(p_val_adj >= 0 & p_val_adj < 0.05) %>%
  ggplot() +
  aes(x = cluster_index, y = avg_log2FC, fill = cluster_index) +
  geom_col(width = 0.3) +
  scale_fill_brewer(palette = "RdYlBu", direction = 1) +
  labs(x = "Cluster", y = "Fold-change", fill = "Cluster") +
  coord_flip() +
  theme_minimal() +
  facet_wrap(
    vars(ID),
    scales = "free_y",
    ncol = 1L,
    nrow = 6L
  )



bar_plot_genes <- final_df %>%
  filter(p_val_adj >= 0 & p_val_adj < 0.05) %>%
  ggplot() +
  aes(x = cluster_index, y = avg_log2FC, fill = cluster_index) +
  geom_col(width =bar_width) + # Adjust the width here
  scale_fill_hue(direction = 1) +
  labs(x = "Cluster", y = "Fold-change", fill = "Cluster") +
  coord_flip() +
  theme_minimal() +
  facet_wrap(vars(ID), 
             scales = "free_y",
             ncol = 1L, 
             nrow = 6L, 
             shrink = TRUE)

bar_plot_genes
plot_name <- "oneplot_bar_plot_genes_pval_filter_001"
tiff_filename <- paste0("data/output/SC/heatmap/", plot_name, ".tiff")
pdf_filename <- paste0("data/output/SC/heatmap/", plot_name, ".pdf")
ggsave(filename = tiff_filename, plot = bar_plot_genes, device = "tiff", width = 4, height = 6, units = "in")
ggsave(filename = pdf_filename, plot = bar_plot_genes, device = "pdf", width = 4, height = 6, units = "in")

