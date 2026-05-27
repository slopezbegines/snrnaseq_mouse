
# libraries ####

source("./code/00_packages.R")
# Global variables
tiff_extension <- ".tiff"
pdf_extension <- ".pdf" #Vectorial format

# Setting thresholds ####
# Setting threshold for p-value and Fold-Change
p_val <- 0.05
FC <- 0.1

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
input_data <- excel_sheet_reader(paste0(output_path,"tables/","cluster_list.xlsx"))


enrichr_libs <- rba_enrichr_libs()


# Get list of annotation dataset. Choose one for annot_dataset option into rba_panther_enrich
annots <- rba_panther_info(what = "datasets")

#species: 9606 for Human, 10090 for mouse, 7955 for zebrafish
perform_enrichR <- function(gene_set) {
  enriched <- rba_enrichr(gene_list = gene_set$ID)
  return(enriched)
}
# Apply the enrichGO function to each gene set in differential_genes
results_enrichR <- lapply(input_data, perform_enrichR)

# Export results ####
save(results_enrichR, file = paste0(output_path,"RData/results_enrichR.RData")) #File for Strings DB

# Load data ####
#Load dep_analysis object
#load("data/output/RData/SC/results_enrichR.RData")


test_mouse_genes <- input_data$cluster0$ID

mouse_enrichr_output <- enrichR::enrichr(test_mouse_genes, databases = "Mouse_Gene_Atlas")

mouse_enrichr_output %>% as.data.frame() %>%View

results_enrichR$cluster0$Mouse_Gene_Atlas %>%View
