

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

id_list <- lapply(split_filtered_list, function(df) data.frame(gene = df$gene))

# Strings ####
# load string database
#9606 for Human, 10090 for mouse, 7955 for zebrafish
string_db <- STRINGdb$new(version = "11.5", species = 10090, score_threshold = 200, input_directory="")
class(string_db)


string_function <- function(x) {
  string_example <- string_db$map(x, "gene", removeUnmappedRows = TRUE)
  dimension <- dim(string_example)[1]
  hits <- string_example$STRING_id[1:dimension]
  link <- as.character(string_db$get_link(hits))
  return(link)
}

STRING_plotnames <- lapply(names(id_list), function(name) {
  link <- string_function(id_list[[name]])
  data.frame(Name = name, Link = link, stringsAsFactors = FALSE)
})


### Save String plot links ####
STRING_plotnames_df <- do.call(rbind, STRING_plotnames)
STRING_plotnames_table <- STRING_plotnames_df %>%
  kable("html") %>%
  kable_styling()

print(STRING_plotnames_table)

STRING_plotnames_df %>% 
  write_xlsx(paste0(output_path, "tables/","STRING_plotnames_table.xlsx"))

