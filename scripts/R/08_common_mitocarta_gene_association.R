# Script to find common overlapping genes between MitoCarta and found common genes
# Call necessary libraries

library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidydr)
library(tidyverse)
library(writexl)

mito_df <- read_excel("D:/rbp-project/data/raw/mitocarta_dataset/human_mitocarta_dataset.xls", sheet = "A Human MitoCarta3.0") #Imported Mitochondrial Protein Data set
com_df <- read_csv("data/processed/common_rbps_with_prot_name.csv") #Imported common RBP genes

# path to the mitocarta datasets D:\rbp-project\data\raw\mitocarta_dataset\Human.MitoCarta3.0
# path to the common rbp dataset D:\rbp-project\data\processed\common_rbps_with_prot_name

# Data Cleaning MitoCarta dataset
mito_df <- rename(mito_df, gene_name = "Symbol")

# Clean IDs: trim, uppercase, drop blanks/NA, de-duplicate ---------------
clean_ids <- function(df) {
    df |>
      mutate(
        gene_name = str_trim(gene_name),
        gene_name = toupper(gene_name)
      ) |>
      filter(!is.na(gene_name), gene_name != "") |>
      distinct(gene_name,description, .keep_all = TRUE)
  }

mito_df <- clean_ids(mito_df)

any(is.na(mito_df$gene_name))

mito_df <- na.omit(mito_df)


# Check if there are NA values in MitoCarta dataset
if(any(is.na(mito_df))) {
  print("There are NA values in the MitoCarta dataframe.")
} else {
  print("There are no NA values in the MitoCarta dataframe.")
}

# Getting common genes from MitoCarta Dataset and genes across databases 
common_genes_with_mito <- mito_df %>%
  inner_join(com_df)

# Get dataframe as Excel file and store in data/processed folder
write_xlsx(common_genes_with_mito, "data/processed/common_genes_with_mito.xlsx")
print("Written common_genes_with_mito dataframe to Excel file in data/processed folder")




