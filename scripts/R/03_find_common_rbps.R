# To find common genes between the four dataframes and to create a new dataframe with those common genes
# Load necessary libraries
library(dplyr)

common_genes <- inner_join(rbpgo, rbpdb, rbpimage, rbpworld, by = "Gene_ID")
