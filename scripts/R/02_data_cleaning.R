# Clean and prepare data for analysis
# Load necessary libraries
library(dplyr)
library(readr)

#1. Data cleaning rbpdb dataframe  
# First to drop columns which are not required 

drop <- c("1226","2010-04-24","2010-04-24","2010-03-18","9606")

rbpdb <- rbpdb[,!(names(rbpdb) %in% drop)]
print('Modified dataframe:-')
rbpdb

# Droping columns from rbpdb dataframe using indexing
# Remove the last 5 columns
n_to_remove <- 5
rbpdb <- rbpdb[, 1:(ncol(rbpdb) - n_to_remove)]

# Add column names to rbpdb dataframe
new_names <- c("ensembl_id", "gene_name", "function", "species")
colnames(rbpdb) <- new_names





