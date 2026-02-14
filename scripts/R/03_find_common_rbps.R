# To find common genes between the four dataframes and to create a new dataframe with those common genes
# Load necessary libraries
library(dplyr)
library(readr)
library(readxl)
library(writexl)

# Read Excel files from data/processed folder
rbpgo <- read_excel("data/processed/rbpgo.xlsx")
rbpimage <- read_excel("data/processed/rbpimage.xlsx")
rbpworld <- read_excel("data/processed/rbpworld.xlsx")

# Find common genes 
common_rbps <- rbpgo %>%
  inner_join(rbpimage) %>%
  inner_join(rbpworld)

# Get rbpdb dataframe as Excel file and store in data/processed folder
write_xlsx(common_rbps, "data/processed/common_rbps.xlsx")
print("Written common_rbps dataframe to Excel file in data/processed folder")
