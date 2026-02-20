# Clean and prepare data for analysis
# Load necessary libraries
library(dplyr)
library(readr)
library(writexl)
library(tidyr)


#1. Data cleaning rbpgo dataframe
# Dropping columns from rbpgo dataframe
drop <- c("Entry_Name","RBP2GO_Score","Protein_Name","Nb_Datasets","Listing_Counts","AVG10_Int_Listing_Count","Mass_kDa","Length_AA","pI","Listing_Count")
rbpgo <- rbpgo[,!(names(rbpgo) %in% drop)]
print('Modified dataframe:-')
rbpgo

# Keep specific columns by index
rbpgo <- rbpgo[, c(1, 2, 3)]
print('Modified dataframe:-')
rbpgo

# Add column names to rbpdb dataframe
new_names <- c("uniprot_id", "gene_name", "alias_names")
colnames(rbpgo) <- new_names

# Remove duplicates from rbpgo dataframe
rbpgo <- rbpgo[!duplicated(rbpgo), ]

# Check if there are any dupicates in rbpgo dataframe
if(any(duplicated(rbpgo))) {
  print("There are duplicate rows in the rbpgo dataframe.")
} else {
  print("There are no duplicate rows in the rbpgo dataframe.")
}

# Drop NA values from rbpgo dataframe
rbpgo <- na.omit(rbpgo)

# Check if there are any NA values in the rbpgo dataframe
if(any(is.na(rbpgo))) {
  print("There are NA values in the rbpgo dataframe.")
} else {
  print("There are no NA values in the rbpgo dataframe.")
}

# Get rbpdb dataframe as Excel file and store in data/processed folder
write_xlsx(rbpgo, "data/processed/rbpgo.xlsx")
print("Written rbpgo dataframe to Excel file in data/processed folder")


#2. Data cleaning rbpimage dataframe
# Keep specific columns by index
rbpimage <- rbpimage[, c(1, 3)]
print('Modified dataframe:-')
rbpimage

# Add column names to rbpimage dataframe
new_names <- c("gene_name", "ensembl_id")
colnames(rbpimage) <- new_names

# Remove duplicates from rbpimage dataframe
rbpimage <- rbpimage[!duplicated(rbpimage), ]

# Check if there are any dupicates in rbpgo dataframe
if(any(duplicated(rbpimage))) {
  print("There are duplicate rows in the rbpimage dataframe.")
} else {
  print("There are no duplicate rows in the rbpimage dataframe.")
}

# Drop NA values from rbpgo dataframe
rbpimage <- na.omit(rbpimage)

# Check if there are any NA values in the rbpgo dataframe
if(any(is.na(rbpimage))) {
  print("There are NA values in the rbpimage dataframe.")
} else {
  print("There are no NA values in the rbpimage dataframe.")
}

# Get rbpimage dataframe as Excel file and store in data/processed folder
write_xlsx(rbpimage, "data/processed/rbpimage.xlsx")
print("Written rbpimage dataframe to Excel file in data/processed folder")


#3. Data cleaning rbpworld dataframe 
# Dropping columns from rbpworld dataframe
drop <- c("No. RBPome")
rbpworld <- rbpworld[,!(names(rbpworld) %in% drop)]
print('Modified dataframe:-')
rbpworld

# Add column names to rbpimage dataframe
new_names <- c("ensembl_id", "gene_name", "rbp_type")
colnames(rbpworld) <- new_names

# Remove duplicates from rbpimage dataframe
rbpimage <- rbpimage[!duplicated(rbpimage), ]

# Check if there are any dupicates in rbpgo dataframe
if(any(duplicated(rbpimage))) {
  print("There are duplicate rows in the rbpimage dataframe.")
} else {
  print("There are no duplicate rows in the rbpimage dataframe.")
}

# Drop NA values from rbpgo dataframe
rbpimage <- na.omit(rbpimage)

# Check if there are any NA values in the rbpgo dataframe
if(any(is.na(rbpimage))) {
  print("There are NA values in the rbpimage dataframe.")
} else {
  print("There are no NA values in the rbpimage dataframe.")
}


# Get rbpworld dataframe as Excel file and store in data/processed folder
write_xlsx(rbpworld, "data/processed/rbpworld.xlsx")
print("Written rbpworld dataframe to Excel file in data/processed folder")


#4. Data cleaning uniprot dataframe
# Dropping columns from uniprot dataframe
uniprot <- uniprot[, c(4, 5, 7)] 
print('Modified dataframe:-')
uniprot

# Add column names to rbpimage dataframe
new_names <- c("protein_name", "gene_name", "aa_length")
colnames(uniprot) <- new_names

# Remove duplicates from uniprot dataframe
uniprot <- uniprot[!duplicated(uniprot), ]

# Check if there are any duplicate rows in the uniprot dataframe
if(any(duplicated(uniprot))) {
  print("There are duplicate rows in the uniprot dataframe.")
} else {
  print("There are no duplicate rows in the uniprot dataframe.")
}

# Drop NA values from uniprot dataframe
uniprot <- uniprot %>% 
  drop_na()      # removes all rows that have at least one NA in any column

# Check if there are any NA values in the rbpdb dataframe
if(any(is.na(uniprot))) {
  print("There are NA values in the uniprot dataframe.")
} else {
  print("There are no NA values in the uniprot dataframe.")
}

# Get uniprot dataframe as Excel file and store in data/processed folder
write_xlsx(uniprot, "data/processed/uniprot.xlsx")


#5. Data cleaning ncbi dataframe
# Dropping columns from ncbi dataframe
ncbi <- ncbi[, c(6, 10, 11, 17)] 
print('Modified dataframe:-')
ncbi

# Add column names to rbpimage dataframe
# Rename the columns by index
colnames(ncbi)[1] <- "gene_name"
colnames(ncbi)[4] <- "omim_id"
print('Modified dataframe:-')
ncbi

# Remove duplicates from uniprot dataframe
ncbi <- ncbi[!duplicated(ncbi), ]

# Check if there are any duplicate rows in the uniprot dataframe
if(any(duplicated(ncbi))) {
  print("There are duplicate rows in the ncbi dataframe.")
} else {
  print("There are no duplicate rows in the ncbi dataframe.")
}

# Drop NA values from ncbi dataframe
ncbi <- ncbi %>% 
  drop_na()      # removes all rows that have at least one NA in any column

# Check if there are any NA values in the ncbi dataframe
if(any(is.na(ncbi))) {
  print("There are NA values in the ncbi dataframe.")
} else {
  print("There are no NA values in the ncbi dataframe.")
}

# Get uniprot dataframe as Excel file and store in data/processed folder
write_xlsx(ncbi, "data/processed/ncbi.xlsx")


#6. Data cleaning postar dataframe
# Dropping columns from postar dataframe
postar <- postar[, c(1, 5)] 
print('Modified dataframe:-')
postar

# Add column names to rbpimage dataframe
# Rename the columns by index
colnames(postar)[1] <- "gene_name"
colnames(postar)[2] <- "function"
print('Modified dataframe:-')
postar

# Remove duplicates from uniprot dataframe
postar <- postar[!duplicated(ncbi), ]

# Check if there are any duplicate rows in the uniprot dataframe
if(any(duplicated(postar))) {
  print("There are duplicate rows in the postar dataframe.")
} else {
  print("There are no duplicate rows in the postar dataframe.")
}

# Drop NA values from ncbi dataframe
postar <- postar %>% 
  drop_na()      # removes all rows that have at least one NA in any column

# Check if there are any NA values in the ncbi dataframe
if(any(is.na(postar))) {
  print("There are NA values in the postar dataframe.")
} else {
  print("There are no NA values in the postar dataframe.")
}

# Get uniprot dataframe as Excel file and store in data/processed folder
write_xlsx(postar, "data/processed/postar.xlsx")







