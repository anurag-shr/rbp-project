# Clean and prepare data for analysis
# Load necessary libraries
library(dplyr)
library(readr)
library(writexl)

#1. Data cleaning rbpdb dataframe
# Drop columns which are not required

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

# Remove duplicate from rbpdb dataframe
rbpdb <- rbpdb[!duplicated(rbpdb), ]

# Check if there are any duplicate rows in the rbpdb dataframe
if(any(duplicated(rbpdb))) {
  print("There are duplicate rows in the rbpdb dataframe.")
} else {
  print("There are no duplicate rows in the rbpdb dataframe.")
}

# Drop NA values from rbpdb dataframe
rbpdb <- na.omit(rbpdb)

# Check if there are any NA values in the rbpdb dataframe
if(any(is.na(rbpdb))) {
  print("There are NA values in the rbpdb dataframe.")
} else {
  print("There are no NA values in the rbpdb dataframe.")
}

# Get rbpdb dataframe as Excel file and store in data/processed folder
write_xlsx(rbpdb, "data/processed/rbpdb.xlsx")


#2. Data cleaning rbpgo dataframe
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


#3. Data cleaning rbpimage dataframe
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


#4. Data cleaning rbpworld dataframe 
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






