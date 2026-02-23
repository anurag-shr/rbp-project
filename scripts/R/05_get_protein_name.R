# Import libraries
library(dplyr)
library(readr)
library(tidyr)
library(biomaRt)
library(uniprotREST)

# Convert common rbp new csv file to dataframe
common_rbps_new <- read_csv("data/processed/common_rbps_new.csv") %>% 
  as_tibble()

# Using uniprotREST package to get the protein names for the genes
# example: your dataframe
# df has a column "uniprot_id" and an empty/NA column "protein_name"
# df <- data.frame(uniprot_id = c("P22682","P47941"), protein_name = NA_character_)

# query UniProt – choose the fields you need
res <- uniprot_map(
  ids    = common_rbps_new$uniprot_id,
  from   = "UniProtKB_AC-ID",
  to     = "UniProtKB",
  format = "tsv",
  fields = c("accession", "protein_name")
)

# Keep specific columns by index
res <- res[, c(1, 3)]
print('Modified dataframe:-')
res


new_names <- c("accession", "accession", "protein_name")
colnames(res) <- new_names

# res will have at least 'accession' and 'protein_name' columns [web:20][web:25]

# join back to your dataframe
df_filled <- common_rbps_new %>%
  left_join(res, by = c("uniprot_id" = "accession")) %>%
  mutate(protein_name = ifelse(is.na(protein_name.x), protein_name.y, protein_name.x)) %>%
  dplyr::select(-protein_name.x, -protein_name.y)

df_filled <- df_filled %>% relocate(protein_name, .before = length_aa)

# Get df_filled as a csv file 
write.csv(df_filled, "data/processed/common_rbps_with_prot_name.csv")



