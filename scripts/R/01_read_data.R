# Import required libraries
library(dplyr)
library(tidyverse)
library(UpSetR)
library(readxl)
library(readr)



# Import raw data into objects using assignment operator

rbpworld <- read_csv("data/raw/rbpworld/Homo_sapiens-RBPs.csv")
rbpgo <- read_tsv("data/raw/rbpgo/Table_HS_RBP.txt")
rbpimage <- read_excel("data/raw/rbpimage/Hela_FeaturesMatrix.xlsx")
uniprot <- read_tsv("data/raw/uniprot/uniprotkb_keyword_KW_0694_AND_organism_2026_02_17.tsv")
ncbi <- read_tsv("data/raw/ncbi/gene_result.txt")
postar <- read.csv("data/raw/postar3/POSTAR3_CLIPdb_module_browse_RBP_sub_info.csv")


