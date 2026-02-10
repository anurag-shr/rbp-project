# Import required libraries
library(dplyr)
library(tidyverse)
library(UpSetR)
library(readxl)
library(readr)



# Import raw data into objects using assignment operator

rbpworld <- read_csv("data/raw/rbpworld/Homo_sapiens-RBPs.csv")
rbpgo <- read_tsv("data/raw/rbpgo/Table_HS_RBP.txt")
rbpdb <- read_csv("data/raw/rbpdb/RBPDB_v1.3.1_proteins_human_2012-11-21.csv")
rbpimage <- read_excel("data/raw/rbpimage/Hela_FeaturesMatrix.xlsx")


