library(ggVennDiagram)
library(ggplot2)
library(dplyr)
library(sf)

# Define sets
sets <- list(
  RBP_GO  = rbpgo    %>% distinct(gene_name) %>% pull(gene_name),
  RBP_IMG = rbpimage %>% distinct(gene_name) %>% pull(gene_name),
  RBP_WRLD = rbpworld %>% distinct(gene_name) %>% pull(gene_name)
)