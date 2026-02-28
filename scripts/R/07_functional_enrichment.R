# Load necessary libraries
library(clusterProfiler)
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggplotify)

# Importing gene dataset

genes <- read_tsv("genes.tsv")

library(clusterProfiler)
library(org.Hs.eg.db)   # Human gene annotation; change for other organisms
library(enrichplot)
library(ggplot2)

# ── 1. Your gene list ────────────────────────────────────────────────────────
# Replace with your actual gene symbols
gene_list <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN", "AKT1", "CDK2", "MDM2")

# ── 2. Convert gene symbols to Entrez IDs ────────────────────────────────────
gene_entrez <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_ids <- gene_entrez$ENTREZID

# ── 3. Gene Ontology (GO) Enrichment ─────────────────────────────────────────
# ont: "BP" (Biological Process), "MF" (Molecular Function),
#       "CC" (Cellular Component), or "ALL"
go_results <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE   # converts Entrez IDs back to gene symbols
)

# View top results
head(as.data.frame(go_results))

# ── 4. KEGG Pathway Enrichment ────────────────────────────────────────────────
kegg_results <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",     # "hsa" for human; e.g., "mmu" for mouse
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Convert KEGG Entrez IDs to symbols for readability
kegg_results <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(as.data.frame(kegg_results))

# ── 5. Visualizations ─────────────────────────────────────────────────────────

# Dot plot (top GO terms)
dotplot(go_results, showCategory = 20, title = "GO Biological Process")

# Bar plot
barplot(go_results, showCategory = 15, title = "GO Enrichment")

# Gene-concept network (links genes to enriched terms)
cnetplot(go_results, showCategory = 5)

# Enrichment map (clusters similar GO terms)
go_results_pair <- pairwise_termsim(go_results)
emapplot(go_results_pair, showCategory = 20)

# KEGG dot plot
dotplot(kegg_results, showCategory = 20, title = "KEGG Pathway Enrichment")

# ── 6. Save results to CSV ────────────────────────────────────────────────────
write.csv(as.data.frame(go_results),   "GO_enrichment_results.csv",   row.names = FALSE)
write.csv(as.data.frame(kegg_results), "KEGG_enrichment_results.csv", row.names = FALSE)