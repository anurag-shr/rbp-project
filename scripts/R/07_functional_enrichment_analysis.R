# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)   # Human gene annotation; change for other organisms
library(enrichplot)
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(scales)
library(ragg)


# ── 1. Genes dataframe  ───────────────────────────────────────────────────────
# Importing data as dataframe
genes <- read_csv("data/processed/common_rbps_with_prot_name.csv")
# Renaming gene_name column in genes dataframe
genes <- rename(genes, gene_symbol = "gene_name")
# Drop first column from the genes dataframe using dplyr
genes <- select(genes, -1)

# ── 2. Convert gene symbols to Entrez IDs ────────────────────────────────────
gene_entrez <- bitr(
  genes$gene_symbol,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_ids <- gene_entrez$ENTREZID

# Merge the Entrez id vector to genes
genes <- genes %>%
  left_join(gene_entrez, 
            by = c("gene_symbol" = "SYMBOL"))

# Rename ENTREZID column from genes
colnames(genes)[10] <- "entrez_id"


# Check if there are any dupicates in genes dataframe
if(any(duplicated(genes))) {
  print("There are duplicate rows in the genes dataframe.")
} else {
  print("There are no duplicate rows in the genes dataframe.")
}

# ── 3. Gene Ontology (GO) Enrichment ─────────────────────────────────────────-

# gene_vector <- unique(na.omit(genes$entrez_id))

# ── 3.1 Ontology: "BP" (Biological Process)  
go_bp_results <- enrichGO(
  gene          = genes$entrez_id,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Removing redundant GO terms
go_bp_simplified <- clusterProfiler::simplify(
  go_bp_results,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# View top results
head(as.data.frame(go_bp_simplified))


# Saving GO results 
write.csv(as.data.frame(go_bp_simplified),   "data/processed/GO_bp_enrichment_results.csv",   row.names = FALSE)


# ── 3.2 Ontology: "MF" (Molecular Function)

go_mf_results <- enrichGO(
  gene          = genes$entrez_id,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Removing redundant GO terms
go_mf_simplified <- clusterProfiler::simplify(
  go_mf_results,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# View top results
head(as.data.frame(go_mf_simplified))


# Saving GO results 
write.csv(as.data.frame(go_mf_simplified),   "data/processed/GO_mf_enrichment_results.csv",   row.names = FALSE)


# ── 3.3  Ontology: "CC" (Cellular Component)

go_cc_results <- enrichGO(
  gene          = genes$entrez_id,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Removing redundant GO terms
go_cc_simplified <- clusterProfiler::simplify(
  go_cc_results,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# View top results
head(as.data.frame(go_cc_simplified))


# Saving GO results 
write.csv(as.data.frame(go_cc_simplified),   "data/processed/GO_cc_enrichment_results.csv",   row.names = FALSE)


# ── 4. Visualizations ─────────────────────────────────────────────────────────
# Bubble plot for GO Biological Function.
# 1) Read your enrichment results (this file has GeneRatio like "76/253") [file:30]
go_bp_df <- readr::read_csv("data/processed/GO_bp_enrichment_results.csv", show_col_types = FALSE)

top_n <- 15

go_bf_plot <- go_bp_df %>%
  mutate(
    # GeneRatio comes as "k/n" -> numeric
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio))
  ) %>%
  filter(!is.na(GeneRatio_num), !is.na(Description), !is.na(Count), !is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  mutate(
    # order terms by significance, most significant at top
    Description = factor(Description, levels = rev(Description))
  )

## 2. Dot plot for GO BP enrichment
p_bf <- ggplot(go_bf_plot,
               aes(y = Description,
                   x = GeneRatio_num,
                   color = p.adjust,
                   size  = Count)) +
  geom_point(alpha = 0.9) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_gradient2(
    low  = "red",
    mid  = "white",
    high = "blue",
    trans = "log10",
    name = "p.adjust"
  ) +
  labs(
    y     = "GO Biological Process",
    x     = "Gene Ratio (%)",
    title = "GO Biological Process Enrichment",
    size  = "Gene count"
  ) +
  guides(
    size  = guide_legend(order = 1),
    color = guide_colourbar(order = 2)
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title    = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.x   = element_text(size = 10, colour = "black"),
    axis.text.y   = element_text(size = 8,  colour = "black"),
    
    legend.title  = element_text(size = 12, face = "bold", colour = "black"),
    legend.text   = element_text(size = 10, colour = "black"),
    
    # Legends outside panel, right side, bottom-aligned, stacked vertically
    legend.position      = "right",        # right margin
    legend.justification = "center",       # vertically centered
    legend.box          = "vertical",      # stack size above color
    legend.box.just     = "bottom",        # align stack to bottom
    legend.margin       = margin(l = 8),   # space between plot and legend
    legend.spacing.y    = unit(0.3, "cm"), # vertical space between legend items
    
    panel.grid.major  = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(colour = "grey92", linewidth = 0.25),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

print(p_bf)



# Bubble plot for GO Molecular Function.
# 1) Read your enrichment results (this file has GeneRatio like "33/262") [file:30]
go_mf_df <- readr::read_csv("data/processed/GO_mf_enrichment_results.csv", show_col_types = FALSE)

top_n <- 15

go_plot_mf <- go_mf_df %>%
  mutate(
    # GeneRatio comes as "k/n" -> numeric
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio))
  ) %>%
  filter(!is.na(GeneRatio_num), !is.na(Description), !is.na(Count), !is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  mutate(
    # order terms by significance, most significant at top
    Description = factor(Description, levels = rev(Description))
  )

## 2. Dot plot for GO MF enrichment
p_mf <- ggplot(go_plot_mf,
               aes(y = Description,
                   x = GeneRatio_num,
                   color = p.adjust,
                   size  = Count)) +
  geom_point(alpha = 0.9) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_gradient2(
    low  = "red",
    mid  = "white",
    high = "blue",
    trans = "log10",
    name = "p.adjust"
  ) +
  labs(
    y     = "GO Molecular Function",
    x     = "Gene Ratio (%)",
    title = "GO Molecular Function Enrichment",
    size  = "Gene count"
  ) +
  guides(
    size  = guide_legend(order = 1),
    color = guide_colourbar(order = 2)
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title    = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.x   = element_text(size = 10, colour = "black"),
    axis.text.y   = element_text(size = 8,  colour = "black"),
    
    legend.title  = element_text(size = 12, face = "bold", colour = "black"),
    legend.text   = element_text(size = 10, colour = "black"),
    
    # Legends outside panel, right side, bottom-aligned, stacked vertically
    legend.position      = "right",        # right margin
    legend.justification = "center",       # vertically centered
    legend.box          = "vertical",      # stack size above color
    legend.box.just     = "bottom",        # align stack to bottom
    legend.margin       = margin(l = 8),   # space between plot and legend
    legend.spacing.y    = unit(0.3, "cm"), # vertical space between legend items
    
    panel.grid.major  = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(colour = "grey92", linewidth = 0.25),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

print(p_mf)



# Bubble plot for GO Cellular Components.
# 1) Read your enrichment results (this file has GeneRatio like "33/262") [file:30]
go_cc_df <- readr::read_csv("data/processed/GO_cc_enrichment_results.csv", show_col_types = FALSE)

top_n <- 15

go_plot_cc <- go_cc_df %>%
  mutate(
    # GeneRatio comes as "k/n" -> numeric
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio))
  ) %>%
  filter(!is.na(GeneRatio_num), !is.na(Description), !is.na(Count), !is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  mutate(
    # order terms by significance, most significant at top
    Description = factor(Description, levels = rev(Description))
  )

## 2. Dot plot for GO CC enrichment
p_cc <- ggplot(go_plot_cc,
               aes(y = Description,
                   x = GeneRatio_num,
                   color = p.adjust,
                   size  = Count)) +
  geom_point(alpha = 0.9) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_gradient2(
    low  = "red",
    mid  = "white",
    high = "blue",
    trans = "log10",
    name = "p.adjust"
  ) +
  labs(
    y     = "GO Cellular Component",
    x     = "Gene Ratio (%)",
    title = "GO Cellular Component Enrichment",
    size  = "Gene count"
  ) +
  guides(
    size  = guide_legend(order = 1),
    color = guide_colourbar(order = 2)
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title    = element_text(size = 12, face = "bold", colour = "black"),
    axis.text.x   = element_text(size = 10, colour = "black"),
    axis.text.y   = element_text(size = 8,  colour = "black"),
    
    legend.title  = element_text(size = 12, face = "bold", colour = "black"),
    legend.text   = element_text(size = 10, colour = "black"),
    
    # Legends outside panel, right side, bottom-aligned, stacked vertically
    legend.position      = "right",        # right margin
    legend.justification = "center",       # vertically centered
    legend.box          = "vertical",      # stack size above color
    legend.box.just     = "bottom",        # align stack to bottom
    legend.margin       = margin(l = 8),   # space between plot and legend
    legend.spacing.y    = unit(0.3, "cm"), # vertical space between legend items
    
    panel.grid.major  = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_line(colour = "grey92", linewidth = 0.25),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

print(p_cc)


# ── 5. Saving the Plots ──────────────────────────────────────────────────────

# Use your existing folder
outdir <- "plots"

# Create it only if missing (won't warn if it already exists)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)  # dir.create details [web:71]

# Set publication-friendly size
w <- 7  # inches
h <- 5  # inches

# A) GO BP Functional Enrichment
# 1) Vector PDF (best for journals; crisp text/lines)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_BP_bubble.pdf"),
  plot = p_bf,
  width = w, height = h, units = "in",
  device = grDevices::cairo_pdf   # Cairo often improves font handling [web:47][web:54]
)

# 2) High-DPI PNG (good for submissions that require raster)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_BP_bubble.png"),
  plot = p_bf,
  width = w, height = h, units = "in",
  dpi = 600,
  device = ragg::agg_png
)

# 3) Optional: TIFF (some journals prefer this)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_BP_bubble.tiff"),
  plot = p_bf,
  width = w, height = h, units = "in",
  dpi = 600,
  compression = "lzw"
)


# A) GO MF Functional Enrichment
# 1) Vector PDF (best for journals; crisp text/lines)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_MF_bubble.pdf"),
  plot = p_mf,
  width = w, height = h, units = "in",
  device = grDevices::cairo_pdf   # Cairo often improves font handling [web:47][web:54]
)

# 2) High-DPI PNG (good for submissions that require raster)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_MF_bubble.png"),
  plot = p_mf,
  width = w, height = h, units = "in",
  dpi = 600,
  device = ragg::agg_png
)

# 3) Optional: TIFF (some journals prefer this)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_MF_bubble.tiff"),
  plot = p_mf,
  width = w, height = h, units = "in",
  dpi = 600,
  compression = "lzw"
)


# A) GO CC Functional Enrichment
# 1) Vector PDF (best for journals; crisp text/lines)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_CC_bubble.pdf"),
  plot = p_cc,
  width = w, height = h, units = "in",
  device = grDevices::cairo_pdf   # Cairo often improves font handling [web:47][web:54]
)

# 2) High-DPI PNG (good for submissions that require raster)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_CC_bubble.png"),
  plot = p_cc,
  width = w, height = h, units = "in",
  dpi = 600,
  device = ragg::agg_png
)

# 3) Optional: TIFF (some journals prefer this)
ggplot2::ggsave(
  filename = file.path(outdir, "GO_CC_bubble.tiff"),
  plot = p_cc,
  width = w, height = h, units = "in",
  dpi = 600,
  compression = "lzw"
)




































# ── 4. KEGG Pathway Enrichment ────────────────────────────────────────────────
kegg_results <- enrichKEGG(
  gene          = genes$entrez_id,
  organism      = "hsa",     # "hsa" for human; 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Convert KEGG Entrez IDs to symbols for readability
kegg_results <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(as.data.frame(kegg_results))



