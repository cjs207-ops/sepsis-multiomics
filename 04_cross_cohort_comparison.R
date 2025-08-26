
##############################################################
# Script: 04_cross_cohort_comparison.R
# Purpose: Compare gene expression fold changes across
#          discovery and validation cohorts; assess correlation;
#          identify shared dysregulated genes; visualize patterns
# Author: CJS
# Date: 2025-08-06
##############################################################

# ========================
# Load Required Libraries
# ========================
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(VennDiagram)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggraph)
library(igraph)
library(cowplot)
library(gridExtra)
library(reshape2)
library(viridis)
library(randomForest)
library(pROC)

# ========================
# Load DEG Lists
# ========================
dfl.chain <- read.csv("data/dfl.chain.csv")
dfl.disco <- read.csv("data/dfl.disco.csv")

conflicted::conflicts_prefer(base::intersect)

# Harmonize by Ensembl ID
rownames(dfl.chain) <- dfl.chain$ensembl_id
rownames(dfl.disco) <- dfl.disco$ensembl_id

common_ids <- intersect(rownames(dfl.chain), rownames(dfl.disco))
dfl.chain <- dfl.chain[common_ids, ]
dfl.disco <- dfl.disco[common_ids, ]

# ========================
# Merge and Annotate
# ========================
dfl.unite <- data.frame(
  gene_symbol = dfl.disco$gene_symbol,
  ensembl_id = dfl.chain$ensembl_id,
  log2foldchange.disco = dfl.disco$log2foldchange,
  log2foldchange.chain = dfl.chain$log2foldchange,
  unadj.pval.disco = dfl.disco$pvalue,
  unadj.pval.chain = dfl.chain$pvalue
)


# Adjust p-values
dfl.unite$adj.pval.disco <- p.adjust(dfl.unite$unadj.pval.disco, method = "fdr")
dfl.unite$adj.pval.chain <- p.adjust(dfl.unite$unadj.pval.chain, method = "fdr")

# Subset genes significant in at least one cohort

dfl.unite.nom.sig <- dfl.unite[dfl.unite$adj.pval.chain<0.05 | dfl.unite$adj.pval.disco<0.05,]

# ========================
# Cross-Cohort Fold Change Correlation
# ========================
ggplot(dfl.unite.nom.sig, aes(x = log2foldchange.chain, y = log2foldchange.disco)) +
  geom_point(alpha = 0.3, color = "gray40") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  xlim(-2, 2) + ylim(-2, 2) +
  stat_cor(
    method = "spearman",
    size = 5,
    color = "red",
    label.x = 0.1,
    label.y = -1.5,
    p.accuracy = 0.001,   
    r.accuracy = 0.1     
  )+
  labs(
    x = expression(paste("Validation cohort – ", log[2], " fold change")),
    y = expression(paste("Discovery cohort – ", log[2], " fold change"))
  ) +
  theme_minimal(base_size = 16)



# =======================================================
# Venn Diagrams: Upregulated and Downregulated Genes
# =======================================================

# Classify direction of regulation
dfl.unite.nom.sig <- dfl.unite.nom.sig %>%
  mutate(
    up.disco = log2foldchange.disco > 0,
    down.disco = log2foldchange.disco < 0,
    up.chain = log2foldchange.chain > 0,
    down.chain = log2foldchange.chain < 0
  )

# ---- DOWNREGULATED GENES ----
disco_down <- unique(dfl.unite.nom.sig[dfl.unite.nom.sig$down.disco, "gene_symbol"])
chain_down <- unique(dfl.unite.nom.sig[dfl.unite.nom.sig$down.chain, "gene_symbol"])
overlap_down <- length(intersect(disco_down, chain_down))
dev.off()
venn_down <- draw.pairwise.venn(
  area1 = length(disco_down),
  area2 = length(chain_down),
  cross.area = overlap_down,
  
  
  category = c("Validation\ncohort", "Discovery\ncohort"),
  #fill = c("#fde725ff", "#21908dff"),
  fill = c(alpha("darkgreen",0.6), alpha("purple",0.6)),
  alpha = 0.5,
  fontfamily = "sans",
  cat.fontfamily="sans",
  lwd=6,
  margin = 0.2,  # Increased margin size
  
  cat.col = c("darkgreen", "purple"),
  cat.pos = c(-30, 30),
  cat.dist = 0.08,
  cex = 0,   
  cat.cex = 1.5,# suppress count font size
  #print.mode = "none",         # suppress printing of counts entirely
  col = "black",
  ext.line.lty = 0
)

grid.text(paste0(length(disco_down), "\ndownregulated\ngenes"), x = 0.14, y = 0.5, gp = gpar(fontsize = 19))
grid.text(paste0(length(chain_down), "\ndownregulated\ngenes"), x = 0.86, y = 0.5, gp = gpar(fontsize = 19))
grid.text(paste0(overlap_down, " overlap"), x = 0.5, y = 0.5, gp = gpar(fontsize = 19))

# ---- UPREGULATED GENES ----
disco_up <- unique(dfl.unite.nom.sig[dfl.unite.nom.sig$up.disco, "gene_symbol"])
chain_up <- unique(dfl.unite.nom.sig[dfl.unite.nom.sig$up.chain, "gene_symbol"])
overlap_up <- length(intersect(disco_up, chain_up))

dev.off()
venn_up <- draw.pairwise.venn(
  area1 = length(disco_down),
  area2 = length(chain_down),
  cross.area = overlap_down,
  
  
  category = c("Validation\ncohort", "Discovery\ncohort"),
  #fill = c("#fde725ff", "#21908dff"),
  fill = c(alpha("darkgreen",0.6), alpha("purple",0.6)),
  alpha = 0.5,
  fontfamily = "sans",
  cat.fontfamily="sans",
  lwd=6,
  margin = 0.2,  # Increased margin size
  
  cat.col = c("darkgreen", "purple"),
  cat.pos = c(-30, 30),
  cat.dist = 0.08,
  cex = 0,   
  cat.cex = 0,# suppress count font size
  #print.mode = "none",         # suppress printing of counts entirely
  col = "black",
  ext.line.lty = 0
)
grid.text(paste0(length(disco_up), " upregulated\ngenes"), x = 0.14, y = 0.5, gp = gpar(fontsize = 19))
grid.text(paste0(length(chain_up), " upregulated\ngenes"), x = 0.86, y = 0.5, gp = gpar(fontsize = 19))
grid.text(paste0(overlap_up, " overlap"), x = 0.5, y = 0.5, gp = gpar(fontsize = 19))

