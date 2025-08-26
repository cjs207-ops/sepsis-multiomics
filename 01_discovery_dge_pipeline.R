##########################################################
# Script: 01_discovery_deg_pipeline.R
# Purpose: Preprocess RNA-seq data from the discovery cohort 
#          and generate differentially expressed gene (DEG) list
# Input: Raw HTSeq count matrix, phenotype metadata
# Output: DEG results and annotated data objects
# Author: CJS
# Date: 2023-08-06
##########################################################

# ========================
# Load Required Libraries
# ========================
library(edgeR)
library(limma)
library(biomaRt)
library(tibble)
library(ggplot2)
library(ggrepel)

# ==============================
# Load Raw Count and Metadata
# ==============================
# Set working directory to project root (optional if using here::here)
# setwd("path/to/project")

countdata <- as.matrix(read.csv("data/globin.depleted.htseq.count.data.csv", row.names = 1))
phenodata <- read.csv("data/phenodata.csv", row.names = 1, stringsAsFactors = FALSE)


# ==========================
# Clean Sample Identifiers
# ==========================
colnames(countdata) <- gsub("_htseq_counts.txt", "", sub("globin_depleted_WTCHG_", "", colnames(countdata)))

# ==========================
# Remove QC Failures / Outliers
# ==========================
outliers <- c('202118','275127','209107','265102','276139','267126','287176', '270162', '277151','266114')
phenodata <- phenodata[!(rownames(phenodata) %in% outliers), ]
countdata <- countdata[, !(colnames(countdata) %in% outliers)]

# ==========================
# Keep Only In-Hospital Deaths + Survivors
# ==========================
phenodata <- subset(phenodata, death_status %in% c("inhospital", "survivors"))

# Align phenotype and count data
countdata <- countdata[, match(rownames(phenodata), colnames(countdata))]
stopifnot(identical(rownames(phenodata), colnames(countdata)))

# ======================
# Create DGEList Object
# ======================
dge <- DGEList(countdata)
dge$samples$group <- as.factor(phenodata$outcome)
dge$samples$batch <- as.factor(phenodata$batch)
dge$samples$qc <- as.factor(phenodata$QC)
dge$samples$sex <- as.factor(phenodata$sex)

# Additional continuous covariates
temp_axilla <- as.numeric(phenodata$temp_axilla)
age <- as.numeric(phenodata$age_months)

# ==============================
# Add Gene Annotations from Ensembl
# ==============================
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "gene_biotype", "description"),
               filters = "ensembl_gene_id",
               values = rownames(dge$counts),
               mart = ensembl,
               uniqueRows = TRUE,
               useCache = FALSE)


# Deduplicate annotations
annot <- annot[!duplicated(annot$ensembl_gene_id), ]
rownames(annot) <- annot$ensembl_gene_id
genes <- na.omit(annot[rownames(dge$counts), ])
genes<-genes[rownames(dge$counts),]
is<-base::intersect(rownames(dge$counts), rownames(genes))
dge$counts<-dge$counts[is,]
genes<-genes[is,]
stopifnot(identical(rownames(dge$counts), genes$ensembl_gene_id))



# Assign to DGEList
dge$genes <- genes

# Remove genes with duplicated external names (optional cleanup)
dge$genes <- dge$genes[!duplicated(dge$genes$external_gene_name), ]
dge$counts <- dge$counts[rownames(dge$genes), ]
stopifnot(identical(rownames(dge$counts), dge$genes$ensembl_gene_id))


# ======================
# Normalize Counts (TMM)
# ======================
dge <- calcNormFactors(dge, method = "TMM")

# =========================
# Differential Expression
# =========================
design <- model.matrix(~0 + dge$samples$group)
colnames(design) <- gsub("dge\\$samples\\$group", "", colnames(design))
colnames(design) <- make.names(colnames(design), allow_ = T)

# Contrast: survivors vs non-survivors
contr.matrix <- makeContrasts(
  survivor.vs.non.survivor = survivor - non.survivor,
  levels = colnames(design)
)

# ====================
# Voom Transformation
# ====================
par(mfrow = c(1,2), mar = c(5,5,5,5))
v <- voom(dge, design, plot = TRUE, save.plot = TRUE)

# Fit model and shrink variances
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main = "Final model: Mean-variance trend")

# ===============================
# Extract DEGs (adjusted p < 0.05)
# ===============================
dt <- decideTests(efit, adjust.method = "BH")
deg <- topTable(efit, coef = "survivor.vs.non.survivor", number = Inf, adjust.method = "BH")
significant <- deg[deg$adj.P.Val < 0.05, ]

# ======================
# Final DEG gene list
# ======================
df_limma <- tibble(

  log2foldchange = efit$coefficients[, 1],
  pvalue = efit$p.value[, 1],
  ensembl_id = efit$genes$ensembl_gene_id,
  gene_symbol = efit$genes$external_gene_name,
)

dfl.disco<-df_limma

write.csv(dfl.disco,"data/dfl.disco.csv")


# ---- Export discovery cohort GEX matrix for downstream scripts ----

# 1) log-CPM matrix
lcpm <- edgeR::cpm(dge, log = TRUE)

# 2) Gene names on rows (handle duplicates safely)
gene_names <- dge$genes$external_gene_name
# if duplicated gene symbols exist, make them unique to avoid column clashes after transpose
gene_names <- make.unique(ifelse(is.na(gene_names) | gene_names == "", rownames(lcpm), gene_names))
rownames(lcpm) <- gene_names

# 3) Samples on rows, genes on columns
data_expr <- as.data.frame(t(lcpm))

# 4) Sanity checks: sample alignment
stopifnot(identical(rownames(data_expr), rownames(phenodata)))

# 5) Bind minimal metadata up front (serial, outcome) — adjust names if needed
stopifnot(all(c("serial","outcome") %in% colnames(phenodata)))
datax <- phenodata[, c("serial", "outcome"), drop = FALSE] |>
  dplyr::bind_cols(data_expr)

# 6) Ensure deterministic column order: meta then genes (optional but tidy)
meta_cols <- c("serial", "outcome")
gene_cols <- setdiff(colnames(datax), meta_cols)
datax <- datax[, c(meta_cols, sort(gene_cols)), drop = FALSE]

# 7) Write to standardized repo path
write.csv(datax, file = "data/datax.disco.csv", row.names = TRUE)



