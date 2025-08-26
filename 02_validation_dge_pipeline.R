
##############################################################
# Script: 02_validation_deg_pipeline.R
# Purpose: Preprocess RNA-seq data from the CHAIN validation cohort
#          and generate a DEG list (alive vs dead at hospital discharge)
# Input: Transcriptome matrix from CHAIN, clinical metadata
# Output: Annotated DEG results and voom-normalized data
# Author: CJS
# Date: 2025-08-06
##############################################################

# ========================
# Load Required Libraries
# ========================
library(edgeR)
library(limma)
library(biomaRt)
library(readxl)
library(tibble)
library(ggrepel)
library(ggpubr)
library(stringr)

# ============================
# Load CHAIN Expression Data
# ============================
chain.read1 <- read_excel("data/Read1_host transcriptome.xlsx")
chain.countdata <- chain.read1
rownames(chain.countdata) <- chain.countdata$record_id
rn <- rownames(chain.countdata)
chain.countdata <- chain.countdata[ , -1:-6]
rownames(chain.countdata) <- rn
colnames(chain.countdata) <- str_split_fixed(colnames(chain.countdata), "\\.", 2)[,1]
chain.countdata <- t(chain.countdata)

# ========================
# Load CHAIN Metadata
# ========================
chain.metadata <- read.csv("data/czh_chain_data.csv", stringsAsFactors = FALSE)
chain.outcome <- chain.metadata[, c("record_id", "died", "hivstatus_adm", "age_adm", 
                                    "muac_enrol", "site", "adm_sex", "temp_adm", 
                                    "categ_enrol", "ref_deathdate")]
rownames(chain.outcome) <- chain.outcome$record_id

# Subset to samples with both expression and metadata
samples <- base::intersect(colnames(chain.countdata), rownames(chain.outcome))
chain.outcome <- chain.outcome[samples, ]
chain.countdata <- chain.countdata[, samples]
stopifnot(identical(rownames(chain.outcome), colnames(chain.countdata)))

# ============================
# Define Outcome Categories
# ============================
chain.outcome$outcome <- ifelse(chain.outcome$died == 1, "dead", "alive")
chain.outcome$outcome <- ifelse(chain.outcome$categ_enrol == "Community", "Community", chain.outcome$outcome)

# Exclude community-enrolled children
chain.outcome <- chain.outcome[chain.outcome$outcome != "Community", ]
chain.countdata <- chain.countdata[, rownames(chain.outcome)]
stopifnot(identical(rownames(chain.outcome), colnames(chain.countdata)))

# ======================
# Create DGEList Object
# ======================
dge <- DGEList(chain.countdata)
dge$samples$group <- as.factor(chain.outcome$outcome)
dge$samples$batch <- as.factor(chain.outcome$site)
dge$samples$sex <- as.factor(chain.outcome$adm_sex)

# Continuous covariates (not modeled here but extracted for future use)
temp_axilla <- as.numeric(chain.outcome$temp_adm)
age <- as.numeric(chain.outcome$age_adm)

# ================================
# Add Gene Annotations via biomaRt
# ================================
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
               values = rownames(dge$counts),
               mart = ensembl,
               uniqueRows = TRUE,
               useCache = FALSE)

# Remove duplicates and align annotation
annot <- annot[!duplicated(annot$ensembl_gene_id), ]
rownames(annot) <- annot$ensembl_gene_id
genes <- annot[base::intersect(rownames(dge), rownames(annot)), ]

# Subset DGEList to matched genes
dge <- dge[rownames(genes), ]
dge$genes <- genes
stopifnot(identical(rownames(dge$counts), dge$genes$ensembl_gene_id))

# Remove duplicated gene symbols
dge$genes <- dge$genes[!duplicated(dge$genes$external_gene_name), ]
dge$counts <- dge$counts[rownames(dge$genes), ]
stopifnot(identical(rownames(dge$counts), dge$genes$ensembl_gene_id))


# ==========================
# Normalize Expression (TMM)
# ==========================
dge <- calcNormFactors(dge, method = "TMM")

# ==========================
# Design Matrix & Contrasts
# ==========================
group <- dge$samples$group
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- make.names(colnames(design), allow_ = TRUE)

# Contrast: alive vs dead
contr.matrix <- makeContrasts(
  survivor.vs.non.survivor = alive - dead,
  levels = colnames(design)
)

# =======================
# Voom Transformation
# =======================
par(mfrow = c(1, 2), mar = c(5, 5, 5, 5))
v <- voom(dge, design, plot = TRUE, save.plot = TRUE)

# =============================
# Fit Linear Model & Shrinkage
# =============================
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main = "Final model: Mean-variance trend")

# ===========================
# Extract DEGs (adjusted p < 0.05)
# ===========================
dt <- decideTests(efit, adjust.method = "BH")
deg <- topTable(efit, coef = "survivor.vs.non.survivor", number = Inf, adjust.method = "BH")
significant <- deg[deg$adj.P.Val < 0.05, ]

# ====================
# Volcano Plot (ggplot)
# ====================
df_limma <- tibble(
  log2foldchange = efit$coefficients[, 1],
  pvalue = efit$p.value[, 1],
  ensembl_id = efit$genes$ensembl_gene_id,
  gene_symbol = efit$genes$external_gene_name,
)

dfl.chain=df_limma
write.csv(dfl.chain,"data/dfl.chain.csv")

