# 07_dmsr_metric.R
# -------------------
# Title: Distance from Median Survivor Response (DMSR) metric from NMDS of plasma proteome
# Author: CJS
# Date: 2025-08-14
#
# Description
# -----------
# This script computes the DMSR metric using an NMDS ordination of selected immune-related
# plasma proteins (SomaScan). DMSR is defined as the Euclidean distance, in NMDS space,
# from each sample to the centroid (median) of the survivor cluster. The script then:
#   1) Visualizes the NMDS with group ellipses and survivor centroid.
#   2) Plots DMSR across outcomes (Community, alive, dead) with CIs and pairwise tests.
#   3) Tests association between DMSR and survival time among non-survivors (early deaths).
#   4) Correlates DMSR with gene expression for a curated immune gene set.
#

# Dependencies
# ------------
# CRAN: data.table, dplyr, tidyr, stringr, ggplot2, ggrepel, vegan, scales,
#       wesanderson, hrbrthemes (optional), ggh4x, ggbeeswarm, broom, AnnotationDbi,
#       org.Hs.eg.db, GO.db
#
# ==============================================================================

# ========================
# CONFIG
# ========================
SOMASCAN_FILE <- "data/chain.somascan.csv"   # e.g. "Somascan_data/chain.somascan.csv"
TBL_FILE      <- "data/chain.tbl.csv"        # e.g. "Somascan_data/chain.tbl.csv"
METADATA_FILE <- "data/czh_chain_data.csv"   # e.g. "Metadata/czh_chain_data.csv"
GEX_FILE      <- "data/data.m.csv"           # Long-format expression table with cols: record_id, gene, outcome2, expression

OUTPUT_DIR    <- "data/output"
SEED_NMDS     <- 1234

# NMDS settings
NMDS_K        <- 2
NMDS_TRYMAX   <- 100

# Thresholds/filters
EARLY_DEATH_WINDOW <- 10  # keep deaths occurring within 10 days of admission for NMDS/DMSR computations

# ========================
# Libraries
# ========================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ggbeeswarm)
  library(vegan)
  library(scales)
  library(wesanderson)
  library(ggh4x)
  library(broom)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

# Optional theme (fallback if not installed)
.use_ipsum <- requireNamespace("hrbrthemes", quietly = TRUE)
if (.use_ipsum) {
  theme_base <- hrbrthemes::theme_ipsum()
} else {
  theme_base <- theme_minimal(base_size = 12)
}

# ========================
# Helpers
# ========================
.dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

.check_cols <- function(df, required, df_name = "data") {
  missing <- setdiff(required, colnames(df))
  if (length(missing)) {
    stop(sprintf("Missing required columns in %s: %s", df_name, paste(missing, collapse = ", ")))
  }
}

# Ellipse for ggplot
veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(circle %*% chol(cov)))
}

# Summaries with CI for a variable by group
summarise_ci <- function(df, value, group) {
  value <- rlang::ensym(value)
  group <- rlang::ensym(group)
  df %>%
    group_by(!!group) %>%
    summarise(
      mean = mean(!!value, na.rm = TRUE),
      sd   = sd(!!value, na.rm = TRUE),
      n    = sum(is.finite(!!value)),
      se   = sd / sqrt(pmax(n, 1)),
      lower = mean - qt(1 - (0.05/2), pmax(n - 1, 1)) * se,
      upper = mean + qt(1 - (0.05/2), pmax(n - 1, 1)) * se,
      .groups = "drop"
    )
}

# Safe rescale to [0,100]
rescale_0_100 <- function(x) {
  scales::rescale(x, to = c(0, 100))
}

# ========================
# Data Ingest & Wrangling
# ========================
.message <- function(...) cat(sprintf("[DMSR] %s\n", sprintf(...)))
.dir_create(OUTPUT_DIR)

.message("Loading data…")
chain_soma <- fread(SOMASCAN_FILE, sep = "auto")
if ("V1" %in% names(chain_soma)) chain_soma$V1 <- NULL

# Extract record_id from SubjectID pattern: "<something> - <record_id> <rest>"
if (!"record_id" %in% names(chain_soma)) {
  .check_cols(chain_soma, c("SubjectID"), "chain_soma")
  chain_soma$record_id <- trimws(
    stringr::str_split_fixed(
      stringr::str_split_fixed(chain_soma$SubjectID, " - ", 2)[, 2],
      "\\ ", 2
    )[, 1]
  )
}

chain_tbl <- fread(TBL_FILE, sep = "auto")
chain_metadata <- fread(METADATA_FILE, sep = "auto")

# Keep explicit cols
keep_cols <- c("hivstatus_adm", "age_adm", "muac_enrol", "site", "temp_adm", "categ_enrol",
               "record_id", "date_adm", "sc_lastvitalstat_date", "died")
chain_metadata <- chain_metadata[, ..keep_cols]

# Enrich metadata
.message("Processing metadata…")
rownames(chain_metadata) <- as.character(chain_metadata$record_id)
chain_metadata$outcome <- ifelse(chain_metadata$died == 1, "dead", "alive")
chain_metadata$outcome2 <- ifelse(chain_metadata$categ_enrol == "Community", "Community", chain_metadata$outcome)
chain_metadata$date_adm <- as.Date(chain_metadata$date_adm, "%Y-%m-%d")
chain_metadata$sc_lastvitalstat_date <- as.Date(chain_metadata$sc_lastvitalstat_date, "%Y-%m-%d")
chain_metadata$time_to_death <- as.numeric(chain_metadata$sc_lastvitalstat_date - chain_metadata$date_adm)
chain_metadata$time_to_death <- ifelse(chain_metadata$outcome == "dead", chain_metadata$time_to_death, NA)
chain_metadata$record_id <- as.character(chain_metadata$record_id)

# Merge SomaScan with metadata (one row per SubjectID measurement)
soma_metadata <- merge(chain_metadata, chain_soma, by = "record_id")

# ========================
# Curated immune-related proteins 
# ========================
sel <- c(
  "IL18BP","IL1A","IL1B","IL1F10",
  "IL1F3","IL1RA","IL1F5","IL1F6","IL1F7",
  "IL1F8","IL1RL2","IL1F9","IL33",
  # IL1 receptors
  "IL18R1","IL18RAP","IL1R1","IL1R2","IL1R3","IL1R8","IL1R9","IL1RL1","SIGIRR",
  # TNF family
  "BAFF","4-1BBL","TNFSF8","CD70","CD95L",
  "CD178","EDA-A1","TNFSF14","LTA","TNFB",
  "LTB","TNFA","TNFSF10","TNFSF11","TNFSF12",
  "TNFSF13","TNFSF15","TNFSF4","TNF",
  # TNF receptors
  "TNFRSF9","BAFFR","TNFRSF7","CD40","CD95","TNFRSF6B","TNFRSF21","EDA2R",
  "EDAR","TNFRSF19L","TNFR1","TNFR2","TNFRSF11A","TNFRSF11B",
  "TNFRSF12A","TNFRSF13B","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF19","TNFRSF25",
  "LTBR","TNFRSF4","TNFRSF8","TRAILR1","TRAILR2","TRAILR3","TRAILR4",
  # Interferon related genes
  "IFNA1","IFNA10","IFNA13","IFNA14","IFNA2","IFNA4","IFNA7",
  "IFNB1","IFNE","IFNG","IFNZ","IFNA8","IFNA5",
  "IFNAG","IFNW1",
  # Interferon receptors
  "IFNAR1","IFNAR2","IFNGR1","IFNGR2",
  # IL6 family
  "CLCF1","CNTF","IL11","IL31","IL6","LEP","LIF","OSM",
  # IL6 receptors
  "CNTFR","IL11RA","IL6R","LEPR","LIFR","OSMR","IL31RA",
  # IL10 family
  "IL10","IL19","IL20","IL22","IL24","IL28B","IL28A","IL29",
  # IL10 family receptors
  "IL10RA","IL10RB","IL20RA","IL20RB","IL22RA2","IL22R",
  # TGF beta family
  "TGFB1","TGFB2","TGFB3",
  # TGF beta family receptors
  "ACVR1C","ATF2","CD105","ENG","TGFBR1","TGFBR2","TGFBR3",
  # integrins
  "ITGAL","ITGAM","ITGAX","ITGAD",
  # T cells and antigen presentation
  "CD38","CCR4","CCR6","CCR10","CXCR3","CXCR5","CCR7","LCK","ZAP70","CD4","CD8A","CD8B","CD3E","TAP1","PSME1","B2M","CALR","PDIA3",
  # interferon stimulated genes
  "RSAD2","IFI35","IFI44L","IFI6","IFIT1","IFIT2","IFIT3","MX1","MX2","MOV10","OAS3","IFI35","MOV10","OAS3","IFI35",
  # innate antimicrobial effectors
  "CAMP","GZMB","GZMA",
  # chemokines
  "CCL1","TCA3","CCL11","CCL12","MCP-5","CCL13",
  "MCP-4","CCL14","CCL15","CCL16","CCL17","TARC",
  "CCL18","CCL19","CCL2","MCP1","CCL20","CCL21","CCL22",
  "MDC","CCL23","CCL24","CCL25","CCL26","CCL27","CCL28",
  "CCL3","CCL3L3","CCL4","CCL4L1","LAG1","CCL5","CCL6","CCL7",
  "CCL8","CCL9","CX3CL1","CXCL1","CXCL10","CXCL11","CXCL12",
  "CXCL13","CXCL14","CXCL15","CXCL16","CXCL17","CXCL2",
  "MIP2","CXCL3","CXCL4","CXCL5","CXCL6","CXCL7",
  "Ppbp","CXCL9","IL8","CXCL8","XCL1","XCL2","FAM19A1",
  "FAM19A2","FAM19A3","FAM19A4","FAM19A5",
  "IL12A","IL12B","IL1R1","IL37","TNF","IL36RN","IL36B",
  "IL18RAP","IL18R1","IL10"
)

.message("Selecting aptamers for curated immune proteins…")
.check_cols(chain_tbl, c("EntrezGeneSymbol", "AptName"), "chain_tbl")

apts_df <- chain_tbl[chain_tbl$EntrezGeneSymbol %in% sel, ]
apts    <- apts_df$AptName

# Reduce soma_metadata to means per record_id for seq.* columns
.message("Computing per-record means for aptamers…")
soma_metadata <- as.data.frame(soma_metadata)
seq_cols <- c("record_id", grep("^seq\\.", colnames(soma_metadata), value = TRUE))
soma_metadata_seqs <- soma_metadata[, seq_cols]

soma_metadata_means <- soma_metadata_seqs %>%
  group_by(record_id) %>%
  summarise(across(starts_with("seq."), mean), .groups = "drop")

soma_metadata_md <- soma_metadata[!duplicated(soma_metadata$record_id),
                                  c("record_id","hivstatus_adm","age_adm","muac_enrol","site","temp_adm",
                                    "categ_enrol","date_adm","sc_lastvitalstat_date","died","outcome","outcome2","time_to_death")]
stopifnot(identical(soma_metadata_means$record_id, soma_metadata_md$record_id))

soma_metadata <- cbind(soma_metadata_md, soma_metadata_means[, -1])

# Build data matrix: rows = record_id, columns = selected aptamers
.message("Building analysis matrix…")
.check_cols(soma_metadata, c("record_id"), "soma_metadata")

seq_keep <- apts
missing_seq <- setdiff(seq_keep, colnames(soma_metadata))
if (length(missing_seq)) {
  warning(sprintf("%d selected aptamer columns not found and will be dropped: %s",
                  length(missing_seq), paste(head(missing_seq, 10), collapse = ", ")))
}
seq_keep <- base::intersect(seq_keep, colnames(soma_metadata))

# Mean per record_id for selected aptamers (already meaned); construct table
mat <- soma_metadata[, c("record_id", seq_keep), drop = FALSE]

# Rename columns to \"AptName_Gene\" using chain_tbl mapping
apts_df$full.name <- paste(apts_df$AptName, apts_df$EntrezGeneSymbol, sep = "_")
cns <- dplyr::tibble(AptName = colnames(mat)[-1]) %>%
  left_join(apts_df[, c("AptName", "full.name")], by = "AptName")
colnames(mat)[-1] <- cns$full.name

# Keep only valid columns
mat <- mat[, !is.na(colnames(mat))]

# Early-death filter (retain non-deaths and deaths outside window; drop deaths beyond window per original logic)
.message("Applying early-death filter (<= %d days kept; later deaths removed)…", EARLY_DEATH_WINDOW)
keep_idx <- with(soma_metadata, ifelse(time_to_death >= EARLY_DEATH_WINDOW & outcome2 == "dead", FALSE, TRUE))
mat <- mat[keep_idx, , drop = FALSE]
mat <- na.omit(mat)

# ========================
# NMDS Ordination
# ========================
.message("Running NMDS (k=%d, trymax=%d)…", NMDS_K, NMDS_TRYMAX)
set.seed(SEED_NMDS)
# rescale features to [0,100]
mat_scaled <- mat
mat_scaled[, -1] <- lapply(mat_scaled[, -1], rescale_0_100)

nmds <- metaMDS(mat_scaled[, -1], k = NMDS_K, trymax = NMDS_TRYMAX)
.message("NMDS stress: %.4f", nmds$stress)

points <- data.frame(record_id = mat_scaled$record_id, nmds$points)
points <- left_join(points, soma_metadata, by = "record_id")
points$outcome2 <- factor(points$outcome2, levels = c("Community", "alive", "dead"))

# Ellipses by outcome group
.message("Computing group ellipses…")
df_ell <- data.frame()
for (g in levels(points$outcome2)) {
  sub <- subset(points, outcome2 == g)
  if (nrow(sub) >= 3) {
    ell <- veganCovEllipse(
      cov.wt(cbind(sub$MDS1, sub$MDS2), wt = rep(1/length(sub$MDS1), length(sub$MDS1)))$cov,
      center = c(mean(sub$MDS1), mean(sub$MDS2))
    )
    df_ell <- rbind(df_ell, cbind(as.data.frame(ell), group = g))
  }
}
colnames(df_ell)[1:2] <- c("MDS1", "MDS2")

# Colors
cols <- wesanderson::wes_palette(n = 5, name = "Zissou1")[c(1, 3, 5)]

# Centroids
points_dist <- points %>% select(record_id, MDS1, MDS2, outcome2)
centroids <- points_dist %>% group_by(outcome2) %>% summarise(mean.MDS1 = mean(MDS1), mean.MDS2 = mean(MDS2), .groups = "drop")

# ========================
# Plot: NMDS ordination with ellipses
# ========================
.message("Plotting NMDS ordination…")
p_nmds <- ggplot(points, aes(x = MDS1, y = MDS2, color = outcome2)) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_point(data = centroids, aes(x = mean.MDS1, y = mean.MDS2, color=outcome2), size = 3, inherit.aes = FALSE) +
  geom_path(data = df_ell, aes(x = MDS1, y = MDS2, colour = group), size = 0.5, linetype = "2121", inherit.aes = FALSE) +
  coord_flip() +
  labs(x = "MDS2", y = "MDS1", color = "Outcome") +
  theme_classic()+
  theme_ipsum() +
  theme(legend.position="bottom",
        axis.title.x  =element_text(size=14),
        axis.title.y  =element_text(size=14),
        axis.text.x = element_text(angle=0,hjust=1,color="black",size=14),
        axis.text.y = element_text(angle=0,hjust=1,color="black",size=14),
        panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        
  ) +
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)

p_nmds
ggsave(file.path(OUTPUT_DIR, "nmds_ordination.png"), p_nmds, width = 5, height = 5, dpi = 300)

# ========================
# DMSR: Distance to survivor centroid
# ========================
.message("Computing DMSR distances to survivor centroid…")
# Use alive centroid as in original code
alive_centroid <- centroids %>% filter(outcome2 == "alive")
if (nrow(alive_centroid) != 1) stop("Alive centroid not available for DMSR.")

points_dist <- points_dist %>%
  mutate(
    alive.mean.MDS1 = alive_centroid$mean.MDS1[1],
    alive.mean.MDS2 = alive_centroid$mean.MDS2[1],
    distance = sqrt((MDS1 - alive.mean.MDS1)^2 + (MDS2 - alive.mean.MDS2)^2)
  )

# ========================
# DMSR by outcome (incl. Community); pairwise comparisons
# ========================
.message("Summarising DMSR by outcome…")
# IL10 panel relied on merging with IL10 GEX, but DMSR is outcome-level; summarise directly
points_dist$outcome2 <- factor(points_dist$outcome2, levels = c("Community", "alive", "dead"))
summary_by_outcome <- summarise_ci(points_dist, distance, outcome2)

p_dmsr <- ggplot(points_dist, aes(y = distance, x = outcome2)) +
  ggbeeswarm::geom_quasirandom(color = "gray80", width = 0.2, size = 0.7) +
  geom_hline(yintercept = mean(points_dist$distance[points_dist$outcome2 == "alive"], na.rm = TRUE),
             linetype = "2121") +
  geom_errorbar(data = summary_by_outcome, aes(y = mean, ymin = lower, ymax = upper), width = 0.2) +
  geom_point(data = summary_by_outcome, aes(y = mean), size = 1.2) +
  labs(y = "Proteome deviation from median\nsurvival response (DMSR)", x = "") +
  theme_base +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, colour = "gray50")) +
  scale_color_manual(values = cols[-1]) +
  scale_fill_manual(values = cols[-1])

p_dmsr
ggsave(file.path(OUTPUT_DIR, "dmsr_by_outcome.png"), p_dmsr, width = 4, height = 5, dpi = 300)

# Optional pairwise tests (alive vs dead, alive vs Community)
try({
  suppressPackageStartupMessages(library(ggpubr))
  p_dmsr_tests <- p_dmsr +
    ggpubr::stat_compare_means(comparisons = list(c("alive", "dead")), label = "p.signif", size = 4, vjust = 0.65) +
    ggpubr::stat_compare_means(comparisons = list(c("alive", "Community")), label = "p.signif", size = 4, vjust = 0.2)
  ggsave(file.path(OUTPUT_DIR, "dmsr_by_outcome_tests.png"), p_dmsr_tests, width = 4, height = 5, dpi = 300)
}, silent = TRUE)

# ========================
# DMSR vs. survival time among deaths (capped at 8-10 days)
# ========================
.message("Analysing DMSR vs survival time among deaths…")
points_dist2 <- left_join(points_dist, chain_metadata[, c("record_id", "time_to_death")], by = "record_id")
points_dist2 <- points_dist2[!is.na(points_dist2$time_to_death), ]
points_dist2$time_to_death <- ifelse(points_dist2$time_to_death > 7, 8, points_dist2$time_to_death)

summary_time <- points_dist2 %>%
  group_by(time_to_death) %>%
  summarise(mean = median(distance, na.rm = TRUE), sd = sd(distance, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(se = sd / sqrt(pmax(n, 1)),
         lower = mean - qt(1 - (0.05/2), pmax(n - 1, 1)) * se,
         upper = mean + qt(1 - (0.05/2), pmax(n - 1, 1)) * se)

p_time <- ggplot(points_dist2, aes(y = distance, x = time_to_death)) +
  ggbeeswarm::geom_quasirandom(width = 0.1, color = "gray40", pch = 21) +
  geom_point(data = summary_time, aes(y = mean), size = 2.2) +
  geom_smooth(data = summary_time, aes(y = mean, x = as.numeric(time_to_death)), method = "loess", se = FALSE, size = 0.7, linetype = "2121") +
  coord_cartesian(ylim = c(0, 0.45)) +
  labs(y = "Proteome deviation from median\nsurvival response (DMSR)", x = "Survival time (days)") +
  theme_base +
  theme(legend.position = "none",
        axis.text.x = element_text(hjust = 0.5, colour = "gray50")) +
  scale_x_continuous(breaks = 0:8, labels = c("0","1","2","3","4","5","6","7","8-10"))

p_time
ggsave(file.path(OUTPUT_DIR, "dmsr_vs_survival_time.png"), p_time, width = 6, height = 4, dpi = 300)

# ========================
# Correlation: DMSR vs gene expression for immune genes
# ========================
.message("Loading gene expression for DMSR correlations…")
if (file.exists(GEX_FILE)) {
  data_m <- fread(GEX_FILE, sep = "auto")
  data_m <- data_m[, !names(data_m) %in% c("V1", "X"), with = FALSE]
  .check_cols(data_m, c("record_id", "gene", "outcome2", "expression"), "GEX long table")
  
  data_m$record_id <- as.character(data_m$record_id)
  soma_gex <- subset(data_m, gene %in% sel)
  soma_gex <- left_join(soma_gex, points_dist[, c("record_id", "outcome2", "distance")], by = c("record_id", "outcome2"))
  soma_gex <- soma_gex[, c("record_id", "gene", "outcome2", "expression", "distance")]
  soma_gex$expression <- ifelse(soma_gex$expression == 0, NA, soma_gex$expression)
  
  # Spearman per gene (>=3 paired observations)
  res_spear <- soma_gex %>%
    group_by(gene) %>%
    filter(sum(is.finite(expression) & is.finite(distance)) >= 3) %>%
    summarise(
      p.value = cor.test(expression, distance, method = "spearman")$p.value,
      rho     = cor.test(expression, distance, method = "spearman")$estimate,
      .groups = "drop"
    ) %>%
    mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>%
    arrange(adj.p.value)
  
  # Annotate with gene names
  gene_names <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = res_spear$gene,
                                      keytype = "SYMBOL",
                                      columns = c("SYMBOL", "GENENAME"))
  res_spear <- left_join(res_spear, gene_names, by = c("gene" = "SYMBOL"))
  
  # GO flagging for cytokine/chemokine/secreted (optional)
  genes_for_go <- unique(res_spear$gene)
  go_annots <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = genes_for_go,
                                     columns = c("SYMBOL", "GO", "ONTOLOGY"),
                                     keytype = "SYMBOL") %>%
    filter(!is.na(GO))
  go_terms <- AnnotationDbi::select(GO.db,
                                    keys = unique(go_annots$GO),
                                    columns = c("TERM"),
                                    keytype = "GOID")
  go_annot_full <- left_join(go_annots, go_terms, by = c("GO" = "GOID"))
  relevant_terms <- c("cytokine activity", "chemokine activity", "extracellular region")
  go_flags <- go_annot_full %>%
    filter(tolower(TERM) %in% relevant_terms) %>%
    mutate(category = dplyr::case_when(
      TERM == "cytokine activity" ~ "Cytokine",
      TERM == "chemokine activity" ~ "Chemokine",
      TERM == "extracellular region" ~ "Secreted",
      TRUE ~ NA_character_
    )) %>%
    distinct(SYMBOL, category)
  go_summary <- go_flags %>%
    group_by(SYMBOL) %>%
    summarise(GO_annotation = paste(unique(category), collapse = "; "), .groups = "drop")
  
  res_spear <- res_spear %>%
    left_join(go_summary, by = c("gene" = "SYMBOL")) %>%
    filter(!grepl("receptor", GENENAME, ignore.case = TRUE)) %>%
    filter(!is.na(GO_annotation))
  
  # Plot volcano-style bar: -log10(FDR), dot size = |rho|
  res_spear$gene <- factor(res_spear$gene, levels = rev(res_spear$gene))
  res_spear$correlation <- ifelse(res_spear$adj.p.value < 0.05, "Significant", "Not significant")
  
  p_corr <- ggplot(res_spear, aes(y = -log10(adj.p.value), x = gene)) +
    geom_col(color = "black", width = 0.05) +
    geom_point(aes(size = abs(rho), color = correlation)) +
    geom_point(aes(size = abs(rho)), color = "black", pch = 21, stroke = 0.3) +
    coord_flip() +
    scale_size(range = c(0.05, 7)) +
    labs(x = "", y = "Adjusted p.value (1/p)") +
    scale_y_continuous(breaks = 0:5,
                       labels = as.character(round(1/10^(seq(from = 0, to = 5, by = 1)), 4))) +
    theme_base +
    theme(legend.position = "right",
          axis.text.x = element_text(colour = "gray50")) +
    scale_color_manual(values = cols[-2]) +
    scale_fill_manual(values = cols[-2])
  
  fwrite(res_spear, file.path(OUTPUT_DIR, "dmsr_gex_spearman.csv"))
  ggsave(file.path(OUTPUT_DIR, "dmsr_gex_spearman.png"), p_corr, width = 6, height = 7, dpi = 300)
  
} else {
  .message("GEX file not found: %s (skipping DMSR~GEX correlations)", GEX_FILE)
}

p_corr
# ========================
# Save key objects
# ========================
.message("Saving key objects…")
saveRDS(list(
  nmds = nmds,
  points = points,
  centroids = centroids,
  points_dist = points_dist
), file = file.path(OUTPUT_DIR, "dmsr_nmds_objects.rds"))

.message("Done. Outputs written to '%s'", OUTPUT_DIR)

