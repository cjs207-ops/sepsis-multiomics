# ================================================================
# Script 6: Correlation Between Plasma IL10 and Gene Expression
# ================================================================
# This script performs the following:
# 1. Loads SomaScan plasma protein abundance data and GEX data.
# 2. Merges IL10 RFU (plasma) with gene expression from the validation cohort.
# 3. Computes gene-wise correlation between IL10 protein and gene expression.
# 4. Performs pathway enrichment (Reactome) on IL10-correlated genes.
# 5. Identifies immune-related pathways with significant association.
# 6. Visualizes regression slopes and R² by pathway.
# ================================================================

# -----------------------------
# Load libraries
# -----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(hrbrthemes)
library(readr)
library(data.table)
library(readxl)
library(stringr)
library(RColorBrewer)
library(scales)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(reactome.db)

# -----------------------------
# Load input data
# -----------------------------
chain.soma <- read.csv("data/chain.somascan.csv", stringsAsFactors = FALSE)
chain.tbl  <- read.csv("data/chain.tbl.csv", stringsAsFactors = FALSE)
chain.metadata <- read.csv("data/czh_chain_data.csv", stringsAsFactors = FALSE)

# -----------------------------
# Preprocess metadata
# -----------------------------
chain.metadata <- chain.metadata %>%
  select(record_id, hivstatus_adm, age_adm, muac_enrol, site, temp_adm, categ_enrol,
         date_adm, sc_lastvitalstat_date, died) %>%
  mutate(
    record_id = as.character(record_id),
    outcome = ifelse(died == 1, "dead", "alive"),
    outcome2 = ifelse(categ_enrol == "Community", "Community", outcome),
    date_adm = as.Date(date_adm),
    sc_lastvitalstat_date = as.Date(sc_lastvitalstat_date),
    time_to_death = ifelse(outcome == "dead",
                           as.numeric(sc_lastvitalstat_date - date_adm), NA)
  )

# -----------------------------
# Preprocess SomaScan data
# -----------------------------
chain.soma$record_id <- trimws(str_split_fixed(str_split_fixed(chain.soma$SubjectID, " - ", 2)[, 2], "\\ ", 2)[, 1])
soma.metadata <- merge(chain.metadata, chain.soma, by = "record_id")

# Long-form SomaScan with gene mapping
dta <- soma.metadata %>%
  select(record_id, outcome2, all_of(chain.tbl$AptName)) %>%
  group_by(record_id) %>%
  summarise(across(-outcome2, median)) %>%
  right_join(
    soma.metadata[!duplicated(soma.metadata$record_id), c("record_id", "outcome2")],
    by = "record_id"
  )

dta.m <- reshape2::melt(dta, id.vars = c("outcome2", "record_id"))
colnames(dta.m) <- c("outcome2", "record_id", "aptamer", "RFU")
dta.m$outcome2 <- factor(dta.m$outcome2, levels = c("Community", "alive", "dead"))

# Map aptamers to gene symbols
chain.tbl$aptamer <- chain.tbl$AptName
conflicted::conflicts_prefer(dplyr::rename)
dta.m <- dta.m %>%
  left_join(chain.tbl[, c("aptamer", "EntrezGeneSymbol")], by = "aptamer") %>%
  rename(gene = EntrezGeneSymbol) %>%
  mutate(apt.gene = paste(gene, aptamer, sep = "_"))

# -----------------------------
# Merge with gene expression
# -----------------------------
data.m <- fread("data/data.m.csv") %>%
  select(-V1) %>%
  mutate(gex.gene = gene)

gex <- data.m %>%
  select(record_id, time_to_death, muac_enrol, gex.gene, expression) %>%
  mutate(record_id = as.character(record_id))

# -----------------------------
# Merge IL10 RFU with expression data
# -----------------------------
il10_aptamer_id <- "seq.13723.6"
plasma_il10 <- dta.m %>%
  filter(aptamer == il10_aptamer_id) %>%
  mutate(record_id = as.character(record_id))

merged <- left_join(plasma_il10, gex, by = "record_id") %>%
  filter(!is.na(expression)) %>%
  mutate(
    muac = ifelse(muac_enrol < 11.5, "low", "high"),
    rfu.cut = cut(RFU, breaks = c(100, 900, 1000, 1200, 1400, 1900, 2200, 10000))
  ) %>%
  filter(outcome2 != "Community")

# -----------------------------
# Gene-wise correlation with IL10
# -----------------------------
med.sm <- merged %>%
  group_by(gex.gene, rfu.cut) %>%
  summarise(med1 = median(expression)) %>%
  mutate(il10_strata = as.numeric(rfu.cut)) %>%
  summarise(
    r2 = summary(lm(med1 ~ il10_strata))$r.squared,
    p.val = summary(lm(med1 ~ il10_strata))$coefficients[2, 4],
    lm.slope = coef(lm(med1 ~ il10_strata))[2],
    rho = cor.test(med1, il10_strata, method = "spearman")$estimate
  )

# -----------------------------
# Map genes to Reactome pathways
# -----------------------------

#define mapping function
gene_to_reactome <- function(gene_symbols,
                             keep_unmapped = FALSE,
                             add_ids = FALSE,
                             add_hierarchy = TRUE,
                             hierarchy_with_ids = FALSE) {
  stopifnot(is.character(gene_symbols))
  
  # Map HGNC symbols → Entrez
  map_df <- clusterProfiler::bitr(gene_symbols,
                                  fromType = "SYMBOL",
                                  toType   = "ENTREZID",
                                  OrgDb    = org.Hs.eg.db)
  
  # Entrez → Reactome (IDs + names)
  rx <- AnnotationDbi::select(reactome.db,
                              keys    = unique(map_df$ENTREZID),
                              keytype = "ENTREZID",
                              columns = c("PATHID", "PATHNAME"))
  out <- map_df %>%
    inner_join(rx, by = "ENTREZID") %>%
    transmute(SYMBOL,
              Pathway    = PATHNAME,
              ReactomeID = PATHID,
              ENTREZID)
  
  # Add hierarchy if possible
  if (add_hierarchy && nrow(out)) {
    uniq_ids <- unique(out$ReactomeID)
    parent_tbl <- .get_parent_map(uniq_ids)
    
    if (is.null(parent_tbl)) {
      warning("This reactome.db build does not expose parent–child relations; returning without hierarchy. Set add_hierarchy = FALSE to silence.")
      out$Hierarchy <- NA_character_
    } else {
      parents_list <- setNames(parent_tbl$parents, parent_tbl$PATHID)
      name_map <- .get_path_names(unique(c(uniq_ids, unlist(parents_list))))
      hier_map <- lapply(uniq_ids, function(pid) {
        hs <- .build_hierarchies_from_map(pid, parents_list, name_map,
                                          include_ids = hierarchy_with_ids)
        paste(unique(hs), collapse = " | ")
      })
      out <- out %>% left_join(tibble(ReactomeID = uniq_ids,
                                      Hierarchy = unlist(hier_map)),
                               by = "ReactomeID")
    }
  }
  
  # Include unmapped symbols if requested
  if (keep_unmapped) {
    missing_syms <- setdiff(gene_symbols, unique(out$SYMBOL))
    if (length(missing_syms)) {
      add_cols <- tibble(SYMBOL = missing_syms,
                         Pathway = NA_character_,
                         ReactomeID = NA_character_,
                         ENTREZID = NA_character_)
      if (add_hierarchy) add_cols$Hierarchy <- NA_character_
      out <- bind_rows(out, add_cols)
    }
  }
  
  out <- distinct(out) %>% arrange(SYMBOL, Pathway)
  if (!add_ids) {
    sel <- c("SYMBOL", "Pathway", if (add_hierarchy) "Hierarchy")
    out <- out[, sel, drop = FALSE]
  }
  out
}

#map genes to reactome table
genes <- unique(na.omit(med.sm$gex.gene))
tbl <- gene_to_reactome(genes, keep_unmapped = TRUE, add_ids = TRUE, add_hierarchy = TRUE)
tbl$Pathway <- str_split_fixed(tbl$Pathway, ": ", 2)[, 2]


# reactome pathways related to the immune response
reactome_immune_pathways <- c(
  "Interleukin-6 signaling", "Platelet degranulation", "Downstream signaling events of B Cell Receptor (BCR)", "Activation of NF-kappaB in B cells", "Activation of RAS in B cells", "ISG15 antiviral mechanism", "Antiviral mechanism by IFN-stimulated genes", "Latent infection - Other responses of Mtb to phagocytosis", "Cross-presentation of particulate exogenous antigens (phagosomes)", "Antigen processing-Cross presentation", "Cross-presentation of soluble exogenous antigens (endosomes)", "Interleukin-7 signaling", "Cytokine Signaling in Immune system", "Adaptive Immune System", "Beta defensins", "Defensins", "Alpha-defensins", "ZBP1(DAI) mediated induction of type I IFNs", "IRF3 mediated activation of type 1 IFN", "HIV Infection",
  "Nef mediated downregulation of CD28 cell surface expression", "Nef mediated downregulation of MHC class I complex cell surface expression", "The role of Nef in HIV-1 replication and disease pathogenesis", "Toll Like Receptor 4 (TLR4) Cascade", "MyD88-independent TLR4 cascade", "Complement cascade", "Lectin pathway of complement activation", "Initial triggering of complement", "Terminal pathway of complement", "RNA Pol II CTD phosphorylation and interaction with CE during HIV infection", "Trafficking and processing of endosomal TLR", "Toll Like Receptor 9 (TLR9) Cascade", "Toll Like Receptor 10 (TLR10) Cascade", "Toll Like Receptor 3 (TLR3) Cascade", "Toll Like Receptor 5 (TLR5) Cascade", "Toll Like Receptor TLR1:TLR2 Cascade", "Toll Like Receptor 7/8 (TLR7/8) Cascade", "Toll Like Receptor TLR6:TLR2 Cascade", "Innate Immune System", "Influenza Infection",
  "Immune System", "Entry of Influenza Virion into Host Cell via Endocytosis", "NOD1/2 Signaling Pathway", "Toll-like Receptor Cascades", "DDX58/IFIH1-mediated induction of interferon-alpha/beta", "Signaling by TGF-beta Receptor Complex", "Classical antibody-mediated complement activation", "Alternative complement activation", "Interactions of Vpr with host cellular proteins", "Interactions of Tat with host cellular proteins", "Interactions of Rev with host cellular proteins", "APOBEC3G mediated resistance to HIV-1 infection", "Toll Like Receptor 2 (TLR2) Cascade", "STING mediated induction of host immune responses", "Cytosolic sensors of pathogen-associated DNA", "Fcgamma receptor (FCGR) dependent phagocytosis", "FCGR activation", "MHC class II antigen presentation", "Downregulation of TGF-beta receptor signaling", "TGF-beta receptor signaling activates SMADs",
  "TGF-beta receptor signaling in EMT (epithelial to mesenchymal transition)", "Ficolins bind to repetitive carbohydrate structures on the target cell surface", "DEx/H-box helicases activate type I IFN and inflammatory cytokines production", "LRR FLII-interacting protein 1 (LRRFIP1) activates type I IFN production", "Regulation of innate immune responses to cytosolic DNA", "STAT6-mediated induction of chemokines", "IRF3-mediated induction of type I IFN", "Signaling by TGF-beta Receptor Complex in Cancer", "Chemokine receptors bind chemokines", "Costimulation by the CD28 family", "CD28 co-stimulation", "CD28 dependent PI3K/Akt signaling", "CD28 dependent Vav1 pathway", "CTLA4 inhibitory signaling", "PD-1 signaling", "Interleukin-1 family signaling", "Interleukin-12 family signaling", "Interleukin-17 signaling", "Interleukin-1 processing", "Signaling by Interleukins",
  "Other interleukin signaling", "Interleukin-2 family signaling", "SUMOylation of immune response proteins", "Interleukin-3, Interleukin-5 and GM-CSF signaling", "Diseases of Immune System", "TNFR1-induced proapoptotic signaling", "Regulation of TNFR1 signaling", "TNFR1-induced NFkappaB signaling pathway", "Diseases associated with the TLR signaling cascade", "TLR3 deficiency - HSE", "MyD88 deficiency (TLR2/4)", "MyD88 deficiency (TLR5)", "IKBKG deficiency causes anhidrotic ectodermal dysplasia with immunodeficiency (EDA-ID) (via TLR)", "IRAK4 deficiency (TLR5)", "IRAK4 deficiency (TLR2/4)", "Defective SLC40A1 causes hemochromatosis 4 (HFE4) (macrophages)", "TNFR1-mediated ceramide production", "CLEC7A/inflammasome pathway", "TNFR2 non-canonical NF-kB pathway", "TNFs bind their physiological receptors",
  "TNF receptor superfamily (TNFSF) members mediating non-canonical NF-kB pathway", "Regulation of TLR by endogenous ligand", "CD22 mediated BCR regulation", "Inflammasomes", "Interleukin-6 family signaling", "Interleukin-10 signaling", "Interleukin-4 and Interleukin-13 signaling", "IL-6-type cytokine receptor ligand interactions", "Antimicrobial peptides", "Ion influx/efflux at host-pathogen interface", "Mitotic Telophase/Cytokinesis", "TNF signaling", "The NLRP1 inflammasome", "The NLRP3 inflammasome", "The AIM2 inflammasome", "The IPAF inflammasome", "Interferon gamma signaling", "Regulation of IFNG signaling", "Interleukin-20 family signaling", "InlB-mediated entry of Listeria monocytogenes into host cell",
  "Listeria monocytogenes entry into host cells", "InlA-mediated entry of Listeria monocytogenes into host cells", "RUNX1 and FOXP3 control the development of regulatory T lymphocytes (Tregs)", "RUNX1 regulates transcription of genes involved in BCR signaling", "RUNX1 regulates transcription of genes involved in interleukin signaling", "RUNX3 Regulates Immune Response and Cell Migration", "Gene and protein expression by JAK-STAT signaling after Interleukin-12 stimulation", "Interleukin-15 signaling", "OAS antiviral response", "Interleukin-35 Signalling", "Interleukin-9 signaling", "Interleukin-38 signaling", "Interleukin-37 signaling", "Interleukin-18 signaling", "TLR3-mediated TICAM1-dependent programmed cell death", "Interleukin-33 signaling", "Interleukin-2 signaling", "Interleukin-12 signaling", "Interleukin-1 signaling", "Interleukin-23 signaling",
  "Interleukin-27 signaling", "Interleukin-21 signaling", "Interferon alpha/beta signaling", "Interleukin receptor SHC signaling", "Regulation of IFNA/IFNB signaling", "Interferon Signaling", "TRIF(TICAM1)-mediated TLR4 signaling", "TRAF6-mediated induction of TAK1 complex within TLR4 complex", "Insertion of tail-anchored proteins into the endoplasmic reticulum membrane", "HCMV Infection", "Infection with Mycobacterium tuberculosis", "Modulation by Mtb of host immune system", "Leishmania infection", "ADORA2B mediated anti-inflammatory cytokines production", "Purinergic signaling in leishmaniasis infection", "CD163 mediating an anti-inflammatory response", "Anti-inflammatory response favouring Leishmania parasite infection", "FCGR3A-mediated IL10 synthesis", "Parasite infection", "FCGR3A-mediated phagocytosis",
  "Cell recruitment (pro-inflammatory response)", "SARS-CoV-1 Infection", "SARS-CoV Infections", "SARS-CoV-1 activates/modulates innate immune responses", "SARS-CoV-2 Infection", "SARS-CoV-2 activates/modulates innate and adaptive immune responses", "TRAF6 mediated IRF7 activation in TLR7/8 or 9 signaling", "TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation", "IRAK1 recruits IKK complex upon TLR7/8 or 9 stimulation", "IRAK2 mediated activation of TAK1 complex upon TLR7/8 or 9 stimulation", "Epithelial-Mesenchymal Transition (EMT) during gastrulation", "Early SARS-CoV-2 Infection Events", "Late SARS-CoV-2 Infection Events", "Regulation of Complement cascade", "Antigen processing: Ubiquitination & Proteasome degradation", "Class I MHC mediated antigen processing & presentation", "Antigen Presentation: Folding, assembly and peptide loading of class I MHC", "Antigen activates B Cell Receptor (BCR) leading to generation of second messengers", "Signaling by the B Cell Receptor (BCR)"
)

# Filter for immune-related pathways (user-defined list)
tbl <- tbl[tbl$Pathway %in% reactome_immune_pathways, ]
tbl <- tbl %>% rename(gex.gene = SYMBOL)

# Merge with correlation results
reactome.full <- tbl %>%
  left_join(med.sm, by = "gex.gene") %>%
  filter(!is.na(r2)) %>%
  rename(pname = Pathway)

# -----------------------------
# Identify significantly associated pathways
# -----------------------------
tt <- reactome.full %>%
  group_by(pname) %>%
  filter(n() > 2) %>%
  summarise(p.val = t.test(lm.slope, mu = 0)$p.value) %>%
  mutate(p.adj = p.adjust(p.val, method = "fdr"))

sig <- tt %>% filter(p.val < 0.05) %>% pull(pname)

ord <- reactome.full %>%
  group_by(pname) %>%
  summarise(md = median(lm.slope)) %>%
  arrange(-md) %>%
  pull(pname)

reactome.full$pname <- factor(reactome.full$pname, levels = ord)
r.sig <- reactome.full[reactome.full$pname %in% sig, ]

# -----------------------------
# Plot correlation by pathway
# -----------------------------
r.sig.err <- r.sig %>%
  group_by(pname) %>%
  summarise(
    mean.lm.slope = mean(lm.slope),
    mean.r2 = mean(r2),
    sd.lm.slope = sd(lm.slope),
    n = n()
  ) %>%
  mutate(
    se = sd.lm.slope / sqrt(n),
    lower = mean.lm.slope - qt(0.975, df = n - 1) * se,
    upper = mean.lm.slope + qt(0.975, df = n - 1) * se
  )

r.sig.err$color <- colorRampPalette(brewer.pal(11, "RdYlGn"))(nrow(r.sig.err))[rank(r.sig.err$mean.r2)]
r.sig.err$pname <- factor(r.sig.err$pname, levels = r.sig.err$pname[order(r.sig.err$mean.r2)])

ggplot(r.sig.err, aes(y = pname, x = mean.lm.slope)) +
  geom_jitter(data = r.sig, aes(x = lm.slope), color = "gray75", size = 0.2, height = 0.2) +
  geom_bar(stat = "identity", aes(fill = mean.r2), width = 0.6, alpha = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.4, linewidth = 0.35) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdYlGn")), name = "Adj R²") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_ipsum() +
  labs(x = "Regression coefficient", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, size = 8),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    panel.grid = element_blank()
  ) +
  xlim(-85, 25)
