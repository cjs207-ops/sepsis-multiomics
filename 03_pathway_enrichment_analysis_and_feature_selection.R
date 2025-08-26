
# ========================
# Script 03: Pathway Enrichment, Heatmap & Feature selection
# ========================
# Objective:
# 1. Use clusterProfiler to identify enriched immune regulatory pathways 
#    from DEGs in the discovery cohort (restricted to innate immune system genes).
# 2. Visualize scaled gene expression of the most enriched pathway 
#    ("negative regulation of immune response") using a heatmap and dendrogram.

# ========================
# ========================
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(enrichplot)
library(tidyr)
library(dplyr)
library(reshape2)
library(viridis)
library(wesanderson)
library(ggdendro)
library(grid)
library(cowplot)
library(data.table)

# ========================
# Load DEG List and Gene Set
# ========================
dfl.disco <- read.csv("data/dfl.disco.csv")     # Discovery cohort DEG list

# Predefined Reactome innate immune system gene list

reactome.innate.immune.system_gene.list= c(
  "DEFA1", "CPN2", "NBEAL2", "CLEC6A", "CAPZA1", "FUCA1", "PSMB1", "NCKAP1", "PPP3CB", "DEFB4A", "ISG15", "GNLY", "CR1", "MASP2", "LILRB2",
  "LILRB3", "PTPN6", "TLR1", "DEFB131A", "RAB44", "NFATC1", "LTF", "CTSA", "CLEC4C", "SARM1", "PRTN3", "GSDME", "CALM1", "DEFB114", "TRIM21",
  "CD68", "NOD2", "PSMB9", "STK11IP", "CASP9", "PSMB4", "ELK1", "CAPN1", "FUCA2", "NCK1", "MLEC", "TICAM1", "DEFB135", "KCMF1", "PSMD13",
  "BAIAP2", "PGM1", "MANBA", "CARD11", "LPCAT1", "DNAJC5", "LTA4H", "RNF125", "FLG2", "PKM", "TRIM4", "SLC11A1", "RASGRP4", "RBSN", "HSPA1B",
  "CD247", "RNASET2", "APP", "MYO9B", "DGAT1", "BRI3", "ABCA13", "GPR84", "MEF2C", "DEFB104A", "ATP6V0E1", "ORM1", "CLEC4E", "MAP3K14", "ELMO1",
  "SERPINA3", "RNASE6", "PGLYRP3", "PIK3R2", "LAMTOR1", "TAX1BP1", "FCN1", "SIRPA", "PIGR", "CYB5R3", "PSEN1", "COMMD9", "ERP44", "SKP1", "GPI",
  "MAP3K7", "HVCN1", "MYD88", "CYBA", "WASF3", "CALML5", "LBP", "PLCG1", "PTK2", "YES1", "PELI3", "DYNC1LI1", "CDC34", "DEFA3", "ITGAL",
  "TXNDC5", "ITK", "PPIE", "BIRC3", "DDX3X", "MAP2K3", "CTSS", "CTSC", "IFIH1", "PSMD1", "ATP6V1C2", "VAV1", "ATP7A", "ADAM8", "UBA52",
  "AGL", "DEFB107B", "SEM1", "TIFA", "C6", "C1QB", "PPP2R1B", "TUBB", "PGAM1", "LIMK1", "FGL2", "MUC12", "RAB5B", "C5AR2", "CFP",
  "AOC1", "KLRD1", "ACAA1", "CFD", "SIGLEC5", "ASAH1", "LILRA3", "MUC5B", "CD59", "CREBBP", "PSMD2", "ARHGAP45", "PSMB3", "CTNNB1", "IFI16",
  "ITGB2", "CD4", "CRK", "PSMA1", "TMBIM1", "DEFB127", "CHI3L1", "ATP6V1C1", "PSMD9", "CHIT1", "MAPKAPK2", "GLB1", "GSTP1", "CD209", "TLR10",
  "MUCL1", "TLR9", "ATAD3B", "CPB2", "RAP1B", "MUC15",  "PIK3CA", "TREX1", "UBE2D3", "ACTR10", "TANK", "JUP", "PSMD4", "PSMA8",
  "EP300", "MGAM", "QPCT", "DNM1", "OSCAR", "ARL8A", "TXK", "RAB6A", "DSN1", "CD36", "ARMC8", "NFAM1", "DHX9", "PGRMC1", "NCKIPSD",
  "VAMP8", "PFKL", "S100A1", "AGA", "MAPK11", "AIM2", "ACP3", "BIRC2", "RPS6KA5", "GLA", "MAPK3", "CYFIP2", "PRG2", "CFHR2", "KCNAB2",
  "ELMO2", "ACTR2", "DEFB106B", "BCL10", "GHDC", "DNAJC13", "CASP4", "MAGT1", "PRSS3", "OSTF1", "DERA", "C2", "MYO10", "IMPDH2", "CYSTM1",
  "PLAU", "DEFB108B", "LYZ", "IFNA13", "DEFA6", "ZBP1", "NCF1", "DOK3", "EEF1A1", "TTR", "TAB3", "RIPK3", "PLD4", "IFNA2", "DHX58",
  "PSMB2", "TRAF3", "BPIFB2", "TLR8", "TUBB4B", "OTUD5", "RETN", "VAV3", "SNAP25", "FTH1", "CPNE3", "TOLLIP", "RAP1A", "B4GALT1", "DDOST",
  "SFTPA1", "STBD1", "MUC21", "WIPF3", "IRAG2", "PRG3", "CCR6", "TCIRG1", "MME", "SLPI", "SERPINA1", "CSTB", "SPTAN1", "CLU", "GMFG",
  "MMP9", "CD180", "PSMA6", "SLC44A2", "GSDMD", "UBB", "XRCC6", "NEU1", "GOLGA7", "RIGI", "SIGLEC15", "CAP1", "TNIP2", "ARSA", "SUGT1",
  "PRKCE", "PSMC6", "CFHR1", "MYH2", "PIN1", "QSOX1", "ATP6V0D2", "DUSP7", "MBL2", "SLC2A5", "COPB1", "TLR4", "ARSB", "XRCC5", "KPNB1",
  "GNS", "NCF4", "RIPK1", "BST2", "ATP6V1E2", "DSP", "ATG5", "PSME4", "VPS35L", "PRCP", "ACTR3", "C4A", "VAPA", "MOSPD2", "NHLRC3",
  "TMEM63A", "PNP", "NFKBIA", "ILF2", "PSMD5", "GRN", "CD93", "ATP6V1G2", "FCN3", "CYBB", "CD33", "VNN1", "CD3G", "GZMM", "UBE2V1",
  "CD47", "PRKACG", "MAP3K8", "KIR2DS5", "S100A7A", "PSMA4", "PPP2CA", "RAB27A", "ATF2", "ARHGAP9", "REG3G", "DUSP3", "DEFB107A", "GAB2", "S100A8",
  "DEFA5", "AGER", "BPI", "ALDOA", "USP18", "MUC16", "HTN1", "CEACAM6", "RAB3A", "DEFB1", "PSMB10", "LRG1", "FBXW11", "RNASE7", "DEFB113",
  "NCF2", "NFATC3", "LAMP2", "IDH1", "NLRX1", "IRAK4", "PSMD12", "CAT", "LAMTOR3", "CTSZ", "UBE2M", "RAF1", "DEFB136", "PSMD8", "TMEM179B",
  "DUSP6", "RELB", "HRNR", "C1S", "UBA7", "PPP2CB", "NFATC2", "FTL", "ATP6V0C", "DEGS1", "BPIFB6", "GRAP2", "IFNA5", "AZU1", "DEFB123",
  "CASP8", "GLIPR1", "DEFB124", "C3", "MS4A2", "TARM1", "AGPAT2", "LCN2", "PRKACA", "CAND1", "NLRP4", "PI3", "MAPK12", "PSMF1", "RPS6KA3",
  "TAB2", "DNM3", "CKAP4", "RAC1", "CPN1", "SIGLEC9", "TSPAN14", "RHOG", "PLA2G6", "ALPK1", "NPC2", "RAB10", "SHC1", "BPIFB4", "KLRC2",
  "C3AR1", "CARD9", "LY96", "ATP6V1D", "CDA", "TIMP2", "DNAJC3", "STAT6", "MVP", "POLR1D", "DEFB132", "PSMD10", "CANT1", "TAB1", "ATF1",
  "SERPINB12", "CNN2", "ATP8A1", "NF2", "RAB37", "CAPZA2", "POLR2E", "BPIFA1", "AP2A2", "PSMB5", "CRP", "YPEL5", "SELL", "TRAPPC1", "ATP6V0A2",
  "NRAS", "N4BP1", "DYNLL1", "CCR2", "RAP2B", "CD300A", "LAT", "CPNE1", "SIRPB1", "TNFAIP6", "SERPING1", "NAPRT", "JUN", "IFNA1", "PSMD7",
  "ACTB", "ATP6AP2", "PTPN4", "BRK1", "MNDA", "EPPIN", "TRPM2", "COTL1", "NOD1", "PSMB6", "CDC42", "LCK", "C5AR1", "VRK3", "SERPINB1",
  "FCGR3B", "DYNLT1", "MEFV", "APOB", "ATP6V0D1", "PPP3R1", "SIGIRR", "BIN2", "AP1M1", "GSN", "MS4A3", "DEFB103B", "SDCBP", "TLR7", "CD46",
  "BST1", "PRKCQ", "DEFB134", "TMC6", "C8A", "ORMDL3", "BPIFA2", "COMMD3", "METTL7A", "NFASC", "FCAR", "ANPEP", "CD300E", "STING1", "CMTM6",
  "FCER1A", "NLRC3", "SFTPD", "SNAP23", "FYN", "SVIP", "HUWE1", "CLEC10A", "TXNIP", "SURF4", "CRCP", "PPIA", "TOMM70", "PRSS2", "PLEKHO2",
  "B2M", "NLRC5", "HSP90B1", "IFNA6", "PGLYRP1", "PLA2G2A", "MAN2B1", "FOLR3", "TYROBP", "PAK1", "ATP6V1E1", "LCP2", "ORM2", "PIK3CB", "PYGB",
  "DEFB130B", "DOCK1", "HRAS", "CXCR2", "PGM2", "CD55", "TKFC", "PIK3R1", "PGLYRP2", "CRISP3", "ADA2", "IKBKB", "ATP6V0B", "SYNGR1", "CFI",
  "RELA", "FGB", "HK3", "TXN", "CD300LB", "S100B", "ACTG1", "IFNA8", "MAPK1", "TRAF2", "ANXA2", "KRAS", "RNF135", "PSMC3", "CHRNB4",
  "KIR2DS1", "TRIM56", "NCSTN", "LY86", "RHOF", "IFNA10", "PDPK1", "MRE11", "PSTPIP1", "PELI1", "POLR2H", "LAIR1", "RAB31", "DEFA4", "POLR3F",
  "ARPC5", "C5", "FGR", "OLR1", "BTK", "SAA1", "NKIRAS2", "MAP2K1", "ADGRE3", "IRAK3", "CLEC4D", "APAF1", "DEFB110", "PRDX4", "ACTR1B",
  "FABP5", "FPR1", "VCL", "COLEC10", "CPPED1", "MAPK14", "TLR6", "ITGAM", "IFNA14", "CAB39", "ARPC4", "LGALS3", "KIR2DS2", "APRT", "KLRK1",
  "IKBKG", "C7", "NKIRAS1", "TBK1", "SYK", "P2RX7", "RAB5C", "POLR3E", "SIGLEC14", "RAC2", "DSG1", "HCK", "SIKE1", "HEXB", "IST1",
  "ROCK1", "HMGB1", "GYG1", "C8G", "GM2A", "S100A11", "DEFB4B", "CUL1", "IRAK2", "MUC5AC", "CLEC5A", "PSMB11", "PSMA2", "AHSG", "RIPK2",
  "SRC", "DEFB128", "ITGAX", "S100A12", "LOC102725035", "CD81", "RNASE2", "TLR2", "PAK3", "PSMB8", "ACLY", "PECAM1", "COLEC11", "ICAM2", "PPP2R5D",
  "AMPD3", "PYGL", "SCAMP1", "LPO", "ATG12", "APEH", "PANX1", "NFKB1", "CGAS", "MPO", "C1orf35", "MUC7", "FOS", "CTSL", "UNC93B1",
  "VAT1", "HLA-E", "PYCARD", "PELI2", "MAPK13", "PTPN11", "DEFB126", "CLEC4A", "EPPIN-WFDC6", "MUC19", "CD44", "PRKDC", "REG3A", "FCGR2A", "HMOX1",
  "TNFAIP3", "ATP6V1G1", "POLR3D", "DEFB115", "KIR2DS4", "S100P", "IFNB1", "IRF3", "MUC20", "MAVS", "DTX4", "DEFB118", "TOM1", "WASF2", "MMP8",
  "MASP1", "UBR4", "CD58", "PA2G4", "C6orf120", "ART1", "PPBP", "PAFAH1B2", "AAMP", "ELANE", "SFTPA2", "CYLD", "CCT2", "LEAP2", "CLEC7A",
  "PLAUR", "UBA3", "MCEMP1", "RPS6KA2", "HTN3", "GCA", "HSP90AB1", "DSC1", "C1QC", "CXCR1", "SLC15A4", "ATP6V1H", "AHCYL1", "CALM3", "MIF",
  "NOS1", "TREM2", "USP14", "MAP3K1", "DEFB116", "FCN2", "PSAP", "CFHR4", "TMED7-TICAM2", "TBC1D10C", "IRF7", "RNF216", "RAB7A", "NFKB2", "PSMC1",
  "DEFA1B", "TRIM32", "DEFB121", "WIPF2", "FGA", "ATP11A", "CFL1", "PIK3R4", "EEF2", "PDZD11", "CALM2", "PTX3", "ALDH3B1", "IFNA17", "ATG7",
  "DEFB125", "CFH", "PDAP1", "FCER1G", "TLR5", "TRAF6", "DEFB129", "WIPF1", "CD63", "GAA", "CNPY3", "CFB", "ADGRE5", "UBE2L6", "STOM",
  "IGF2R", "CASP2", "ATOX1", "MGST1", "LOC107987462", "NIT2", "DCD", "IKBIP", "RAB14", "HPSE", "ICAM3", "RPS6KA1", "ATP11B", "LAMP1", "FCGR1A",
  "DDX41", "FGG", "C1QA", "CTSB", "RAB3D", "NME2", "PLD1", "HLA-B", "NLRC4", "SERPINB3", "POLR3G", "CTSH", "C4B_2", "ITPR2", "MUC3A",
  "FAF2", "ARPC1A", "DHX36", "RASGRP1", "LRRFIP1", "PTPRB", "KRT1", "TRIM25", "CEACAM1", "SEMG1", "HBB", "DUSP4", "PLD2", "ITLN1", "PKP1",
  "RAB24", "NDUFC2", "PGLYRP4", "SNAP29", "FCGR3A", "A1BG", "PPP2R1A", "CREG1", "PSME3", "ABI1", "POLR3B", "LAT2", "HMOX2", "GALNS", "ATP6V1B1",
  "CD19", "TNFRSF1B", "GGH", "PTPRJ", "LRRC14", "POLR3K", "MEF2A", "ANO6", "RAB18", "RASGRP2", "C4B", "UBC", "PLCG2", "PIK3C3", "PRKCSH",
  "PRKCD", "STK10", "RAP2C", "FPR2", "RNASE8", "ARPC2", "PLAC8", "HGSNAT", "SOS1", "HSPA6", "NOS2", "ENPP4", "MAPK9", "IFNA16", "IL1B",
  "RAB9B", "CEACAM8", "VCP", "HSPA8", "RHOA", "MYH9", "PSMA3", "MAPK8", "KLRC4-KLRK1", "ITPR1", "MAP2K7", "C9", "TCN1", "CRISPLD2", "CYFIP1",
  "TIRAP", "ATP6V1F", "FRK", "FRMPD3", "PROS1", "FADD", "PCBP2", "SLC27A2", "ATP6V0A1", "PSMA5", "DBNL", "CFHR5", "MUC6", "UBE2D1", "SERPINB6",
  "VAV2", "POLR3GL", "PSMC2", "CTSV", "ATP8B4", "DYNC1H1", "TMEM30A", "BTRC", "TEC", "SLCO4C1", "LAMTOR2", "NCR2", "PSME1", "PSMD14", "MUC1",
  "IFNA7", "IQGAP1", "POLR2L", "GUSB", "S100A9", "OLFM4", "MYO5A", "PAK2", "MAP2K4", "CLEC12A", "DPP7", "HLA-A", "DOCK2", "NOS3", "ITCH",
  "DEFB103A", "NLRP1", "ECSIT", "ITGAV", "NME1-NME2", "POLR1C", "ATP6V1B2", "CFHR3", "DNM2", "IQGAP2", "ABI2", "CR2", "MUC4", "CHUK", "ALOX5",
  "MUC17", "CCT8", "DEFB105B", "PSMD11", "WAS", "MAP2K6", "ALDOC", "PTPRN2", "CTSD", "PSMC5", "POLR3A", "WASF1", "MAPKAPK3", "PSME2", "CDK13",
  "P2RX1", "HP", "PSMC4", "CD53", "GDI2", "HSP90AA1", "TREM1", "C1R", "UNC13D", "CSNK2B", "DEFB130A", "ALAD",  "PLPP5", "NLRP3",
  "PLPP4", "PTAFR", "PSMA7", "PPP3CA", "PDXK", "ABL1", "GRB2", "BCL2L1", "POLR3H", "ITPR3", "PSMD6", "SRP14", "CHGA", "TLR3", "DEFB119",
  "ARPC3", "EEA1", "MUC13", "LYN", "UBE2D2", "PEDS1-UBE2V1", "POLR3C", "WASL", "NFKBIB", "POLR2F", "CXCL1", "CASP10", "IMPDH1", "UBE2N", "CASP1",
  "NCKAP1L", "SERPINB10", "CD14", "C4BPA", "HSPA1A", "F2", "RPS27A", "IFNA4", "ATP6V0A4", "SLC2A3", "MYO1C", "BPIFB1", "PTPRC", "PTGES2", 
  "CAMP", "KIR3DS1", "PSMB7", "CEP290", "CREB1", "VTN", "UBE2K", "ARG1", "MAPK10", "EPX", "BCL2", "DEFB105A", "CTSG", "LGMN", "IKBKE",
  "DEFB106A", "PRDX6", "PLD3", "SOCS1", "HEBP2", "ADAM10", "IRAK1", "CEACAM3", "TICAM2", "ARPC1B", "C8B", "TP53", "ATP6V1G3", "PADI2", "RAB4B",
  "DIAPH1", "HERC5", "S100A7", "PSMD3", "ATP6V1A", "DEFB104B", "IFNA21", "CTSK", "MMP25", "PRKACB", "ATP6V0E2", "HLA-C", "MAPK7", "POLR2K", "CD177",
  "CST3", "MALT1", "RNASE3", "DNASE1L1", "TMEM176A", "TNFRSF12A", "SNAI2", "HGF", "UBR2", "IFNGR1", "NEDD4L", "PPP2R5B", "ADGRF5", "LNX1", "GRAMD4",
  "UBE2A", "RAPGEF4",  "JAG1", "WFDC2", "BMX", "CCL22", "CD276", "NDRG1", "FGL1", "ICAM4", "CCL24", "TMEM176B", "AMBP", 
  "FANCE", "TFEB", "VEGFA", "RNF130", "BTNL8", "IL5", "LIFR", "BCL6", "IL1R2", "IL1R1", "IL1RL1", "IL18R1", "UBE2B",
  "TCP1", "PTK2B", "FLT3", "TWIST1", "CNTFR", "WFDC3",  "CDH26", "BMP4", "GPR108", "TNFSF14", "IL22", "DLL4", "SEC14L1", "GDF15",
  "SMPDL3B", "THEMIS2", "IL13RA1", "RARA", "SPINK5", "VTCN1", "GPNMB", "IL10", "IL36B", "DAB2IP", "RNF144B", "HECW2", "ANXA3", "SOD1", "LYST",
  "GPR17", "HSPD1", "SLIT2", "HERC4", "TRIM48", "ADAM17", "RAB3C", "MERTK", "FBXL18", "VSIG4", "FZD7", "CXCL13", "DYNC1I1", "APOA2", "UBE2J2",
  "FBXL13", "CD300LG", "IL22RA2", "RAET1E", "SHH", "NFIL3", "SVEP1", "LRFN5", "INPPL1", "HPRT1", "KCTD6", "ONECUT1", "SCG2", "KLHL6", "TRIB1",
  "CLCF1", "SFN", "CACNB4", "PLCB1", "WFDC10B", "PRKG1"
)


# ========================
# Subset DEGs and Run GO Enrichment
# ========================
ig <- dfl.disco[dfl.disco$pvalue < 0.05 & 
                  dfl.disco$gene_symbol %in% reactome.innate.immune.system_gene.list, ]
immune_genes <- ig$gene_symbol

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(immune_genes, fromType = "SYMBOL", 
                    toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO Biological Process enrichment
ego <- enrichGO(gene          = gene_entrez$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# ========================
# Plot Top 20 Enriched GO Terms
# ========================
ego_result <- ego@result
ego_result$GeneRatio <- sapply(strsplit(ego_result$GeneRatio, "/"), 
                               function(x) as.numeric(x[1]) / as.numeric(x[2]))
ego_top20 <- ego_result[order(ego_result$p.adjust), ][1:20, ]

ggplot(ego_top20, aes(x = -log10(p.adjust), y = reorder(Description, -p.adjust))) +
  geom_col(fill = "gray45") +
  geom_text(aes(label = ID), hjust = 1, nudge_x = -0.5, size = 3, color = "white") +
  xlab("-log10 (q-value)") + ylab("") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 13),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

# ========================
# Load Expression and Phenodata
# ========================
datax.disco <- as.data.frame(fread("data/datax.disco.csv", sep = "auto"))
phenodata   <- read.csv("data/phenodata.csv", stringsAsFactors = FALSE)

datax.disco$V1 <- NULL
phenodata$X    <- NULL
rownames(datax.disco) <- datax.disco$serial
rownames(phenodata)   <- phenodata$serial
is=base::intersect(rownames(phenodata),rownames(datax.disco))
phenodata=phenodata[is,]
datax.disco=datax.disco[is,]
stopifnot(identical(rownames(phenodata), rownames(datax.disco)))

# ========================
# Heatmap of selected genes in "Negative Regulation of Immune Response"
# ========================
path1 <- c("PGLYRP1", "LTF", "ARG1", "GPNMB", "BCL6", "HMOX1", "IL10", 
           "RARA", "PLCB1", "GPR17", "STAT6", "TMEM176A", "TMEM176B", 
           "VSIG4", "IRAK3", "MERTK", "TRIB1", "DUSP3", "PADI2")

# Combine expression and metadata
dtx <- data.frame(
  outcome     = datax.disco$outcome,
  serial      = datax.disco$serial,
  date_death  = as.Date(phenodata$date_of_death, "%Y-%m-%d"),
  date_adm    = as.Date(phenodata$date_admitted, "%Y-%m-%d"),
  ttd         = NA,
  datax.disco[, path1]
)
dtx$ttd <- as.numeric(dtx$date_death - dtx$date_adm)

# Melt and scale
ord <- dtx %>% arrange(outcome) %>% pull(serial)
dtx$serial <- factor(dtx$serial, levels = ord)

dtx.m <- reshape2::melt(dtx, id.vars = c("outcome", "serial", "date_adm", "date_death", "ttd"))
dtx.m <- dtx.m %>% 
  group_by(variable) %>% 
  mutate(value.s = scales::rescale(1.2^value, c(0, 100000)))  # for clustering

# ========================
# Hierarchical Clustering
# ========================
wide_mat <- dtx.m %>% 
  select(serial, variable, value.s) %>%
  pivot_wider(names_from = serial, values_from = value.s, values_fill = 0) %>%
  column_to_rownames("variable")

dist_mat <- dist(wide_mat, method = "euclidean")
clust <- hclust(dist_mat, method = "ward.D")
protein_order <- clust$labels[clust$order]

# Update factor levels for plotting
dtx.m$variable <- factor(dtx.m$variable, levels = protein_order)
dtx.m$serial   <- factor(as.character(dtx.m$serial), levels = ord)

# ========================
# Heatmap Construction
# ========================
# Final scaling for heatmap
dtx.m <- dtx.m %>%
  group_by(variable) %>%
  mutate(value.s = scales::rescale(value, c(0, 100)))

# Color palettes
pal1 <- wes_palette("Zissou1", 100000, type = "continuous")

# Dendrogram
dendro_data <- as.dendrogram(clust)
dend_data <- dendro_data(dendro_data)

segment_data <- with(segment(dend_data), 
                     data.frame(x = y, y = x, xend = yend, yend = xend))
gene_pos_table <- with(dend_data$labels, 
                       data.frame(y_center = x, gene = as.character(label), height = 1))
gene_axis_limits <- with(gene_pos_table, 
                         c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))) + c(-0.1, 0.1)

p1 <- ggplot(segment_data) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_y_continuous(breaks = gene_pos_table$y_center, labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, expand = c(0, 0)) +
  theme_classic() +
  theme_void()

# Main heatmap
heatmap_gg <- ggplot(dtx.m, aes(x = as.numeric(serial), y = variable, fill = value.s)) +
  scale_fill_viridis(option = "magma", labels = c("min", "", "", "", "max")) +
  geom_tile(color = "gray60", size = 0.08) +
  ylab("GO:0002683") + xlab("Paediatric inpatients") +
  theme_classic(base_family = "Helvetica") +
  theme(
    legend.position = "left",
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    plot.background = element_rect(fill = "white")
  )

# Combine plots
hm <- plot_grid(heatmap_gg, p1, align = "h", rel_widths = c(1, 0.2))

hm



# ========================
# Script 05: Random Forest Analysis
# ========================
# Objective:
# - Train a Random Forest model to predict mortality using expression of 
#   genes from the most enriched GO pathway in the discovery cohort.
# - Evaluate performance via AUROC and extract top predictors.

# ========================
# Load Required Libraries
# ========================
library(randomForest)
library(pROC)
library(dplyr)
library(ggplot2)
library(readr)
library(stringi)
library(stringr)
library(data.table)

# ========================
# Load DEG List and Expression Data
# ========================
dfl.disco <- read.csv("data/dfl.disco.csv")  # DEG table from discovery cohort

# Extract gene list from top enriched pathway (assumes `ego` object is in environment)
sel <- stringi::stri_remove_empty(str_split_fixed(ego@result[1, 8], "/", 300))

# Filter DEGs for this pathway and order by log2FC
sel <- dfl.disco %>%
  filter(gene_symbol %in% sel, pvalue < 0.05) %>%
  arrange(log2foldchange) %>%
  pull(gene_symbol)

# Load gene expression matrix
datax.disco <- as.data.frame(fread("data/datax.disco.csv", sep = "auto"))
datax.disco$V1 <- NULL

# Subset to genes of interest (and keep outcome)
data_rf <- data.frame(outcome = datax.disco$outcome,
                      datax.disco[, colnames(datax.disco) %in% sel[1:35]])

# ========================
# Random Forest Classifier
# ========================

# Format outcome as binary factor
data_rf$outcome <- as.factor(data_rf$outcome)  # 0 = alive, 1 = dead

# Separate predictors and labels
predictors <- data_rf[, !(names(data_rf) %in% c("outcome", "id"))]
outcome <- data_rf$outcome

# Fit model
set.seed(126573)
rf_model <- randomForest(x = predictors, y = outcome, 
                         importance = TRUE, ntree = 10000)

# ========================
# Model Evaluation
# ========================

# Predict probability for class 'dead'
rf_probs <- predict(rf_model, type = "prob")[, "non-survivor"]

# ROC curve and AUROC
conflicted::conflicts_prefer(pROC::roc)
roc_obj <- roc(outcome, rf_probs)

# Plot ROC
conflicted::conflicts_prefer(pROC::auc)
plot(roc_obj, col = "#1f77b4", 
     main = paste("AUROC =", round(auc(roc_obj), 3)))

# Get optimal threshold (Youden’s index)
opt_thresh <- coords(roc_obj, "best", best.method = "youden",
                     ret = c("threshold", "sensitivity", "specificity"))
print(opt_thresh)

# ggplot-style ROC curve
roc_df <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(
    title = paste("AUROC =", round(auc(roc_obj), 3)),
    x = "1 - Specificity",
    y = "Sensitivity"
  )

# ========================
# Variable Importance Plot
# ========================

importance_df <- as.data.frame(importance(rf_model, type = 1))  # MeanDecreaseAccuracy
importance_df$gene <- rownames(importance_df)
colnames(importance_df)[1] <- "MeanDecreaseAccuracy"

# Order genes by importance
importance_df <- importance_df %>%
  arrange(MeanDecreaseAccuracy) %>%
  mutate(gene = factor(gene, levels = gene))

# Final ggplot2 importance plot
ggplot(importance_df, aes(x = MeanDecreaseAccuracy, y = gene)) +
  geom_point(color = "red", size = 5, alpha = 0.5) +
  geom_point(pch = 21, size = 5) +
  theme_minimal() +
  labs(
    title = "Top Predictors of Mortality",
    x = "Mean decrease in accuracy",
    y = NULL
  ) +
  scale_x_continuous(breaks = NULL) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  ) +
  annotate("text", x = 0.6, y = 0.6, 
           label = paste("Specificity:", round(opt_thresh[1, 3] * 100, 1), "%"),
           hjust = -0.5, vjust = -11, size = 4) +
  annotate("text", x = 0.6, y = 0.6, 
           label = paste("Sensitivity:", round(opt_thresh[1, 2] * 100, 1), "%"),
           hjust = -0.5, vjust = -13, size = 4)
