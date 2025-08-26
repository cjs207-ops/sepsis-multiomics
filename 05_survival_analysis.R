
# ========================
# Script 05: Survival Analysis
# ========================
# Objectives:
# 1. Generate Kaplan–Meier survival curves comparing high vs. low IL10 expression.
#    - Full validation cohort
#    - Subgroup: children with MUAC ≥ 11.5 cm (non-malnourished)
# 2. Fit Cox proportional hazards models for IL10 and ETS1 (adjusted for age).
# 3. Visualize hazard ratios with confidence intervals.

# ========================
# Load Required Libraries
# ========================
library(survival)
library(survminer)
library(readxl)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(stringr)

# ========================
# Load Expression Data
# ========================
chain.read1 <- read_excel("data/Read1_host transcriptome.xlsx")
chain.countdata <- chain.read1
rownames(chain.countdata) <- chain.countdata$record_id
chain.countdata <- chain.countdata[ , -1:-6]
colnames(chain.countdata) <- str_split_fixed(colnames(chain.countdata), "\\.", 2)[,1]
chain.countdata <- t(chain.countdata)

# ========================
# Load Metadata
# ========================
data.m <- read.csv("data/data.m.csv", stringsAsFactors = FALSE)
data.m$X <- NULL

# ========================
# Survival Dataset Prep
# ========================
gene_of_interest <- "IL10"
dta.plt <- subset(data.m, gene == gene_of_interest)
dta.plt$muac_low <- ifelse(dta.plt$muac_enrol < 11.5, "muac < 11.5", "muac ≥ 11.5")

surv.set <- dta.plt %>%
  filter(outcome2 != "Community") %>%
  transmute(
    outcome2,
    date_adm = as.Date(date_adm),
    date_last = as.Date(sc_lastvitalstat_date),
    time.to.event = as.numeric(as.Date(sc_lastvitalstat_date) - as.Date(date_adm)),
    death = ifelse(outcome2 == "alive", 0, 1),
    log.expression = log10(expression + 0.01),
    IL10 = ifelse(log10(expression + 0.01) < median(log10(expression + 0.01)), "Low", "High"),
    muac_low,
    age_adm
  )

# Subsets
surv.set.normnut <- subset(surv.set, muac_low == "muac ≥ 11.5")
surv.set.malnourished <- subset(surv.set, muac_low == "muac < 11.5")

# ========================
# Fit Kaplan–Meier Curves
# ========================
fit_all <- survfit(Surv(time.to.event, death) ~ IL10, data = surv.set)
fit_normnut <- survfit(Surv(time.to.event, death) ~ IL10, data = surv.set.normnut)

# ========================
# Plot Kaplan–Meier: Full Cohort
# ========================
ggsurvplot(
  fit_all,
  data = surv.set,
  risk.table = TRUE,
  conf.int = FALSE,
  risk.table.pos = "in",
  break.time.by = 60,
  xlab = "Time from admission (days)",
  font.tickslab = c(15, color = "black"),
  tables.theme = theme_ipsum(
    grid=F,
    plot_title_face="plain",
    plot_margin = ggplot2::margin(1, 1, 5, 1),
    plot_title_size = 10
  ),
  palette = c("black", "darkred"),
  ggtheme = theme_ipsum(axis = TRUE, grid = FALSE, axis_title_size = 14, axis_text_size = 14)
)

# ========================
# Plot Kaplan–Meier: Non-malnourished
# ========================
ggsurvplot(
  fit_normnut,
  data = surv.set.normnut,
  risk.table = TRUE,
  conf.int = FALSE,
  risk.table.pos = "in",
  break.time.by = 60,
  xlab = "Time from admission (days)",
  font.tickslab = c(15, color = "black"),
  tables.theme = theme_ipsum(
    grid=F,
    plot_title_face="plain",
    plot_margin = ggplot2::margin(1, 1, 5, 1),
    plot_title_size = 10
  ),
  palette = c("black", "darkred"),
  ggtheme = theme_ipsum(axis = TRUE, grid = FALSE, axis_title_size = 14, axis_text_size = 14)
)

# ========================
# Cox Models
# ========================
cox_all <- coxph(Surv(time.to.event, death) ~ IL10 + muac_low + age_adm, data = surv.set)
cox_normnut <- coxph(Surv(time.to.event, death) ~ IL10 + age_adm, data = surv.set.normnut)

summary(cox_all)
summary(cox_normnut)

# ========================
# Cox Models: Gene-wise (IL10, ETS1)
# ========================

gene.list <- c("IL10", "ETS1")

# Unified survival dataset -> choose this or the one below
surv.set <- data.m %>%
  filter(outcome2 != "Community") %>%
  mutate(
    muac_low = ifelse(muac_enrol < 11.5, "muac < 11.5", "muac ≥ 11.5"),
    time.to.event = as.numeric(as.Date(sc_lastvitalstat_date) - as.Date(date_adm)),
    death = ifelse(outcome2 == "alive", 0, 1),
    date_adm = as.Date(date_adm),
    sc_lastvitalstat_date = as.Date(sc_lastvitalstat_date),
    log.expression = log10(expression + 0.01)
  )

# Non-malnourished survival dataset
surv.set <- data.m %>%
  filter(outcome2 != "Community" & muac_enrol>11) %>%
  mutate(
    muac_low = ifelse(muac_enrol < 11.5, "muac < 11.5", "muac ≥ 11.5"),
    time.to.event = as.numeric(as.Date(sc_lastvitalstat_date) - as.Date(date_adm)),
    death = ifelse(outcome2 == "alive", 0, 1),
    date_adm = as.Date(date_adm),
    sc_lastvitalstat_date = as.Date(sc_lastvitalstat_date),
    log.expression = log10(expression + 0.01)
  )

# Function to compute HR for each gene
hr.calc <- function(gene_symbol) {
  cat("\r", "\033[K", paste0("Processing ", gene_symbol, "..."))
  df <- subset(surv.set, gene == gene_symbol)
  df$high.low <- ifelse(df$log.expression < median(df$log.expression), "high", "low")
  
  # Skip if insufficient data
  if (length(unique(df$high.low)) < 2 || any(table(df$high.low) < 6)) return(NA)
  
  # Fit Cox model
  fit <- withCallingHandlers(
    coxph(Surv(time.to.event, death) ~ age_adm + high.low, data = df),
    warning = function(w) {
      if (grepl("coefficient may be infinite", w$message) |
          grepl("did not converge", w$message)) invokeRestart("muffleWarning")
    }
  )
  
  # Return HR and CI
  data.frame(
    gene = gene_symbol,
    hr = exp(coef(fit))[2],
    lower.ci = exp(confint(fit))[2, 1],
    upper.ci = exp(confint(fit))[2, 2],
    p.value = summary(fit)$coefficients[2, 5]
  )
}

# Run for gene list
res <- do.call(rbind, lapply(gene.list, hr.calc))
res$adj.p.value <- p.adjust(res$p.value, method = "fdr")
res <- na.omit(res)
res <- res[!is.infinite(res$upper.ci), ]
res$sig <- res$adj.p.value < 0.05

# Subset to IL10 and ETS1
il10_regs <- c("ETS1", "IL10")  # Order for plotting
res_sig <- subset(res, gene %in% il10_regs)
res_sig <- res_sig[match(il10_regs, res_sig$gene), ]
res_sig$gene <- factor(res_sig$gene, levels = il10_regs)

# Labels
res_sig$lab <- c(
  paste0("HR: ", round(res_sig$hr[1], 1), "\n(95% CI: ", round(res_sig$lower.ci[1], 1), "-", round(res_sig$upper.ci[1], 1), ")"),
  paste0("HR: ", round(res_sig$hr[2], 1), "\n(95% CI: ", round(res_sig$lower.ci[2], 1), "-", round(res_sig$upper.ci[2], 1), ")")
)

# ========================
# Plot HRs (IL10 & ETS1)
# ========================
ggplot(res_sig, aes(x = gene, y = hr)) +
  geom_point(color = "red", size = 4) +
  geom_text(aes(label = lab), vjust = -0.6, size = 3) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.01) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans = 'log10') +
  coord_flip() +
  theme_classic() +
  xlab("Gene") +
  ylab("Hazard of Death (log scale)") +
  theme_ipsum() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(color = "gray50")
  )

