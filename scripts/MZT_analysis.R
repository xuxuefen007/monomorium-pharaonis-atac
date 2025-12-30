rm(list = ls()); gc()
options(stringsAsFactors = FALSE, warn = -1)

# ==================== 0) Working directory and packages ====================
setwd("C:/Users/Xuefen Xu/Documents/RNA/SRR")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(grid)
})

# ==================== 1) Developmental stage order ====================
stage_order <- c(
  "Embryo_0-12h", "Embryo_12-24h", "Embryo_36-48h", "Embryo_60-72h",
  "Embryo_84-96h", "Embryo_108-120h", "Embryo_132-144h",
  "Embryo_156-168h", "Embryo_180-192h"
)

# ==================== 2) Locate all quant.sf files ====================
quant_paths <- list.files(
  pattern   = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

cat("Total samples found:", length(quant_paths), "\n")
if(length(quant_paths) == 0) {
  stop("No quant.sf files found. Please check the working directory or the folder structure.")
}

# ==================== 3) Extract gene expression from quant.sf ====================
# This function reads each Salmon quant.sf, extracts TPM/NumReads for a target transcript/gene ID,
# and attaches sample metadata inferred from the file path.
extract_gene_expression <- function(quant_paths, gene_id) {

  all_data <- data.frame()

  for(i in seq_along(quant_paths)) {

    path <- quant_paths[i]

    # Sample name is defined as the directory two levels above quant.sf
    sample_name <- basename(dirname(dirname(path)))

    # Extract developmental stage from the file path (e.g., Embryo_12-24h)
    developmental_stage <- ifelse(
      grepl("Embryo_[0-9]+-[0-9]+h", path),
      gsub(".*(Embryo_[0-9]+-[0-9]+h).*", "\\1", path),
      "Unknown"
    )

    # Infer strain based on sample naming convention
    strain <- ifelse(grepl("_4030", sample_name), "4030", "D03")

    # Read quant.sf
    dat <- read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "")

    # Extract the row matching the target gene/transcript ID
    gene_row <- dat[dat$Name == gene_id, , drop = FALSE]

    # If the gene/transcript is not present, set TPM/NumReads to 0
    if(nrow(gene_row) > 0) {
      all_data <- rbind(all_data, data.frame(
        Sample = sample_name,
        DevelopmentalStage = developmental_stage,
        Strain = strain,
        TPM = gene_row$TPM[1],
        NumReads = gene_row$NumReads[1],
        GenePresent = TRUE,
        stringsAsFactors = FALSE
      ))
    } else {
      all_data <- rbind(all_data, data.frame(
        Sample = sample_name,
        DevelopmentalStage = developmental_stage,
        Strain = strain,
        TPM = 0,
        NumReads = 0,
        GenePresent = FALSE,
        stringsAsFactors = FALSE
      ))
    }
  }

  return(all_data)
}

# ==================== 4) Unified maternal clearance analysis (0–12h vs 12–24h) ====================
# This function quantifies maternal transcript clearance by:
# - comparing expression between Embryo_0-12h and Embryo_12-24h,
# - computing fold-change, log2 fold-change, percent change,
# - performing a two-sample t-test,
# - computing effect size (Cohen's d).
analyze_maternal_clearance_unified <- function(gene_expression, gene_name) {

  early_data <- gene_expression %>%
    filter(DevelopmentalStage %in% c("Embryo_0-12h", "Embryo_12-24h"))

  # Require at least 2 samples per stage (total >= 4)
  if(nrow(early_data) < 4) return(NULL)

  stage_stats <- early_data %>%
    group_by(DevelopmentalStage) %>%
    summarise(
      n = n(),
      Mean_TPM = mean(TPM),
      SD_TPM = sd(TPM),
      SE_TPM = sd(TPM) / sqrt(n()),
      .groups = "drop"
    )

  mean_0_12h  <- stage_stats$Mean_TPM[stage_stats$DevelopmentalStage == "Embryo_0-12h"]
  mean_12_24h <- stage_stats$Mean_TPM[stage_stats$DevelopmentalStage == "Embryo_12-24h"]

  fold_change <- mean_12_24h / mean_0_12h
  log2_fc     <- log2(fold_change)

  # Percent change (12–24h relative to 0–12h; negative value indicates a decrease)
  percent_change <- ((mean_12_24h - mean_0_12h) / mean_0_12h) * 100

  # Two-sample t-test
  t_test_result <- t.test(TPM ~ DevelopmentalStage, data = early_data)

  # Cohen's d (positive if 0–12h > 12–24h)
  sd_0 <- stage_stats$SD_TPM[stage_stats$DevelopmentalStage == "Embryo_0-12h"]
  sd_1 <- stage_stats$SD_TPM[stage_stats$DevelopmentalStage == "Embryo_12-24h"]
  cohens_d <- (mean_0_12h - mean_12_24h) / sqrt((sd_0^2 + sd_1^2) / 2)

  out <- data.frame(
    gene_name = gene_name,
    mean_0_12h = mean_0_12h,
    mean_12_24h = mean_12_24h,
    fold_change = fold_change,
    log2_fold_change = log2_fc,
    percent_change = percent_change,
    p_value = t_test_result$p.value,
    cohens_d = cohens_d,
    n_0_12h = stage_stats$n[stage_stats$DevelopmentalStage == "Embryo_0-12h"],
    n_12_24h = stage_stats$n[stage_stats$DevelopmentalStage == "Embryo_12-24h"],
    stringsAsFactors = FALSE
  )

  return(out)
}

# ==================== 5) Target genes for MZT analysis ====================
target_genes <- data.frame(
  GeneID = c("MphaG08722.3", "MphaG09302.3", "MphaG10013.1", "MphaG03798.1",
             "MphaG08516.1", "MphaG00620.3", "MphaG04660.1", "MphaG13938.1",
             "MphaG01178.1", "MphaG01177.1"),
  GeneName = c("Trl", "Clamp", "Mad", "opa", "brk", "ewg", "hry", "vfl", "mirr", "ara"),
  stringsAsFactors = FALSE
)

cat("Starting comprehensive MZT analysis for", nrow(target_genes), "genes...\n\n")

# Containers for outputs
maternal_results <- list()
all_expression_data <- list()

# ==================== 6) Main loop: extract expression, run tests ====================
for(i in 1:nrow(target_genes)) {

  gene_id   <- target_genes$GeneID[i]
  gene_name <- target_genes$GeneName[i]

  cat("Analyzing", gene_name, "(", gene_id, ")...\n")

  # Extract expression for this gene across all samples
  gene_expression <- extract_gene_expression(quant_paths, gene_id)

  # Keep only valid embryo stages and enforce ordering
  gene_expression <- gene_expression[gene_expression$DevelopmentalStage %in% stage_order, ]
  gene_expression$DevelopmentalStage <- factor(gene_expression$DevelopmentalStage, levels = stage_order)

  all_expression_data[[gene_name]] <- gene_expression

  # Maternal clearance analysis (0–12h vs 12–24h)
  mat_res <- analyze_maternal_clearance_unified(gene_expression, gene_name)
  if(!is.null(mat_res)) maternal_results[[gene_name]] <- mat_res
}

# ==================== 7) Multiple testing correction and result annotation ====================

# ---- Maternal clearance results ----
maternal_df <- do.call(rbind, maternal_results)

# Benjamini–Hochberg FDR correction across all tested genes
maternal_df$FDR <- p.adjust(maternal_df$p_value, method = "fdr")

# Significance label
maternal_df$Clearance_Significance <- ifelse(
  maternal_df$FDR < 0.001, "***",
  ifelse(maternal_df$FDR < 0.01, "**",
         ifelse(maternal_df$FDR < 0.05, "*", "ns"))
)

# Clearance pattern classification (based on FDR and fold change thresholds)
maternal_df$Clearance_Pattern <- ifelse(
  maternal_df$FDR < 0.05 & maternal_df$fold_change < 0.5, "Strong Clearance",
  ifelse(maternal_df$FDR < 0.05 & maternal_df$fold_change < 0.7, "Moderate Clearance",
         "No Significant Clearance")
)

# ==================== 8) Output: maternal clearance ====================
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("MATERNAL CLEARANCE ANALYSIS RESULTS\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Significant maternal clearance
significant_clearance <- maternal_df %>%
  filter(FDR < 0.05 & fold_change < 0.7) %>%
  arrange(fold_change)

cat("Significantly Cleared Maternal Transcripts (FDR < 0.05 & fold_change < 0.7):\n")
print(significant_clearance %>%
        select(gene_name, mean_0_12h, mean_12_24h, percent_change,
               fold_change, log2_fold_change, FDR, Clearance_Pattern, Clearance_Significance))

# ==================== 9) Focus genes: formatted maternal clearance table ====================
genes_of_interest <- c("Mad", "ewg", "Clamp", "vfl", "opa")

focus_clearance <- maternal_df %>%
  filter(gene_name %in% genes_of_interest) %>%
  mutate(
    mean_0_12h = round(mean_0_12h, 2),
    mean_12_24h = round(mean_12_24h, 2),
    percent_change = round(percent_change, 1),
    fold_change = round(fold_change, 3),
    FDR = signif(FDR, 4),
    p_value = signif(p_value, 4),
    Significance = Clearance_Significance,
    Expression_Change = sprintf("%.2f → %.2f TPM", mean_0_12h, mean_12_24h),
    Percent_Decrease = sprintf("%.1f%%", abs(percent_change)),
    Statistical_Significance = sprintf("FDR=%g %s", FDR, Significance)
  ) %>%
  select(gene_name, Expression_Change, Percent_Decrease, fold_change, p_value, FDR, Statistical_Significance) %>%
  arrange(match(gene_name, genes_of_interest))

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("FOCUS GENES: MATERNAL CLEARANCE (FORMATTED TABLE)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")
print(focus_clearance)

# ==================== 10) save outputs ====================
# write.csv(maternal_df, file = "maternal_clearance_results.csv", row.names = FALSE)
# write.csv(focus_clearance, file = "focus_genes_maternal_clearance_table.csv", row.names = FALSE)

cat("\nAll done!\n")

