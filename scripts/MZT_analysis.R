rm(list = ls()); gc()
options(stringsAsFactors = FALSE, warn = -1)

# ==================== 0) Working directory ====================
setwd("C:/Users/Xuefen Xu/Documents/RNA/SRR")

# ==================== 0.1) Packages ====================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(DESeq2)
  library(ggplot2)
})

# ==================== 1) Developmental stage order ====================
stage_order <- c(
  "Embryo_0-12h", "Embryo_12-24h", "Embryo_36-48h", "Embryo_60-72h",
  "Embryo_84-96h", "Embryo_108-120h", "Embryo_132-144h",
  "Embryo_156-168h", "Embryo_180-192h"
)

BASELINE_STAGE <- "Embryo_0-12h"
MZT_STAGE      <- "Embryo_12-24h"
POST_START     <- "Embryo_36-48h"

OUTDIR <- "MZT_RESULTS_OPTIMIZED"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ==================== 2) Locate all quant.sf files ====================
quant_paths <- list.files(pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)
cat("Total samples found:", length(quant_paths), "\n")
if (length(quant_paths) == 0) stop("No quant.sf files found. Check folder structure.")

# ==================== 3) Build sample_info (Sample / Stage / Strain / Path) ====================
sample_info <- tibble(quant_path = quant_paths) %>%
  mutate(
    Sample = basename(dirname(dirname(quant_path))),  # two levels above quant.sf (as you used)
    Stage  = ifelse(
      grepl("Embryo_[0-9]+-[0-9]+h", quant_path),
      str_extract(quant_path, "Embryo_[0-9]+-[0-9]+h"),
      NA_character_
    ),
    Strain = ifelse(grepl("_4030", Sample), "4030", "D03")
  )

if (any(is.na(sample_info$Stage))) {
  bad <- sample_info$quant_path[is.na(sample_info$Stage)][1:min(10, sum(is.na(sample_info$Stage)))]
  stop(paste0("[ERROR] Cannot parse Stage from some paths. Examples:\n", paste(bad, collapse="\n")))
}

sample_info$Stage <- factor(sample_info$Stage, levels = stage_order)
if (any(is.na(sample_info$Stage))) {
  bad2 <- unique(as.character(sample_info$Stage[is.na(sample_info$Stage)]))
  stop("[ERROR] Some parsed stages are not in stage_order. Please update stage_order to match your names.")
}

write.csv(sample_info, file.path(OUTDIR, "sample_info.csv"), row.names = FALSE)

# ==================== 4) Read quant.sf -> build gene-level matrices (NumReads + TPM) ====================
read_quant_min <- function(path) {
  dat <- fread(path)
  dat[, .(Name, TPM, NumReads)]
}

cat("Reading quant.sf files and building matrices...\n")
quant_list <- lapply(sample_info$quant_path, read_quant_min)
names(quant_list) <- sample_info$Sample

# union of gene IDs across samples
all_gene_ids <- Reduce(union, lapply(quant_list, function(x) x$Name))

# build counts matrix
count_mat <- sapply(quant_list, function(df) {
  v <- df$NumReads
  names(v) <- df$Name
  out <- v[all_gene_ids]
  out[is.na(out)] <- 0
  as.numeric(out)
})
rownames(count_mat) <- all_gene_ids
colnames(count_mat) <- names(quant_list)

# build TPM matrix
tpm_mat <- sapply(quant_list, function(df) {
  v <- df$TPM
  names(v) <- df$Name
  out <- v[all_gene_ids]
  out[is.na(out)] <- 0
  as.numeric(out)
})
rownames(tpm_mat) <- all_gene_ids
colnames(tpm_mat) <- names(quant_list)

write.csv(count_mat, file.path(OUTDIR, "gene_counts_NumReads.csv"))
write.csv(tpm_mat,   file.path(OUTDIR, "gene_TPM.csv"))

# ==================== 5) Target genes (GeneID -> GeneName) ====================
target_genes <- data.frame(
  GeneID = c(
    "MphaG08722.3", # Trl
    "MphaG09302.3", # Clamp
    "MphaG10013.1", # Mad
    "MphaG08516.1", # opa
    "MphaG00620.3", # ewg
    "MphaG13938.1", # zld
    "MphaG01178.1", # mirr
    "MphaG01177.1"  # ara
  ),
  GeneName = c("Trl", "Clamp", "Mad", "opa", "ewg", "zld", "mirr", "ara"),
  stringsAsFactors = FALSE
)

# ==================== 6) DESeq2 for 12–24h vs 0–12h (count-based, genome-wide) ====================
cat("Running DESeq2 (counts) for 12–24h vs 0–12h...\n")

coldata <- sample_info %>%
  select(Sample, Stage, Strain) %>%
  as.data.frame()
rownames(coldata) <- coldata$Sample

# ensure sample order aligned
count_mat2 <- count_mat[, rownames(coldata), drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat2),
  colData   = coldata,
  design    = ~ Strain + Stage
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("Stage", MZT_STAGE, BASELINE_STAGE)) %>% as.data.frame()
res$GeneID <- rownames(res)
res$pct_change <- (2^res$log2FoldChange - 1) * 100

res_out <- res %>%
  select(GeneID, log2FoldChange, pct_change, pvalue, padj, baseMean) %>%
  arrange(padj)

write.csv(res_out, file.path(OUTDIR, "DESeq2_12_24_vs_0_12_full.csv"), row.names = FALSE)

# Focus summary for target genes
res_focus <- res_out %>%
  filter(GeneID %in% target_genes$GeneID) %>%
  left_join(target_genes, by = "GeneID") %>%
  transmute(
    GeneName, GeneID,
    pct_change = round(pct_change, 1),
    log2FoldChange = log2FoldChange,
    FDR = padj
  ) %>%
  arrange(match(GeneName, target_genes$GeneName))

write.csv(res_focus, file.path(OUTDIR, "DESeq2_focus_genes_summary.csv"), row.names = FALSE)
print(res_focus)

# Write text-ready lines
fmt_sci <- function(x, digits=2) ifelse(is.na(x), NA, format(x, scientific=TRUE, digits=digits))
text_lines <- res_focus %>%
  mutate(
    pct_txt = ifelse(pct_change >= 0, paste0("+", pct_change, "%"), paste0(pct_change, "%")),
    fdr_txt = fmt_sci(FDR, digits = 2),
    line = paste0(GeneName, " ", pct_txt, " (FDR = ", fdr_txt, ")")
  ) %>% pull(line)

writeLines(text_lines, con = file.path(OUTDIR, "WRITE_READY_clearance_lines.txt"))

# ==================== 7) TPM time-course summaries for GOI + plotting ====================
cat("Summarizing TPM time-course and plotting...\n")

tpm_long <- as.data.frame(tpm_mat)
tpm_long$GeneID <- rownames(tpm_long)

tpm_long <- tpm_long %>%
  pivot_longer(-GeneID, names_to = "Sample", values_to = "TPM") %>%
  left_join(sample_info %>% select(Sample, Stage, Strain), by = c("Sample" = "Sample"))

stage_mean <- tpm_long %>%
  group_by(GeneID, Stage) %>%
  summarise(mean_TPM = mean(TPM, na.rm=TRUE), .groups="drop")

write.csv(stage_mean, file.path(OUTDIR, "stage_mean_TPM_all_genes.csv"), row.names = FALSE)

# Plot only target genes
plot_df <- stage_mean %>%
  filter(GeneID %in% target_genes$GeneID) %>%
  left_join(target_genes, by = "GeneID") %>%
  mutate(Stage = factor(Stage, levels = stage_order))

if (nrow(plot_df) == 0) {
  message("[WARNING] plot_df is empty. Check whether GeneID in target_genes exists in tpm_mat rownames.")
} else {
  p <- ggplot(plot_df, aes(x = Stage, y = mean_TPM, group = GeneName)) +
    geom_line() + geom_point() +
    facet_wrap(~ GeneName, scales = "free_y") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x=NULL, y="Mean TPM", title="Time-course expression (mean TPM per stage)")

  ggsave(file.path(OUTDIR, "GOI_timecourse_TPM_faceted.pdf"), p, width = 12, height = 6)
}

# ==================== 8) ara/mirr strictly zygotic-like summary (pre low, post rise, peak stage) ====================
cat("Computing zygotic-like summaries for ara/mirr...\n")

summarize_zygotic_pattern <- function(gene_id, gene_name,
                                     pre_stages = c("Embryo_0-12h","Embryo_12-24h"),
                                     post_start = "Embryo_36-48h",
                                     stage_order,
                                     sustained_frac = 0.5) {

  df <- tpm_long %>%
    filter(GeneID == gene_id, as.character(Stage) %in% stage_order) %>%
    group_by(Stage) %>%
    summarise(meanTPM = mean(TPM, na.rm=TRUE), .groups="drop") %>%
    mutate(Stage = factor(Stage, levels=stage_order)) %>%
    arrange(Stage)

  pre_mean <- df %>%
    filter(as.character(Stage) %in% pre_stages) %>%
    summarise(pre = mean(meanTPM)) %>% pull(pre)

  post_stages <- stage_order[match(post_start, stage_order):length(stage_order)]
  post_df <- df %>% filter(as.character(Stage) %in% post_stages)

  peak_stage <- as.character(post_df$Stage[which.max(post_df$meanTPM)])
  peak_val   <- max(post_df$meanTPM)

  # “remain elevated”: all stages after peak >= sustained_frac * peak
  idx_peak <- match(peak_stage, stage_order)
  later_stages <- stage_order[idx_peak:length(stage_order)]
  later_df <- df %>% filter(as.character(Stage) %in% later_stages)

  sustained_high <- all(later_df$meanTPM >= sustained_frac * peak_val)

  data.frame(
    GeneName = gene_name,
    GeneID = gene_id,
    pre_mean_TPM = pre_mean,
    peak_stage = peak_stage,
    peak_TPM = peak_val,
    sustained_high = sustained_high,
    stringsAsFactors = FALSE
  )
}

ara_id  <- target_genes$GeneID[target_genes$GeneName == "ara"]
mirr_id <- target_genes$GeneID[target_genes$GeneName == "mirr"]

zygo_summary <- bind_rows(
  summarize_zygotic_pattern(ara_id, "ara", stage_order=stage_order),
  summarize_zygotic_pattern(mirr_id, "mirr", stage_order=stage_order)
)

write.csv(zygo_summary, file.path(OUTDIR, "ara_mirr_zygotic_summary.csv"), row.names = FALSE)
print(zygo_summary)

cat("\nDONE. All outputs saved in: ", OUTDIR, "\n")
