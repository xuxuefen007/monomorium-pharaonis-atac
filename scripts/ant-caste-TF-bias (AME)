# ============================== 1) TF importance ranking (stage-weighted) ==============================
#
# Input files:
#   *_with_log2.txt files are derived from the original AME (MEME Suite) motif enrichment outputs.
#   Each file corresponds to one developmental stage and contains motif-level enrichment statistics
#   comparing a positive set (e.g. caste-biased peaks) against a negative/background set.
#
#   Compared to the raw AME output, these *_with_log2.txt files include one additional column:
#     - log2_TP_over_FP_p05
#
#   This column represents the log2 enrichment of motif hits in the positive set relative to the
#   negative set, and is computed as:
#
#       log2_TP_over_FP_p05 = log2( ((TP + 0.5) / N_pos) / ((FP + 0.5) / N_neg) )
#
#   where:
#     TP     = number of motif hits in the positive set,
#     FP     = number of motif hits in the negative set,
#     N_pos  = total number of sequences in the positive set,
#     N_neg  = total number of sequences in the negative set.
#
#   The pseudocount of 0.5 is added to both TP and FP to avoid division by zero and stabilize the
#   estimate for rare motifs (hence the suffix "_p05").
#
# Note on caste-specific processing:
#   Gyne and worker are processed separately. This block computes TF importance for the gyne dataset
#   only (input_dir points to gyne-specific AME-derived *_with_log2.txt files). An analogous pipeline
#   should be run for worker in a separate directory (e.g. worker_with_log2), and integration is done
#   in later sections of this script.

# Set input/output directories
input_dir <- "/data/work/caste_biased_list/gyne_with_log2"
output_dir <- "/data/work/caste_biased_list/gyne_with_log2/integrated_analysis-a"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# List all AME result files
ame_files <- list.files(path = input_dir, pattern = "*_with_log2\\.txt$", full.names = TRUE)

# Extract developmental stage names
stage_names <- gsub("_with_log2\\.txt$", "", basename(ame_files))
stage_names <- gsub("_filtered$", "", stage_names)

cat("Detected developmental stages:\n")
print(stage_names)

# Initialize an empty list to store all data
all_data <- list()

# Read and process each file
for (i in seq_along(ame_files)) {
  stage <- stage_names[i]
  file <- ame_files[i]

  cat("Reading:", stage, "\n")

  # Read data
  df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  # Add stage annotations
  df$stage <- stage
  df$stage_rank <- df$rank

  # Select key columns
  key_cols <- c("motif_ID", "motif_alt_ID", "consensus", "p-value", "adj_p-value",
                "E-value", "%TP", "%FP", "log2_TP_over_FP_p05", "stage", "stage_rank")

  # Keep only available columns
  available_cols <- key_cols[key_cols %in% colnames(df)]
  df_sub <- df[, available_cols]

  # Store into list
  all_data[[stage]] <- df_sub
}

# Merge all stages
combined_data <- do.call(rbind, all_data)
rownames(combined_data) <- NULL

# Load required library
library(dplyr)

# De-duplication at the TF level:
#   AME outputs contain one row per motif per stage, so the same TF (motif_ID/motif_alt_ID) appears
#   repeatedly across developmental stages. Here we collapse these repeated entries by grouping on
#   (motif_ID, motif_alt_ID) and summarizing across stages, resulting in ONE row per TF.
tf_scores_raw <- combined_data %>%
  group_by(motif_ID, motif_alt_ID) %>%
  summarise(
    # Cross-stage reproducibility
    n_stages = n_distinct(stage),
    stages_present = paste(sort(unique(stage)), collapse = ", "),

    # Statistical significance
    mean_neg_log10_adj_p = mean(-log10(`adj_p-value` + 1e-323)),

    # Enrichment strength
    mean_log2_enrich = mean(log2_TP_over_FP_p05, na.rm = TRUE),

    # Specificity proxy
    mean_specificity = mean(`%TP` - `%FP`, na.rm = TRUE),

    # Mean rank across stages
    mean_rank = mean(stage_rank, na.rm = TRUE),

    .groups = 'drop'
  )

# Safe min-max normalization function
safe_normalize <- function(x) {
  if (length(unique(x)) == 1 || all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    if (all(is.na(x))) return(rep(NA, length(x)))
    return(rep(0.5, length(x)))
  } else {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
}

# Compute normalized scores (stage count is weighted higher)
tf_scores <- tf_scores_raw %>%
  mutate(
    score_n_stages = safe_normalize(n_stages),
    score_significance = safe_normalize(mean_neg_log10_adj_p),
    score_enrichment = safe_normalize(mean_log2_enrich),
    score_specificity = safe_normalize(mean_specificity),

    composite_score = (score_n_stages * 0.4) +
                     (score_significance * 0.25) +
                     (score_enrichment * 0.25) +
                     (score_specificity * 0.10)
  ) %>%
  filter(!is.na(mean_log2_enrich), !is.na(mean_neg_log10_adj_p)) %>%
  arrange(desc(composite_score)) %>%
  mutate(overall_rank = 1:n())

# Reorder columns
final_result <- tf_scores %>%
  select(overall_rank, motif_ID, motif_alt_ID, composite_score, n_stages, stages_present,
         mean_neg_log10_adj_p, mean_log2_enrich, mean_specificity, mean_rank,
         score_n_stages, score_significance, score_enrichment, score_specificity)

# Save results
output_file <- file.path(output_dir, "TF_importance_ranking_stage_weighted.csv")
write.csv(final_result, output_file, row.names = FALSE)

# Note: worker uses the same procedure (not shown here)
# worker_ame_dir <- "/data/work/caste_biased_list/worker_with_log2"

# ============================== 2) TF bias analysis between gyne and worker ==============================

# Set paths
gyne_ame_dir <- "/data/work/caste_biased_list/gyne_with_log2"
worker_ame_dir <- "/data/work/caste_biased_list/worker_with_log2"
output_dir <- "/data/work/caste_biased_list/duplicate"

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)

# Read AME raw tables from a directory
read_ame_files <- function(dir_path) {
  ame_files <- list.files(dir_path, pattern = "*_with_log2\\.txt$", full.names = TRUE)
  if (length(ame_files) == 0) {
    stop("No AME files found. Please check path: ", dir_path)
  }

  all_data <- list()
  for (file in ame_files) {
    stage <- gsub("_with_log2\\.txt$", "", basename(file))
    df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    df$stage <- stage
    all_data[[stage]] <- df
  }

  combined <- do.call(rbind, all_data)
  rownames(combined) <- NULL
  return(combined)
}

cat("Reading gyne data...\n")
gyne_ame_raw <- read_ame_files(gyne_ame_dir)
cat("Reading worker data...\n")
worker_ame_raw <- read_ame_files(worker_ame_dir)

cat("Gyne data: ", nrow(gyne_ame_raw), " rows; ", n_distinct(gyne_ame_raw$motif_ID), " TFs\n")
cat("Worker data: ", nrow(worker_ame_raw), " rows; ", n_distinct(worker_ame_raw$motif_ID), " TFs\n")

# Compute mean metrics per TF across stages
calculate_tf_stats <- function(ame_data, caste_name) {
  tf_stats <- ame_data %>%
    group_by(motif_ID, motif_alt_ID) %>%
    summarise(
      n_stages = n_distinct(stage),
      stages_present = paste(sort(unique(stage)), collapse = ", "),
      mean_log2_enrich = mean(log2_TP_over_FP_p05, na.rm = TRUE),
      mean_neg_log10_adj_p = mean(-log10(`adj_p-value` + 1e-323), na.rm = TRUE),
      mean_specificity = mean(`%TP` - `%FP`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(caste = caste_name)
  return(tf_stats)
}

cat("Computing TF-level statistics...\n")
gyne_stats <- calculate_tf_stats(gyne_ame_raw, "gyne")
worker_stats <- calculate_tf_stats(worker_ame_raw, "worker")

# Identify shared TFs (appear in both gyne and worker)
common_tfs <- intersect(gyne_stats$motif_ID, worker_stats$motif_ID)
cat("Number of shared TFs:", length(common_tfs), "\n")

if (length(common_tfs) == 0) {
  stop("No shared TFs found; bias analysis cannot proceed.")
}

gyne_common <- gyne_stats %>% filter(motif_ID %in% common_tfs)
worker_common <- worker_stats %>% filter(motif_ID %in% common_tfs)

# Compute bias metrics for shared TFs
cat("Computing bias metrics...\n")
bias_analysis <- gyne_common %>%
  inner_join(worker_common,
            by = c("motif_ID", "motif_alt_ID"),
            suffix = c("_gyne", "_worker")) %>%
  mutate(
    log2_enrich_diff = mean_log2_enrich_gyne - mean_log2_enrich_worker,
    significance_diff = mean_neg_log10_adj_p_gyne - mean_neg_log10_adj_p_worker,
    composite_bias_score = 0.6 * log2_enrich_diff + 0.4 * significance_diff,

    enrich_ratio = ifelse(mean_log2_enrich_worker == 0 | mean_log2_enrich_worker < 0.1,
                         mean_log2_enrich_gyne / 0.1,
                         mean_log2_enrich_gyne / mean_log2_enrich_worker),
    log2_enrich_ratio = log2(pmax(enrich_ratio, 0.001)),

    bias_direction = ifelse(composite_bias_score > 0, "Gyne", "Worker"),
    abs_bias_score = abs(composite_bias_score)
  )

# ================== Scientific threshold setting ==================
cat("Setting thresholds...\n")

bias_abs_scores <- abs(bias_analysis$composite_bias_score)
log2_ratios <- abs(bias_analysis$log2_enrich_ratio)

# Method 1: quantile-based thresholds
threshold_strong_quantile <- quantile(bias_abs_scores, 0.85, na.rm = TRUE)   # top 15%
threshold_moderate_quantile <- quantile(bias_abs_scores, 0.65, na.rm = TRUE) # top 35%

# Method 2: MAD-based thresholds
score_stats <- list(
  median = median(bias_abs_scores, na.rm = TRUE),
  mad = mad(bias_abs_scores, na.rm = TRUE)
)
threshold_strong_mad <- score_stats$median + 2.5 * score_stats$mad
threshold_moderate_mad <- score_stats$median + 1.5 * score_stats$mad

# Method 3: log2 ratio-based thresholds
log2_threshold_strong <- quantile(log2_ratios, 0.85, na.rm = TRUE)
log2_threshold_moderate <- quantile(log2_ratios, 0.65, na.rm = TRUE)

# Combine multiple methods
final_threshold_strong <- mean(c(threshold_strong_quantile, threshold_strong_mad,
                                log2_threshold_strong), na.rm = TRUE)
final_threshold_moderate <- mean(c(threshold_moderate_quantile, threshold_moderate_mad,
                                  log2_threshold_moderate), na.rm = TRUE)

cat("Threshold summary:\n")
cat("Strong threshold:", round(final_threshold_strong, 3), "\n")
cat("Moderate threshold:", round(final_threshold_moderate, 3), "\n")
cat("Log2 ratio strong threshold:", round(log2_threshold_strong, 3), "\n")
cat("Log2 ratio moderate threshold:", round(log2_threshold_moderate, 3), "\n")

# Apply thresholds to classify TFs
bias_analysis <- bias_analysis %>%
  mutate(
    bias_strength = case_when(
      abs_bias_score >= final_threshold_strong |
      abs(log2_enrich_ratio) >= log2_threshold_strong ~ "Strong",

      abs_bias_score >= final_threshold_moderate |
      abs(log2_enrich_ratio) >= log2_threshold_moderate ~ "Moderate",

      TRUE ~ "Weak"
    ),

    absolute_bias_category = paste(bias_strength, bias_direction, "biased", sep = "_"),

    confidence_score = abs_bias_score * (0.5 + 0.3 * (n_stages_gyne + n_stages_worker) / 10),

    stage_coverage = paste0("G", n_stages_gyne, "/W", n_stages_worker)
  )

# Summarize classification results
cat("\n=== Classification summary ===\n")
classification_summary <- bias_analysis %>%
  group_by(absolute_bias_category) %>%
  summarise(
    count = n(),
    percent = round(n() / nrow(bias_analysis) * 100, 1),
    avg_bias_score = mean(abs_bias_score),
    avg_log2_ratio = mean(abs(log2_enrich_ratio)),
    avg_confidence = mean(confidence_score),
    .groups = "drop"
  ) %>%
  arrange(desc(avg_bias_score))

print(classification_summary)

# Visualize threshold distribution
threshold_plot <- ggplot(bias_analysis, aes(x = abs_bias_score, fill = bias_strength)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  geom_vline(xintercept = c(final_threshold_moderate, final_threshold_strong),
             linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = final_threshold_moderate, y = 5,
           label = paste("Moderate:", round(final_threshold_moderate, 2)),
           hjust = -0.1, color = "red") +
  annotate("text", x = final_threshold_strong, y = 3,
           label = paste("Strong:", round(final_threshold_strong, 2)),
           hjust = -0.1, color = "red") +
  labs(title = "Bias score distribution and thresholds",
       x = "|Composite Bias Score|", y = "Count", fill = "Strength") +
  theme_minimal()

ggsave(file.path(output_dir, "scientific_threshold_distribution.pdf"),
       threshold_plot, width = 10, height = 6)

# Sort and export results by category
cat("\n=== Final results sorted by category ===\n")

category_order <- c("Strong_Gyne_biased", "Moderate_Gyne_biased", "Weak_Gyne_biased",
                   "Strong_Worker_biased", "Moderate_Worker_biased", "Weak_Worker_biased")

bias_analysis_sorted <- bias_analysis %>%
  mutate(absolute_bias_category = factor(absolute_bias_category, levels = category_order)) %>%
  group_by(absolute_bias_category) %>%
  arrange(desc(confidence_score), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    overall_rank = row_number(),
    category_rank = ave(desc(confidence_score), absolute_bias_category,
                       FUN = function(x) rank(x, ties.method = "first"))
  )

# Print top 3 TFs per category
for (category in category_order) {
  if (category %in% bias_analysis_sorted$absolute_bias_category) {
    cat("\n", category, "(Top 3):\n")
    top_tfs <- bias_analysis_sorted %>%
      filter(absolute_bias_category == category) %>%
      head(3) %>%
      select(category_rank, motif_alt_ID, confidence_score,
             composite_bias_score, log2_enrich_ratio, stage_coverage)
    print(top_tfs)
  }
}

# Save bias analysis results
output_file <- file.path(output_dir, "TF_scientific_bias_classification.csv")
write.csv(bias_analysis_sorted, output_file, row.names = FALSE)

cat("\nBias classification completed.\n")
cat("Shared TFs:", length(common_tfs), "\n")
cat("Output file:", output_file, "\n")
cat("Threshold plot: scientific_threshold_distribution.pdf\n")

# Final summary
cat("\n=== Final summary ===\n")
cat("Strong TFs:", sum(grepl("Strong", bias_analysis$absolute_bias_category)), "\n")
cat("Moderate TFs:", sum(grepl("Moderate", bias_analysis$absolute_bias_category)), "\n")
cat("Weak TFs:", sum(grepl("Weak", bias_analysis$absolute_bias_category)), "\n")
cat("Gyne-biased TFs:", sum(grepl("Gyne", bias_analysis$absolute_bias_category)), "\n")
cat("Worker-biased TFs:", sum(grepl("Worker", bias_analysis$absolute_bias_category)), "\n")


##-------------------------- 3) Integrate into final unique TF lists -------------------#
# ================== Set paths ==================
gyne_ranking_file <- "/data/work/caste_biased_list/duplicate/gyne_TF_importance_ranking_stage_weighted.csv"
worker_ranking_file <- "/data/work/caste_biased_list/duplicate/worker_TF_importance_ranking_stage_weighted.csv"
bias_analysis_file <- "/data/work/caste_biased_list/duplicate/TF_scientific_bias_classification.csv"
output_dir <- "/data/work/caste_biased_list/duplicate/final_results"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

library(dplyr)

cat("Reading input tables...\n")
gyne_ranking <- read.csv(gyne_ranking_file, stringsAsFactors = FALSE)
worker_ranking <- read.csv(worker_ranking_file, stringsAsFactors = FALSE)
bias_analysis <- read.csv(bias_analysis_file, stringsAsFactors = FALSE)

# Rename columns to avoid conflicts
gyne_ranking <- gyne_ranking %>%
  rename_with(~ paste0(., "_gyne"), -c(motif_ID, motif_alt_ID))

worker_ranking <- worker_ranking %>%
  rename_with(~ paste0(., "_worker"), -c(motif_ID, motif_alt_ID))

# Identify unique and shared TFs across gyne/worker rankings
all_gyne_tfs <- unique(gyne_ranking$motif_ID)
all_worker_tfs <- unique(worker_ranking$motif_ID)
common_tfs <- intersect(all_gyne_tfs, all_worker_tfs)

cat("TF counts:\n")
cat("Gyne TFs (total):", length(all_gyne_tfs), "\n")
cat("Worker TFs (total):", length(all_worker_tfs), "\n")
cat("Shared TFs:", length(common_tfs), "\n")

# De-duplication by assignment:
#   Shared TFs (common_tfs) must NOT appear in both final outputs.
#   We assign each shared TF to exactly one caste (gyne or worker) using the bias classification.
bias_classification <- bias_analysis %>%
  select(motif_ID, motif_alt_ID, absolute_bias_category, confidence_score)

gyne_biased_tfs <- bias_classification %>%
  filter(grepl("Gyne", absolute_bias_category)) %>%
  pull(motif_ID)

worker_biased_tfs <- bias_classification %>%
  filter(grepl("Worker", absolute_bias_category)) %>%
  pull(motif_ID)

cat("Gyne-biased shared TFs:", length(gyne_biased_tfs), "\n")
cat("Worker-biased shared TFs:", length(worker_biased_tfs), "\n")

# Build final gyne TF list (de-duplicated):
#   Keep (i) gyne-unique TFs, and (ii) shared TFs classified as Gyne-biased.
#   Shared TFs classified as Worker-biased are excluded here to prevent duplication.
final_gyne_tfs <- gyne_ranking %>%
  filter(!motif_ID %in% common_tfs | motif_ID %in% gyne_biased_tfs) %>%
  left_join(bias_classification %>% select(motif_ID, absolute_bias_category, confidence_score),
            by = "motif_ID") %>%
  mutate(
    tf_type = ifelse(motif_ID %in% common_tfs, "Shared_Gyne_biased", "Gyne_unique"),
    final_rank = row_number()
  ) %>%
  select(final_rank, motif_ID, motif_alt_ID, tf_type, absolute_bias_category,
         composite_score_gyne, confidence_score, everything())

# Build final worker TF list (de-duplicated):
#   Keep (i) worker-unique TFs, and (ii) shared TFs classified as Worker-biased.
#   Shared TFs classified as Gyne-biased are excluded here to prevent duplication.
final_worker_tfs <- worker_ranking %>%
  filter(!motif_ID %in% common_tfs | motif_ID %in% worker_biased_tfs) %>%
  left_join(bias_classification %>% select(motif_ID, absolute_bias_category, confidence_score),
            by = "motif_ID") %>%
  mutate(
    tf_type = ifelse(motif_ID %in% common_tfs, "Shared_Worker_biased", "Worker_unique"),
    final_rank = row_number()
  ) %>%
  select(final_rank, motif_ID, motif_alt_ID, tf_type, absolute_bias_category,
         composite_score_worker, confidence_score, everything())

# Save final outputs
gyne_output <- file.path(output_dir, "final_gyne_unique_tfs.csv")
write.csv(final_gyne_tfs, gyne_output, row.names = FALSE)

worker_output <- file.path(output_dir, "final_worker_unique_tfs.csv")
write.csv(final_worker_tfs, worker_output, row.names = FALSE)

# Save shared TF assignment details
shared_tfs_output <- file.path(output_dir, "shared_tfs_assignment.csv")
shared_tfs_info <- bias_analysis %>%
  select(motif_ID, motif_alt_ID, absolute_bias_category, confidence_score,
         composite_bias_score, log2_enrich_ratio,
         mean_log2_enrich_gyne, mean_log2_enrich_worker) %>%
  arrange(desc(confidence_score))
write.csv(shared_tfs_info, shared_tfs_output, row.names = FALSE)

# Report summary
cat("\n=== Final TF assignment summary ===\n")
cat("Final gyne TFs:", nrow(final_gyne_tfs), "\n")
cat("  - Gyne-unique:", sum(final_gyne_tfs$tf_type == "Gyne_unique"), "\n")
cat("  - Shared (gyne-biased):", sum(final_gyne_tfs$tf_type == "Shared_Gyne_biased"), "\n")

cat("Final worker TFs:", nrow(final_worker_tfs), "\n")
cat("  - Worker-unique:", sum(final_worker_tfs$tf_type == "Worker_unique"), "\n")
cat("  - Shared (worker-biased):", sum(final_worker_tfs$tf_type == "Shared_Worker_biased"), "\n")

cat("\nTop 5 gyne TFs:\n")
print(final_gyne_tfs %>% head(5) %>% select(final_rank, motif_alt_ID, tf_type, composite_score_gyne))

cat("\nTop 5 worker TFs:\n")
print(final_worker_tfs %>% head(5) %>% select(final_rank, motif_alt_ID, tf_type, composite_score_worker))

# Summary table
summary_stats <- data.frame(
  Category = c("Total_Gyne_TFs", "Total_Worker_TFs",
               "Unique_Gyne_TFs", "Unique_Worker_TFs",
               "Shared_Gyne_Biased", "Shared_Worker_Biased"),
  Count = c(nrow(final_gyne_tfs), nrow(final_worker_tfs),
            sum(final_gyne_tfs$tf_type == "Gyne_unique"),
            sum(final_worker_tfs$tf_type == "Worker_unique"),
            sum(final_gyne_tfs$tf_type == "Shared_Gyne_biased"),
            sum(final_worker_tfs$tf_type == "Shared_Worker_biased"))
)

write.csv(summary_stats, file.path(output_dir, "final_assignment_summary.csv"), row.names = FALSE)

cat("\nIntegration completed.\n")
cat("Output directory:", output_dir, "\n")
cat("final_gyne_unique_tfs.csv - Final gyne TF list (ranked)\n")
cat("final_worker_unique_tfs.csv - Final worker TF list (ranked)\n")
cat("shared_tfs_assignment.csv - Shared TF assignment details\n")
cat("final_assignment_summary.csv - Summary statistics\n")


