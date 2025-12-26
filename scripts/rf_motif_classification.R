# ============================================================
setwd("C:/Users/Xuefen Xu/Documents/RF")  # Set working directory
rm(list = ls()); gc()  # Clear workspace and run garbage collection
options(stringsAsFactors = FALSE)  # Disable automatic conversion of strings to factors
set.seed(42)  # Set random seed for reproducibility

# -------------------- Basic Configuration --------------------
INPUT_PFM  <- "gyne96_worker159_motif.txt"   
# This file contains sequence data for 96 high-confidence gyne motifs and 159 worker motifs (see Table S3). 
OUTDIR     <- file.path(getwd(), "rf_gw_eval")  # Output directory for results
POS_CLASS  <- "gyne"   # Define the positive class label (case-insensitive)

# -------------------- Boundary Thresholds --------------------
MARGIN_CUTOFF <- 0.20   # |P(gyne)-P(worker)| < 0.20 → Probability boundary for uncertainty
SEP_CUTOFF    <- 0.05   # Proximity separation < 0.05 → Similarity boundary (negative values indicate stronger similarity)

# Create output directory if it doesn't exist
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -------------------- Install Required Packages --------------------
need_pkgs <- c("tidyverse", "randomForest", "pheatmap", "RColorBrewer")
# Install missing packages from CRAN repository
for (p in need_pkgs) if (!requireNamespace(p, quietly = TRUE)) {
  install.packages(p, repos = "https://cloud.r-project.org")
}
# Load necessary libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(pheatmap)
  library(RColorBrewer)
})

# -------------------- 1) Parse PFM File --------------------
# This function parses the PFM file. Each motif should have a header like:
# >gyne_xxx\tTFNAME or >worker_xxx\tTFNAME, followed by four rows of A/C/G/T base frequency vectors.
parse_pfm_file <- function(file_path) {
  lines <- readLines(file_path)  # Read all lines from the PFM file
  motifs <- list(); labels <- character(); full_names <- character()  # Initialize empty lists to store motifs and labels
  current_matrix <- NULL; current_name <- ""; current_full_name <- ""  # Temporary variables for the current motif

  i <- 1
  while (i <= length(lines)) {
    line <- trimws(lines[i])  # Remove leading/trailing whitespace

    if (startsWith(line, ">")) {
      # Finalize previous motif before starting a new one
      if (!is.null(current_matrix) && ncol(current_matrix) > 0) {
        # Ensure all A/C/G/T bases are included
        if (!all(c("A","C","G","T") %in% rownames(current_matrix))) {
          miss <- setdiff(c("A","C","G","T"), rownames(current_matrix))
          for (b in miss) {
            current_matrix <- rbind(
              current_matrix,
              setNames(data.frame(t(rep(0, ncol(current_matrix)))), b)
            )
          }
          current_matrix <- current_matrix[c("A","C","G","T"), , drop = FALSE]
        }
        motifs <- c(motifs, list(as.matrix(current_matrix)))  # Add the current motif to the list
        labels <- c(labels, current_name)  # Store the label for this motif
        full_names <- c(full_names, current_full_name)  # Store the full name of this motif
      }
      head_line <- sub("^>", "", line)  # Remove ">" from the header line
      current_full_name <- head_line  # Save the full name
      current_name <- strsplit(head_line, "\t")[[1]][1]  # Extract the motif name (e.g., gyne_MAxxxx or worker_MAxxxx)
      current_matrix <- matrix(nrow = 0, ncol = 0)  # Reset the matrix for the new motif
      i <- i + 1; next  # Move to the next line
    }

    if (grepl("^[ACGT]\\s*\\[", line)) {
      # Parse the A/C/G/T base frequency rows
      parts <- strsplit(line, "\\s*\\[\\s*|\\s*\\]")[[1]]
      base  <- trimws(parts[1])  # The base (A, C, G, or T)
      nums  <- as.numeric(strsplit(trimws(parts[2]), "\\s+")[[1]])  # The frequency values for the base
      if (is.null(current_matrix) || ncol(current_matrix) == 0)
        current_matrix <- matrix(0, nrow = 0, ncol = length(nums))  # Initialize the matrix if it's empty
      if (!base %in% rownames(current_matrix)) {
        current_matrix <- rbind(current_matrix, nums)  # Add the new base to the matrix
        rownames(current_matrix)[nrow(current_matrix)] <- base
      } else {
        current_matrix[base, ] <- nums  # Update the existing base row with new values
      }
    }
    i <- i + 1
  }

  # Finalize the last motif after finishing the loop
  if (!is.null(current_matrix) && ncol(current_matrix) > 0) {
    if (!all(c("A","C","G","T") %in% rownames(current_matrix))) {
      miss <- setdiff(c("A","C","G","T"), rownames(current_matrix))
      for (b in miss) {
        current_matrix <- rbind(
          current_matrix,
          setNames(data.frame(t(rep(0, ncol(current_matrix)))), b)
        )
      }
      current_matrix <- current_matrix[c("A","C","G","T"), , drop = FALSE]
    }
    motifs <- c(motifs, list(as.matrix(current_matrix)))
    labels <- c(labels, current_name)
    full_names <- c(full_names, current_full_name)
  }

  message(sprintf("✅ Successfully parsed %d motifs", length(motifs)))
  list(motifs = motifs, labels = labels, full_names = full_names)
}

# -------------------- 2) Feature Engineering --------------------
# This function extracts various features from a PFM matrix, such as information content, motif length, GC content, etc.
extract_features <- function(pfm_matrix) {
  col_sums <- colSums(pfm_matrix) + 1e-8  # Sum of each column (with a small offset)
  prob_matrix <- sweep(pfm_matrix, 2, col_sums, "/")  # Normalize the matrix by column
  prob_matrix <- pmin(pmax(prob_matrix, 1e-8), 1 - 1e-8)  # Clip values to avoid numerical issues

  ic_per_position <- 2 + colSums(prob_matrix * log2(prob_matrix), na.rm = TRUE)  # Information content per position
  H <- -colSums(prob_matrix * log(prob_matrix), na.rm = TRUE)  # Entropy
  conservation <- 1 - H / log(4)  # Conservation score

  total <- sum(pfm_matrix) + 1e-8  # Total sum of all base frequencies
  tibble(
    ic_mean = mean(ic_per_position, na.rm = TRUE),
    ic_sd   = sd(ic_per_position, na.rm = TRUE),
    ic_max  = max(ic_per_position, na.rm = TRUE),
    motif_length = ncol(pfm_matrix),
    total_count  = sum(pfm_matrix),
    pref_mean = mean(apply(prob_matrix, 2, max)),
    pref_sd   = sd(apply(prob_matrix, 2, max)),
    cons_mean = mean(conservation, na.rm = TRUE),
    cons_sd   = sd(conservation, na.rm = TRUE),
    gc_content = (sum(pfm_matrix[c("C","G"), , drop=FALSE]) / total),
    A_content  = (sum(pfm_matrix["A", , drop=FALSE]) / total),
    C_content  = (sum(pfm_matrix["C", , drop=FALSE]) / total),
    G_content  = (sum(pfm_matrix["G", , drop=FALSE]) / total),
    T_content  = (sum(pfm_matrix["T", , drop=FALSE]) / total)
  )
}

# -------------------- 3) Random Forest Training --------------------
# Prepare the dataset by extracting features and labels from the parsed motifs
parsed  <- parse_pfm_file(INPUT_PFM)
dataset <- prepare_dataset(parsed$motifs, parsed$labels, parsed$full_names)

# Select features and labels for training
x_all <- dplyr::select(dataset, -label, -tf_full_info, -tf_id)
y_all <- dataset$label
rownames(x_all) <- dataset$tf_id

# Set random seed and train the Random Forest model
set.seed(42)
rf_full <- randomForest(
  x = x_all,
  y = y_all,
  ntree = 500,  # Number of trees in the forest
  importance = TRUE,  # Calculate feature importance
  proximity = TRUE,  # Enable proximity calculation
  oob.prox = TRUE   # Use out-of-bag proximity for more conservative results
)

# -------------------- 4) OOB Quantitative Results --------------------
# 4.1 Built-in OOB confusion matrix and error rates
conf_oob_builtin <- rf_full$confusion[1:2, 1:2, drop = FALSE]
oob_rates        <- rf_full$err.rate[nrow(rf_full$err.rate), ]

# Save confusion matrix and error rates to CSV files
write.csv(as.data.frame(conf_oob_builtin),
          file.path(OUTDIR, "rf_OOB_confusion_matrix.csv"))
write.csv(as.data.frame(t(oob_rates)),
          file.path(OUTDIR, "rf_OOB_error_rates.csv"), row.names = TRUE)

# Print confusion matrix and error rates
cat("\n=== OOB Confusion Matrix (builtin) ===\n"); print(conf_oob_builtin)
cat("\n=== OOB Error Rates (builtin) ===\n");  print(oob_rates)

# 4.2 Derive metrics from the confusion matrix (Accuracy, Sensitivity, Specificity, Balanced Accuracy)
oob_acc  <- sum(diag(conf_oob_builtin)) / sum(conf_oob_builtin)
oob_sens <- conf_oob_builtin["gyne","gyne"]     / sum(conf_oob_builtin["gyne",])
oob_spec <- conf_oob_builtin["worker","worker"] / sum(conf_oob_builtin["worker",])
oob_bal  <- mean(c(oob_sens, oob_spec))

# Create a summary table for the metrics
metrics_tbl <- tibble(
  metric = c("OOB_Accuracy","OOB_Sensitivity_gyne","OOB_Specificity_worker","OOB_Balanced_Accuracy"),
  value  = c(oob_acc, oob_sens, oob_spec, oob_bal)
)

# Save metrics to a CSV file
write.csv(metrics_tbl, file.path(OUTDIR, "rf_OOB_metrics.csv"), row.names = FALSE)
cat(sprintf("\nOOB: Acc=%.3f, Sens(gyne)=%.3f, Spec(worker)=%.3f, BalAcc=%.3f\n",
            oob_acc, oob_sens, oob_spec, oob_bal))

# -------------------- 5) In-sample Confusion Matrix --------------------
# Predict on the training set (in-sample) and create confusion matrix
pred_in <- predict(rf_full, newdata = x_all, type = "response")
conf_in <- table(True = y_all, Pred = pred_in)

# Ensure the row/column order is consistent (gyne, worker)
conf_in <- conf_in[c("gyne", "worker"), c("gyne", "worker")]

# Save confusion matrix to CSV
write.csv(as.data.frame(conf_in),
          file.path(OUTDIR, "rf_in_sample_confusion_matrix.csv"))
cat("\n=== In-sample Confusion Matrix (Training Data) ===\n"); print(conf_in)

# -------------------- 6) Boundary Degree Calculation --------------------
# OOB probability matrix (two columns: gyne / worker)
prob_oob <- rf_full$votes
colnames(prob_oob) <- tolower(colnames(prob_oob))
if (!"gyne" %in% colnames(prob_oob)) stop("Column 'gyne' not found in votes")
prob_oob <- prob_oob[, c("gyne","worker")]  # Ensure column order

# Compute margin and proximity separation for each sample
p_gy <- prob_oob[, "gyne"]; p_wk <- prob_oob[, "worker"]
pred_oob <- rf_full$predicted
correct  <- pred_oob == y_all

# Calculate probability margin (lower margin indicates uncertainty)
oob_margin <- abs(p_gy - p_wk)

# Compute OOB-based proximity (from oob.prox=TRUE)
prox <- rf_full$proximity; prox[is.na(prox)] <- 0
labs <- y_all
sep_vec <- numeric(nrow(prox))
for (i in seq_len(nrow(prox))) {
  same <- prox[i, labs == labs[i] & seq_len(nrow(prox)) != i]
  othr <- prox[i, labs != labs[i]]
  m_same <- if (length(same)) mean(same) else 0
  m_othr <- if (length(othr)) mean(othr) else 0
  sep_vec[i] <- m_same - m_othr
}
sep_vec[!is.finite(sep_vec)] <- 0

# Summarize each sample's information (include per-sample CSV for review)
per_sample <- tibble(
  tf_id        = rownames(x_all),
  tf_full_info = dataset$tf_full_info[match(rownames(x_all), dataset$tf_id)],
  true_label   = as.character(y_all),
  pred_oob     = as.character(pred_oob),
  correct      = as.logical(correct),
  oob_p_gyne   = as.numeric(p_gy),
  oob_p_worker = as.numeric(p_wk),
  oob_margin   = as.numeric(oob_margin),
  prox_sep     = as.numeric(sep_vec)
) |>
  mutate(
    OOB_Boundary  = (oob_margin < MARGIN_CUTOFF | !correct),
    Prox_Boundary = (prox_sep   < SEP_CUTOFF),
    boundary_union    = OOB_Boundary | Prox_Boundary,
    boundary_intersec = OOB_Boundary & Prox_Boundary
  ) |>
  arrange(oob_margin, prox_sep)

# Save boundary summary data to CSV files
write.csv(per_sample, file.path(OUTDIR, "rf_per_sample_summary.csv"), row.names = FALSE)
write.csv(per_sample %>% filter(boundary_union),
          file.path(OUTDIR, "rf_boundary_TFs_union.csv"), row.names = FALSE)
write.csv(per_sample %>% filter(boundary_intersec),
          file.path(OUTDIR, "rf_boundary_TFs_intersection.csv"), row.names = FALSE)
write.csv(per_sample, file.path(OUTDIR, "rf_boundary_TFs_full_table.csv"), row.names = FALSE)

# -------------------- 7) Feature Importance --------------------
# Calculate and save feature importance from Random Forest
imp <- as.data.frame(rf_full$importance)
imp <- tibble(feature = rownames(imp), MeanDecreaseGini = imp$MeanDecreaseGini) |>
  arrange(desc(MeanDecreaseGini))
write.csv(imp, file.path(OUTDIR, "rf_feature_importance.csv"), row.names = FALSE)

# -------------------- 8) Heatmap Visualization --------------------
# Generate row annotations for heatmap visualization
row_anno <- data.frame(
  TrueLabel     = y_all,
  PredOOB       = pred_oob,
  Correct       = ifelse(correct, "Yes", "No"),
  OOB_Boundary  = ifelse(per_sample$OOB_Boundary,  "Yes", "No"),
  Prox_Boundary = ifelse(per_sample$Prox_Boundary, "Yes", "No")
)
rownames(row_anno) <- rownames(x_all)

# Define annotation colors
ann_colors <- list(
  TrueLabel     = c(gyne = "#d95f02", worker = "#1b9e77"),
  PredOOB       = c(gyne = "#fc8d62", worker = "#66c2a5"),
  Correct       = c(Yes = "#4daf4a", No = "#e41a1c"),
  OOB_Boundary  = c(Yes = "#984ea3", No = "#cccccc"),
  Prox_Boundary = c(Yes = "#377eb8", No = "#cccccc")
)

# (1) Proximity Heatmap (Sample × Sample)
dist_mat <- as.dist(1 - prox)
pdf(file.path(OUTDIR, "rf_proximity_heatmap.pdf"), width = 9.5, height = 9.5)
pheatmap(
  mat = prox,
  color = colorRampPalette(brewer.pal(9, "YlGnBu"))(200),
  clustering_distance_rows = dist_mat,
  clustering_distance_cols = dist_mat,
  clustering_method = "average",
  annotation_row = row_anno,
  annotation_colors = ann_colors,
  show_rownames = FALSE, show_colnames = FALSE,
  main = "Random Forest Proximity Heatmap (OOB accumulation)"
)
dev.off()

# (2) OOB Probability Heatmap (Sample × Class)
prob_mat <- as.matrix(prob_oob)
rownames(prob_mat) <- rownames(x_all)
row_dist <- dist(prob_mat, method = "euclidean")

pdf(file.path(OUTDIR, "rf_probability_heatmap.pdf"), width = 7, height = 9.5)
pheatmap(
  mat = prob_mat,
  color = colorRampPalette(brewer.pal(9, "PuBuGn"))(200),
  clustering_distance_rows = row_dist,
  clustering_method = "average",
  annotation_row = row_anno,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  main = "OOB Predicted Probability Heatmap (sample × class)"
)
dev.off()
