## =========================================================

rm(list = ls()); gc()  # Clear the workspace and run garbage collection
options(stringsAsFactors = FALSE)  # Disable automatic conversion of strings to factors
set.seed(42)  # Set random seed for reproducibility

## -------------------------
## Path Configuration
## -------------------------
setwd("C:/Users/Xuefen Xu/Documents/")  # Set working directory
INPUT_PFM  <- "gyne96_worker159_motif.txt"  # Path to the input PFM file containing sequence data for 96 high-confidence gyne motifs and 159 worker motifs (see Table S3)
OUTDIR    <- file.path(getwd(), "rf_shap_perm")  # Output directory to save results
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)  # Create output directory if it doesn't exist

## -------------------------
## Install Required Packages
## -------------------------
pkgs <- c(
  "tidyverse", "randomForest", "fastshap", "ggbeeswarm",  # Data manipulation and model building
  "RColorBrewer", "gridExtra", "grid", "vip", "yardstick", "ggplot2"  # Visualization and performance metrics
)
# Install any missing packages from CRAN repository
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
# Suppress package startup messages
suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(fastshap)
  library(ggbeeswarm)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(vip)
  library(yardstick)
  library(ggplot2)
})

## -------------------------
# Function to calculate standard deviation, ignoring non-finite values
sd0 <- function(x) {
  x <- x[is.finite(x)]  # Remove non-finite values
  if (length(x) < 2) return(0)  # Return 0 if less than 2 valid values
  sd(x)  # Standard deviation
}

## -------------------------
## Parse PFM File
## -------------------------
# Function to parse PFM (Position Frequency Matrix) files
parse_pfm_file <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)  # Read file lines
  motifs <- list(); labels <- character(); full_names <- character()  # Initialize lists to store motifs and labels
  current_matrix <- NULL; current_name <- ""; current_full_name <- ""  # Temporary variables for motif parsing
  i <- 1
  while (i <= length(lines)) {
    line <- trimws(lines[i])  # Trim whitespace from the line
    if (startsWith(line, ">")) {  # If the line starts with ">", it's a motif header
      if (!is.null(current_matrix) && ncol(current_matrix) > 0) {  # If a motif is already parsed
        # Ensure all A/C/G/T bases are included
        if (!all(c("A", "C", "G", "T") %in% rownames(current_matrix))) {
          miss <- setdiff(c("A", "C", "G", "T"), rownames(current_matrix))  # Find missing bases
          for (b in miss) {
            current_matrix <- rbind(
              current_matrix,
              setNames(data.frame(t(rep(0, ncol(current_matrix)))), b)
            )
          }
        }
        current_matrix <- current_matrix[c("A", "C", "G", "T"), , drop = FALSE]  # Reorder bases to A/C/G/T
        motifs     <- c(motifs, list(as.matrix(current_matrix)))  # Add parsed motif to list
        labels     <- c(labels, current_name)  # Add label
        full_names <- c(full_names, current_full_name)  # Add full name
      }
      head_line <- sub("^>", "", line)  # Remove ">" from header line
      current_full_name <- head_line  # Save full name
      current_name <- strsplit(head_line, "\t")[[1]][1]  # Extract motif name (e.g., gyne_MAxxxx or worker_MAxxxx)
      current_matrix <- matrix(nrow = 0, ncol = 0)  # Reset matrix for the next motif
      i <- i + 1; next  # Move to the next line
    }
    if (grepl("^[ACGT]\\s*\\[", line)) {  # If the line represents base frequencies
      parts <- strsplit(line, "\\s*\\[\\s*|\\s*\\]")[[1]]  # Split base and frequencies
      base  <- parts[1]  # Extract base (A/C/G/T)
      nums  <- suppressWarnings(as.numeric(strsplit(trimws(parts[2]), "\\s+")[[1]]))  # Convert frequency values to numeric
      if (is.null(current_matrix) || ncol(current_matrix) == 0)  # If matrix is not initialized
        current_matrix <- matrix(0, nrow = 0, ncol = length(nums))  # Initialize matrix
      if (!base %in% rownames(current_matrix)) {  # If the base is not present in the matrix
        current_matrix <- rbind(current_matrix, nums)  # Add new base row
        rownames(current_matrix)[nrow(current_matrix)] <- base  # Set row name
      } else {
        current_matrix[base, ] <- nums  # Update base row
      }
    }
    i <- i + 1  # Move to next line
  }
  # Final check at the end of the file
  if (!is.null(current_matrix) && ncol(current_matrix) > 0) {
    if (!all(c("A", "C", "G", "T") %in% rownames(current_matrix))) {
      miss <- setdiff(c("A", "C", "G", "T"), rownames(current_matrix))  # Ensure all bases are included
      for (b in miss) {
        current_matrix <- rbind(
          current_matrix,
          setNames(data.frame(t(rep(0, ncol(current_matrix)))), b)
        )
      }
    }
    current_matrix <- current_matrix[c("A", "C", "G", "T"), , drop = FALSE]  # Reorder bases to A/C/G/T
    motifs     <- c(motifs, list(as.matrix(current_matrix)))  # Add final motif to list
    labels     <- c(labels, current_name)  # Add label
    full_names <- c(full_names, current_full_name)  # Add full name
  }
  message(sprintf("✅ Successfully parsed %d motifs", length(motifs)))  # Message indicating the number of parsed motifs
  list(motifs = motifs, labels = labels, full_names = full_names)  # Return parsed motifs, labels, and full names
}


## -------------------------
# Function to extract features from a PFM matrix, including information content, GC content, AT content, etc.
extract_features <- function(mat) {
  col_sums <- colSums(mat)  # Sum of each column
  col_sums[col_sums <= 0] <- 1e-8  # Avoid division by zero
  prob <- sweep(mat, 2, col_sums, "/")  # Normalize by column sums
  prob <- pmin(pmax(prob, 1e-8), 1 - 1e-8)  # Ensure values stay within the range [1e-8, 1 - 1e-8]

  ic_pos <- 2 + colSums(prob * log2(prob))  # Information content for each position
  H2   <- -colSums(prob * log2(prob))  # Entropy
  cons <- 1 - H2 / 2  # Conservation score
  max_pref <- apply(prob, 2, max)  # Maximum base preference per position

  total <- sum(mat)  # Total sum of all base frequencies
  if (total <= 0) total <- 1e-8  # Prevent division by zero
  A <- sum(mat["A", ], na.rm = TRUE)  # Sum of A base frequencies
  C <- sum(mat["C", ], na.rm = TRUE)  # Sum of C base frequencies
  G <- sum(mat["G", ], na.rm = TRUE)  # Sum of G base frequencies
  T <- sum(mat["T", ], na.rm = TRUE)  # Sum of T base frequencies

  tibble(
    ic_mean      = mean(ic_pos, na.rm = TRUE),
    ic_sd        = sd0(ic_pos),  # Standard deviation
    ic_max       = max(ic_pos, na.rm = TRUE),
    motif_length = ncol(mat),
    pref_mean    = mean(max_pref, na.rm = TRUE),
    pref_sd      = sd0(max_pref),  # Standard deviation of base preference
    cons_mean    = mean(cons, na.rm = TRUE),
    cons_sd      = sd0(cons),  # Standard deviation of conservation score
    GC_content   = (G + C) / total,  # GC content
    A_content    = A / total,  # A content
    C_content    = C / total,  # C content
    G_content    = G / total,  # G content
    T_content    = T / total,  # T content
    AT_content   = (A + T) / total  # AT content
  ) |> mutate(across(everything(), ~replace_na(., 0)))  # Replace NA with 0
}

## -------------------------
## Build Dataset & Train Random Forest
## -------------------------
# Parse PFM file and prepare dataset
parsed  <- parse_pfm_file(INPUT_PFM)
dataset <- prepare_dataset(parsed$motifs, parsed$labels)
X <- dataset %>% select(-label, -tf_id)  # Features (excluding labels and tf_id)
y <- dataset$label  # Labels (gyne or worker)

# Handle missing values in features
if (anyNA(X)) X <- randomForest::na.roughfix(X)

# Train Random Forest classifier
set.seed(42)
rf_fit <- randomForest(x = X, y = y, ntree = 500, importance = TRUE)  # Train the model with 500 trees

## -------------------------
## RF Gini Importance
## -------------------------
# Extract and visualize feature importance (Gini)
gini <- as.data.frame(randomForest::importance(rf_fit))
gini$feature <- rownames(gini)
col <- intersect(names(gini), c("MeanDecreaseGini", "IncNodePurity"))[1]
gini <- gini[, c("feature", col)]
names(gini)[2] <- "MeanDecreaseGini"
gini <- gini[order(-gini$MeanDecreaseGini), ]

# Plot feature importance (Gini)
p_rf <- ggplot(gini, aes(x = reorder(feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "#D53E4F") +
  coord_flip() +
  labs(title = "RF Feature Importance", x = "Feature", y = "MeanDecreaseGini") +
  theme_bw(base_size = 12) +
  theme(plot.margin = ggplot2::margin(10, 20, 10, 10, unit = "pt"))

# Save the feature importance plot
ggsave(file.path(OUTDIR, "rf_gini_importance.pdf"),
       p_rf, width = 4.5, height = 4.2)

## -------------------------
## Permutation Importance (AUC)
## -------------------------
# Function to wrap predictions for ROC AUC calculation
pred_wrapper <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "gyne"]
}
metric_auc <- function(truth, estimate) {
  yardstick::roc_auc_vec(truth = truth, estimate = estimate, event_level = "first")
}

# Permutation importance calculation with AUC metric
set.seed(42)
vi_perm <- vip::vi_permute(
  object = rf_fit,
  feature_names = colnames(X),
  train = X,
  target = y,
  metric = metric_auc,
  pred_wrapper = pred_wrapper,
  nsim = 20,
  smaller_is_better = FALSE
)

# Plot and save permutation importance (AUC)
perm_df <- vi_perm %>%
  as.data.frame() %>%
  rename(feature = Variable, Importance = Importance) %>%
  arrange(desc(Importance))

p_perm <- ggplot(perm_df, aes(x = reorder(feature, Importance), y = Importance)) +
  geom_col(fill = "#D53E4F") +
  coord_flip() +
  labs(title = "Permutation Importance (ROC AUC)", x = "Feature", y = "Importance") +
  theme_bw(base_size = 12) +
  theme(plot.margin = ggplot2::margin(10, 20, 10, 10, unit = "pt"))

# Save the permutation importance plot
ggsave(file.path(OUTDIR, "perm_importance.pdf"),
       p_perm, width = 4.5, height = 4.2)

## -------------------------
## SHAP Calculation & Visualization
## -------------------------
# Compute SHAP values using the fastshap package
set.seed(42)
shap_vals <- fastshap::explain(
  rf_fit, X = X,
  pred_wrapper = pred_wrapper,
  nsim = 200, adjust = TRUE
)

# SHAP Feature Importance Plot
shap_imp <- tibble(
  feature = colnames(shap_vals),
  mean_abs_shap = colMeans(abs(shap_vals))
) %>% arrange(desc(mean_abs_shap))

p_shap_imp <- ggplot(shap_imp,
                     aes(x = reorder(feature, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(fill = "#D53E4F") +
  coord_flip() +
  labs(title = "SHAP Feature Importance", x = "Feature", y = "mean(|SHAP|)") +
  theme_bw(base_size = 12) +
  theme(plot.margin = ggplot2::margin(10, 20, 10, 10, unit = "pt"))

# Save SHAP feature importance plot
ggsave(file.path(OUTDIR, "shap_importance.pdf"),
       p_shap_imp, width = 4.5, height = 4.2)

# SHAP Bee Swarm Plot for Top 10 Features
shap_long <- shap_vals %>%
  as.data.frame() %>%
  mutate(id = 1:nrow(.)) %>%
  tidyr::pivot_longer(
    cols = -id,
    names_to = "feature",
    values_to = "shap_value"
  ) %>%
  left_join(
    X %>% mutate(id = 1:nrow(.)) %>%
      tidyr::pivot_longer(
        cols = -id,
        names_to = "feature",
        values_to = "feature_value"
      ),
    by = c("id", "feature")
  ) %>%
  mutate(
    true_label = y[id],   
    true_label = factor(true_label, levels = c("gyne", "worker"))
  )

p_shap_swarm_labelled <- ggplot(
  shap_long %>% filter(feature %in% top10_feats),
  aes(
    x = shap_value,
    y = reorder(feature, abs(shap_value), mean),
    color = true_label
  )) +
  geom_quasirandom(size = 1.5, alpha = 0.8, width = 0.3) +
  scale_color_manual(
    values = c("gyne" = "#D53E4F",
               "worker" = "#6bafda")
  ) +
  labs(
    title = "SHAP Bee Swarm (Top 10 Features)",
    x = "SHAP value (impact on P(gyne))",
    y = "Feature",
    color = "caste"
  ) +
  theme_bw(base_size = 12)

# Save SHAP Bee Swarm plot
ggsave(file.path(OUTDIR, "shap_beeswarm.pdf"),
       p_shap_swarm_labelled, width = 7, height = 5)
