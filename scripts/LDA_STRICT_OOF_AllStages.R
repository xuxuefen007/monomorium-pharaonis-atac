rm(list=ls()); gc()
options(stringsAsFactors = FALSE, warn = -1)

# ============================================================
# Purpose:
#   Stage-wise motif classification (Gyne vs Worker) using
#   strictly out-of-fold (OOF) Linear Discriminant Analysis (LDA).
#
# Input motif files (FILES):
#   These files were generated from a multi-step ATAC-seq pipeline:
#   1) DiffBind identifies caste-differential ATAC-seq peaks (Gyne vs Worker)
#      separately for each developmental stage.
#   2) AME (MEME Suite) performs motif enrichment on stage-specific differential peaks,
#      producing two enriched motif sets: Gyne-biased and Worker-biased.
#   3) For each stage, Gyne- and Worker-enriched motifs are merged into one file for
#      joint downstream analyses.
#
# File naming convention:
#   <stage>_G<Ng>W<Nw>_merged.txt
#     Ng = number of Gyne-biased enriched motifs
#     Nw = number of Worker-biased enriched motifs
#
# Motif format:
#   Each file is in standard JASPAR position frequency matrix (PFM) format.
#   Each motif starts with a header line, for example:
#     >worker_MA0013.1  br
#   followed by four rows (A, C, G, T) containing PFM counts/frequencies.
# ============================================================

# =========================
# Parameters
# =========================
BASE_DIR <- "C:/Users/Xuefen Xu/Documents/"
STAGES   <- c("larva_2nd", "larva_3rd", "prepupa", "young_pupa", "old_pupa", "adult")
FILES    <- c(
  larva_2nd   = file.path(BASE_DIR, "larva_2nd_G91W131_merged.txt"),
  larva_3rd   = file.path(BASE_DIR, "larva_3rd_G85W119_merged.txt"),
  prepupa     = file.path(BASE_DIR, "prepupa_G137W65_merged.txt"),
  young_pupa  = file.path(BASE_DIR, "young_pupa_G81W136_merged.txt"),
  old_pupa    = file.path(BASE_DIR, "old_pupa_G94W45_merged.txt"),
  adult       = file.path(BASE_DIR, "adult_G80W50_merged.txt")
)

SEED       <- 123
K_FOLDS    <- 5
N_PERM     <- 5000     
PARALLEL   <- TRUE
N_CORES    <- max(1, parallel::detectCores() - 1)

OUTDIR_ALL <- file.path(BASE_DIR, "LDA_STRICT_OOF_AllStages")
dir.create(OUTDIR_ALL, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

suppressPackageStartupMessages({
  library(MASS)
  library(tidyverse)
  library(caret)
  library(pROC)
  library(patchwork)
  library(doParallel)
})

# =========================
# parallel backend
# =========================
if (PARALLEL) {
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  message("Parallel backend registered: ", N_CORES, " cores.")
} else {
  registerDoSEQ()
}

on.exit({
  if (PARALLEL) try(stopCluster(cl), silent = TRUE)
}, add = TRUE)

# ============================================================
# Parse JASPAR PFM file into a list of 4xL matrices
# ============================================================
parse_pfm <- function(file_path){
  if (!file.exists(file_path)) stop("Input file not found: ", file_path)

  lines <- readLines(file_path, warn = FALSE)
  motifs <- list(); labels <- character()

  cur <- NULL; name <- ""

  for (ln in lines){
    s <- trimws(ln)
    if (s == "") next

    # New motif header
    if (startsWith(s, ">")){
      if (!is.null(cur) && ncol(cur) > 0){
        motifs <- c(motifs, list(cur))
        labels <- c(labels, name)
      }
      name <- sub("^>", "", s)
      cur  <- matrix(0, nrow = 4, ncol = 0)
      rownames(cur) <- c("A","C","G","T")
      next
    }

    # Typical JASPAR PFM line: A [ 1 2 3 ... ]
    if (grepl("^[ACGT]\\s*\\[", s)){
      parts <- strsplit(s, "\\s*\\[\\s*|\\s*\\]")[[1]]
      base  <- trimws(parts[1])
      nums  <- as.numeric(strsplit(trimws(parts[2]), "\\s+")[[1]])

    # Alternative line: A 1 2 3 ...
    } else if (grepl("^[ACGT]\\s+[-0-9\\.eE\\s]+$", s)){
      parts <- strsplit(s, "\\s+")[[1]]
      base  <- parts[1]
      nums  <- as.numeric(parts[-1])

    } else {
      next
    }

    if (ncol(cur) == 0){
      cur <- matrix(0, nrow = 4, ncol = length(nums))
      rownames(cur) <- c("A","C","G","T")
    }

    if (base %in% rownames(cur)){
      cur[base, 1:length(nums)] <- nums
    }
  }

  # Append last motif
  if (!is.null(cur) && ncol(cur) > 0){
    motifs <- c(motifs, list(cur))
    labels <- c(labels, name)
  }

  if (length(motifs) == 0) stop("No motifs parsed from file: ", file_path)
  list(motifs = motifs, labels = labels)
}

# ============================================================
# Feature extraction from a PFM matrix
#   - converts raw counts to probabilities per column
#   - computes information content, entropy-based conservation, base composition, etc.
# ============================================================
extract_features <- function(pfm_matrix){
  req <- c("A","C","G","T")

  # Ensure A/C/G/T rows exist
  if (!all(req %in% rownames(pfm_matrix))){
    for (b in setdiff(req, rownames(pfm_matrix))){
      pfm_matrix <- rbind(pfm_matrix, rep(0.001, ncol(pfm_matrix)))
      rownames(pfm_matrix)[nrow(pfm_matrix)] <- b
    }
    pfm_matrix <- pfm_matrix[req, , drop = FALSE]
  }

  pfm_matrix <- pfm_matrix + 1e-8
  col_sums   <- colSums(pfm_matrix)
  prob_mat   <- sweep(pfm_matrix, 2, col_sums, "/")
  prob_mat   <- pmin(pmax(prob_mat, 1e-8), 1 - 1e-8)

  # Per-position information content in bits: 2 + sum(p*log2 p)
  ic_pos <- 2 + colSums(prob_mat * log2(prob_mat), na.rm = TRUE)

  # Entropy-based "conservation" (scaled to [0,1])
  H_ln <- -colSums(prob_mat * log(prob_mat), na.rm = TRUE)
  cons <- 1 - H_ln / log(4)

  total <- sum(pfm_matrix) + 1e-8
  gc_w  <- sum(pfm_matrix[c("C","G"), , drop=FALSE]) / total

  A_w <- sum(pfm_matrix["A", , drop=FALSE]) / total
  C_w <- sum(pfm_matrix["C", , drop=FALSE]) / total
  G_w <- sum(pfm_matrix["G", , drop=FALSE]) / total
  T_w <- sum(pfm_matrix["T", , drop=FALSE]) / total

  c(
    ic_mean = mean(ic_pos),
    ic_sd   = sd(ic_pos),
    ic_max  = max(ic_pos),
    motif_length = ncol(pfm_matrix),
    total_count  = total,
    pref_mean = mean(apply(prob_mat, 2, max)),
    pref_sd   = sd(apply(prob_mat, 2, max)),
    cons_mean = mean(cons),
    cons_sd   = sd(cons),
    gc_content = gc_w,
    A_content = A_w, C_content = C_w, G_content = G_w, T_content = T_w
  )
}

# ============================================================
# Train/test preprocessing applied within each fold
#   - remove zero-variance features (training set)
#   - z-score standardize using training statistics
#   - remove highly correlated features (training set)
# ============================================================
prep_apply <- function(Xtr, Xte){
  sds   <- apply(Xtr, 2, sd)
  keep1 <- which(is.finite(sds) & sds > 0)
  Xtr   <- Xtr[, keep1, drop = FALSE]
  Xte   <- Xte[, keep1, drop = FALSE]

  mu  <- apply(Xtr, 2, mean)
  sdv <- apply(Xtr, 2, sd); sdv[sdv == 0] <- 1

  Xtr_s <- sweep(sweep(Xtr, 2, mu, "-"), 2, sdv, "/")
  Xte_s <- sweep(sweep(Xte, 2, mu, "-"), 2, sdv, "/")

  if (ncol(Xtr_s) > 1){
    cor_mat <- suppressWarnings(cor(Xtr_s))
    rm_idx  <- tryCatch(findCorrelation(cor_mat, cutoff = 0.95, exact = TRUE),
                        error = function(e) integer(0))
    if (length(rm_idx) > 0){
      keep2 <- setdiff(seq_len(ncol(Xtr_s)), rm_idx)
      Xtr_s <- Xtr_s[, keep2, drop = FALSE]
      Xte_s <- Xte_s[, keep2, drop = FALSE]
    }
  }
  list(Xtr = Xtr_s, Xte = Xte_s)
}

# ============================================================
# One fold training:
#   - Fit LDA on training fold
#   - Predict class/posterior on test fold (OOF)
#   - Compute LD1 and enforce a consistent orientation:
#       mean(LD1_gyne) > mean(LD1_worker)
#   - Return OOF predictions + absolute loadings as feature importance proxy
# ============================================================
train_one_fold <- function(tr_idx, te_idx, X, y){
  Xtr <- X[tr_idx, , drop = FALSE]
  Xte <- X[te_idx, , drop = FALSE]
  ytr <- y[tr_idx]
  yte <- y[te_idx]

  pp <- prep_apply(Xtr, Xte)

  df_tr <- data.frame(pp$Xtr, caste = ytr)
  fit   <- lda(caste ~ ., data = df_tr)

  pr  <- predict(fit, newdata = as.data.frame(pp$Xte))
  LD1 <- as.numeric(pr$x[, 1])

  # Force a consistent sign across folds (for density plot comparability)
  if (mean(LD1[yte == "gyne"]) < mean(LD1[yte == "worker"])) LD1 <- -LD1

  coef_vec <- rep(NA_real_, ncol(X)); names(coef_vec) <- colnames(X)
  sc <- fit$scaling[, 1]; names(sc) <- colnames(pp$Xtr)
  coef_vec[names(sc)] <- abs(sc)

  list(
    oof = tibble(row_id = te_idx,
                 caste = yte,
                 pred  = pr$class,
                 prob_gyne = pr$posterior[, "gyne"],
                 LD1 = LD1),
    coef = coef_vec
  )
}

# =========================
# Main loop across stages
# =========================
density_plots <- list()
summary_rows  <- list()

for (st in STAGES){
  message("▶ Stage: ", st)
  PFM_FILE <- FILES[[st]]
  OUTDIR   <- file.path(OUTDIR_ALL, st)
  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

  parsed   <- parse_pfm(PFM_FILE)
  feat_lst <- lapply(parsed$motifs, extract_features)

  X <- as.data.frame(do.call(rbind, feat_lst))
  tf_id <- parsed$labels

  # Label rule: motif headers contain "worker" or "gyne" prefix
  y <- ifelse(grepl("worker", tf_id, ignore.case = TRUE), "worker", "gyne")
  y <- factor(y, levels = c("gyne","worker"))

  n_gyne   <- sum(y == "gyne")
  n_worker <- sum(y == "worker")

  # Create stratified folds
  folds <- caret::createFolds(y, k = K_FOLDS, list = TRUE)

  # Run folds (sequential; safe & reproducible). 
  oof_list  <- vector("list", K_FOLDS)
  coef_list <- vector("list", K_FOLDS)
  fold_acc  <- numeric(K_FOLDS)

  for (k in seq_len(K_FOLDS)){
    te_idx <- folds[[k]]
    tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    res <- train_one_fold(tr_idx, te_idx, X, y)
    oof_list[[k]]  <- res$oof
    coef_list[[k]] <- res$coef
    fold_acc[k]    <- mean(res$oof$pred == res$oof$caste)
  }

  oof <- bind_rows(oof_list) %>% arrange(row_id)
  coef_mat <- do.call(rbind, coef_list); colnames(coef_mat) <- colnames(X)

  # ROC/AUC on OOF posterior
  roc_obj <- pROC::roc(response = oof$caste,
                       predictor = oof$prob_gyne,
                       levels = c("worker","gyne"),
                       direction = ">")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  if (auc_val < 0.5) auc_val <- 1 - auc_val

  # Cohen's d on LD1 separation
  gy <- oof$LD1[oof$caste == "gyne"]
  wk <- oof$LD1[oof$caste == "worker"]
  cohen_d <- (mean(gy) - mean(wk)) / sqrt((var(gy) + var(wk)) / 2)

  acc_all <- mean(oof$pred == oof$caste)

  # ---- Plot: OOF LD1 density
  p_density <- ggplot(oof, aes(x = LD1, fill = caste)) +
    geom_density(alpha = 0.6, adjust = 1.3, trim = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_fill_manual(values = c(gyne = "#E41A1C", worker = "#377EB8"),
                      name = "Caste",
                      labels = c("Gyne","Worker")) +
    labs(
      title = paste0("LDA OOF Density (", st, ")"),
      subtitle = sprintf("n(G)= %d | n(W)= %d | Acc= %.1f%% | AUC= %.3f | d= %.2f",
                         n_gyne, n_worker, 100 * acc_all, auc_val, cohen_d),
      x = "OOF LD1",
      y = "Density"
    ) +
    coord_cartesian(xlim = range(oof$LD1) + c(-1, 1)) +  # safer than xlim()
    theme_minimal(base_size = 13)

  density_plots[[st]] <- p_density

  ggsave(file.path(OUTDIR, "density_oof.pdf"), p_density, width = 8, height = 5)

  # ---- Save OOF table and fold accuracies
  write.csv(oof, file.path(OUTDIR, "oof_predictions.csv"), row.names = FALSE)
  write.csv(data.frame(fold = seq_len(K_FOLDS), acc = fold_acc),
            file.path(OUTDIR, "fold_accuracy.csv"), row.names = FALSE)

  # ---- Save mean absolute loadings across folds (feature importance proxy)
  coef_mean <- colMeans(abs(coef_mat), na.rm = TRUE)
  imp_df <- tibble(feature = names(coef_mean), mean_abs_loading = coef_mean) %>%
    arrange(desc(mean_abs_loading))
  write.csv(imp_df, file.path(OUTDIR, "feature_importance_mean_abs_loading.csv"),
            row.names = FALSE)

  # ---- Collect stage summary
  summary_rows[[st]] <- tibble(
    stage = st,
    n_gyne = n_gyne,
    n_worker = n_worker,
    acc = acc_all,
    auc = auc_val,
    cohen_d = cohen_d
  )
}

# =========================
# Combine all stages: one PDF + summary CSV
# =========================
summary_df <- bind_rows(summary_rows) %>% arrange(match(stage, STAGES))
write.csv(summary_df, file.path(OUTDIR_ALL, "summary_all_stages.csv"), row.names = FALSE)

p_density_all <- wrap_plots(density_plots[STAGES], ncol = 2) +
  plot_annotation(
    title = "LDA Density Distributions Across Developmental Stages",
    subtitle = "Gyne (red) vs Worker (blue); strictly out-of-fold (OOF) predictions"
  )

ggsave(file.path(OUTDIR_ALL, "density_all_stages_A4.pdf"),
       p_density_all, width = 11.7, height = 8.3)

message("Done. Outputs written to: ", OUTDIR_ALL)
