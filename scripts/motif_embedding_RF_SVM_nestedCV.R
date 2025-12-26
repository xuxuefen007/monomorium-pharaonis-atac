rm(list = ls()); gc()
options(stringsAsFactors = FALSE, warn = -1, error = NULL)

# ============================================================
# 0) Paths and global configuration
# ============================================================
setwd("C:/Users/Xuefen Xu/Documents/tSNE")

FILE_PFM <- "C:/Users/Xuefen Xu/Documents/tSNE/gyne_worker_G75W140.txt"  # High-confidence motif set (75 gyne, 140 worker; boundary motifs removed)
OUTDIR   <- file.path(getwd(), "SUPERVISED_NESTED_CV_ALL")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

SEED    <- 123
K_OUT   <- 5    # Outer folds
K_IN    <- 3    # Inner folds
RPT_OUT <- 3    # Increase to 5 for the final run
EMBEDDER_METHOD <- "glove"   # Try GloVe first; fallback to PPMI+SVD if unstable

# Reproducibility: lock RNG and restrict BLAS/OpenMP threads
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(SEED)
Sys.setenv(
  "OMP_NUM_THREADS" = "1",
  "MKL_NUM_THREADS" = "1",
  "OPENBLAS_NUM_THREADS" = "1",
  "VECLIB_MAXIMUM_THREADS" = "1",
  "NUMEXPR_NUM_THREADS" = "1"
)

# ============================================================
# 1) Dependencies
# ============================================================
pkgs <- c(
  "motifStack","TFBSTools","text2vec","Matrix","RSpectra",
  "readr","e1071","randomForest","pROC",
  "tibble","dplyr","purrr","yaml"
)
new_packages <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_packages)) install.packages(new_packages, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({ lapply(pkgs, library, character.only = TRUE) })

# ============================================================
# 2) Helper functions
# ============================================================
col_IC <- function(p4){
  eps <- 1e-12
  p4  <- pmax(p4, eps)
  sum(p4 * log2(p4 / 0.25))
}

pcm_to_pfm <- function(pcm){
  mat <- as.matrix(pcm)
  cs  <- colSums(mat)
  cs[cs == 0] <- 1
  sweep(mat, 2, cs, "/")
}

trim_lowIC <- function(pfm, thr = 0.2){
  keep <- apply(pfm, 2, function(p4) col_IC(p4) >= thr)
  if (!any(keep)) pfm else pfm[, keep, drop = FALSE]
}

revcomp_pfm <- function(pfm){
  rc <- pfm[c("T","G","C","A"), , drop = FALSE]
  rc[, ncol(rc):1, drop = FALSE]
}

infer_label <- function(nm){
  if (grepl("^gyne", nm, ignore.case = TRUE)) "gyne"
  else if (grepl("^worker", nm, ignore.case = TRUE)) "worker"
  else "other"
}

prob_to_iupac <- function(p4, major_thr = 0.85, pair_thr = 0.25){
  names(p4) <- c("A","C","G","T")
  if (max(p4) >= major_thr) return(names(which.max(p4)))

  sel <- names(p4[p4 >= pair_thr])
  if (length(sel) == 0) sel <- names(sort(p4, decreasing = TRUE))[1:2]

  key <- paste(sort(sel), collapse = "")
  map <- c(
    "A"="A","C"="C","G"="G","T"="T",
    "AG"="R","CT"="Y","CG"="S","AT"="W","GT"="K","AC"="M",
    "CGT"="B","AGT"="D","ACT"="H","ACG"="V",
    "ACGT"="N"
  )
  if (!is.na(map[key])) return(unname(map[key]))
  "N"
}

# Robustly extract the raw matrix from motifStack/TFBSTools objects
get_mat_from_pf <- function(x){
  if (is.matrix(x)) return(x)
  m <- tryCatch(slot(x, "mat"), error = function(e) NULL)
  if (is.null(m)) m <- tryCatch(slot(x, "profileMatrix"), error = function(e) NULL)
  if (is.null(m)) m <- tryCatch(as.matrix(x), error = function(e) NULL)
  m
}

# Choose a canonical strand by lexicographic comparison of IUPAC strings
canon_orient <- function(pfm){
  rc <- revcomp_pfm(pfm)
  c1 <- paste(apply(pfm, 2, prob_to_iupac), collapse = "")
  c2 <- paste(apply(rc,  2, prob_to_iupac), collapse = "")
  if (c2 < c1) rc else pfm
}

# Stratified K-fold split returning indices per fold
make_stratified_folds <- function(y, K = 5, seed = 123){
  set.seed(seed)
  y <- factor(y)
  idx_list <- split(seq_along(y), y)

  folds <- vector("list", K)
  for (k in 1:K) folds[[k]] <- integer(0)

  for (cls in names(idx_list)){
    idx <- sample(idx_list[[cls]])
    parts <- split(idx, cut(seq_along(idx), breaks = K, labels = FALSE))
    for (k in 1:K) folds[[k]] <- c(folds[[k]], parts[[k]])
  }
  folds
}

# Safe extraction of positive-class probabilities from a prob matrix
safe_prob_col <- function(prob_mat, pos = "gyne"){
  if (!is.null(colnames(prob_mat)) && pos %in% colnames(prob_mat)) return(prob_mat[, pos])
  other <- setdiff(colnames(prob_mat), pos)
  if (length(other) >= 1) return(1 - prob_mat[, other[1]])
  as.numeric(prob_mat)
}

# -------------------- Embedding (GloVe preferred; fallback to PPMI+SVD) --------------------
fit_embedder <- function(tokens_list, window = 2, embed_dim_max = 16, method = EMBEDDER_METHOD, seed = SEED){
  set.seed(seed)

  it <- text2vec::itoken(tokens_list, progressbar = FALSE)
  vocab <- text2vec::create_vocabulary(it)
  vectorizer <- text2vec::vocab_vectorizer(vocab)
  tcm <- text2vec::create_tcm(it, vectorizer, skip_grams_window = window)

  V <- ncol(tcm)
  embed_dim <- min(embed_dim_max, max(2, V - 1))

  if (tolower(method) == "glove"){
    for (lr in c(0.05, 0.02, 0.01, 0.005, 0.002)){
      glove <- text2vec::GlobalVectors$new(rank = embed_dim, x_max = 5, learning_rate = lr)
      tcm_eps <- tcm
      if (length(tcm_eps@x) > 0) tcm_eps@x <- tcm_eps@x + 1e-8

      W <- tryCatch({
        wv_main <- glove$fit_transform(tcm_eps, n_iter = 100, convergence_tol = 1e-4, n_threads = 1)
        wv_ctx  <- glove$components
        W <- wv_main + t(wv_ctx)
        rownames(W) <- rownames(tcm)
        W
      }, error = function(e) NULL)

      if (!is.null(W)){
        return(list(token_vectors = W, token_index = rownames(W), embed_dim = ncol(W)))
      }
    }
  }

  # Fallback: PPMI + truncated SVD
  X  <- as(tcm, "dgCMatrix")
  rs <- Matrix::rowSums(X)
  cs <- Matrix::colSums(X)
  S  <- sum(X)

  PMI  <- as.matrix(log((X * S + 1e-12) / ((rs %o% cs) + 1e-12)))
  PPMI <- pmax(PMI, 0)

  k_eff <- min(embed_dim, max(1, min(nrow(PPMI), ncol(PPMI)) - 1))
  sv <- RSpectra::svds(PPMI, k = k_eff)

  W <- sv$u %*% diag(sqrt(pmax(sv$d, 1e-12)))
  rownames(W) <- rownames(PPMI)
  list(token_vectors = W, token_index = rownames(W), embed_dim = ncol(W))
}

# Transform motif token sequences into motif-level vectors (IC-weighted mean)
transform_embedder <- function(embedder, tokens_list, ic_list = NULL){
  token_index <- embedder$token_index
  D <- embedder$embed_dim

  out <- lapply(seq_along(tokens_list), function(i){
    toks <- tokens_list[[i]]
    ics  <- if (!is.null(ic_list)) ic_list[[i]] else rep(1, length(toks))

    keep <- toks %in% token_index
    toks <- toks[keep]
    ics  <- ics[keep]

    if (length(toks) == 0) return(rep(0, D))

    M <- embedder$token_vectors[toks, , drop = FALSE]
    w <- ics / sum(ics)
    vec <- as.numeric(colSums(M * w))
    vec[!is.finite(vec)] <- 0
    vec
  })

  do.call(rbind, out)
}

# ============================================================
# 3) Read motifs and build tokens/IC tracks
# ============================================================
message("Reading motifs (JASPAR-like format): ", FILE_PFM)

pf_list <- motifStack::importMatrix(FILE_PFM, format = "jaspar")

motif_list <- lapply(seq_along(pf_list), function(i){
  nm  <- names(pf_list)[i]
  mat <- get_mat_from_pf(pf_list[[i]])

  if (is.null(rownames(mat)) && nrow(mat) == 4) rownames(mat) <- c("A","C","G","T")

  pfm <- pcm_to_pfm(mat)
  pfm <- trim_lowIC(pfm, 0.2)
  pfm <- canon_orient(pfm)

  L <- ncol(pfm)
  toks <- character(L)
  icv  <- numeric(L)

  for (j in seq_len(L)){
    p4 <- pfm[, j]
    toks[j] <- prob_to_iupac(p4)
    icv[j]  <- col_IC(p4)
  }

  list(name = nm, label = infer_label(nm), tokens = toks, ic = icv)
})

motif_tbl <- do.call(rbind, lapply(motif_list, function(x) data.frame(name = x$name, label = x$label)))
motif_tokens <- lapply(motif_list, `[[`, "tokens")
motif_ic     <- lapply(motif_list, `[[`, "ic")

# IMPORTANT: factor level order is kept exactly as in your original script
y_all <- factor(motif_tbl$label, levels = c("worker","gyne"))

# ============================================================
# 4) SVM nested CV
# ============================================================
run_svm_nestedCV <- function(){
  pred_rows_all <- list()
  rpt_seeds <- sample(1:1e6, RPT_OUT)

  for (r in seq_len(RPT_OUT)){
    outer_splits <- make_stratified_folds(y_all, K = K_OUT, seed = rpt_seeds[r])

    for (k in seq_len(K_OUT)){
      message(sprintf("[SVM Outer] rep=%d/%d, fold=%d/%d", r, RPT_OUT, k, K_OUT))

      test_idx  <- outer_splits[[k]]
      train_idx <- setdiff(seq_along(y_all), test_idx)

      embedder <- fit_embedder(motif_tokens[train_idx], method = EMBEDDER_METHOD)
      X_tr <- transform_embedder(embedder, motif_tokens[train_idx], motif_ic[train_idx])
      X_te <- transform_embedder(embedder, motif_tokens[test_idx],  motif_ic[test_idx])

      y_tr <- y_all[train_idx]
      y_te <- y_all[test_idx]

      # Inner CV for tuning C (linear SVM)
      inner_splits <- make_stratified_folds(y_tr, K = K_IN, seed = SEED)
      gridC <- c(0.1, 1, 10)

      cv_acc <- sapply(gridC, function(Cval){
        accs <- c()
        for (j in 1:K_IN){
          val_idx <- inner_splits[[j]]
          tr_idx  <- setdiff(seq_along(y_tr), val_idx)

          fit <- svm(
            x = X_tr[tr_idx, ], y = y_tr[tr_idx],
            type = "C-classification", kernel = "linear", cost = Cval,
            probability = TRUE, scale = TRUE
          )
          pred <- predict(fit, X_tr[val_idx, ])
          accs <- c(accs, mean(pred == y_tr[val_idx]))
        }
        mean(accs)
      })

      bestC <- gridC[which.max(cv_acc)]

      model <- svm(
        x = X_tr, y = y_tr,
        type = "C-classification", kernel = "linear",
        cost = bestC, probability = TRUE
      )
      pred <- predict(model, X_te, probability = TRUE)
      prob <- attr(pred, "probabilities")[, "gyne"]

      pred_rows_all[[length(pred_rows_all) + 1]] <- data.frame(
        model = "SVM", rep_id = r, fold_id = k,
        name = motif_tbl$name[test_idx],
        true_label = as.character(y_te),
        Prob = prob, Pred = as.character(pred),
        bestC = bestC
      )
    }
  }

  dplyr::bind_rows(pred_rows_all)
}

# ============================================================
# 5) RF nested CV
# ============================================================
run_rf_nestedCV <- function(){
  param_grid <- expand.grid(
    ntree = c(200, 500),
    mtry_factor = c(1/3, 1/2, NA),
    stringsAsFactors = FALSE
  )

  pred_rows_all <- list()
  rpt_seeds <- sample(1:1e6, RPT_OUT)

  for (r in seq_len(RPT_OUT)){
    outer_splits <- make_stratified_folds(y_all, K = K_OUT, seed = rpt_seeds[r])

    for (k in seq_len(K_OUT)){
      message(sprintf("[RF Outer] rep=%d/%d, fold=%d/%d", r, RPT_OUT, k, K_OUT))

      test_idx  <- outer_splits[[k]]
      train_idx <- setdiff(seq_along(y_all), test_idx)

      embedder <- fit_embedder(motif_tokens[train_idx], method = EMBEDDER_METHOD)

      X_tr_all <- transform_embedder(embedder, motif_tokens[train_idx], motif_ic[train_idx])
      X_te     <- transform_embedder(embedder, motif_tokens[test_idx],  motif_ic[test_idx])

      y_tr_all <- y_all[train_idx]
      y_te     <- y_all[test_idx]

      # Inner CV for selecting (ntree, mtry)
      inner_splits <- make_stratified_folds(y_tr_all, K = K_IN, seed = SEED)
      perf <- list()

      for (pg in 1:nrow(param_grid)){
        nt <- param_grid$ntree[pg]
        mf <- param_grid$mtry_factor[pg]

        for (j in seq_along(inner_splits)){
          val_idx <- inner_splits[[j]]
          tr_idx  <- setdiff(seq_along(y_tr_all), val_idx)

          X_tr  <- X_tr_all[tr_idx,  , drop = FALSE]
          y_tr  <- y_tr_all[tr_idx]
          X_val <- X_tr_all[val_idx, , drop = FALSE]
          y_val <- y_tr_all[val_idx]

          p <- ncol(X_tr)
          mtry_val <- if (is.na(mf)) max(1, round(sqrt(p))) else max(1, round(p * mf))

          rf_model <- randomForest(x = X_tr, y = y_tr, ntree = nt, mtry = mtry_val)
          rf_prob  <- predict(rf_model, X_val, type = "prob")[, "gyne"]
          auc_val  <- as.numeric(pROC::auc(y_val, rf_prob))

          perf[[length(perf) + 1]] <- data.frame(ntree = nt, mtry = mtry_val, auc = auc_val)
        }
      }

      best <- dplyr::bind_rows(perf) %>%
        dplyr::group_by(ntree, mtry) %>%
        dplyr::summarise(mean_auc = mean(auc), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(mean_auc)) %>%
        dplyr::slice(1)

      rf_model <- randomForest(x = X_tr_all, y = y_tr_all, ntree = best$ntree, mtry = best$mtry)
      rf_prob  <- predict(rf_model, X_te, type = "prob")[, "gyne"]
      rf_pred  <- predict(rf_model, X_te)

      pred_rows_all[[length(pred_rows_all) + 1]] <- data.frame(
        model = "RF", rep_id = r, fold_id = k,
        name = motif_tbl$name[test_idx],
        true_label = as.character(y_te),
        Prob = rf_prob, Pred = as.character(rf_pred),
        best_ntree = best$ntree, best_mtry = best$mtry
      )
    }
  }

  dplyr::bind_rows(pred_rows_all)
}

# ============================================================
# 6) Run and save predictions
# ============================================================
res_svm <- run_svm_nestedCV()
res_rf  <- run_rf_nestedCV()

readr::write_csv(res_svm, file.path(OUTDIR, "SVM_nested_predictions.csv"))
readr::write_csv(res_rf,  file.path(OUTDIR, "RF_nested_predictions.csv"))

message("Completed. Results saved to: ", OUTDIR)

# ============================================================
# 7) Build per-sample misclassification list (aggregated over repeats/folds)
# ============================================================
make_misclassified_list <- function(pred_df, model_name){
  agg <- pred_df %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      true_label = dplyr::first(true_label),
      prob_mean  = mean(Prob, na.rm = TRUE),
      n_repeats  = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(final_pred = ifelse(prob_mean >= 0.5, "gyne", "worker"))

  mis <- agg %>% dplyr::filter(final_pred != true_label)
  out_file <- file.path(OUTDIR, paste0("PerSample_misclassified_", model_name, ".csv"))
  readr::write_csv(mis, out_file)

  message(sprintf("[%s] Misclassified motifs: %d (out of %d)", model_name, nrow(mis), nrow(agg)))
  mis
}

mis_svm <- make_misclassified_list(res_svm, "SVM")
mis_rf  <- make_misclassified_list(res_rf,  "RF")

# ============================================================
# 8) Console summary and per-motif summary export
# ============================================================
summarise_results <- function(pred_df, model_name){
  cat("\n==================== ", model_name, " ====================\n"); flush.console()

  # Counts
  N_total  <- length(unique(pred_df$name))
  N_gyne   <- sum(motif_tbl$label == "gyne")
  N_worker <- sum(motif_tbl$label == "worker")

  # Aggregate per motif by mean probability
  agg <- pred_df %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(
      true_label = dplyr::first(true_label),
      prob_mean  = mean(Prob, na.rm = TRUE),
      vote_pred  = ifelse(mean(Prob, na.rm = TRUE) >= 0.5, "gyne", "worker"),
      .groups = "drop"
    )

  n_correct <- sum(agg$vote_pred == agg$true_label)
  n_wrong   <- sum(agg$vote_pred != agg$true_label)

  # Confusion matrix
  y_true <- factor(agg$true_label, levels = c("worker","gyne"))
  y_pred <- factor(agg$vote_pred,  levels = c("worker","gyne"))
  cm <- table(Predicted = y_pred, Actual = y_true)

  TN <- cm["worker","worker"]; TP <- cm["gyne","gyne"]
  FN <- cm["worker","gyne"];   FP <- cm["gyne","worker"]

  acc  <- (TP + TN) / (TP + TN + FP + FN)
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  prec <- TP / (TP + FP)
  f1   <- ifelse((prec + sens) > 0, 2 * prec * sens / (prec + sens), NA)
  auc_val <- as.numeric(pROC::auc(y_true, agg$prob_mean))

  cat(sprintf("Motifs: %d | Gyne: %d | Worker: %d\n", N_total, N_gyne, N_worker)); flush.console()
  cat(sprintf("Accuracy = %.3f | AUC = %.3f\n", acc, auc_val)); flush.console()
  cat(sprintf("Sensitivity = %.3f | Specificity = %.3f\n", sens, spec)); flush.console()
  cat(sprintf("Precision = %.3f | F1 = %.3f\n", prec, f1)); flush.console()
  cat(sprintf("Correct motifs: %d | Wrong motifs: %d\n", n_correct, n_wrong)); flush.console()
  print(cm); flush.console()

  # Save per-motif summary
  readr::write_csv(agg, file.path(OUTDIR, paste0("PerMotif_summary_", model_name, ".csv")))

  list(
    model = model_name, N_total = N_total, N_gyne = N_gyne, N_worker = N_worker,
    ACC = round(acc, 3), AUC = round(auc_val, 3),
    Sensitivity = round(sens, 3), Specificity = round(spec, 3),
    Precision = round(prec, 3), F1 = round(f1, 3),
    n_correct = n_correct, n_wrong = n_wrong
  )
}

summary_svm <- summarise_results(res_svm, "SVM")
summary_rf  <- summarise_results(res_rf,  "RF")

# ============================================================
# 9) Probability density plots
# ============================================================
library(ggplot2)

plot_density <- function(pred_df, model_name){
  ggplot(pred_df, aes(x = Prob, fill = true_label)) +
    geom_density(alpha = 0.5, adjust = 1.2) +
    scale_fill_manual(values = c("worker" = "#79a2f2", "gyne" = "#E74C53")) +
    labs(
      title = paste("Predicted Probability Density -", model_name),
      x = "Predicted Probability (Gyne)",
      y = "Density"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

p1 <- plot_density(res_svm, "SVM")
p2 <- plot_density(res_rf,  "RF")

ggsave(file.path(OUTDIR, "Density_SVM.pdf"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(OUTDIR, "Density_RF.pdf"),  p2, width = 6, height = 5, dpi = 300)

message("Saved density plots to: ", OUTDIR)
