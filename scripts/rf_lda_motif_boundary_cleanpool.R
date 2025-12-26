rm(list = ls()); gc()
options(stringsAsFactors = FALSE)

# -------------------- Configuration --------------------
config <- list(
  INPUT_PFM   = "gyne96_worker159_motif.txt",  # Input PFM file containing gyne and worker motifs (see Table S3)
  OUTDIR      = file.path(getwd(), sprintf("rf_gw_eval_%s", format(Sys.time(), "%Y%m%d_%H%M%S"))),
  POS_CLASS   = "gyne",

  # Random Forest (RF) boundary criteria based on out-of-bag (OOB) estimates
  MARGIN_CUTOFF = 0.20,   # Probability-margin boundary: |P(gyne) - P(worker)| < 0.20
  SEP_CUTOFF    = 0.05,   # Proximity-separation boundary: (mean within-class prox - mean cross-class prox) < 0.05

  # LDA boundary criterion based on cross-validated posteriors
  LDA_CUTOFF    = 0.20,   # LDA boundary: |Post(gyne) - Post(worker)| < 0.20

  # TF-level aggregation consistency threshold:
  # if >= this fraction of motifs for a TF point to the same caste, label the TF as that caste
  TF_CONSISTENCY = 0.75,

  TOP_N_LDA    = 10,      # Number of top LDA features to plot
  cols_label   = c(gyne = "#d95f02", worker = "#1b9e77")
)

# -------------------- Packages --------------------
need_pkgs <- c("tidyverse","randomForest","MASS","ggplot2","ggdist","ggbeeswarm","readr")
invisible(lapply(need_pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

# Use explicit prefixes (e.g., dplyr::) where helpful to avoid namespace conflicts
RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(42)

# -------------------- Utilities --------------------
safe_write_csv <- function(df, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    utils::write.csv(df, path, row.names = row.names),
    error = function(e) {
      alt <- sub("\\.csv$", paste0("_", format(Sys.time(), "%H%M%S"), ".csv"), path)
      utils::write.csv(df, alt, row.names = row.names)
      message("Write failed; saved to: ", alt)
    }
  )
}

# Select a PDF device: prefer cairo_pdf when available; otherwise fall back to pdf
get_pdf_device <- function() {
  if (isTRUE(capabilities("cairo")) && "cairo_pdf" %in% getNamespaceExports("grDevices")) {
    return(grDevices::cairo_pdf)
  } else {
    message("cairo_pdf not available; using standard pdf device.")
    return(grDevices::pdf)
  }
}

# Ensure the PFM has rows A/C/G/T in that order; add missing bases as all-zero rows
normalize_pfm_matrix <- function(pfm_matrix) {
  if (!all(c("A","C","G","T") %in% rownames(pfm_matrix))) {
    miss <- setdiff(c("A","C","G","T"), rownames(pfm_matrix))
    for (b in miss) {
      pfm_matrix <- rbind(pfm_matrix, setNames(data.frame(t(rep(0, ncol(pfm_matrix)))), b))
    }
    pfm_matrix <- pfm_matrix[c("A","C","G","T"), , drop = FALSE]
  }
  as.matrix(pfm_matrix)
}

# Parse a PFM text file:
# - motif header lines start with ">"
# - matrix rows are "A [ ... ]", "C [ ... ]", etc.
parse_pfm_file <- function(file_path) {
  lines <- readLines(file_path)
  motifs <- list(); labels <- character(); full_names <- character()

  current_matrix <- NULL
  current_name <- ""
  current_full_name <- ""

  i <- 1
  while (i <= length(lines)) {
    line <- trimws(lines[i])

    if (startsWith(line, ">")) {
      # Flush previous motif
      if (!is.null(current_matrix) && ncol(current_matrix) > 0) {
        motifs <- c(motifs, list(normalize_pfm_matrix(current_matrix)))
        labels <- c(labels, current_name)
        full_names <- c(full_names, current_full_name)
      }

      head_line <- sub("^>", "", line)
      current_full_name <- head_line
      current_name <- strsplit(head_line, "\t")[[1]][1]
      current_matrix <- matrix(nrow = 0, ncol = 0)
      i <- i + 1
      next
    }

    if (grepl("^[ACGT]\\s*\\[", line)) {
      parts <- strsplit(line, "\\s*\\[\\s*|\\s*\\]")[[1]]
      base  <- trimws(parts[1])
      nums  <- as.numeric(strsplit(trimws(parts[2]), "\\s+")[[1]])

      if (is.null(current_matrix) || ncol(current_matrix) == 0) {
        current_matrix <- matrix(0, nrow = 0, ncol = length(nums))
      }

      if (!base %in% rownames(current_matrix)) {
        current_matrix <- rbind(current_matrix, nums)
        rownames(current_matrix)[nrow(current_matrix)] <- base
      } else {
        current_matrix[base, ] <- nums
      }
    }

    i <- i + 1
  }

  # Flush last motif
  if (!is.null(current_matrix) && ncol(current_matrix) > 0) {
    motifs <- c(motifs, list(normalize_pfm_matrix(current_matrix)))
    labels <- c(labels, current_name)
    full_names <- c(full_names, current_full_name)
  }

  message(sprintf("Parsed %d motifs successfully.", length(motifs)))
  list(motifs = motifs, labels = labels, full_names = full_names)
}

# -------------------- Feature extraction (weighted/global definitions) --------------------
extract_features <- function(pfm_matrix) {
  col_sums <- colSums(pfm_matrix) + 1e-8
  prob_matrix <- sweep(pfm_matrix, 2, col_sums, "/")
  prob_matrix <- pmin(pmax(prob_matrix, 1e-8), 1 - 1e-8)

  # Information content / entropy / conservation per position
  ic_per_position <- 2 + colSums(prob_matrix * log2(prob_matrix), na.rm = TRUE)  # 2 - H2
  H_ln <- -colSums(prob_matrix * log(prob_matrix), na.rm = TRUE)
  conservation <- 1 - H_ln / log(4)

  # Weighted/global base composition (weighted by total counts across all positions)
  total <- sum(pfm_matrix) + 1e-8
  gc_weighted <- (sum(pfm_matrix[c("C","G"), , drop = FALSE]) / total)
  at_weighted <- (sum(pfm_matrix[c("A","T"), , drop = FALSE]) / total)
  A_weighted  <- (sum(pfm_matrix["A", , drop = FALSE]) / total)
  C_weighted  <- (sum(pfm_matrix["C", , drop = FALSE]) / total)
  G_weighted  <- (sum(pfm_matrix["G", , drop = FALSE]) / total)
  T_weighted  <- (sum(pfm_matrix["T", , drop = FALSE]) / total)

  tibble::tibble(
    ic_mean = mean(ic_per_position),
    ic_sd   = stats::sd(ic_per_position),
    ic_max  = max(ic_per_position),
    motif_length = ncol(pfm_matrix),
    total_count  = sum(pfm_matrix),
    pref_mean = mean(apply(prob_matrix, 2, max)),
    pref_sd   = stats::sd(apply(prob_matrix, 2, max)),
    cons_mean = mean(conservation),
    cons_sd   = stats::sd(conservation),
    gc_content = gc_weighted,
    at_content = at_weighted,
    A_content  = A_weighted,
    C_content  = C_weighted,
    G_content  = G_weighted,
    T_content  = T_weighted
  )
}

prepare_dataset <- function(motifs, labels, full_names) {
  feats <- purrr::map(motifs, extract_features)

  target <- purrr::map_chr(labels, ~ {
    if (grepl("^worker", .x, ignore.case = TRUE)) "worker"
    else if (grepl("^gyne", .x, ignore.case = TRUE)) "gyne"
    else { message("Unrecognized label prefix: ", .x); NA_character_ }
  })

  df <- dplyr::bind_rows(feats)
  df$label <- factor(target, levels = c("gyne","worker"))
  df$tf_full_info <- full_names
  df$tf_id <- labels

  df <- df[!is.na(df$label), ]
  df <- df[, colSums(is.na(df)) == 0]

  message("Class distribution:")
  print(table(df$label))

  df
}

# -------------------- Random Forest: training + true OOB probability/boundaries --------------------
perform_rf_analysis <- function(x, y, outdir) {
  RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(42)

  rf_full <- randomForest::randomForest(
    x = x, y = y, ntree = 500,
    importance = TRUE,
    proximity  = TRUE,
    oob.prox   = TRUE,
    keep.inbag = TRUE
  )

  conf_oob <- rf_full$confusion[1:2, 1:2, drop = FALSE]
  safe_write_csv(as.data.frame(conf_oob), file.path(outdir, "rf_OOB_confusion_matrix.csv"))

  oob_acc  <- sum(diag(conf_oob)) / sum(conf_oob)
  oob_sens <- conf_oob["gyne","gyne"]     / sum(conf_oob["gyne",])
  oob_spec <- conf_oob["worker","worker"] / sum(conf_oob["worker",])
  oob_bal  <- mean(c(oob_sens, oob_spec))

  metrics_tbl <- tibble::tibble(
    metric = c("OOB_Accuracy","OOB_Sensitivity_gyne","OOB_Specificity_worker","OOB_Balanced_Accuracy"),
    value  = c(oob_acc, oob_sens, oob_spec, oob_bal)
  )
  safe_write_csv(metrics_tbl, file.path(outdir, "rf_OOB_metrics.csv"))

  imp <- as.data.frame(rf_full$importance)
  imp <- tibble::tibble(feature = rownames(imp), MeanDecreaseGini = imp$MeanDecreaseGini) |>
    dplyr::arrange(dplyr::desc(MeanDecreaseGini))
  safe_write_csv(imp, file.path(outdir, "rf_feature_importance.csv"))

  rf_full
}

calculate_rf_boundaries <- function(rf_model, x, y, dataset, outdir) {
  RNGkind("Mersenne-Twister","Inversion","Rejection"); set.seed(42)

  pa <- predict(rf_model, x, predict.all = TRUE, type = "response")
  indiv <- pa$individual
  inbag <- rf_model$inbag
  if (is.null(inbag)) stop("rf_model$inbag is NULL. Please train with keep.inbag=TRUE.")

  n <- nrow(indiv); Tn <- ncol(indiv)
  oob_p_gyne <- oob_p_worker <- numeric(n)
  oob_n_trees <- integer(n)
  pred_oob_cls <- character(n)

  # For each sample, compute OOB vote proportions only across trees where the sample was OOB
  for (i in 1:n) {
    oob_trees <- which(inbag[i, ] == 0)
    oob_n_trees[i] <- length(oob_trees)

    # Fallback: if no OOB trees (rare), use all trees
    if (length(oob_trees) == 0) oob_trees <- 1:Tn

    votes <- indiv[i, oob_trees]
    g <- sum(votes == "gyne")
    w <- sum(votes == "worker")
    s <- g + w; if (s == 0) s <- 1

    oob_p_gyne[i]   <- g / s
    oob_p_worker[i] <- w / s

    # Tie-breaking: assign to worker when g == w (i.e., "strict" behavior via g > w)
    pred_oob_cls[i] <- ifelse(g > w, "gyne", "worker")
  }

  correct    <- pred_oob_cls == as.character(y)
  oob_margin <- abs(oob_p_gyne - oob_p_worker)

  # Proximity-based separation: mean proximity to same-class minus mean proximity to other-class
  prox <- rf_model$proximity
  prox[is.na(prox)] <- 0

  sep_vec <- vapply(1:nrow(prox), function(ii) {
    same <- prox[ii, y == y[ii] & (1:nrow(prox)) != ii]
    othr <- prox[ii, y != y[ii]]
    m_same <- if (length(same) > 0) mean(same) else 0
    m_othr <- if (length(othr) > 0) mean(othr) else 0
    m_same - m_othr
  }, numeric(1))

  per_sample <- tibble::tibble(
    tf_id        = rownames(x),
    tf_full_info = dataset$tf_full_info[match(rownames(x), dataset$tf_id)],
    true_label   = as.character(y),
    pred_oob     = pred_oob_cls,
    correct      = correct,
    oob_p_gyne   = oob_p_gyne,
    oob_p_worker = oob_p_worker,
    oob_margin   = oob_margin,
    prox_sep     = sep_vec,
    oob_n_trees  = oob_n_trees
  ) |>
    dplyr::mutate(
      OOB_Boundary      = (oob_margin < config$MARGIN_CUTOFF | !correct),
      Prox_Boundary     = (prox_sep   < config$SEP_CUTOFF),
      boundary_union    = OOB_Boundary | Prox_Boundary,
      boundary_intersec = OOB_Boundary & Prox_Boundary
    ) |>
    dplyr::arrange(oob_margin, prox_sep)

  safe_write_csv(per_sample, file.path(outdir, "rf_per_sample_summary.csv"))
  per_sample
}

# -------------------- Cross-validated LDA (jackknife via MASS::lda CV=TRUE) --------------------
perform_lda_analysis_cv <- function(x, y, dataset, outdir) {
  x_scaled <- as.data.frame(scale(x))
  df_lda <- x_scaled
  df_lda$label <- y

  lda_fit <- MASS::lda(label ~ ., data = df_lda)             # For LD1 direction and coefficients
  lda_cv  <- MASS::lda(label ~ ., data = df_lda, CV = TRUE)  # Cross-validated prediction and posteriors

  ld1_in <- as.numeric(predict(lda_fit, newdata = df_lda)$x[, 1])

  # Enforce a consistent LD1 direction (worker -> gyne) using medians
  med_g <- stats::median(ld1_in[y == "gyne"])
  med_w <- stats::median(ld1_in[y == "worker"])
  flip  <- (med_g < med_w)

  lda_scores <- tibble::tibble(
    tf_id           = dataset$tf_id,
    tf_full_info    = dataset$tf_full_info,
    true_label      = as.character(y),
    lda_class       = as.character(lda_cv$class),
    lda_correct     = lda_cv$class == y,
    LD1             = if (flip) -ld1_in else ld1_in,
    posterior_gyne  = lda_cv$posterior[, "gyne"],
    posterior_worker= lda_cv$posterior[, "worker"]
  ) |>
    dplyr::mutate(
      posterior_diff = abs(posterior_gyne - posterior_worker),
      LDA_Boundary   = posterior_diff < config$LDA_CUTOFF
    )

  coef_tbl <- tibble::tibble(
    feature = colnames(x),
    LD1_coef = as.numeric(lda_fit$scaling[, 1]) * (if (flip) -1 else 1)
  ) |>
    dplyr::arrange(dplyr::desc(abs(LD1_coef))) |>
    dplyr::mutate(direction = ifelse(LD1_coef > 0, "gyne", "worker"))

  safe_write_csv(lda_scores, file.path(outdir, "lda_per_sample_summary.csv"))
  safe_write_csv(coef_tbl,  file.path(outdir, "lda_LD1_coefficients.csv"))

  lda_scores
}

# -------------------- Merge RF + LDA boundaries and derive a "clean" motif pool --------------------
merge_boundary_analysis <- function(rf_results, lda_results, outdir) {
  merged_tbl <- rf_results |>
    dplyr::rename(
      RF_Boundary_Union    = boundary_union,
      RF_Boundary_Intersec = boundary_intersec
    ) |>
    dplyr::left_join(
      dplyr::select(
        lda_results,
        tf_id, lda_class, lda_correct, LD1,
        posterior_gyne, posterior_worker, posterior_diff, LDA_Boundary
      ),
      by = "tf_id"
    ) |>
    dplyr::mutate(
      LDA_Boundary = dplyr::coalesce(LDA_Boundary, FALSE),
      boundary_union    = RF_Boundary_Union | LDA_Boundary,
      boundary_intersec = RF_Boundary_Intersec & LDA_Boundary,
      Boundary_Source = dplyr::case_when(
        RF_Boundary_Union & LDA_Boundary  ~ "Both",
        RF_Boundary_Union & !LDA_Boundary ~ "RF_only",
        !RF_Boundary_Union & LDA_Boundary ~ "LDA_only",
        TRUE                              ~ "None"
      )
    )

  safe_write_csv(merged_tbl, file.path(outdir, "rf_lda_merged_summary.csv"))

  # Clean pool definition:
  # - RF OOB prediction correct
  # - LDA CV prediction correct
  # - Not in the union of boundary calls
  clean_pool <- merged_tbl %>%
    dplyr::filter(correct == TRUE, lda_correct == TRUE, !boundary_union)

  gyne_final_motif   <- dplyr::filter(clean_pool, true_label == "gyne")
  worker_final_motif <- dplyr::filter(clean_pool, true_label == "worker")

  safe_write_csv(gyne_final_motif,   file.path(outdir, "gyne_final_motif_pool.csv"))
  safe_write_csv(worker_final_motif, file.path(outdir, "worker_final_motif_pool.csv"))

  list(
    merged = merged_tbl,
    clean_pool = clean_pool,
    gyne_final_motif = gyne_final_motif,
    worker_final_motif = worker_final_motif
  )
}

# -------------------- TF-level aggregation (direction based on RF OOB probabilities) --------------------
parse_tf_name <- function(x) {
  parts <- strsplit(x, "\t")[[1]]
  if (length(parts) >= 2) trimws(parts[2]) else ""
}

aggregate_to_TF <- function(clean_pool, tf_consistency = 0.75, outdir = getwd()) {
  mt <- clean_pool %>%
    dplyr::mutate(
      TF_name = vapply(tf_full_info, parse_tf_name, character(1)),
      pred_dir = ifelse(oob_p_gyne >= oob_p_worker, "gyne", "worker")
    )

  agg <- mt %>%
    dplyr::group_by(TF_name) %>%
    dplyr::summarise(
      n_motifs      = dplyr::n(),
      n_gyne_pred   = sum(pred_dir == "gyne"),
      n_worker_pred = sum(pred_dir == "worker"),
      gyne_frac     = n_gyne_pred / n_motifs,
      worker_frac   = n_worker_pred / n_motifs,
      mean_margin   = mean(oob_margin),
      min_margin    = min(oob_margin),
      mean_prox_sep = mean(prox_sep),
      lda_cv_acc    = mean(lda_correct),
      rf_oob_acc    = mean(correct),
      tf_ids        = paste(tf_id, collapse = ";"),
      .groups = "drop"
    )

  gy_tbl <- agg %>%
    dplyr::filter(gyne_frac >= tf_consistency) %>%
    dplyr::mutate(final_class = "gyne") %>%
    dplyr::arrange(dplyr::desc(gyne_frac), dplyr::desc(mean_margin), dplyr::desc(mean_prox_sep))

  wk_tbl <- agg %>%
    dplyr::filter(worker_frac >= tf_consistency) %>%
    dplyr::mutate(final_class = "worker") %>%
    dplyr::arrange(dplyr::desc(worker_frac), dplyr::desc(mean_margin), dplyr::desc(mean_prox_sep))

  safe_write_csv(gy_tbl, file.path(outdir, "gyne_final_TF_list_TFlevel.csv"))
  safe_write_csv(wk_tbl, file.path(outdir, "worker_final_TF_list_TFlevel.csv"))
  safe_write_csv(agg,    file.path(outdir, "TFlevel_all_candidates.csv"))

  list(gy = gy_tbl, wk = wk_tbl, agg = agg)
}

# -------------------- Visualizations (ECDF and raincloud) --------------------
create_visualizations <- function(lda_scores, outdir) {
  ld1_g <- lda_scores$LD1[lda_scores$true_label == "gyne"]
  ld1_w <- lda_scores$LD1[lda_scores$true_label == "worker"]
  ks_res <- suppressWarnings(stats::ks.test(ld1_g, ld1_w))

  p_ecdf <- ggplot2::ggplot(lda_scores, ggplot2::aes(LD1, color = true_label)) +
    ggplot2::stat_ecdf(size = 1.2, alpha = 0.8) +
    ggplot2::scale_color_manual(values = config$cols_label, name = "Label") +
    ggplot2::labs(
      title = "Empirical CDF of LD1",
      subtitle = sprintf("KS test: D = %.3f, p = %s",
                         ks_res$statistic, formatC(ks_res$p.value, format = "e", digits = 2)),
      x = "LD1 (worker -> gyne)",
      y = "Cumulative Probability"
    ) +
    ggplot2::theme_classic(base_size = 12)

  p_rain <- ggplot2::ggplot(lda_scores, ggplot2::aes(x = true_label, y = LD1, fill = true_label)) +
    ggdist::stat_halfeye(
      adjust = 0.7, width = 0.6, justification = -0.25,
      .width = c(0.5, 0.8, 0.95),
      point_colour = NA,
      slab_alpha = 0.6
    ) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95, fatten = 1.2) +
    ggbeeswarm::geom_quasirandom(ggplot2::aes(color = true_label), width = 0.10, alpha = 0.75, size = 1.6) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = config$cols_label, guide = "none") +
    ggplot2::scale_color_manual(values = config$cols_label, guide = "none") +
    ggplot2::labs(title = "LDA LD1 Distribution (Raincloud)", x = NULL, y = "LD1 (worker -> gyne)") +
    ggplot2::theme_classic(base_size = 12)

  devfun <- get_pdf_device()
  ggplot2::ggsave(file.path(outdir, "LD1_ECDF_with_KS.pdf"), p_ecdf, width = 7.0, height = 5.5, device = devfun)
  ggplot2::ggsave(file.path(outdir, "LDA_LD1_Raincloud.pdf"), p_rain, width = 7.2, height = 4.6, device = devfun)
}

# -------------------- Top-N LDA feature barplot --------------------
find_coef_file <- function() {
  if (exists("config") && is.list(config)) {
    f <- file.path(config$OUTDIR, "lda_LD1_coefficients.csv")
    if (file.exists(f)) return(f)
  }

  cand <- Sys.glob(file.path(getwd(), "rf_gw_eval_*", "lda_LD1_coefficients.csv"))
  if (length(cand) > 0) {
    cand <- cand[order(file.info(cand)$mtime, decreasing = TRUE)][1]
    return(cand)
  }

  if (file.exists("lda_LD1_coefficients.csv")) return("lda_LD1_coefficients.csv")

  stop("Cannot find lda_LD1_coefficients.csv. Please run the analysis first.")
}

plot_lda_top_features <- function(top_n = config$TOP_N_LDA, outdir = config$OUTDIR) {
  coef_file <- find_coef_file()
  message("Using LDA coefficient file: ", normalizePath(coef_file, winslash = "/"))

  coef_tbl <- tryCatch(
    readr::read_csv(coef_file, show_col_types = FALSE),
    error = function(e) utils::read.csv(coef_file, stringsAsFactors = FALSE, check.names = FALSE)
  )

  need_cols <- c("feature","LD1_coef","direction")
  if (!all(need_cols %in% colnames(coef_tbl))) {
    stop("Missing required columns in coefficient table: ",
         paste(setdiff(need_cols, colnames(coef_tbl)), collapse = ", "))
  }

  # Normalize direction to match palette keys
  coef_tbl$direction <- tolower(trimws(coef_tbl$direction))
  coef_tbl$direction[coef_tbl$direction %in% c("g","gy","gyn")] <- "gyne"
  coef_tbl$direction[coef_tbl$direction %in% c("w","wk","work")] <- "worker"
  bad <- setdiff(unique(coef_tbl$direction), c("gyne","worker"))
  if (length(bad) > 0) stop("Unexpected values in 'direction': ", paste(bad, collapse = ", "))
  coef_tbl$direction <- factor(coef_tbl$direction, levels = c("gyne","worker"))

  coef_tbl$abs_coef <- abs(coef_tbl$LD1_coef)
  top_feats <- coef_tbl[order(-coef_tbl$abs_coef), ][seq_len(min(top_n, nrow(coef_tbl))), ]

  fill_vals <- c(
    gyne   = unname(config$cols_label["gyne"]),
    worker = unname(config$cols_label["worker"])
  )

  p <- ggplot2::ggplot(
    top_feats,
    ggplot2::aes(
      x = stats::reorder(feature, abs_coef),
      y = LD1_coef,
      fill = direction
    )
  ) +
    ggplot2::geom_col(width = 0.7, alpha = 0.95) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fill_vals, breaks = c("gyne","worker"), drop = FALSE) +
    ggplot2::labs(
      title = paste("Top", nrow(top_feats), "LDA Features (|coef|)"),
      x = "Feature",
      y = "LD1 coefficient",
      fill = "Direction"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  devfun <- get_pdf_device()
  out_pdf <- file.path(outdir, sprintf("LDA_TopFeatures_top%d.pdf", top_n))
  ggplot2::ggsave(filename = out_pdf, plot = p, width = 7, height = 5, device = devfun)

  message("Saved Top-N LDA feature plot to: ", normalizePath(out_pdf, winslash = "/"))
}

# -------------------- Base composition visualizations --------------------
plot_base_composition <- function(dataset, outdir) {
  # Long format for base composition summaries
  base_data <- dataset %>%
    dplyr::select(tf_id, label, A_content, C_content, G_content, T_content, gc_content, at_content) %>%
    tidyr::pivot_longer(
      cols = c(A_content, C_content, G_content, T_content, gc_content, at_content),
      names_to = "base_type",
      values_to = "content"
    )

  # Fixed colors for base categories 
  base_colors <- c(
    A_content  = "#E41A1C",
    C_content  = "#377EB8",
    G_content  = "#4DAF4A",
    T_content  = "#984EA3",
    gc_content = "#FF7F00",
    at_content = "#A65628"
  )

  # Boxplot + jitter for base composition by caste
  p_base <- ggplot2::ggplot(base_data, ggplot2::aes(x = base_type, y = content, fill = base_type)) +
    ggplot2::geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(ggplot2::aes(color = base_type), width = 0.2, alpha = 0.6, size = 1) +
    ggplot2::facet_wrap(~label, ncol = 2) +
    ggplot2::scale_fill_manual(values = base_colors) +
    ggplot2::scale_color_manual(values = base_colors) +
    ggplot2::labs(
      title = "Base Composition by Caste",
      x = "Base Type",
      y = "Weighted Content"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Scatter plot: GC vs AT
  p_scatter <- ggplot2::ggplot(dataset, ggplot2::aes(x = gc_content, y = at_content, color = label)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::scale_color_manual(values = config$cols_label) +
    ggplot2::labs(
      title = "GC vs AT Content",
      x = "GC Content",
      y = "AT Content",
      color = "Caste"
    ) +
    ggplot2::theme_classic(base_size = 12)

  devfun <- get_pdf_device()
  ggplot2::ggsave(file.path(outdir, "base_composition_boxplot.pdf"), p_base, width = 10, height = 6, device = devfun)
  ggplot2::ggsave(file.path(outdir, "GC_vs_AT_scatter.pdf"), p_scatter, width = 7, height = 5, device = devfun)
}

# -------------------- Threshold sensitivity scan --------------------
threshold_sensitivity_scan <- function(
  rf_results, lda_results, outdir,
  margins = c(0.15, 0.20, 0.25),
  seps    = c(0.05, 0.07, 0.10),
  ldas    = c(0.15, 0.20, 0.25)
) {
  ans <- list()
  k <- 1

  for (m in margins) for (s in seps) for (l in ldas) {
    tmp_rf  <- rf_results %>% dplyr::mutate(
      OOB_Boundary      = (oob_margin < m | !correct),
      Prox_Boundary     = (prox_sep   < s),
      boundary_union    = OOB_Boundary | Prox_Boundary,
      boundary_intersec = OOB_Boundary & Prox_Boundary
    )

    tmp_lda <- lda_results %>%
      dplyr::mutate(LDA_Boundary = abs(posterior_gyne - posterior_worker) < l)

    tmp_merged <- tmp_rf %>%
      dplyr::rename(
        RF_Boundary_Union    = boundary_union,
        RF_Boundary_Intersec = boundary_intersec
      ) %>%
      dplyr::left_join(
        dplyr::select(
          tmp_lda,
          tf_id, lda_class, lda_correct, LD1,
          posterior_gyne, posterior_worker,
          posterior_diff, LDA_Boundary
        ),
        by = "tf_id"
      ) %>%
      dplyr::mutate(boundary_union = RF_Boundary_Union | LDA_Boundary)

    gy_n <- sum(tmp_merged$true_label == "gyne"   & !tmp_merged$boundary_union & tmp_merged$correct & tmp_merged$lda_correct)
    wk_n <- sum(tmp_merged$true_label == "worker" & !tmp_merged$boundary_union & tmp_merged$correct & tmp_merged$lda_correct)

    ans[[k]] <- data.frame(MARGIN = m, SEP = s, LDA = l, gyne = gy_n, worker = wk_n, total = gy_n + wk_n)
    k <- k + 1
  }

  out <- dplyr::bind_rows(ans)
  safe_write_csv(out, file.path(outdir, "threshold_sensitivity_grid.csv"))
  out
}

# -------------------- Main pipeline --------------------
main <- function() {
  if (!file.exists(config$INPUT_PFM)) stop("Input file not found: ", config$INPUT_PFM)

  dir.create(config$OUTDIR, recursive = TRUE, showWarnings = FALSE)
  message("Output directory: ", normalizePath(config$OUTDIR))

  # 1) Parse PFM and extract features
  parsed  <- parse_pfm_file(config$INPUT_PFM)
  dataset <- prepare_dataset(parsed$motifs, parsed$labels, parsed$full_names)

  x_all <- as.data.frame(dplyr::select(dataset, -label, -tf_full_info, -tf_id))
  y_all <- dataset$label
  rownames(x_all) <- dataset$tf_id

  # 2) Random Forest: true OOB probabilities and OOB proximity
  message("Running Random Forest ...")
  rf_model   <- perform_rf_analysis(x_all, y_all, config$OUTDIR)
  rf_results <- calculate_rf_boundaries(rf_model, x_all, y_all, dataset, config$OUTDIR)

  # 3) Cross-validated LDA: conservative posteriors
  message("Running cross-validated LDA ...")
  lda_results <- perform_lda_analysis_cv(x_all, y_all, dataset, config$OUTDIR)

  # 4) Merge RF and LDA boundaries; define clean motif pool
  message("Merging RF and LDA boundary calls ...")
  merged_out <- merge_boundary_analysis(rf_results, lda_results, config$OUTDIR)

  # 5) TF-level aggregation
  message("Aggregating to TF level ...")
  tf_agg <- aggregate_to_TF(
    merged_out$clean_pool,
    tf_consistency = config$TF_CONSISTENCY,
    outdir = config$OUTDIR
  )

  # 6) Visualizations
  message("Generating plots ...")
  create_visualizations(lda_results, config$OUTDIR)
  plot_lda_top_features(top_n = config$TOP_N_LDA, outdir = config$OUTDIR)
  plot_base_composition(dataset, config$OUTDIR)

  # 7) Sensitivity scan for thresholds
  message("Running threshold sensitivity scan ...")
  sensitivity <- threshold_sensitivity_scan(rf_results, lda_results, config$OUTDIR)

  # 8) Summary files
  stats_line <- sprintf(
    paste0(
      "Summary:\n",
      "- Clean motif pool: gyne=%d, worker=%d\n",
      "- TF-level (consistency >= %.2f): gyne=%d, worker=%d\n"
    ),
    nrow(merged_out$gyne_final_motif),
    nrow(merged_out$worker_final_motif),
    config$TF_CONSISTENCY,
    nrow(tf_agg$gy), nrow(tf_agg$wk)
  )

  writeLines(stats_line, con = file.path(config$OUTDIR, "summary_counts.txt"))
  writeLines(c(capture.output(sessionInfo())), con = file.path(config$OUTDIR, "sessionInfo.txt"))

  message("Done. Output directory: ", normalizePath(config$OUTDIR))

  invisible(list(
    rf_model       = rf_model,
    dataset        = dataset,
    x_all          = x_all,
    y_all          = y_all,
    rf_results     = rf_results,
    lda_results    = lda_results,
    merged_results = merged_out$merged,
    clean_pool     = merged_out$clean_pool,
    tf_agg         = tf_agg,
    sensitivity    = sensitivity
  ))
}

# -------------------- Run --------------------
res <- main()

# quick check of boundary counts under different LDA cutoffs
for (thr in c(0.15, 0.20, 0.25)) {
  tmp <- res$lda_results %>% dplyr::mutate(LDA_Boundary = abs(posterior_gyne - posterior_worker) < thr)
  cat("LDA cutoff =", thr, " -> number of boundary motifs =", sum(tmp$LDA_Boundary), "\n")
}
