rm(list = ls()); gc()
options(stringsAsFactors = FALSE, warn = -1)

# ============================================================
# 0) Paths and global settings
# ============================================================
setwd("C:/Users/Xuefen Xu/Documents/tSNE")

FILE_PFM <- "C:/Users/Xuefen Xu/Documents/tSNE/gyne_worker_G75W140.txt"  # High-confidence motif set (75 gyne, 140 worker; boundary motifs removed)
OUTDIR   <- file.path(getwd(), "FULL_PIPELINE_OUT")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

SEED <- 123
set.seed(SEED)

# ============================================================
# 1) Packages
# ============================================================
suppressPackageStartupMessages({
  library(motifStack)
  library(TFBSTools)
  library(universalmotif)
  library(text2vec)
  library(Matrix)
  library(RSpectra)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(vegan)
  library(proxy)
})

# ============================================================
# 2) Utility functions
# ============================================================
col_IC <- function(p4) {
  eps <- 1e-9
  p4  <- pmax(p4, eps)
  sum(p4 * log2(p4 / 0.25))
}

pcm_to_pfm <- function(pcm) {
  mat <- as.matrix(pcm)
  cs  <- colSums(mat)
  cs[cs == 0] <- 1
  pfm <- sweep(mat, 2, cs, "/")
  rownames(pfm) <- c("A", "C", "G", "T")
  pfm
}

trim_lowIC <- function(pfm, thr = 0.2) {
  # Keep columns with sufficient information content; if none, keep at least one (the middle column)
  keep <- apply(pfm, 2, function(p4) col_IC(p4) >= thr)
  if (!any(keep)) {
    mid <- ceiling(ncol(pfm) / 2)
    pfm[, mid, drop = FALSE]
  } else {
    pfm[, keep, drop = FALSE]
  }
}

revcomp_pfm <- function(pfm) {
  rc <- pfm[c("T", "G", "C", "A"), , drop = FALSE]
  rc[, ncol(rc):1, drop = FALSE]
}

total_IC <- function(pfm) sum(apply(pfm, 2, col_IC))

infer_label <- function(nm) {
  if (grepl("^gyne", nm, ignore.case = TRUE)) "gyne"
  else if (grepl("^worker", nm, ignore.case = TRUE)) "worker"
  else "other"
}

prob_to_iupac <- function(p4, major_thr = 0.85, pair_thr = 0.25) {
  names(p4) <- c("A", "C", "G", "T")
  if (max(p4) >= major_thr) return(names(which.max(p4)))

  sel <- names(p4[p4 >= pair_thr])
  if (length(sel) == 0) sel <- names(sort(p4, decreasing = TRUE))[1:2]
  sel <- paste(sort(sel), collapse = "")

  map <- c(
    "A" = "A", "C" = "C", "G" = "G", "T" = "T",
    "AG" = "R", "CT" = "Y", "CG" = "S", "AT" = "W",
    "GT" = "K", "AC" = "M",
    "CGT" = "B", "AGT" = "D", "ACT" = "H", "ACG" = "V",
    "ACGT" = "N"
  )
  if (!is.na(map[sel])) map[sel] else "N"
}

# ============================================================
# 3) Read and preprocess motifs
# ============================================================
message("Reading motifs (JASPAR format): ", FILE_PFM)
motifs_raw <- motifStack::importMatrix(FILE_PFM, format = "jaspar")

motif_list <- lapply(seq_along(motifs_raw), function(i) {
  nm  <- names(motifs_raw)[i]
  pcm <- motifs_raw[[i]]@mat
  pfm <- pcm_to_pfm(pcm)

  # (1) Trim low-IC columns
  pfm <- trim_lowIC(pfm, 0.2)

  # (2) Choose strand with larger total information content
  pfm_rc <- revcomp_pfm(pfm)
  if (total_IC(pfm_rc) > total_IC(pfm)) pfm <- pfm_rc

  # (3) Convert each position to an IUPAC token and record IC
  L <- ncol(pfm)
  toks <- character(L)
  icv  <- numeric(L)
  for (j in seq_len(L)) {
    toks[j] <- prob_to_iupac(pfm[, j])
    icv[j]  <- col_IC(pfm[, j])
  }

  tibble(
    name   = nm,
    label  = infer_label(nm),
    len    = L,
    tokens = list(toks),
    ic     = list(icv),
    pfm    = list(pfm)
  )
})

motif_tbl <- bind_rows(motif_list) %>% filter(label %in% c("gyne", "worker"))

# Filter invalid motifs and enforce matrix format
valid_idx <- vapply(motif_tbl$pfm, function(m) {
  is.matrix(m) && nrow(m) == 4 && ncol(m) > 0
}, logical(1))

if (any(!valid_idx)) {
  message("Removing invalid motifs: ", sum(!valid_idx))
  motif_tbl <- motif_tbl[valid_idx, ]
}

motif_tbl$pfm <- lapply(motif_tbl$pfm, function(m) {
  m <- as.matrix(m)
  rownames(m) <- c("A", "C", "G", "T")
  m
})

stopifnot(nrow(motif_tbl) > 1)

# ============================================================
# 4) Build motif embeddings (GloVe; fallback to PPMI + SVD)
# ============================================================
it <- itoken(motif_tbl$tokens, progressbar = FALSE)
vocab <- create_vocabulary(it)
vectorizer <- vocab_vectorizer(vocab)
tcm <- create_tcm(it, vectorizer, skip_grams_window = 2, skip_grams_window_context = "symmetric")

V <- length(vocab$term)
EMBED_DIM <- min(16, max(2, V - 1))

word_vectors <- NULL
glove_ok <- FALSE

try({
  glove <- GlobalVectors$new(rank = EMBED_DIM, x_max = 5, learning_rate = 0.03)
  tcm_eps <- tcm
  if (length(tcm_eps@x) > 0) tcm_eps@x <- tcm_eps@x + 1e-8

  wv_main <- glove$fit_transform(tcm_eps, n_iter = 150, convergence_tol = 1e-4, n_threads = 0)
  wv_ctx  <- glove$components
  word_vectors <- wv_main + t(wv_ctx)
  glove_ok <- TRUE
}, silent = TRUE)

if (!glove_ok) {
  message("GloVe was unstable; using PPMI + truncated SVD fallback.")
  X  <- as(tcm, "dgCMatrix")
  rs <- rowSums(X)
  cs <- colSums(X)
  S  <- sum(X)

  PMI  <- log((X * S) / (rs %o% cs + 1e-9) + 1e-9)
  PPMI <- pmax(PMI, 0)

  sv <- RSpectra::svds(PPMI, k = EMBED_DIM)
  word_vectors <- sv$u %*% diag(sqrt(pmax(sv$d, 1e-9)))
  rownames(word_vectors) <- rownames(PPMI)
}

token_index <- rownames(word_vectors)
EMBED_DIM   <- ncol(word_vectors)

# IC-weighted average of token vectors per motif
motif_embed <- lapply(seq_len(nrow(motif_tbl)), function(i) {
  toks <- motif_tbl$tokens[[i]]
  ics  <- motif_tbl$ic[[i]]

  keep <- toks %in% token_index
  toks <- toks[keep]
  ics  <- ics[keep]

  if (length(toks) == 0) return(rep(0, EMBED_DIM))

  M <- word_vectors[toks, , drop = FALSE]
  w <- ics
  if (sum(w) <= 0) w <- rep(1, length(ics))
  w <- w / sum(w)
  colSums(M * w)
})

embed_mat <- do.call(rbind, motif_embed)
rownames(embed_mat) <- motif_tbl$name

# ============================================================
# 5) PWM alignment similarity and Mantel analysis
# ============================================================
revcomp_mat <- function(pfm) pfm[c("T", "G", "C", "A"), ncol(pfm):1, drop = FALSE]

colsim_pcc <- function(A, B, eps = 1e-8) {
  if (!is.matrix(A) || !is.matrix(B) || ncol(A) == 0 || ncol(B) == 0) {
    return(matrix(0, nrow = 1, ncol = 1))
  }
  n <- 4
  muA <- matrix(colMeans(A), n, ncol(A), byrow = TRUE)
  muB <- matrix(colMeans(B), n, ncol(B), byrow = TRUE)
  ZA <- A - muA
  ZB <- B - muB
  sdA <- sqrt(pmax(colSums(ZA^2) / (n - 1), eps))
  sdB <- sqrt(pmax(colSums(ZB^2) / (n - 1), eps))
  num <- t(ZA) %*% ZB
  den <- (n - 1) * (sdA %o% sdB)
  S <- num / den
  S[!is.finite(S)] <- 0
  S
}

sw_pwm_score <- function(A, B, gap_pen = 0.3, min_overlap = 5) {
  if (!is.matrix(A) || !is.matrix(B) || ncol(A) == 0 || ncol(B) == 0) {
    return(list(score = 0, aln_len = 0))
  }

  score_one <- function(A, B) {
    li <- ncol(A); lj <- ncol(B)
    S <- colsim_pcc(A, B)
    H <- matrix(0, li, lj)
    L <- matrix(0, li, lj)
    best <- 0; bestL <- 0

    for (i in 1:li) for (j in 1:lj) {
      d <- if (i > 1 && j > 1) H[i - 1, j - 1] + S[i, j] else max(0, S[i, j])
      u <- if (i > 1) H[i - 1, j] - gap_pen else -Inf
      l <- if (j > 1) H[i, j - 1] - gap_pen else -Inf
      h <- max(0, d, u, l)
      H[i, j] <- h

      L[i, j] <- if (h == 0) 0 else if (h == d) (if (i > 1 && j > 1) L[i - 1, j - 1] else 0) + 1
      else if (h == u) (if (i > 1) L[i - 1, j] else 0) else (if (j > 1) L[i, j - 1] else 0)

      if (h > best) { best <- h; bestL <- L[i, j] }
    }

    norm <- max(1, min(li, lj))
    if (bestL < min_overlap) best <- 0
    list(score = max(0, best) / norm, aln_len = ifelse(best == 0, 0, bestL))
  }

  r1 <- score_one(A, B)
  r2 <- score_one(A, revcomp_mat(B))
  if (r2$score > r1$score) r2 else r1
}

N <- nrow(motif_tbl)
pfm_list <- motif_tbl$pfm
nm <- motif_tbl$name

sim_align <- matrix(0, N, N, dimnames = list(nm, nm))
diag(sim_align) <- 1

for (i in 2:N) {
  for (j in 1:(i - 1)) {
    sc <- sw_pwm_score(pfm_list[[i]], pfm_list[[j]])$score
    sim_align[i, j] <- sim_align[j, i] <- sc
  }
}

D_align <- 1 - sim_align
diag(D_align) <- 0

D_emb <- as.matrix(proxy::dist(embed_mat, method = "cosine"))

# ============================================================
# 6) Mantel test (embedding distance vs alignment distance)
# ============================================================
set.seed(123)
mant <- mantel(as.dist(D_emb), as.dist(D_align), permutations = 9999)

summary_align <- sprintf(
  "Mean alignment distance = %.3f | Mantel r = %.3f, p = %.3g",
  mean(D_align[upper.tri(D_align)]),
  mant$statistic,
  mant$signif
)

writeLines(summary_align)
writeLines(summary_align, con = file.path(OUTDIR, "summary.txt"))

# ============================================================
# 7) PERMANOVA + ANOSIM + beta-dispersion (alignment distances)
# ============================================================
df_meta <- data.frame(label = motif_tbl$label, row.names = motif_tbl$name)

set.seed(123)
perma   <- adonis2(as.dist(D_align) ~ label, data = df_meta, permutations = 9999, by = "margin")
ano     <- anosim(D_align, grouping = df_meta$label, permutations = 9999)
bd      <- betadisper(as.dist(D_align), group = df_meta$label)
bd_perm <- permutest(bd, permutations = 9999)

sink(file.path(OUTDIR, "permanova_anosim_results.txt"))
cat("### PERMANOVA ###\n"); print(perma); cat("\n")
cat("### ANOSIM ###\n");    print(ano);   cat("\n")
cat("### Beta-dispersion ###\n"); print(bd_perm); cat("\n")
sink(NULL)

# ============================================================
# 8) Visualizations (Mantel scatter removed on request)
# ============================================================

# 8.1 NMDS on alignment distances
set.seed(123)
nmds <- metaMDS(D_align, k = 2, trymax = 50, autotransform = FALSE, trace = FALSE)

df_nmds <- as.data.frame(nmds$points)
colnames(df_nmds) <- c("NMDS1", "NMDS2")
df_nmds$label <- motif_tbl$label

p_nmds <- ggplot(df_nmds, aes(NMDS1, NMDS2, color = label)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(aes(fill = label), geom = "polygon", alpha = 0.15, show.legend = FALSE) +
  theme_classic() +
  labs(
    title = "NMDS on alignment distance",
    subtitle = paste("Stress =", round(nmds$stress, 3))
  )

ggsave(file.path(OUTDIR, "NMDS.pdf"), p_nmds, width = 6, height = 5, dpi = 600)

# 8.2 Distance distribution (within vs between groups)
wb <- function(Dmat, labels) {
  lab <- as.factor(labels)
  ut <- upper.tri(Dmat, diag = FALSE)
  same <- lab[row(Dmat)[ut]] == lab[col(Dmat)[ut]]
  list(within = Dmat[ut][same], between = Dmat[ut][!same])
}

wb_al <- wb(D_align, motif_tbl$label)

df_dist <- rbind(
  data.frame(Distance = wb_al$within,  Type = "Within"),
  data.frame(Distance = wb_al$between, Type = "Between")
)

p_den <- ggplot(df_dist, aes(Distance, fill = Type)) +
  geom_density(alpha = 0.6, color = NA) +
  geom_vline(
    data = df_dist %>% group_by(Type) %>% summarize(m = mean(Distance), .groups = "drop"),
    aes(xintercept = m, color = Type),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("Within" = "#27AE60", "Between" = "#E74C3C")) +
  theme_classic() +
  labs(title = "Distribution of alignment distances", subtitle = summary_align)

ggsave(file.path(OUTDIR, "Distance_density.pdf"), p_den, width = 6, height = 5, dpi = 600)

# 8.3 ECDF
p_ecdf <- ggplot(df_dist, aes(Distance, color = Type)) +
  stat_ecdf(linewidth = 1.2) +
  scale_color_manual(values = c("Within" = "#27AE60", "Between" = "#E74C3C")) +
  theme_classic() +
  labs(title = "ECDF of alignment distances", subtitle = summary_align)

ggsave(file.path(OUTDIR, "Distance_ecdf.pdf"), p_ecdf, width = 6, height = 5, dpi = 600)
