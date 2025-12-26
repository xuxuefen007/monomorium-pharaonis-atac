# -------------------- Motif NLP embedding + tSNE/UMAP  --------------------
rm(list=ls()); gc()
options(stringsAsFactors = FALSE)

# ===== 0) Settings =====
setwd("C:/Users/Xuefen Xu/Documents/tSNE")

FILE_PFM <- "C:/Users/Xuefen Xu/Documents/tSNE/gyne_worker_G75W140.txt"  # High-confidence motif set after RF+LDA integration (75 gyne, 140 worker; boundary motifs removed)
SEED <- 123
set.seed(SEED)

OUTDIR_TSNE <- file.path(getwd(), "tSNE_Figures")
OUTDIR_UMAP <- file.path(getwd(), "UMAP_Figures")
dir.create(OUTDIR_TSNE, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTDIR_UMAP, showWarnings = FALSE, recursive = TRUE)

# ===== 1) Packages =====
suppressPackageStartupMessages({
  library(motifStack)
  library(TFBSTools)
  library(universalmotif)

  library(text2vec)
  library(Matrix)
  library(RSpectra)

  library(uwot)
  library(Rtsne)

  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(stringr)

  library(cluster)
  library(mclust)
  library(readr)

  library(vegan)
  library(proxy)
  library(ape)
})

# ===== 2) Utility functions =====
col_IC <- function(p4) {
  eps <- 1e-9
  p4 <- pmax(p4, eps)
  sum(p4 * log2(p4 / 0.25))
}

pcm_to_pfm <- function(pcm) {
  mat <- as.matrix(pcm)
  cs <- colSums(mat); cs[cs == 0] <- 1
  pfm <- sweep(mat, 2, cs, "/")
  rownames(pfm) <- c("A","C","G","T")
  pfm
}

trim_lowIC <- function(pfm, thr = 0.2) {
  keep <- apply(pfm, 2, function(p4) col_IC(p4) >= thr)
  if (!any(keep)) return(pfm)
  pfm[, keep, drop=FALSE]
}

revcomp_pfm <- function(pfm) {
  rc_rows <- pfm[c("T","G","C","A"), , drop=FALSE]
  rc_rows[, ncol(rc_rows):1, drop=FALSE]
}

total_IC <- function(pfm) sum(apply(pfm, 2, col_IC))

infer_label <- function(nm) {
  if (grepl("^gyne", nm, ignore.case = TRUE)) return("gyne")
  if (grepl("^worker", nm, ignore.case = TRUE)) return("worker")
  "other"
}

prob_to_iupac <- function(p4, major_thr = 0.85, pair_thr = 0.25) {
  names(p4) <- c("A","C","G","T")
  if (max(p4) >= major_thr) return(names(which.max(p4)))
  sel <- names(p4[p4 >= pair_thr])
  if (length(sel) == 0) sel <- names(sort(p4, decreasing = TRUE))[1:2]
  sel_sorted <- paste(sort(sel), collapse = "")

  map <- c(
    "A"="A","C"="C","G"="G","T"="T",
    "AG"="R","CT"="Y","CG"="S","AT"="W","GT"="K","AC"="M",
    "CGT"="B","AGT"="D","ACT"="H","ACG"="V",
    "ACGT"="N"
  )
  if (!is.na(map[sel_sorted])) return(unname(map[sel_sorted]))

  # Fallbacks
  if (all(c("A","C","G","T") %in% sel)) return("N")
  if (all(c("C","G","T") %in% sel)) return("B")
  if (all(c("A","G","T") %in% sel)) return("D")
  if (all(c("A","C","T") %in% sel)) return("H")
  if (all(c("A","C","G") %in% sel)) return("V")
  if (all(c("A","G") %in% sel)) return("R")
  if (all(c("C","T") %in% sel)) return("Y")
  if (all(c("C","G") %in% sel)) return("S")
  if (all(c("A","T") %in% sel)) return("W")
  if (all(c("G","T") %in% sel)) return("K")
  if (all(c("A","C") %in% sel)) return("M")
  names(which.max(p4))
}

purity <- function(true, cluster) {
  tab <- table(cluster, true)
  sum(apply(tab, 1, max)) / length(true)
}

theme_sci <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, color = "grey40"),
      axis.title = element_text(face = "bold", size = base_size + 1),
      axis.text = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      legend.key.size = unit(0.6, "lines"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.margin = unit(c(2, 2, 2, 2), "mm")
    )
}

# Create motif table with tokens + IC weights after PFM processing
build_motif_table <- function(FILE_PFM) {
  message("Reading motifs from: ", FILE_PFM)
  motifs_raw <- motifStack::importMatrix(FILE_PFM, format = "jaspar")
  stopifnot(length(motifs_raw) > 0)

  motif_list <- lapply(seq_along(motifs_raw), function(i) {
    nm  <- names(motifs_raw)[i]
    pcm <- motifs_raw[[i]]@mat
    pfm <- pcm_to_pfm(pcm)

    # Keep the same processing as your original script:
    # (1) trim low-IC columns; (2) choose strand with larger total IC
    pfm <- trim_lowIC(pfm, thr = 0.2)
    pfm_rc <- revcomp_pfm(pfm)
    if (total_IC(pfm_rc) > total_IC(pfm)) pfm <- pfm_rc

    L <- ncol(pfm)
    tokens  <- character(L)
    ic_vals <- numeric(L)

    for (j in seq_len(L)) {
      p4 <- pfm[, j]
      tokens[j]  <- prob_to_iupac(p4, major_thr = 0.85, pair_thr = 0.25)
      ic_vals[j] <- col_IC(p4)
    }

    tibble(
      name   = nm,
      label  = infer_label(nm),
      len    = L,
      tokens = list(tokens),
      ic     = list(ic_vals),
      pfm    = list(pfm)
    )
  })

  motif_tbl <- bind_rows(motif_list) %>% filter(label %in% c("gyne","worker"))
  stopifnot(nrow(motif_tbl) > 0)
  motif_tbl
}

# Build embeddings (GloVe; fallback PPMI+SVD)
build_embeddings <- function(motif_tbl, SEED = 123, window_size = 2, ic_weighted_mean = TRUE) {
  set.seed(SEED)

  it <- itoken(motif_tbl$tokens, progressbar = FALSE)
  vocab <- create_vocabulary(it)
  vectorizer <- vocab_vectorizer(vocab)

  tcm <- create_tcm(
    it, vectorizer,
    skip_grams_window = window_size,
    skip_grams_window_context = "symmetric"
  )

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
    message("GloVe fit was unstable; using PPMI + SVD fallback.")
    X  <- as(tcm, "dgCMatrix")
    rs <- Matrix::rowSums(X); cs <- Matrix::colSums(X); S <- sum(X)
    PMI  <- as.matrix(log((X * S) / (rs %o% cs + 1e-9) + 1e-9))
    PPMI <- pmax(PMI, 0)
    sv   <- RSpectra::svds(PPMI, k = EMBED_DIM)
    word_vectors <- sv$u %*% diag(sqrt(pmax(sv$d, 1e-9)))
    rownames(word_vectors) <- rownames(PPMI)
  }

  token_index <- rownames(word_vectors)
  EMBED_DIM   <- ncol(word_vectors)

  motif_embed <- lapply(seq_len(nrow(motif_tbl)), function(i) {
    toks <- motif_tbl$tokens[[i]]
    ics  <- motif_tbl$ic[[i]]

    keep <- toks %in% token_index
    toks <- toks[keep]; ics <- ics[keep]
    if (length(toks) == 0) return(rep(0, EMBED_DIM))

    M <- word_vectors[toks, , drop = FALSE]
    if (ic_weighted_mean) {
      w <- ics; if (sum(w) <= 0) w <- rep(1, length(ics))
      w <- w / sum(w)
      colSums(M * w)
    } else {
      colMeans(M)
    }
  })

  embed_mat <- do.call(rbind, motif_embed)
  rownames(embed_mat) <- motif_tbl$name
  embed_mat
}

compute_clustering_metrics <- function(embed_mat, motif_tbl, cluster_id) {
  D_for_sil <- dist(embed_mat)
  sil <- cluster::silhouette(cluster_id, D_for_sil)
  sil_score <- mean(sil[, 3])

  true_y <- as.integer(motif_tbl$label == "gyne") + 1
  ari1 <- mclust::adjustedRandIndex(true_y, cluster_id)
  ari2 <- mclust::adjustedRandIndex(true_y, ifelse(cluster_id == 1, 2, 1))
  ARI <- max(ari1, ari2)

  pur1 <- purity(motif_tbl$label, factor(cluster_id))
  pur2 <- purity(motif_tbl$label, factor(ifelse(cluster_id == 1, 2, 1)))
  PUR <- max(pur1, pur2)

  correct_count <- sum((cluster_id == 1 & motif_tbl$label == "gyne") |
                         (cluster_id == 2 & motif_tbl$label == "worker"))
  misclassified_count <- nrow(motif_tbl) - correct_count
  correct_percent <- correct_count / nrow(motif_tbl) * 100

  list(
    sil_score = sil_score,
    ARI = ARI,
    PUR = PUR,
    correct_count = correct_count,
    misclassified_count = misclassified_count,
    correct_percent = correct_percent
  )
}

make_plot_df <- function(motif_tbl, embed_xy, cluster_id, axis_names = c("Dim1","Dim2")) {
  embed_xy <- as.data.frame(embed_xy)
  colnames(embed_xy) <- axis_names

  plot_df <- dplyr::bind_cols(
    motif_tbl %>% dplyr::select(name, label),
    as_tibble(embed_xy),
    tibble(cluster = factor(cluster_id))
  )

  plot_df %>%
    mutate(
      cluster = factor(cluster, levels = c(1, 2), labels = c("Cluster 1", "Cluster 2")),
      point_color = ifelse(cluster == "Cluster 1", "#66C2A5", "#FC8D62"),
      point_shape = ifelse(label == "gyne", 16, 17)
    )
}

add_misclass_annotation <- function(p, plot_df, xcol, ycol, label_col = "name") {
  cluster2_gyne <- plot_df %>% filter(cluster == "Cluster 2" & label == "gyne")
  if (nrow(cluster2_gyne) == 0) return(p)

  # place labels to the left to avoid overlap 
  x_range <- diff(range(plot_df[[xcol]]))
  label_positions <- cluster2_gyne %>%
    mutate(
      label_x = .data[[xcol]] - x_range * 0.12,
      label_y = .data[[ycol]]
    )

  p +
    geom_point(
      data = cluster2_gyne,
      shape = 1,
      color = "red",
      size = 4,
      stroke = 1.0,
      alpha = 0.9
    ) +
    geom_segment(
      data = label_positions,
      aes(x = .data[[xcol]], xend = label_x + 0.01, y = .data[[ycol]], yend = label_y),
      color = "red", size = 0.3, alpha = 0.6,
      arrow = arrow(length = unit(0.01, "npc"), type = "closed")
    ) +
    geom_label(
      data = label_positions,
      aes(x = label_x, y = label_y, label = .data[[label_col]]),
      size = 2.0,
      color = "red",
      fill = alpha("white", 0.95),
      label.size = 0.2,
      label.padding = unit(0.1, "lines"),
      hjust = 1
    )
}

# ===== 3) Read motifs -> build embeddings -> clustering  =====
motif_tbl <- build_motif_table(FILE_PFM)

set.seed(SEED)
embed_mat <- build_embeddings(motif_tbl, SEED = SEED, window_size = 2, ic_weighted_mean = TRUE)

set.seed(SEED)
km <- kmeans(embed_mat, centers = 2, nstart = 50)
cluster_id <- km$cluster

metrics <- compute_clustering_metrics(embed_mat, motif_tbl, cluster_id)

# For reporting
cluster2_gyne_count <- sum(motif_tbl$label == "gyne" & cluster_id == 2)

# ===== 4) tSNE  =====
set.seed(SEED)
tsne_result <- Rtsne(embed_mat, dims = 2, perplexity = 30,
                     check_duplicates = FALSE, pca = TRUE, verbose = FALSE)
tsne2 <- tsne_result$Y
plot_df_tsne <- make_plot_df(motif_tbl, tsne2, cluster_id, axis_names = c("tSNE1","tSNE2"))

tsne_range_x <- diff(range(plot_df_tsne$tSNE1))
tsne_range_y <- diff(range(plot_df_tsne$tSNE2))
aspect_ratio <- tsne_range_y / tsne_range_x

p_tsne <- ggplot(plot_df_tsne, aes(x = tSNE1, y = tSNE2)) +
  stat_ellipse(
    data = subset(plot_df_tsne, cluster == "Cluster 1"),
    color = "#66C2A5", geom = "path", alpha = 0.7, level = 0.8,
    linetype = "dashed", size = 0.6, show.legend = FALSE
  ) +
  stat_ellipse(
    data = subset(plot_df_tsne, cluster == "Cluster 2"),
    color = "#FC8D62", geom = "path", alpha = 0.7, level = 0.8,
    linetype = "dashed", size = 0.6, show.legend = FALSE
  ) +
  geom_point(
    aes(color = point_color, shape = factor(point_shape)),
    size = 2.5, alpha = 0.8
  ) +
  scale_color_identity(
    guide = "legend", name = "Cluster",
    breaks = c("#66C2A5", "#FC8D62"),
    labels = c("Cluster 1", "Cluster 2")
  ) +
  scale_shape_manual(
    name = "Caste",
    values = c("16" = 16, "17" = 17),
    labels = c("16" = "Gyne", "17" = "Worker")
  ) +
  labs(
    title = "t-SNE Visualization of Motif Sequence Clustering",
    subtitle = sprintf(
      "ARI: %.3f | Purity: %.3f | Correct: %d/%d (%.1f%%) | Errors: %d | Cluster2-Gyne: %d",
      metrics$ARI, metrics$PUR, metrics$correct_count, nrow(plot_df_tsne),
      metrics$correct_percent, metrics$misclassified_count, cluster2_gyne_count
    ),
    x = "t-SNE1",
    y = "t-SNE2",
    caption = "Dashed ellipses: Cluster boundaries | Red circle: Misclassified Gyne motifs in Cluster 2"
  ) +
  coord_fixed(ratio = aspect_ratio * 0.8) +
  theme_sci() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.1, "cm"),
    legend.key = element_rect(fill = "white"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0.5),
    plot.caption = element_text(size = 8, color = "grey50", hjust = 0.5),
    aspect.ratio = 0.8
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 15, size = 3)),
    shape = guide_legend(order = 2, override.aes = list(color = "black", size = 2.5))
  )

# Add misclassification annotation 
p_tsne <- add_misclass_annotation(p_tsne, plot_df_tsne, xcol = "tSNE1", ycol = "tSNE2")

out_tsne <- file.path(OUTDIR_TSNE, "Final_tSNE_Clustering.pdf")
ggsave(out_tsne, p_tsne, width = 6, height = 5, dpi = 600, bg = "white")
message("Saved t-SNE figure: ", out_tsne)
print(p_tsne)

# Save tSNE clustering statistics
stats_df_tsne <- data.frame(
  Metric = c("ARI", "Purity", "Correct_Count", "Total_Count", "Correct_Percent",
             "Misclassified_Count", "Cluster2_Gyne_Count"),
  Value  = c(metrics$ARI, metrics$PUR, metrics$correct_count, nrow(plot_df_tsne),
             metrics$correct_percent, metrics$misclassified_count, cluster2_gyne_count)
)
write.csv(stats_df_tsne, file.path(OUTDIR_TSNE, "clustering_statistics.csv"), row.names = FALSE)
message("Saved t-SNE statistics: ", file.path(OUTDIR_TSNE, "clustering_statistics.csv"))

# Console report 
cat("\nDetailed clustering summary (t-SNE):\n")
cat(sprintf("Total motifs: %d\n", nrow(plot_df_tsne)))
cat(sprintf("Correct: %d (%.1f%%)\n", metrics$correct_count, metrics$correct_percent))
cat(sprintf("Errors: %d\n", metrics$misclassified_count))
cat(sprintf("ARI: %.3f\n", metrics$ARI))
cat(sprintf("Purity: %.3f\n", metrics$PUR))
cat(sprintf("Cluster 2 gyne count: %d\n", cluster2_gyne_count))

# ===== 5) UMAP =====
set.seed(SEED)
umap_result <- uwot::umap(
  embed_mat,
  n_components = 2,
  n_neighbors  = min(15, nrow(embed_mat) - 1),
  min_dist     = 0.1,
  metric       = "euclidean",
  verbose      = FALSE
)
plot_df_umap <- make_plot_df(motif_tbl, umap_result, cluster_id, axis_names = c("UMAP1","UMAP2"))

cluster2_gyne_umap <- plot_df_umap %>% filter(cluster == "Cluster 2" & label == "gyne")

p_umap <- ggplot(plot_df_umap, aes(x = UMAP1, y = UMAP2)) +
  stat_ellipse(
    data = subset(plot_df_umap, cluster == "Cluster 1"),
    color = "#66C2A5", geom = "path", alpha = 0.7, level = 0.8,
    linetype = "dashed", size = 0.5, show.legend = FALSE
  ) +
  stat_ellipse(
    data = subset(plot_df_umap, cluster == "Cluster 2"),
    color = "#FC8D62", geom = "path", alpha = 0.7, level = 0.8,
    linetype = "dashed", size = 0.5, show.legend = FALSE
  ) +
  geom_point(
    aes(color = point_color, shape = factor(point_shape)),
    size = 2.5, alpha = 0.8
  ) +
  geom_point(
    data = cluster2_gyne_umap,
    shape = 1, color = "red",
    size = 4, stroke = 0.6, alpha = 0.8
  ) +
  scale_color_identity(
    guide = "legend", name = "Cluster",
    breaks = c("#66C2A5", "#FC8D62"),
    labels = c("Cluster 1", "Cluster 2")
  ) +
  scale_shape_manual(
    name = "Caste",
    values = c("16" = 16, "17" = 17),
    labels = c("16" = "Gyne", "17" = "Worker")
  ) +
  labs(
    title = "UMAP Visualization of Motif Sequence Clustering",
    subtitle = sprintf(
      "ARI: %.3f | Purity: %.3f | Correct: %d/%d (%.1f%%) | Errors: %d | Cluster2-Gyne: %d",
      metrics$ARI, metrics$PUR, metrics$correct_count, nrow(plot_df_umap),
      metrics$correct_percent, metrics$misclassified_count, nrow(cluster2_gyne_umap)
    ),
    x = "UMAP1",
    y = "UMAP2",
    caption = "Dashed ellipses: Cluster boundaries | Red circle: Gyne motifs in Cluster 2"
  ) +
  coord_fixed(ratio = 1) +
  theme_sci() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.05, "cm"),
    legend.key = element_rect(fill = "white"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5),
    plot.caption = element_text(size = 7, color = "grey50", hjust = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(shape = 15, size = 3)),
    shape = guide_legend(order = 2, override.aes = list(color = "black", size = 2.5))
  )

# Add labels for misclassified gyne in cluster 2 
if (nrow(cluster2_gyne_umap) > 0) {
  label_positions <- cluster2_gyne_umap %>%
    mutate(
      label_x = UMAP1 - diff(range(plot_df_umap$UMAP1)) * 0.12,
      label_y = UMAP2
    )

  p_umap <- p_umap +
    geom_segment(
      data = label_positions,
      aes(x = UMAP1, xend = label_x + 0.01, y = UMAP2, yend = label_y),
      color = "red", size = 0.25, alpha = 0.6,
      arrow = arrow(length = unit(0.01, "npc"), type = "closed")
    ) +
    geom_label(
      data = label_positions,
      aes(x = label_x, y = label_y, label = name),
      size = 2.0,
      color = "red",
      fill = alpha("white", 0.95),
      label.size = 0.2,
      label.padding = unit(0.1, "lines"),
      hjust = 1
    )
}

out_umap <- file.path(OUTDIR_UMAP, "Final_UMAP_Clustering.pdf")
ggsave(out_umap, p_umap, width = 6, height = 4.5, dpi = 600, bg = "white")
message("Saved UMAP figure: ", out_umap)
print(p_umap)

# Save UMAP clustering statistics
stats_df_umap <- data.frame(
  Metric = c("ARI", "Purity", "Correct_Count", "Total_Count", "Correct_Percent",
             "Misclassified_Count", "Cluster2_Gyne_Count"),
  Value  = c(metrics$ARI, metrics$PUR, metrics$correct_count, nrow(plot_df_umap),
             metrics$correct_percent, metrics$misclassified_count, nrow(cluster2_gyne_umap))
)
write.csv(stats_df_umap, file.path(OUTDIR_UMAP, "clustering_statistics.csv"), row.names = FALSE)
message("Saved UMAP statistics: ", file.path(OUTDIR_UMAP, "clustering_statistics.csv"))

# Console report 
cat("\nDetailed clustering summary (UMAP):\n")
cat(sprintf("Total motifs: %d\n", nrow(plot_df_umap)))
cat(sprintf("Correct: %d (%.1f%%)\n", metrics$correct_count, metrics$correct_percent))
cat(sprintf("Errors: %d\n", metrics$misclassified_count))
cat(sprintf("ARI: %.3f\n", metrics$ARI))
cat(sprintf("Purity: %.3f\n", metrics$PUR))
cat(sprintf("Cluster 2 gyne count: %d\n", nrow(cluster2_gyne_umap)))
