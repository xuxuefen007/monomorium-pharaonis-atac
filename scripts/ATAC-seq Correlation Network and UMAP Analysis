# Load required libraries
library(ggplot2)
library(dplyr)
library(igraph)
library(tidygraph)
library(ggraph)
library(DiffBind)
library(vegan)
library(umap)
library(purrr)

# Set working directory
setwd("/data/work/embryo_larva_pupa_adult_new/")

# ==================== 1. DATA PROCESSING ====================
cat("=== Step 1: Data Processing ===\n")
monomorium <- dba(sampleSheet = "monomorium_118.csv")
cat("Loaded", length(monomorium$samples), "samples\n")

consensus_by_condition <- dba.peakset(monomorium, consensus = DBA_CONDITION, minOverlap = 0.66)
all_consensus <- dba(consensus_by_condition, mask = consensus_by_condition$masks$Consensus, minOverlap = 1)
final_peaks <- dba.peakset(all_consensus, bRetrieve = TRUE)
cat("Generated", length(final_peaks), "consensus peaks\n")

monomorium <- dba.count(monomorium, peaks = final_peaks)
save(monomorium, file = "monomorium_count_118.RData")
cat("Saved processed data: monomorium_count_118.RData\n")

# ==================== 2. DATA NORMALIZATION ====================
cat("\n=== Step 2: Data Normalization ===\n")
load("monomorium_count_118.RData")
monomorium <- dba.normalize(monomorium)
norm_counts <- dba.peakset(monomorium, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
count_matrix <- as.matrix(norm_counts[, -c(1:3)])
rownames(count_matrix) <- paste(norm_counts[,1], norm_counts[,2], norm_counts[,3], sep=":")
count_matrix_t <- t(count_matrix)
cat("Data dimensions:", nrow(count_matrix_t), "samples ×", ncol(count_matrix_t), "peaks\n")

# ==================== 3. PEAK FILTERING ====================
cat("\n=== Step 3: Peak Quality Filtering ===\n")
calculate_variation <- function(mat) {
  cv <- apply(mat, 2, function(x) if (mean(x) == 0) return(0) else sd(x) / mean(x))
  mean_counts <- colMeans(mat)
  var_counts <- apply(mat, 2, var)
  dispersion <- var_counts / mean_counts
  list(cv = cv, dispersion = dispersion, mean_expr = mean_counts)
}

var_stats <- calculate_variation(count_matrix_t)
selected <- which(var_stats$cv > 0.398 & var_stats$mean_expr > 178.680)
filtered_matrix <- count_matrix_t[, selected, drop=FALSE]
cat("Retained", length(selected), "peaks (", 
    round(length(selected)/ncol(count_matrix_t)*100, 1), "% of total)\n", sep="")

# ==================== 4. PCA ANALYSIS ====================
cat("\n=== Step 4: PCA Analysis ===\n")
pca_result <- prcomp(count_matrix_t, scale. = TRUE)
pca_var <- summary(pca_result)$importance[2, ]
cum_var <- cumsum(pca_var)
n_pcs_optimal <- which(cum_var >= 0.8)[1]
cat("Using first", n_pcs_optimal, "PCs, cumulative variance:", 
    round(cum_var[n_pcs_optimal]*100, 2), "%\n")

# ==================== 5. UMAP ANALYSIS ====================
cat("\n=== Step 5: UMAP Analysis ===\n")
set.seed(123)
umap_result <- umap(
  pca_result$x[, 1:n_pcs_optimal],
  n_neighbors = 15,
  min_dist = 0.1,
  n_components = 2
)

# Prepare sample metadata
sample_info <- dba.show(monomorium)

umap_df <- data.frame(
  Sample = rownames(count_matrix_t),
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  Condition = sample_info$Condition[match(rownames(count_matrix_t), sample_info$SampleID)]
)

# Define developmental stages
extract_stage <- function(condition) {
  parts <- strsplit(condition, "_")[[1]]
  if(parts[1] == "embryo" && length(parts) >= 3) {
    return(paste(parts[1:3], collapse="_"))
  } else if(parts[1] == "adult") {
    return("adult")
  } else if(parts[1] %in% c("larva", "pupa") && length(parts) >= 2) {
    return(paste(parts[1:2], collapse="_"))
  } else {
    return(paste(parts[1:min(2, length(parts))], collapse="_"))
  }
}

umap_df$Stage <- sapply(umap_df$Condition, extract_stage)
umap_df$Caste <- ifelse(
  grepl("embryo|larva_1st", umap_df$Condition), "unknown",
  ifelse(grepl("gyne", umap_df$Condition), "gyne", "worker")
)

# Define stage order
stage_order <- c(
  "embryo_0_12", "embryo_12_24", "embryo_36_48", "embryo_60_72",
  "embryo_84_96", "embryo_108_120", "embryo_132_144", "embryo_156_168",
  "embryo_180_192", "larva_1st", "larva_2nd", "larva_3rd",
  "pre_pupa", "young_pupa", "old_pupa", "adult"
)

umap_df$Stage <- factor(umap_df$Stage, levels = stage_order)

# ==================== 6. UMAP VISUALIZATION ====================
cat("\n=== Step 6: UMAP Visualization ===\n")

# Define color scheme for developmental stages
stage_colors <- c(
  "embryo_0_12" = "#E6E0B0", "embryo_12_24" = "#36600E", "embryo_36_48" = "#D96558",
  "embryo_60_72" = "#B43970", "embryo_84_96" = "#692F7C", "embryo_108_120" = "#282A62",
  "embryo_132_144" = "#6d8bc3", "embryo_156_168" = "#A3C9D5", "embryo_180_192" = "#AED185",
  "larva_1st" = "#F6C63C", "larva_2nd" = "#e3716e", "larva_3rd" = "#7ac7e2",
  "pre_pupa" = "#fae6e9", "young_pupa" = "#bdb5e1", "old_pupa" = "#E4B7BC", "adult" = "#E9687A"
)

caste_shapes <- c("gyne" = 16, "worker" = 17, "unknown" = 18)

# Create UMAP plot
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Stage, shape = Caste), size = 3, alpha = 0.8) +
  scale_color_manual(values = stage_colors) +
  scale_shape_manual(values = caste_shapes) +
  labs(
    title = "M. pharaonis Chromatin Accessibility UMAP",
    subtitle = paste0("Based on first ", n_pcs_optimal, 
                     " PCs (", round(cum_var[n_pcs_optimal]*100, 1), "% variance)"),
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Developmental Stage",
    shape = "Caste"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    legend.box = "vertical"
  )

# Save UMAP plot
pdf("ATAC_UMAP_visualization.pdf", width = 10, height = 8)
print(umap_plot)
dev.off()
cat("UMAP visualization saved to: ATAC_UMAP_visualization.pdf\n")

# ==================== 7. CORRELATION NETWORK ANALYSIS ====================
cat("\n=== Step 7: Correlation Network Analysis ===\n")
cor_mat <- cor(t(filtered_matrix), method = "spearman")
net <- graph_from_adjacency_matrix(cor_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
net <- delete_edges(net, E(net)[weight < 0.8])
cat("Network: ", vcount(net), "nodes, ", ecount(net), "edges (rho ≥ 0.8)\n")

# Add metadata to network
meta <- data.frame(
  node = V(net)$name,
  stage = umap_df$Stage[match(V(net)$name, umap_df$Sample)],
  caste = umap_df$Caste[match(V(net)$name, umap_df$Sample)]
)

V(net)$stage <- meta$stage
V(net)$caste <- meta$caste

# Network visualization
tidy_net <- as_tbl_graph(net) %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree(), community = as.factor(group_louvain()))

set.seed(2024)
pdf("ATAC_correlation_network.pdf", width = 9, height = 10, useDingbats = FALSE)
ggraph(tidy_net, layout = "graphopt") +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.6, lineend = "round") +
  scale_edge_width(range = c(0.5, 3), name = "Edge Weight") +
  scale_edge_color_gradientn(
    colors = colorRampPalette(c("#FDE725", "#21908C", "#440154"))(100),
    limits = c(0.8, 1.0), breaks = c(0.8, 0.9, 1.0), name = "Spearman rho"
  ) +
  geom_node_point(aes(fill = stage, shape = caste), size = 5, color = "gray30", stroke = 0.3) +
  scale_fill_manual(values = stage_colors, name = "Developmental Stage", drop = FALSE, breaks = stage_order,
    guide = guide_legend(override.aes = list(size = 5, shape = 21), ncol = 1, title.position = "top")) +
  scale_shape_manual(values = c(gyne = 21, worker = 22, unknown = 24), name = "Caste",
    guide = guide_legend(override.aes = list(size = 5, fill = "gray50"), ncol = 1)) +
  theme_graph(base_family = "sans") +
  theme(legend.position = "right", plot.margin = margin(10, 10, 10, 10),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10), legend.spacing.y = unit(0.5, "cm")) +
  labs(title = "ATAC-seq Sample Correlation Network",
    subtitle = "Spearman rho > 0.8 | Layout: graphopt | Top 15% CV & Top 35% Mean")
dev.off()
cat("Correlation network saved to: ATAC_correlation_network.pdf\n")

# ==================== 8. PERMANOVA ANALYSIS ====================
cat("\n=== Step 8: PERMANOVA Analysis ===\n")
target_conditions <- c(
  "larva_2nd_gyne", "larva_2nd_worker",
  "larva_3rd_gyne", "larva_3rd_worker", 
  "pre_pupa_gyne", "pre_pupa_worker",
  "young_pupa_gyne", "young_pupa_worker",
  "old_pupa_gyne", "old_pupa_worker",
  "adult_gyne", "adult_worker"
)

metadata <- data.frame(
  Sample = rownames(count_matrix_t),
  Condition = gsub("\\..*", "", rownames(count_matrix_t)),
  stringsAsFactors = FALSE
) %>% 
  filter(Condition %in% target_conditions) %>%
  mutate(
    Stage = gsub("_gyne|_worker", "", Condition),
    Caste = ifelse(grepl("_gyne", Condition), "gyne", "worker"),
    Stage = factor(Stage, levels = c("larva_2nd", "larva_3rd", "pre_pupa", 
                                   "young_pupa", "old_pupa", "adult"))
  )

count_matrix_filtered <- count_matrix_t[metadata$Sample, ]

# PCA-based PERMANOVA
pca_result_filtered <- prcomp(count_matrix_filtered, scale. = TRUE)
pc_scores <- pca_result_filtered$x[, 1:n_pcs_optimal]

set.seed(123)
adonis_pc <- adonis2(
  dist(pc_scores) ~ Stage * Caste,
  data = metadata,
  permutations = 9999,
  by = "terms"
)

# UMAP-based PERMANOVA
umap_config <- umap.defaults
umap_config$n_neighbors <- min(15, nrow(count_matrix_filtered)-1)
umap_config$min_dist = 0.2

set.seed(123)
umap_pca_result <- umap(pc_scores, config = umap_config)
umap_pca_dist <- dist(umap_pca_result$layout)

set.seed(123)
adonis_umap_pca <- adonis2(
  umap_pca_dist ~ Stage * Caste,
  data = metadata,
  permutations = 9999,
  by = "terms"
)

# ==================== 9. PER-STAGE PERMANOVA ANALYSIS ====================
cat("\n=== Step 9: Per-Stage PERMANOVA Analysis ===\n")

# Filter stages with both castes
stages_with_both <- metadata %>%
  group_by(Stage) %>%
  summarise(
    n_gyne = sum(Caste == "gyne"),
    n_worker = sum(Caste == "worker")
  ) %>%
  filter(n_gyne > 0 & n_worker > 0) %>%
  pull(Stage)

cat("Stages analyzed (both castes present):", paste(stages_with_both, collapse = ", "), "\n")

per_stage_results <- map_df(stages_with_both, function(stage) {
  # Filter samples for current stage
  stage_samples <- metadata %>% filter(Stage == stage) %>% pull(Sample)
  stage_pc <- pc_scores[rownames(pc_scores) %in% stage_samples, ]
  stage_dist <- dist(stage_pc, method = "euclidean")
  
  # Execute PERMANOVA
  set.seed(123)
  res <- adonis2(
    stage_dist ~ Caste,
    data = metadata %>% filter(Sample %in% stage_samples),
    permutations = 9999
  )
  
  # Format results
  data.frame(
    Stage = stage,
    Caste_R2 = res$R2[1],
    Caste_F = res$F[1],
    Caste_p = res$`Pr(>F)`[1],
    SampleSize = paste(
      sum(metadata$Stage == stage & metadata$Caste == "gyne"),
      sum(metadata$Stage == stage & metadata$Caste == "worker"),
      sep = "/"
    )
  )
})

# Add significance markers
per_stage_results <- per_stage_results %>%
  mutate(
    Significance = case_when(
      Caste_p < 0.001 ~ "***",
      Caste_p < 0.01 ~ "**",
      Caste_p < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Save per-stage results
write.csv(per_stage_results, "per_stage_permanova_results.csv", row.names = FALSE)

# ==================== 10. NGI ANALYSIS ====================
calculate_absolute_NGI <- function(pca_res, umap_pca_res) {
  R2_PCA <- pca_res$R2[1:3]
  R2_UMAP_PCA <- umap_pca_res$R2[1:3]
  names(R2_PCA) <- names(R2_UMAP_PCA) <- c("Stage", "Caste", "Stage × Caste")
  
  NGI_Stage <- (R2_UMAP_PCA["Stage"] / R2_PCA["Stage"]) - 1
  NGI_Caste <- (R2_UMAP_PCA["Caste"] / R2_PCA["Caste"]) - 1
  NGI_Interaction <- (R2_UMAP_PCA["Stage × Caste"] / R2_PCA["Stage × Caste"]) - 1
  
  interpret_NGI <- function(val) {
    if (val > 0.3) "Nonlinear"
    else if (val < -0.3) "Linear"
    else "Mixed"
  }
  
  data.frame(
    Term = c("Stage", "Caste", "Stage × Caste"),
    R2_PCA = round(R2_PCA, 4),
    R2_UMAP = round(R2_UMAP_PCA, 4),
    NGI = c(round(NGI_Stage, 4), round(NGI_Caste, 4), round(NGI_Interaction, 4)),
    Pattern_Type = c(interpret_NGI(NGI_Stage),
                     interpret_NGI(NGI_Caste),
                     interpret_NGI(NGI_Interaction)),
    stringsAsFactors = FALSE
  )
}

pattern_results <- calculate_absolute_NGI(adonis_pc, adonis_umap_pca)

# Save PERMANOVA results
permanova_results <- list(
  pca = adonis_pc,
  umap = adonis_umap_pca,
  ngi = pattern_results,
  per_stage = per_stage_results
)

save(permanova_results, file = "permanova_results.RData")

# ==================== 11. MANTEL TEST ANALYSIS ====================
cat("\n=== Step 11: Mantel Test Analysis ===\n")
time_mapping <- c(
  "embryo_0_12" = 1, "embryo_12_24" = 2,
  "embryo_36_48" = 3, "embryo_60_72" = 4,
  "embryo_84_96" = 5, "embryo_108_120" = 6,
  "embryo_132_144" = 7, "embryo_156_168" = 8,
  "embryo_180_192" = 9,
  "larva_1st" = 10, "larva_2nd" = 11,
  "larva_3rd" = 12,
  "pre_pupa" = 13, "young_pupa" = 14,
  "old_pupa" = 15,
  "adult" = 16
)

# Check for undefined stages
undefined_stages <- setdiff(unique(meta$stage), names(time_mapping))
if(length(undefined_stages) > 0) {
  warning("Undefined stages (skipped in Mantel test): ", paste(undefined_stages, collapse = ", "))
}

# Calculate distance matrices
full_cor_mat <- cor(t(count_matrix_t), method = "spearman")
net_dist <- as.dist(1 - full_cor_mat)

time_values <- time_mapping[as.character(meta$stage)]
time_dist <- dist(time_values)

# Perform Mantel test
set.seed(123)
mantel_result <- mantel(time_dist, net_dist, method = "spearman", permutations = 9999)

mantel_results <- list(
  statistic = mantel_result$statistic,
  p_value = mantel_result$signif,
  permutations = mantel_result$permutations
)

save(mantel_results, file = "mantel_test_results.RData")
