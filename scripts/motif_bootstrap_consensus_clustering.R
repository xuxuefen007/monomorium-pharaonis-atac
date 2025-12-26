rm(list = ls()); gc()
options(stringsAsFactors = FALSE, warn = -1)

# ============================================================
# 0) Basic Settings
# ============================================================
setwd("C:/Users/Xuefen Xu/Documents/tSNE")

FILE_PFM <- "C:/Users/Xuefen Xu/Documents/tSNE/gyne_worker_G75W140.txt"  # High-confidence motif set (75 gyne, 140 worker; boundary motifs removed)
OUTDIR   <- file.path(getwd(), "CONSENSUS_CLUSTER_06")  # Output directory for results
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

SEED <- 123
set.seed(SEED)
N_BOOT <- 500        # Number of bootstrap iterations
STAB_THRESH <- 0.6   # Stability threshold for motif inclusion

# ============================================================
# 1) Dependencies
# ============================================================
suppressPackageStartupMessages({
  library(motifStack)      # For motif processing and visualization
  library(universalmotif)  # For motif comparison and similarity analysis
  library(dplyr)           # For data manipulation
  library(Cairo)           # For high-quality plot export
  library(grid)            # For adding text annotations on plots
})

# ============================================================
# 2) Read motifs and convert to PFM
# ============================================================
motifs <- importMatrix(FILE_PFM, format = "jaspar")

# Convert to PFM (Position Frequency Matrix) and trim low information content motifs
pfms <- lapply(seq_along(motifs), function(i) {
  pcm <- motifs[[i]]@mat
  pfm_mat <- apply(pcm, 2, function(x) x / sum(x))  # Normalize columns to get PFM
  rownames(pfm_mat) <- c("A", "C", "G", "T")  # Set row names for nucleotides
  new("pfm", mat = pfm_mat, name = names(motifs)[i])
})
names(pfms) <- names(motifs)

# Trim motifs with low information content (IC < 0.4)
pfms <- lapply(pfms, trimMotif, t = 0.4)

# ============================================================
# 3) PWM (Position Weight Matrix) distance calculation
# ============================================================
# Create a list of universalmotif objects from PFMs
um_list <- lapply(pfms, function(p) universalmotif::create_motif(p@mat, name = p@name))

# Compare motifs using Pearson correlation coefficient (PCC) and reverse complement if needed
sim_pwm <- universalmotif::compare_motifs(um_list, method = "PCC", tryRC = TRUE, use.type = "ICM")

# Convert the similarity matrix into a distance matrix
D_pwm <- 1 - as.matrix(sim_pwm)
diag(D_pwm) <- 0  # Set diagonal to 0 to avoid self-comparison
rownames(D_pwm) <- colnames(D_pwm) <- names(pfms)

# Assign labels based on motif names (gyne vs. worker)
labels <- ifelse(grepl("^gyne", names(pfms), ignore.case = TRUE), "gyne", "worker")

# ============================================================
# 4) Bootstrap Consensus Clustering
# ============================================================
# Matrix to store cluster assignments for each bootstrap iteration
boot_clusters <- matrix(NA, nrow = length(pfms), ncol = N_BOOT,
                        dimnames = list(names(pfms), paste0("B", 1:N_BOOT)))

# Perform bootstrap iterations
for (b in 1:N_BOOT) {
  samp_idx <- sample(1:ncol(D_pwm), replace = TRUE)  # Sample with replacement
  D_boot <- D_pwm[samp_idx, samp_idx]  # Bootstrap distance matrix
  
  hc <- hclust(as.dist(D_boot), method = "ward.D2")  # Hierarchical clustering using Ward's method
  cl <- cutree(hc, k = 2)  # Cut the dendrogram into 2 clusters (gyne and worker)
  
  boot_clusters[names(cl), b] <- cl  # Store the cluster assignments
}

# ============================================================
# 5) Stability Statistics
# ============================================================
stab_df <- data.frame(motif = rownames(boot_clusters),
                      label = labels,
                      prop_agree = NA)  # Store the proportion of agreement for each motif

# Calculate the stability for each motif across bootstrap iterations
for (i in 1:nrow(boot_clusters)) {
  cl_assign <- boot_clusters[i, ]
  majority <- ifelse(mean(cl_assign == "1", na.rm = TRUE) > 0.5, 1, 2)  # Majority vote for stability
  prop <- mean(cl_assign == majority, na.rm = TRUE)  # Proportion of agreement with the majority
  stab_df$prop_agree[i] <- prop
}

# Determine if the motif is stable based on the threshold
stab_df$keep <- stab_df$prop_agree >= STAB_THRESH

# ============================================================
# 6) Export Results
# ============================================================
write.csv(stab_df, file.path(OUTDIR, "motif_stability_thresh06.csv"), row.names = FALSE)

cat("Total motifs: ", nrow(stab_df), "\n")
cat("Stable motifs: ", sum(stab_df$keep), "\n")
cat("Unstable motifs: ", sum(!stab_df$keep), "\n")

# ============================================================
# 7) Create Motif Tree (only stable motifs)
# ============================================================
good_motifs <- stab_df$motif[stab_df$keep]  # Filter stable motifs
pfms_filtered <- pfms[names(pfms) %in% good_motifs]  # Filter PFMs to include only stable motifs

# Split motifs into gyne and worker groups
is_gyne <- grepl("^gyne", names(pfms_filtered), ignore.case = TRUE)
is_worker <- grepl("^worker", names(pfms_filtered), ignore.case = TRUE)

pfms_gyne <- pfms_filtered[is_gyne]
pfms_worker <- pfms_filtered[is_worker]
pfms_ordered <- c(pfms_gyne, pfms_worker)  # Combine the motifs in order

n_red <- length(pfms_gyne)
n_blue <- length(pfms_worker)
col_red <- "#FF4D4D"  # Color for gyne motifs
col_blue <- "#3D7CFF"  # Color for worker motifs
group_colors <- c(rep(col_red, n_red), rep(col_blue, n_blue))

# ============================================================
# 8) Plotting
# ============================================================
while (!is.null(dev.list())) dev.off()  # Close any existing plots
CairoPDF(file.path(OUTDIR, "motif_tree_consensus_thresh06.pdf"), width = 14, height = 14)

# Plot the consensus motif tree using motifStack
motifStack(
  pfms_ordered,
  layout = "radialPhylog",  # Radial phylogeny layout
  circle = 0.75,  # Size of the overall tree
  cleaves = 0.005,  # Distance from leaves to center
  clabel.leaves = 0.38,  # Distance of motif IDs from the leaves
  col.bg = group_colors,  # Background color for each group
  col.bg.alpha = 0.15,  # Transparency of background color
  col.leaves = group_colors,  # Color of leaves
  circle.motif = 0,  # Remove outer circle motif logos
  angle = 360,
  edge.lwd = 0.9,
  edge.col = "#333333"
)

# Add labels to the plot (gyne and worker counts)
grid::grid.text(sprintf("gyne (n=%d)", n_red), 
                x = unit(0.05, "npc"), y = unit(0.95, "npc"),
                gp = gpar(col = "#FF4D4D", fontsize = 14))
grid::grid.text(sprintf("worker (n=%d)", n_blue), 
                x = unit(0.05, "npc"), y = unit(0.90, "npc"),
                gp = gpar(col = "#3D7CFF", fontsize = 14))

dev.off()  # Save the plot
