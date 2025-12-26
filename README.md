# monomorium-pharaonis-atac
ATAC-seq analysis of chromatin accessibility underlying caste differentiation in Monomorium pharaonis

# Outline

### `Developmental_trajectory.R` 
**Developmental Trajectory Analysis**: Comprehensive analysis of chromatin accessibility dynamics across development, including data normalization, UMAP visualization, correlation network construction, PERMANOVA statistical testing, and Mantel test for trajectory analysis.

### `developmental_heatmap.R` 
**Main Analysis Script**: Generates chromatin accessibility heatmaps showing dynamic changes across developmental stages from embryo_0_12h to larva_1st. This script performs differential peak analysis, DESeq2 normalization, and visualizes top variable peaks in a publication-quality heatmap.

### `zga_caste_co_motif_analysis.R` 
**ZGA-Caste Co-motif Interaction Analysis**: Analyzes the co-occurrence of zygotic genome activation (ZGA) transcription factor motifs with caste-specific motifs in chromatin accessibility regions. The script performs 6 different comparative analyses and generates comprehensive statistical reports with visualizations of motif interactions.

### `ant-caste-TF-bias_AME.R`
**Caste-biased Transcription Factor Identification from Motif Enrichment**：This script identifies and ranks caste-biased transcription factors based on motif enrichment analyses.

### `rf_motif_classification.R` 
**Motif Classification Using Random Forest**: This script performs motif classification using a Random Forest (RF) model to distinguish between gyne- and worker-biased motifs in chromatin accessibility regions. The analysis includes parsing position frequency matrices (PFMs), extracting key features from the motifs, training the RF classifier, and evaluating model performance with confusion matrices, error rates, and feature importance analysis. The output also includes visualizations such as heatmaps of motif probabilities and proximity between samples.

### `rf_shap_perm_analysis.R`
**Random Forest Feature Importance with Permutation and SHAP**: This script trains a Random Forest (RF) classifier on PFM-derived motif features to distinguish gyne- and worker-biased motifs, and provides complementary feature-importance interpretations using Gini importance, permutation importance (ROC AUC), and SHAP (fastshap). The workflow includes parsing PFMs, extracting IC/GC/AT and base-composition features, training the RF model, computing permutation importance (nsim = 20) and SHAP values (nsim = 200), and exporting publication-ready plots (Gini barplot, permutation importance barplot, SHAP importance barplot, and SHAP beeswarm) to the output directory.

### `rf_lda_motif_boundary_cleanpool.R`
**RF/LDA-Based Motif Boundary Detection and Clean Pool Construction**：This script parses motif PFMs, extracts sequence-composition and information-content features, and trains a Random Forest (true OOB probabilities + OOB proximity) and cross-validated LDA to identify “boundary” motifs with ambiguous caste assignment. It then merges RF and LDA boundary calls to derive a high-confidence “clean motif pool”, aggregates motif-level predictions to TF-level candidates using a consistency threshold, and generates diagnostic plots (LD1 ECDF, raincloud, Top-N LDA features, and base-composition distributions) plus threshold-sensitivity scans.

### `motif_embedding_tsne_umap.R`
**Motif embedding–based t-SNE and UMAP clustering of high-confidence caste-biased motifs**：This script performs embedding-based clustering and visualization of high-confidence gyne- and worker-biased transcription factor motifs obtained after RF/LDA-based boundary filtering. Motifs are embedded using an information-content–weighted distributional representation (PPMI with truncated SVD, as a stable alternative to GloVe) and visualized with t-SNE and UMAP. Clustering quality is assessed using adjusted Rand index, purity, and classification accuracy, with misclassified motifs explicitly highlighted for diagnostic purposes.

### `motif_embedding_RF_SVM_nestedCV.R`
**Nested cross-validated RF and SVM classification of high-confidence caste-biased motifs**：This script applies an information-content–weighted PPMI–SVD motif embedding followed by nested cross-validation with Random Forest and linear SVM classifiers to evaluate and classify high-confidence gyne- and worker-biased transcription factor motifs. Model performance is assessed using standard metrics, and motifs with inconsistent predictions are highlighted as potential boundary or ambiguous regulators.


