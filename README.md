# monomorium-pharaonis-atac
ATAC-seq analysis of chromatin accessibility underlying caste differentiation in Monomorium pharaonis

# Outline

### `Developmental_trajectory.R` 
**Developmental Trajectory Analysis**: Comprehensive analysis of chromatin accessibility dynamics across development, including data normalization, UMAP visualization, correlation network construction, PERMANOVA statistical testing, and Mantel test for trajectory analysis.

### `developmental_heatmap.R` 
**Main Analysis Script**: Generates chromatin accessibility heatmaps showing dynamic changes across developmental stages from embryo_0_12h to larva_1st. This script performs differential peak analysis, DESeq2 normalization, and visualizes top variable peaks in a publication-quality heatmap.

### `zga_caste_co_motif_analysis.R` 
**ZGA-Caste Co-motif Interaction Analysis**: Analyzes the co-occurrence of zygotic genome activation (ZGA) transcription factor motifs with caste-specific motifs in chromatin accessibility regions. The script performs 6 different comparative analyses and generates comprehensive statistical reports with visualizations of motif interactions.
