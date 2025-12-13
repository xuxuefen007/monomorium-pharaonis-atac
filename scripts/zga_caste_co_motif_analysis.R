# ==============================================================================
# Co-motif Analysis: ZGA and Caste-Specific Transcription Factor Interactions
# 
# This script analyzes the co-occurrence of zygotic genome activation (ZGA) motifs
# and caste-specific transcription factor motifs in chromatin accessibility regions
# across different developmental conditions in Monomorium pharaonis.
# 
# Key Analyses:
# 1. Tests co-occurrence enrichment between ZGA motifs and caste-specific motifs
# 2. Compares motif interactions across different comparisons (gyne vs worker, etc.)
# 3. Generates comprehensive statistical summaries and visualizations
# 
# Input Requirements:
# 1. BED files for different peak sets:
#    - larva_2nd_gyne_higher_peaks_FDR0.01.homer.bed
#    - larva_2nd_worker_higher_peaks_FDR0.01.homer.bed
#    - larva_2nd_non_significant_peaks.bed
# 2. Genome assembly FASTA file
# 3. Pre-defined TF motifs for ZGA and caste-specific factors
# ==============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer); library(Rsamtools)
  library(Biostrings);    library(TFBSTools);  library(motifmatchr)
  library(IRanges);       library(dplyr);      library(tidyr)
  library(tibble);        library(readr);      library(data.table)
  library(ggplot2);       library(patchwork)
})

# ---- Directory Configuration ----
BED_GYNE    <- "/data/work/hello_ZGA/raw_data/larva_2nd_gyne_higher_peaks_FDR0.01.homer.bed"
BED_WORKER  <- "/data/work/hello_ZGA/raw_data/larva_2nd_worker_higher_peaks_FDR0.01.homer.bed"
BED_NONDA   <- "/data/work/hello_ZGA/raw_data/larva_2nd_non_significant_peaks.bed"
GENOME_FA   <- "/data/work/hello_ZGA/raw_data/Mpha.genome-mito.fa"
OUTDIR      <- "/data/work/hello_ZGA/larva_2nd_CO_MOTIF_OUT_FULL-update"
SUMMARY_DIR <- file.path(OUTDIR, "SUMMARY")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SUMMARY_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Analysis Parameters ----
NEAR_BP  <- 100      # Maximum distance (bp) for motif co-occurrence
P_CUTOFF <- 5e-4     # P-value cutoff for motif matching
N_LENBIN <- 6        # Number of bins for length matching
N_GCBIN  <- 6        # Number of bins for GC content matching
SEED     <- 123      # Random seed for reproducibility
RUN_CASE_VS_CASE <- TRUE  # Whether to run case vs case comparisons
set.seed(SEED)

# ---- Utility Functions ----
read_bed <- function(p){ gr <- import(p, "BED"); gr[width(gr)>0] }
make_peak_id <- function(gr) paste0(as.character(seqnames(gr)), ":", start(gr), "-", end(gr))

# Add GC content and length information to peak table
add_gc_len_tbl <- function(gr, fa){
  s <- getSeq(fa, gr)
  tibble(idx=seq_along(gr), width=width(gr), 
         gc=as.numeric(letterFrequency(s, "GC", as.prob=TRUE)))
}

# Strict 1:1 matching based on length and GC content
strict_pair_11 <- function(case_tbl, pool_tbl, n_lenbin=6, n_gcbin=6){
  len_all <- c(case_tbl$width, pool_tbl$width)
  gc_all  <- c(case_tbl$gc,    pool_tbl$gc)
  len_brk <- unique(quantile(len_all, probs=seq(0,1,length.out=n_lenbin+1), na.rm=TRUE))
  gc_brk  <- unique(quantile(gc_all,  probs=seq(0,1,length.out=n_gcbin+1), na.rm=TRUE))
  
  case_tbl <- mutate(case_tbl,
                     len_bin=cut(width,len_brk,include.lowest=TRUE),
                     gc_bin=cut(gc,gc_brk,include.lowest=TRUE),
                     bin_key=interaction(len_bin,gc_bin,drop=TRUE))
  pool_tbl <- mutate(pool_tbl,
                     len_bin=cut(width,len_brk,include.lowest=TRUE),
                     gc_bin=cut(gc,gc_brk,include.lowest=TRUE),
                     bin_key=interaction(len_bin,gc_bin,drop=TRUE))
  
  bins <- union(levels(case_tbl$bin_key), levels(pool_tbl$bin_key))
  pc <- list(); pp <- list()
  
  for (b in bins){
    c_bin <- dplyr::filter(case_tbl, bin_key==b)
    p_bin <- dplyr::filter(pool_tbl, bin_key==b)
    n <- min(nrow(c_bin), nrow(p_bin)); if (n==0) next
    pc[[b]] <- tibble(idx=sample(c_bin$idx, n))
    pp[[b]] <- tibble(idx=sample(p_bin$idx, n))
  }
  
  list(case_idx=bind_rows(pc)$idx, pool_idx=bind_rows(pp)$idx)
}

# Convert position frequency matrix to PWMatrix object
pcm_to_PWMatrix <- function(A,C,G,T,name="x",id="x"){
  m <- rbind(A=A,C=C,G=G,T=T); ppm <- sweep(m,2,pmax(colSums(m),1),"/")
  rownames(ppm) <- c("A","C","G","T"); TFBSTools::PWMatrix(ID=id,name=name,profileMatrix=ppm)
}

# Check if motifs co-occur within specified distance
co_near_from_positions <- function(mmA, mmB, near_bp=100){
  n <- length(mmA[[1]]); flag <- logical(n); mind <- rep(Inf,n)
  
  for (i in seq_len(n)){
    rA <- IRanges(); for(m in mmA) if(length(m[[i]])>0) rA <- c(rA, ranges(m[[i]]))
    rB <- IRanges(); for(m in mmB) if(length(m[[i]])>0) rB <- c(rB, ranges(m[[i]]))
    if (length(rA)==0 || length(rB)==0) next
    
    cA <- (start(rA)+end(rA))%/%2; cB <- (start(rB)+end(rB))%/%2
    d  <- min(abs(outer(cA,cB,"-"))); mind[i] <- d; flag[i] <- (d<=near_bp)
  }
  
  list(flag=flag, min_dist=mind)
}

# Fisher's exact test wrapper
fisher_from_counts <- function(a_co,a_no,b_co,b_no){
  ft <- fisher.test(matrix(c(a_co,a_no,b_co,b_no),nrow=2,byrow=TRUE))
  list(OR=unname(ft$estimate),CI_lo=ft$conf.int[1],CI_hi=ft$conf.int[2],P=ft$p.value)
}

# Extract motif hits as data.table
extract_hits_dt <- function(mm, peak_ids){
  out <- data.table()
  for (i in seq_along(mm)){
    grl <- mm[[i]]; if (length(grl)==0) next
    nm <- names(mm)[i]
    for (j in seq_along(grl)){
      gr <- grl[[j]]; if (length(gr)==0) next
      out <- rbind(out, data.table(peak_id=peak_ids[j], motif=nm,
                                   center = (start(gr)+end(gr))%/%2))
    }
  }
  out
}

# Count motif pairs within specified distance
count_pairs <- function(hA, hB, near_bp=100){
  if (nrow(hA)==0 || nrow(hB)==0) return(tibble(motif_A=character(),motif_B=character(),n_peaks=integer()))
  
  grA <- GRanges(seqnames=hA$peak_id, ranges=IRanges(hA$center,hA$center), motif=hA$motif)
  grB <- GRanges(seqnames=hB$peak_id, ranges=IRanges(hB$center,hB$center), motif=hB$motif)
  ov  <- findOverlaps(grA, grB, maxgap=near_bp, type="any", ignore.strand=TRUE)
  
  if (length(ov)==0) return(tibble(motif_A=character(),motif_B=character(),n_peaks=integer()))
  
  tibble(motif_A=mcols(grA)$motif[queryHits(ov)],
         motif_B=mcols(grB)$motif[subjectHits(ov)],
         peak_id=as.character(seqnames(grA))[queryHits(ov)]) %>%
    distinct(motif_A,motif_B,peak_id) %>% count(motif_A,motif_B,name="n_peaks")
}

# ---- ZGA Motif Definitions ----
# Trl (Trithorax-like) motif
Trl_PWM <- pcm_to_PWMatrix(
  c(4107,6613,91,6743,10,6603,141,5105,1888),
  c(597,47,54,21,38,182,29,845,621),
  c(1902,44,6535,44,6783,45,6638,495,3800),
  c(252,154,178,50,27,28,50,413,549),
  name="Trl", id="MA0205.3"
)

# Clamp motif
Clamp_PWM <- pcm_to_PWMatrix(
  c(332,1071,189,713,74,2510,62,325,31,2451,71,2038,234,1251),
  c(250,455,175,1322,29,32,11,1998,12,28,36,332,502,327),
  c(1888,378,2178,354,2532,43,2562,174,2600,155,2517,118,1515,960),
  c(181,747,109,262,16,66,16,154,8,17,27,163,400,113),
  name="Clamp", id="MA1700.1"
)

# Opa (Opaque) motif
Opa_PWM <- pcm_to_PWMatrix(
  c(0,10,1,1,0,0,0,0,4,0,3,2),
  c(1,7,16,15,16,18,18,14,0,14,1,0),
  c(15,1,0,1,1,0,0,0,13,2,4,16),
  c(2,0,1,1,1,0,0,4,1,2,10,0),
  name="Opa", id="MA0456.1"
)

# Vfl motif
Vfl_PWM <- pcm_to_PWMatrix(
  c(52,11435,43,91,93,11460,2143),
  c(10668,28,58,53,689,125,1746),
  c(180,163,11602,11545,37,94,6863),
  c(831,105,28,42,10912,52,979),
  name="Vfl", id="MA1462.2"
)

# ZGA motif collection
ZGA_pwm <- TFBSTools::PWMatrixList(Trl_PWM, Clamp_PWM, Opa_PWM, Vfl_PWM)
names(ZGA_pwm) <- c("Trl", "Clamp", "Opa", "Vfl")

# ---- Gyne-Specific Motif Definitions ----
Brk_PWM   <- pcm_to_PWMatrix(
  c(1,0,0,0,0,1,1,0), c(5,4,0,0,10,0,8,5),
  c(4,0,10,10,0,9,0,1), c(0,6,0,0,0,0,1,4),
  name="brk", id="MA0213.1"
)

Hairy_PWM <- pcm_to_PWMatrix(
  c(1,3,1,22,1,2,1,1,7,1), c(2,7,31,1,30,1,10,1,17,30),
  c(30,17,1,10,1,30,1,31,7,2), c(1,7,1,1,2,1,22,1,3,1),
  name="hairy", id="MA0449.1"
)

Ewg_PWM   <- pcm_to_PWMatrix(
  c(514,91,146,48,2255,79,50,37,52,260),
  c(172,2559,21,2648,297,141,26,2568,70,1826),
  c(1826,70,2568,26,141,297,2648,21,2559,172),
  c(260,52,37,50,79,2255,48,146,91,514),
  name="ewg", id="MA2309.1"
)

Mad_PWM   <- pcm_to_PWMatrix(
  c(0,35,0,35,0,0,23,0,29,29,11,29,12,32,11),
  c(48,26,14,0,102,0,65,71,0,39,53,0,48,34,28),
  c(47,13,88,56,0,102,0,31,73,34,25,73,27,16,63),
  c(7,28,0,11,0,0,14,0,0,0,13,0,15,20,0),
  name="Mad", id="MA0535.1"
)

GYNE_sig_pwm <- PWMatrixList(Brk_PWM, Hairy_PWM, Ewg_PWM, Mad_PWM)
names(GYNE_sig_pwm) <- c("brk", "hairy", "ewg", "Mad")

# ---- Worker-Specific Motif Definitions ----
Ara_PWM   <- pcm_to_PWMatrix(
  c(14,23,34,0,34), c(0,0,0,34,0), c(5,0,0,0,0), c(15,11,0,0,0),
  name="ara", id="MA0210.1"
)

Mirr_PWM  <- pcm_to_PWMatrix(
  c(20,29,41,0,41), c(1,0,0,41,0), c(1,2,0,0,0), c(19,10,0,0,0),
  name="mirr", id="MA0233.1"
)

Caup_PWM  <- pcm_to_PWMatrix(
  c(4,14,19,0,19), c(1,0,0,19,0), c(2,2,0,0,0), c(12,3,0,0,0),
  name="caup", id="MA0217.1"
)

WORKER_sig_pwm <- PWMatrixList(Ara_PWM, Mirr_PWM, Caup_PWM)
names(WORKER_sig_pwm) <- c("ara", "mirr", "caup")

# ---- Main Analysis Function ----
run_one_with_pairwise <- function(label, gr_case, gr_bg, fa, pwmA, pwmB,
                                  near_bp=NEAR_BP, p_cutoff=P_CUTOFF,
                                  n_lenbin=N_LENBIN, n_gcbin=N_GCBIN, outdir=OUTDIR){
  message(">>> Running experiment: ", label)
  lab_dir <- file.path(outdir, paste0("EXP_", label))
  dir.create(lab_dir, TRUE, TRUE)
  
  # Calculate GC content and length for matching
  case_tbl <- add_gc_len_tbl(gr_case, fa)
  bg_tbl <- add_gc_len_tbl(gr_bg, fa)
  
  # Perform 1:1 matching based on length and GC content
  pr <- strict_pair_11(case_tbl, bg_tbl, n_lenbin, n_gcbin)
  gr_c <- gr_case[pr$case_idx]
  gr_b <- gr_bg[pr$pool_idx]
  
  # Extract DNA sequences
  seq_c <- getSeq(fa, gr_c); names(seq_c) <- make_peak_id(gr_c)
  seq_b <- getSeq(fa, gr_b); names(seq_b) <- make_peak_id(gr_b)
  
  # Motif matching for both sets
  mmA_c <- matchMotifs(pwmA, seq_c, p.cutoff=p_cutoff, out="positions")
  mmB_c <- matchMotifs(pwmB, seq_c, p.cutoff=p_cutoff, out="positions")
  mmA_b <- matchMotifs(pwmA, seq_b, p.cutoff=p_cutoff, out="positions")
  mmB_b <- matchMotifs(pwmB, seq_b, p.cutoff=p_cutoff, out="positions")
  
  # Check co-occurrence within specified distance
  st_c <- co_near_from_positions(mmA_c, mmB_c, near_bp=near_bp)
  st_b <- co_near_from_positions(mmA_b, mmB_b, near_bp=near_bp)
  
  # Overall enrichment test
  a_co <- sum(st_c$flag); a_no <- length(st_c$flag)-a_co
  b_co <- sum(st_b$flag); b_no <- length(st_b$flag)-b_no
  ft <- fisher_from_counts(a_co, a_no, b_co, b_no)
  
  # Pairwise motif interactions
  hitsA_c <- extract_hits_dt(mmA_c, names(seq_c))
  hitsB_c <- extract_hits_dt(mmB_c, names(seq_c))
  hitsA_b <- extract_hits_dt(mmA_b, names(seq_b))
  hitsB_b <- extract_hits_dt(mmB_b, names(seq_b))
  
  pairs_c <- count_pairs(hitsA_c, hitsB_c, near_bp=near_bp)
  pairs_b <- count_pairs(hitsA_b, hitsB_b, near_bp=near_bp)
  
  # Calculate pairwise enrichment
  if (nrow(pairs_c)==0 && nrow(pairs_b)==0){
    pairs_enrich <- tibble(motif_A=character(), motif_B=character(), 
                          nA=integer(), nB=integer(), OR=NA_real_, P=NA_real_, FDR=NA_real_)
  } else {
    all_pairs <- full_join(rename(pairs_c, nA=n_peaks),
                          rename(pairs_b, nB=n_peaks),
                          by=c("motif_A","motif_B")) %>%
      mutate(nA=replace_na(nA,0L), nB=replace_na(nB,0L))
    
    totalA <- length(gr_c); totalB <- length(gr_b)
    pairs_enrich <- all_pairs %>%
      rowwise() %>%
      mutate(
        OR = fisher.test(matrix(c(nA,totalA-nA,nB,totalB-nB),nrow=2,byrow=TRUE))$estimate,
        P  = fisher.test(matrix(c(nA,totalA-nA,nB,totalB-nB),nrow=2,byrow=TRUE))$p.value
      ) %>%
      ungroup() %>% 
      mutate(FDR = p.adjust(P, "BH")) %>% 
      arrange(FDR, desc(OR))
  }
  
  # Save results
  write_csv(
    tibble(
      label=label, case_n=length(st_c$flag), bg_n=length(st_b$flag),
      case_co=a_co, case_no=a_no, bg_co=b_co, bg_no=b_no,
      OR=ft$OR, CI_lo=ft$CI_lo, CI_hi=ft$CI_hi, P=ft$P
    ),
    file.path(lab_dir, "overall_fisher.csv")
  )
  
  write_csv(pairs_enrich, file.path(lab_dir, "pairwise_enrichment.csv"))
  
  invisible(list(overall=ft, pairwise=pairs_enrich))
}

# ---- Summary Report Generation ----
generate_summary_report <- function() {
  message("Generating comprehensive summary report...")
  
  # Read all experiment results
  read_and_label <- function(exp_dir, pattern) {
    files <- list.files(exp_dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
    if (length(files) == 0) return(NULL)
    
    result <- lapply(files, function(f) {
      exp_name <- basename(dirname(f))
      df <- read_csv(f, show_col_types = FALSE)
      df$experiment <- exp_name
      return(df)
    }) %>% bind_rows()
    
    return(result)
  }
  
  overall_results <- read_and_label(OUTDIR, "overall_fisher.csv")
  pairwise_results <- read_and_label(OUTDIR, "pairwise_enrichment.csv")
  
  # Generate overall statistical summary
  generate_overall_summary <- function(overall_results) {
    if (is.null(overall_results)) return(NULL)
    
    col_names <- colnames(overall_results)
    p_col <- ifelse("P" %in% col_names, "P", "p")
    
    summary_table <- overall_results %>%
      mutate(
        significance = case_when(
          .data[[p_col]] < 0.001 ~ "***",
          .data[[p_col]] < 0.01 ~ "**", 
          .data[[p_col]] < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        enrichment_status = ifelse(OR > 1, "Enriched", "Depleted")
      ) %>%
      select(experiment, case_n, bg_n, case_co, bg_co, OR, 
             P = .data[[p_col]], significance, enrichment_status)
    
    return(summary_table)
  }
  
  # Generate pairwise statistical summary
  generate_pairwise_summary <- function(pairwise_results) {
    if (is.null(pairwise_results)) return(NULL)
    
    significant_pairs <- pairwise_results %>%
      filter(FDR < 0.1) %>%
      mutate(
        significance = case_when(
          FDR < 0.001 ~ "***",
          FDR < 0.01 ~ "**",
          FDR < 0.05 ~ "*",
          FDR < 0.1 ~ ".",
          TRUE ~ "ns"
        ),
        interaction_type = paste(motif_A, motif_B, sep = "::")
      ) %>%
      arrange(FDR, desc(OR))
    
    return(significant_pairs)
  }
  
  # Generate visualization plots
  generate_plots <- function(overall_results, pairwise_results) {
    plots <- list()
    
    # Overall enrichment plot
    if (!is.null(overall_results)) {
      p1 <- overall_results %>%
        mutate(experiment = factor(experiment)) %>%
        ggplot(aes(x = experiment, y = OR, fill = ifelse(OR > 1, "Enriched", "Depleted"))) +
        geom_col() +
        geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi), width = 0.2) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
        geom_text(aes(label = ifelse(P < 0.05, 
                                    ifelse(P < 0.001, "***", 
                                           ifelse(P < 0.01, "**", "*")), "")), 
                  vjust = -0.5, size = 5) +
        scale_fill_manual(values = c("Enriched" = "steelblue", "Depleted" = "firebrick")) +
        labs(title = "Overall Co-motif Enrichment", 
             x = "Experiment", y = "Odds Ratio", fill = "Status") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      plots$overall_or <- p1
    }
    
    # Pairwise interaction dot plot
    if (!is.null(pairwise_results) && nrow(pairwise_results) > 0) {
      sig_pairs <- pairwise_results %>% filter(FDR < 0.1)
      if (nrow(sig_pairs) > 0) {
        p2 <- sig_pairs %>%
          mutate(
            log10FDR = -log10(FDR),
            interaction = paste(motif_A, motif_B, sep = "::")
          ) %>%
          ggplot(aes(x = experiment, y = interaction, size = OR, color = log10FDR)) +
          geom_point(alpha = 0.7) +
          scale_size_continuous(range = c(2, 8), name = "Odds Ratio") +
          scale_color_viridis_c(name = "-log10(FDR)") +
          labs(title = "Significant Pairwise Interactions", 
               x = "Experiment", y = "Motif Pairs") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                axis.text.y = element_text(size = 8))
        
        plots$pairwise_dot <- p2
      }
    }
    
    return(plots)
  }
  
  # Generate HTML report
  generate_html_report <- function(overall_summary, pairwise_summary, plots, output_dir) {
    
    get_experiment_descriptions <- function() {
      tribble(
        ~experiment, ~description,
        "EXP_E1_gyne_vs_nonDA_ZGAxGYNEsig", "Gyne vs Non-DA: ZGA motifs × Gyne-specific motifs",
        "EXP_E2_worker_vs_nonDA_ZGAxWORKERsig", "Worker vs Non-DA: ZGA motifs × Worker-specific motifs", 
        "EXP_E3_gyne_vs_nonDA_ZGAxWORKERsig", "Gyne vs Non-DA: ZGA motifs × Worker-specific motifs",
        "EXP_E4_worker_vs_nonDA_ZGAxGYNEsig", "Worker vs Non-DA: ZGA motifs × Gyne-specific motifs",
        "EXP_E5_gyne_vs_worker_ZGAxGYNEsig", "Gyne vs Worker: ZGA motifs × Gyne-specific motifs",
        "EXP_E6_worker_vs_gyne_ZGAxWORKERsig", "Worker vs Gyne: ZGA motifs × Worker-specific motifs"
      )
    }
    
    exp_descriptions <- get_experiment_descriptions()
    
    html_content <- paste0(
      "<!DOCTYPE html>
      <html>
      <head>
        <title>Co-motif Analysis Summary</title>
        <style>
          body { font-family: Arial, sans-serif; margin: 40px; }
          h1 { color: #2c3e50; }
          h2 { color: #34495e; border-bottom: 2px solid #bdc3c7; padding-bottom: 10px; }
          table { border-collapse: collapse; width: 100%; margin: 20px 0; }
          th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
          th { background-color: #f2f2f2; }
          .enriched { background-color: #d4edda; }
          .depleted { background-color: #f8d7da; }
          .plot { margin: 30px 0; }
        </style>
      </head>
      <body>
        <h1>Co-motif Analysis Summary Report</h1>
        <p>Generated on: ", Sys.Date(), "</p>
        
        <h2>Experiment Descriptions</h2>
        <table>
          <tr><th>Experiment</th><th>Description</th></tr>",
          paste(sapply(1:nrow(exp_descriptions), function(i) {
            paste0("<tr><td>", exp_descriptions$experiment[i], "</td><td>", 
                   exp_descriptions$description[i], "</td></tr>")
          }), collapse = ""),
        "</table>"
    )
    
    # Add overall results table
    if (!is.null(overall_summary)) {
      html_content <- paste0(html_content, "
        <h2>Overall Co-motif Enrichment Results</h2>
        <table>
          <tr>
            <th>Experiment</th><th>Case Peaks</th><th>Background Peaks</th>
            <th>Case Co-occur</th><th>BG Co-occur</th><th>OR</th><th>P-value</th><th>Status</th>
          </tr>",
          paste(sapply(1:nrow(overall_summary), function(i) {
            row_class <- ifelse(overall_summary$enrichment_status[i] == "Enriched", "enriched", "depleted")
            paste0("<tr class='", row_class, "'><td>", overall_summary$experiment[i], "</td>",
                   "<td>", overall_summary$case_n[i], "</td>",
                   "<td>", overall_summary$bg_n[i], "</td>",
                   "<td>", overall_summary$case_co[i], "</td>",
                   "<td>", overall_summary$bg_co[i], "</td>",
                   "<td>", round(overall_summary$OR[i], 3), "</td>",
                   "<td>", format.pval(overall_summary$P[i], digits = 3), " ", overall_summary$significance[i], "</td>",
                   "<td>", overall_summary$enrichment_status[i], "</td></tr>")
          }), collapse = ""),
        "</table>"
      )
    }
    
    # Add pairwise results table
    if (!is.null(pairwise_summary) && nrow(pairwise_summary) > 0) {
      html_content <- paste0(html_content, "
        <h2>Significant Pairwise Interactions (FDR < 0.1)</h2>
        <table>
          <tr>
            <th>Experiment</th><th>Motif A</th><th>Motif B</th>
            <th>Case Count</th><th>BG Count</th><th>OR</th><th>FDR</th><th>Significance</th>
          </tr>",
          paste(sapply(1:nrow(pairwise_summary), function(i) {
            paste0("<tr><td>", pairwise_summary$experiment[i], "</td>",
                   "<td>", pairwise_summary$motif_A[i], "</td>",
                   "<td>", pairwise_summary$motif_B[i], "</td>",
                   "<td>", pairwise_summary$nA[i], "</td>",
                   "<td>", pairwise_summary$nB[i], "</td>",
                   "<td>", round(pairwise_summary$OR[i], 3), "</td>",
                   "<td>", format.pval(pairwise_summary$FDR[i], digits = 3), "</td>",
                   "<td>", pairwise_summary$significance[i], "</td></tr>")
          }), collapse = ""),
        "</table>"
      )
    }
    
    html_content <- paste0(html_content, "
      </body>
      </html>"
    )
    
    writeLines(html_content, file.path(output_dir, "co_motif_summary_report.html"))
  }
  
  # Execute summary generation
  overall_summary <- generate_overall_summary(overall_results)
  pairwise_summary <- generate_pairwise_summary(pairwise_results)
  plots <- generate_plots(overall_results, pairwise_results)
  
  # Save summary tables
  if (!is.null(overall_summary)) {
    write_csv(overall_summary, file.path(SUMMARY_DIR, "overall_summary.csv"))
    message("✅ Overall summary table saved")
  }
  
  if (!is.null(pairwise_summary) && nrow(pairwise_summary) > 0) {
    write_csv(pairwise_summary, file.path(SUMMARY_DIR, "significant_pairwise_interactions.csv"))
    message("✅ Significant pairwise interactions table saved")
  }
  
  # Save plots
  if (length(plots) > 0) {
    if (!is.null(plots$overall_or)) {
      ggsave(file.path(SUMMARY_DIR, "overall_enrichment_plot.pdf"), 
             plots$overall_or, width = 10, height = 6)
    }
    if (!is.null(plots$pairwise_dot)) {
      ggsave(file.path(SUMMARY_DIR, "significant_pairwise_interactions.pdf"), 
             plots$pairwise_dot, width = 12, height = 8)
    }
    message("✅ Visualization plots saved")
  }
  
  # Generate HTML report
  generate_html_report(overall_summary, pairwise_summary, plots, SUMMARY_DIR)
  message("✅ HTML summary report generated")
  
  # Console output
  message(paste0("\n", paste(rep("=", 60), collapse = ""), "\n"))
  message("📊 CO-MOTIF ANALYSIS SUMMARY")
  message(paste(rep("=", 60), collapse = ""))
  
  if (!is.null(overall_summary)) {
    message("\n📈 OVERALL ENRICHMENT RESULTS:")
    for (i in 1:nrow(overall_summary)) {
      exp <- overall_summary$experiment[i]
      or <- round(overall_summary$OR[i], 3)
      pval <- format.pval(overall_summary$P[i], digits = 3)
      status <- overall_summary$enrichment_status[i]
      sig <- overall_summary$significance[i]
      
      message(sprintf("  %-40s: OR = %5.3f, P = %7s %2s (%s)", 
                     exp, or, pval, sig, status))
    }
  }
  
  if (!is.null(pairwise_summary) && nrow(pairwise_summary) > 0) {
    message("\n🔗 SIGNIFICANT PAIRWISE INTERACTIONS (FDR < 0.1):")
    top_pairs <- pairwise_summary %>% head(10)
    for (i in 1:nrow(top_pairs)) {
      exp <- top_pairs$experiment[i]
      pair <- paste(top_pairs$motif_A[i], top_pairs$motif_B[i], sep = "::")
      or <- round(top_pairs$OR[i], 3)
      fdr <- format.pval(top_pairs$FDR[i], digits = 3)
      
      message(sprintf("  %-30s %-20s: OR = %5.3f, FDR = %7s", 
                     exp, pair, or, fdr))
    }
    if (nrow(pairwise_summary) > 10) {
      message(sprintf("  ... and %d more significant pairs", nrow(pairwise_summary) - 10))
    }
  } else {
    message("\n🔗 No significant pairwise interactions found (FDR < 0.1)")
  }
  
  message("\n📁 Results saved in: ", SUMMARY_DIR)
  message("🌐 Open the HTML report: ", file.path(SUMMARY_DIR, "co_motif_summary_report.html"))
  message(paste(rep("=", 60), collapse = ""))
}

# ---- Main Execution Pipeline ----
message("🚀 Starting Co-motif analysis...")

# Read peaks and genome
gr_gy <- read_bed(BED_GYNE)
gr_wk <- read_bed(BED_WORKER)
gr_nd <- read_bed(BED_NONDA)
fa <- FaFile(GENOME_FA)
open.FaFile(fa)

message(sprintf("Loaded %d gyne peaks, %d worker peaks, %d non-DA peaks", 
                length(gr_gy), length(gr_wk), length(gr_nd)))

# Run comparative analyses
message("\nRunning comparative analyses...")
run_one_with_pairwise("E1_gyne_vs_nonDA_ZGAxGYNEsig",    gr_gy, gr_nd, fa, ZGA_pwm, GYNE_sig_pwm)
run_one_with_pairwise("E2_worker_vs_nonDA_ZGAxWORKERsig", gr_wk, gr_nd, fa, ZGA_pwm, WORKER_sig_pwm)
run_one_with_pairwise("E3_gyne_vs_nonDA_ZGAxWORKERsig",   gr_gy, gr_nd, fa, ZGA_pwm, WORKER_sig_pwm)
run_one_with_pairwise("E4_worker_vs_nonDA_ZGAxGYNEsig",   gr_wk, gr_nd, fa, ZGA_pwm, GYNE_sig_pwm)

if (RUN_CASE_VS_CASE){
  message("\nRunning case vs case comparisons...")
  run_one_with_pairwise("E5_gyne_vs_worker_ZGAxGYNEsig",   gr_gy, gr_wk, fa, ZGA_pwm, GYNE_sig_pwm)
  run_one_with_pairwise("E6_worker_vs_gyne_ZGAxWORKERsig", gr_wk, gr_gy, fa, ZGA_pwm, WORKER_sig_pwm)
}

message("✅ Co-motif analysis completed")
message("📊 Generating summary report...")

# Generate comprehensive summary report
generate_summary_report()

message("🎉 Analysis pipeline completed successfully!")
