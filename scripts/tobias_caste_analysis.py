import pandas as pd
from pathlib import Path
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration Section
# ============================================================================

# Set up file paths
DATA_DIR = Path("/data/work/tobias_calculate/")
OUTPUT_DIR = Path("/data/work/tobias_calculate/score")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Define all developmental stage files
# Each file corresponds to TOBIAS BINDetect analysis results for a specific developmental stage
STAGE_FILES = {
    "larva_2nd": DATA_DIR / "larva_2nd.txt",      # 2nd instar larval stage
    "larva_3rd": DATA_DIR / "larva_3rd.txt",      # 3rd instar larval stage
    "young_pupa": DATA_DIR / "young_pupa.txt",    # Young pupal stage (Note: prepupal stage excluded due to binding bias reversal)
    "old_pupa": DATA_DIR / "old_pupa.txt",        # Old pupal stage
    "adult": DATA_DIR / "adult.txt"               # Adult stage
}

# File format explanation (using larva_2nd.txt as example):
# - Each file contains TOBIAS BINDetect comparison of TF binding activity between gyne and worker samples
# - Key columns:
#   * motif_id: Unique identifier for transcription factor binding motif (e.g., MA2306.1)
#   * name: Abbreviation or name of transcription factor (e.g., D1, br, etc.)
#   * total_tfbs: Total number of predicted binding sites for this motif in the genome
#   * larva_2nd_gyne_mean_score: Average footprint score in gyne samples (reflects binding strength)
#   * larva_2nd_gyne_bound: Number of predicted bound sites in gyne samples
#   * larva_2nd_worker_mean_score: Average footprint score in worker samples
#   * larva_2nd_worker_bound: Number of predicted bound sites in worker samples
#   * larva_2nd_gyne_larva_2nd_worker_change: Differential binding activity (positive = gyne-biased, negative = worker-biased)
#   * larva_2nd_gyne_larva_2nd_worker_pvalue: Statistical significance p-value for differential binding
# Note: All stage files share the same column structure, differing only in the stage prefix

# Stage weights - larval > pupal > adult, reflecting gradient of caste differentiation potential
# Higher weights assigned to early developmental stages with greater plasticity in caste fate determination

STAGE_WEIGHTS = {
    'larva_2nd': 2.5,     
    'larva_3rd': 2.5,    
    'young_pupa': 1.8,   
    'old_pupa': 1.8,     
    'adult': 1.5          
}

# Scoring weights for different criteria (sum = 1.0)
SCORING_WEIGHTS = {
    'tfbs': 0.2,           # Weight for binding site abundance
    'mean_score': 0.1,     # Weight for binding strength
    'bound_sites': 0.1,    # Weight for proportion of bound sites
    'change': 0.1,         # Weight for differential activity magnitude
    'pvalue': 0.1,         # Weight for statistical significance
    'coverage': 0.4        # Weight for temporal persistence (developmental coverage)
}

# Adjustment factors for motif patterns
ADJUSTMENT_FACTORS = {
    'both': 0.7,           # Penalty factor for motifs detected in both lineages
    'specific': 1.1        # Bonus factor for lineage-specific motifs
}

# Bias score thresholds for classification
BIAS_THRESHOLDS = {
    'strong_gyne': 0.3,
    'gyne': 0.1,
    'strong_worker': -0.3,
    'worker': -0.1
}

# ============================================================================
# Core Functions
# ============================================================================

def calculate_bias_score(gyne_weighted_score: float, 
                        worker_weighted_score: float, 
                        total_weighted_score: float) -> tuple:
    """
    Calculate bias score for motif classification.
    
    Args:
        gyne_weighted_score: Weighted score for gyne lineage
        worker_weighted_score: Weighted score for worker lineage
        total_weighted_score: Total weighted score across stages
    
    Returns:
        tuple: (bias_score, bias_direction)
        - bias_score: Range from -1 to 1, positive values indicate gyne bias
        - bias_direction: Classification string
    """
    if total_weighted_score == 0:
        return 0.0, 'neutral'
    
    # Calculate normalized bias score (-1 to 1)
    bias_score = (gyne_weighted_score - worker_weighted_score) / total_weighted_score
    
    # Classify based on thresholds
    if bias_score > BIAS_THRESHOLDS['strong_gyne']:
        bias_direction = 'strong_gyne'
    elif bias_score > BIAS_THRESHOLDS['gyne']:
        bias_direction = 'gyne'
    elif bias_score < BIAS_THRESHOLDS['strong_worker']:
        bias_direction = 'strong_worker'
    elif bias_score < BIAS_THRESHOLDS['worker']:
        bias_direction = 'worker'
    else:
        bias_direction = 'neutral'
    
    return bias_score, bias_direction


def normalize_score(value: float, max_value: float, log_transform: bool = False) -> float:
    """
    Normalize a score to range 0-1.
    
    Args:
        value: Raw score value
        max_value: Maximum value for normalization
        log_transform: Whether to apply log transformation
    
    Returns:
        Normalized score between 0 and 1
    """
    if value <= 0:
        return 0.0
    
    if log_transform:
        value = np.log1p(value)
        max_value = np.log1p(max_value)
    
    return min(value / max_value, 1.0)


def calculate_motif_score(row: dict, pattern_type: str, total_weighted_stages: float) -> float:
    """
    Calculate composite score for a motif based on multiple criteria.
    
    Args:
        row: Dictionary containing motif data
        pattern_type: Classification of motif pattern
        total_weighted_stages: Sum of all stage weights for normalization
    
    Returns:
        Composite score for the motif
    """
    # Extract scoring components
    total_tfbs = row['total_tfbs']
    mean_score = row['mean_score']
    bound_sites = row['bound_sites']
    change_value = abs(row['change_value'])
    pvalue = -np.log10(row['pvalue'] + 1e-300)
    
    # Calculate weighted stage coverage
    weighted_coverage = row['weighted_stage_score'] / total_weighted_stages
    
    # Normalize each component
    tfbs_score = normalize_score(total_tfbs, 100000, log_transform=True)
    mean_score_norm = normalize_score(mean_score, 0.3)
    bound_score = bound_sites / total_tfbs if total_tfbs > 0 else 0.0
    change_score = normalize_score(change_value, 0.5)
    pvalue_score = normalize_score(pvalue, 50)
    coverage_score = weighted_coverage
    
    # Calculate base score using weighted sum
    base_score = (
        SCORING_WEIGHTS['tfbs'] * tfbs_score +
        SCORING_WEIGHTS['mean_score'] * mean_score_norm +
        SCORING_WEIGHTS['bound_sites'] * bound_score +
        SCORING_WEIGHTS['change'] * change_score +
        SCORING_WEIGHTS['pvalue'] * pvalue_score +
        SCORING_WEIGHTS['coverage'] * coverage_score
    )
    
    # Apply adjustment based on pattern type
    adjustment = ADJUSTMENT_FACTORS['specific'] if pattern_type != 'both' else ADJUSTMENT_FACTORS['both']
    final_score = base_score * adjustment
    
    return final_score


def load_and_validate_stage_data(file_path: Path, stage_name: str) -> pd.DataFrame:
    """
    Load and validate stage-specific data file.
    
    Args:
        file_path: Path to the data file
        stage_name: Name of the developmental stage
    
    Returns:
        DataFrame with validated data, or None if loading fails
    """
    if not file_path.exists():
        print(f"Warning: File {file_path} does not exist, skipping {stage_name} stage")
        return None
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        
        # Validate required columns
        required_columns = ['motif_id', 'name', 'total_tfbs']
        if not all(col in df.columns for col in required_columns):
            print(f"Warning: Missing required columns in {file_path}")
            return None
        
        # Find change column
        change_cols = [col for col in df.columns if '_change' in col]
        if not change_cols:
            print(f"Warning: No _change column found in {file_path}")
            return None
        
        return df
    
    except Exception as e:
        print(f"Error loading file {file_path}: {str(e)}")
        return None


def aggregate_motif_data(stage_files: dict, stage_weights: dict) -> dict:
    """
    Aggregate motif data across all developmental stages.
    
    Args:
        stage_files: Dictionary mapping stage names to file paths
        stage_weights: Dictionary of stage-specific weights
    
    Returns:
        Dictionary containing aggregated data for all motifs
    """
    all_motifs_data = {}
    
    print("=" * 80)
    print("Aggregating TOBIAS BINDetect results across developmental stages")
    print("=" * 80)
    
    for stage_name, file_path in stage_files.items():
        stage_weight = stage_weights.get(stage_name, 1.0)
        
        df = load_and_validate_stage_data(file_path, stage_name)
        if df is None:
            continue
        
        # Find relevant columns
        change_col = next((col for col in df.columns if '_change' in col), None)
        mean_score_col = next((col for col in df.columns if '_mean_score' in col), None)
        bound_col = next((col for col in df.columns if '_bound' in col), None)
        pvalue_col = next((col for col in df.columns if '_pvalue' in col), None)
        
        if not all([change_col, mean_score_col, bound_col, pvalue_col]):
            print(f"Warning: Missing required measurement columns in {file_path}")
            continue
        
        print(f"Processing {stage_name}: {len(df)} motifs")
        
        for _, row in df.iterrows():
            motif_id = row['motif_id']
            motif_name = row['name']
            
            # Initialize data structure for new motifs
            if motif_id not in all_motifs_data:
                all_motifs_data[motif_id] = {
                    'name': motif_name,
                    'total_tfbs': 0.0,
                    'weighted_mean_score': 0.0,
                    'weighted_bound_sites': 0.0,
                    'weighted_change': 0.0,
                    'weighted_pvalue': 0.0,
                    'weighted_stage_score': 0.0,
                    'gyne_weighted_score': 0.0,
                    'worker_weighted_score': 0.0,
                    'stages': set(),
                    'gyne_stages': set(),
                    'worker_stages': set()
                }
            
            # Extract values with safe defaults
            total_tfbs = float(row.get('total_tfbs', 0))
            mean_score_val = float(row.get(mean_score_col, 0))
            bound_sites_val = float(row.get(bound_col, 0))
            change_val = float(row.get(change_col, 0))
            pvalue_val = float(row.get(pvalue_col, 1))
            
            # Update aggregated statistics with stage weighting
            all_motifs_data[motif_id]['total_tfbs'] += total_tfbs * stage_weight
            all_motifs_data[motif_id]['weighted_mean_score'] += mean_score_val * stage_weight
            all_motifs_data[motif_id]['weighted_bound_sites'] += bound_sites_val * stage_weight
            all_motifs_data[motif_id]['weighted_change'] += abs(change_val) * stage_weight
            all_motifs_data[motif_id]['weighted_pvalue'] += -np.log10(pvalue_val + 1e-300) * stage_weight
            all_motifs_data[motif_id]['weighted_stage_score'] += stage_weight
            
            # Track stage presence
            all_motifs_data[motif_id]['stages'].add(stage_name)
            
            # Classify lineage bias and track
            if change_val > 0:
                all_motifs_data[motif_id]['gyne_stages'].add(stage_name)
                all_motifs_data[motif_id]['gyne_weighted_score'] += stage_weight
            elif change_val < 0:
                all_motifs_data[motif_id]['worker_stages'].add(stage_name)
                all_motifs_data[motif_id]['worker_weighted_score'] += stage_weight
    
    return all_motifs_data


def classify_motif_pattern(gyne_count: int, worker_count: int, 
                          gyne_weighted: float, worker_weighted: float,
                          total_weighted: float) -> tuple:
    """
    Classify motif pattern based on lineage specificity and bias.
    
    Args:
        gyne_count: Number of stages with gyne bias
        worker_count: Number of stages with worker bias
        gyne_weighted: Weighted gyne score
        worker_weighted: Weighted worker score
        total_weighted: Total weighted score
    
    Returns:
        tuple: (pattern_type, final_pattern)
    """
    # Exclusive patterns
    if gyne_count > 0 and worker_count == 0:
        return 'gyne_only', 'gyne_only'
    elif worker_count > 0 and gyne_count == 0:
        return 'worker_only', 'worker_only'
    
    # Mixed patterns - calculate bias
    pattern_type = 'both'
    bias_score, bias_direction = calculate_bias_score(
        gyne_weighted, worker_weighted, total_weighted
    )
    
    # Determine final pattern based on bias direction
    if bias_direction in ['strong_gyne', 'gyne']:
        final_pattern = 'gyne_biased'
    elif bias_direction in ['strong_worker', 'worker']:
        final_pattern = 'worker_biased'
    else:
        final_pattern = 'neutral'
    
    return pattern_type, final_pattern


def integrate_and_score_motifs(stage_files: dict, stage_weights: dict) -> tuple:
    """
    Main function to integrate TOBIAS results across stages and score motifs.
    
    Args:
        stage_files: Dictionary mapping stage names to file paths
        stage_weights: Dictionary of stage-specific weights
    
    Returns:
        tuple: (all_motifs_df, gyne_motifs_df, worker_motifs_df)
    """
    # Aggregate data across all stages
    all_motifs_data = aggregate_motif_data(stage_files, stage_weights)
    
    if not all_motifs_data:
        print("No motif data available for analysis")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    total_weighted_stages = sum(stage_weights.values())
    
    # Prepare data structures for scoring
    scored_motifs = []
    gyne_motifs = []
    worker_motifs = []
    
    print(f"\nScoring {len(all_motifs_data)} motifs...")
    
    for motif_id, data in all_motifs_data.items():
        weighted_stage_score = data['weighted_stage_score']
        if weighted_stage_score == 0:
            continue
        
        # Calculate weighted averages
        avg_tfbs = data['total_tfbs'] / weighted_stage_score
        avg_mean_score = data['weighted_mean_score'] / weighted_stage_score
        avg_bound_sites = data['weighted_bound_sites'] / weighted_stage_score
        avg_change = data['weighted_change'] / weighted_stage_score
        avg_pvalue_log = data['weighted_pvalue'] / weighted_stage_score
        
        # Classify motif pattern
        gyne_count = len(data['gyne_stages'])
        worker_count = len(data['worker_stages'])
        
        pattern_type, final_pattern = classify_motif_pattern(
            gyne_count, worker_count,
            data['gyne_weighted_score'], data['worker_weighted_score'],
            weighted_stage_score
        )
        
        # Prepare scoring data
        score_data = {
            'motif_id': motif_id,
            'name': data['name'],
            'pattern_type': pattern_type,
            'final_pattern': final_pattern,
            'n_stages': len(data['stages']),
            'weighted_stage_score': weighted_stage_score,
            'gyne_stages': gyne_count,
            'worker_stages': worker_count,
            'gyne_weighted': data['gyne_weighted_score'],
            'worker_weighted': data['worker_weighted_score'],
            'total_tfbs': avg_tfbs,
            'mean_score': avg_mean_score,
            'bound_sites': avg_bound_sites,
            'change_value': avg_change,
            'pvalue': np.power(10, -avg_pvalue_log) if avg_pvalue_log > 0 else 1.0,
            'gyne_stage_list': ', '.join(sorted(data['gyne_stages'])),
            'worker_stage_list': ', '.join(sorted(data['worker_stages'])),
            'all_stages': ', '.join(sorted(data['stages']))
        }
        
        # Calculate composite score
        composite_score = calculate_motif_score(score_data, pattern_type, total_weighted_stages)
        score_data['composite_score'] = composite_score
        
        scored_motifs.append(score_data)
        
        # Sort into lineage-specific lists
        if final_pattern in ['gyne_only', 'gyne_biased']:
            gyne_motifs.append(score_data.copy())
        elif final_pattern in ['worker_only', 'worker_biased']:
            worker_motifs.append(score_data.copy())
    
    # Create DataFrames
    scored_df = pd.DataFrame(scored_motifs)
    gyne_df = pd.DataFrame(gyne_motifs)
    worker_df = pd.DataFrame(worker_motifs)
    
    # Sort by composite score
    scored_df = scored_df.sort_values('composite_score', ascending=False).reset_index(drop=True)
    gyne_df = gyne_df.sort_values('composite_score', ascending=False).reset_index(drop=True)
    worker_df = worker_df.sort_values('composite_score', ascending=False).reset_index(drop=True)
    
    # Save results
    scored_df.to_csv(OUTPUT_DIR / "all_motifs_scored.csv", index=False)
    gyne_df.to_csv(OUTPUT_DIR / "gyne_unique_motifs_ranked.csv", index=False)
    worker_df.to_csv(OUTPUT_DIR / "worker_unique_motifs_ranked.csv", index=False)
    
    # Print summary statistics
    print_summary(scored_df, gyne_df, worker_df)
    
    return scored_df, gyne_df, worker_df


def print_summary(all_df: pd.DataFrame, gyne_df: pd.DataFrame, worker_df: pd.DataFrame) -> None:
    """
    Print comprehensive analysis summary.
    
    Args:
        all_df: DataFrame with all scored motifs
        gyne_df: DataFrame with gyne-specific motifs
        worker_df: DataFrame with worker-specific motifs
    """
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)
    
    print(f"\nTotal motifs analyzed: {len(all_df):,}")
    print(f"Gyne-specific/biased motifs: {len(gyne_df):,}")
    print(f"Worker-specific/biased motifs: {len(worker_df):,}")
    
    # Pattern distribution
    pattern_counts = all_df['final_pattern'].value_counts()
    print("\nPattern Distribution:")
    for pattern, count in pattern_counts.items():
        percentage = (count / len(all_df)) * 100
        print(f"  {pattern:15s}: {count:5d} motifs ({percentage:5.1f}%)")
    
    # Top candidates
    print("\n" + "-" * 40)
    print("TOP GYNE CANDIDATES (Top 10)")
    print("-" * 40)
    if len(gyne_df) > 0:
        for i, (_, row) in enumerate(gyne_df.head(10).iterrows(), 1):
            print(f"{i:2d}. {row['name']:12s} ({row['motif_id']:10s}): "
                  f"Score={row['composite_score']:.4f} "
                  f"[{row['final_pattern']}, Stages={row['n_stages']}]")
    else:
        print("No gyne-specific motifs found")
    
    print("\n" + "-" * 40)
    print("TOP WORKER CANDIDATES (Top 10)")
    print("-" * 40)
    if len(worker_df) > 0:
        for i, (_, row) in enumerate(worker_df.head(10).iterrows(), 1):
            print(f"{i:2d}. {row['name']:12s} ({row['motif_id']:10s}): "
                  f"Score={row['composite_score']:.4f} "
                  f"[{row['final_pattern']}, Stages={row['n_stages']}]")
    else:
        print("No worker-specific motifs found")
    
    # Score statistics
    print("\n" + "-" * 40)
    print("SCORE STATISTICS")
    print("-" * 40)
    print(f"Mean composite score: {all_df['composite_score'].mean():.4f}")
    print(f"Median composite score: {all_df['composite_score'].median():.4f}")
    print(f"Score range: [{all_df['composite_score'].min():.4f}, "
          f"{all_df['composite_score'].max():.4f}]")
    
    # Stage coverage statistics
    print(f"\nAverage stages per motif: {all_df['n_stages'].mean():.1f}")
    print(f"Motifs present in all 5 stages: {(all_df['n_stages'] == 5).sum()}")
    
    print("\nAnalysis complete! Results saved to:")
    print(f"  All motifs: {OUTPUT_DIR / 'all_motifs_scored.csv'}")
    print(f"  Gyne motifs: {OUTPUT_DIR / 'gyne_unique_motifs_ranked.csv'}")
    print(f"  Worker motifs: {OUTPUT_DIR / 'worker_unique_motifs_ranked.csv'}")


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("TOBIAS MOTIF SCORING AND CANALIZATION ANALYSIS")
    print("=" * 80)
    print(f"Data directory: {DATA_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Developmental stages analyzed: {', '.join(STAGE_FILES.keys())}")
    
    # Execute analysis
    scored_df, gyne_df, worker_df = integrate_and_score_motifs(STAGE_FILES, STAGE_WEIGHTS)
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETED SUCCESSFULLY")
    print("=" * 80)
