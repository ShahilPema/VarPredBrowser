#!/usr/bin/env python3
"""
Generate axis tables for ALL chromosomes with two filter options:
1. any_count_gt0: All positions where rgc_any_count > 0
2. mis_count_gt0: Positions where rgc_mis_count > 0 (missense possible)

This creates compressed coordinate mappings for both viewing modes.

Includes:
- ClinVar annotations
- Training labels
- dbNSFP scores (AlphaMissense, MTR, CCR, ESM1b, AlphaSync)
- gnomAD annotations (exomes and genomes over 20)
- Constraint predictions
- Percentiled constraint metrics
- Domain annotations (array format)

Usage:
    python scripts/preprocess_browser_data.py [--input PATH] [--output PATH]
"""

import argparse
import polars as pl
from pathlib import Path
import sys

# Try to import from browser config
try:
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from browser.backend.config import get_config, get_data_dir
    _config = get_config()
    DEFAULT_INPUT = _config.get('data_sources', {}).get('browser_data', '../rgc_browser_data_merged.parquet')
    DEFAULT_OUTPUT = get_data_dir()
except ImportError:
    DEFAULT_INPUT = '../rgc_browser_data_merged.parquet'
    DEFAULT_OUTPUT = Path('data')

# Filter definitions
FILTERS = {
    'all_sites': {
        'description': 'All positions where any variant is possible',
        'condition': pl.col('rgc_any_count') > 0
    },
    'missense_only': {
        'description': 'Positions where missense variants are possible',
        'condition': pl.col('rgc_mis_count') > 0
    },
    'synonymous_only': {
        'description': 'Positions where synonymous variants are possible',
        'condition': pl.col('rgc_syn_count') > 0
    }
}

# ClinVar columns (new format: array of structs)
CLINVAR_VARIANTS_COLUMN = 'clinvar_variants'  # array<struct{alt, significance, status, mol_csq, variation_id}>

# Legacy ClinVar columns (for backwards compatibility)
CLINVAR_LEGACY_COLUMNS = [
    'clinvar.clinvar_count',
    'clinvar.clinvar_label_list',
    'clinvar.clinvar_status_list',
    'clinvar.clinvar_var_type_list',
]

# Training label columns (flattened from training.train_counts.*)
TRAINING_COLUMNS = [
    'train_unlabelled',
    'train_labelled',
    'train_unlabelled_high_qual',
    'train_labelled_high_qual',
]

# dbNSFP score columns (flattened from dbnsfp.max_*)
# Percentiles are calculated in Polars, not stored here
DBNSFP_SCORE_COLUMNS = [
    'max_AlphaMissense_am_pathogenicity',
    'max_ESM1b_score',
    'max_RGC_MTR_MTR',
    'max_Non_Neuro_CCR_resid_pctile',
    'max_AlphaSync_plddt',
    'max_AlphaSync_plddt10',
    'max_AlphaSync_relasa',
    'max_AlphaSync_relasa10',
]

# gnomAD coverage columns (expanded per browser-data-refactor.md Task 2)
# 12 total: 6 thresholds x 2 cohorts
GNOMAD_COVERAGE_COLUMNS = [
    'gnomad_exomes_over_10',
    'gnomad_exomes_over_15',
    'gnomad_exomes_over_20',
    'gnomad_exomes_over_25',
    'gnomad_exomes_over_30',
    'gnomad_exomes_over_50',
    'gnomad_genomes_over_10',
    'gnomad_genomes_over_15',
    'gnomad_genomes_over_20',
    'gnomad_genomes_over_25',
    'gnomad_genomes_over_30',
    'gnomad_genomes_over_50',
]

# Variant frequency array columns (browser-data-refactor.md Task 4)
VARIANT_FREQUENCY_COLUMNS = [
    'rgc_variants',           # array<struct{alt, af, ac, an, filters}>
    'gnomad_exomes_variants', # array<struct{alt, af, ac, an, filters}>
    'gnomad_genomes_variants', # array<struct{alt, af, ac, an, filters}>
]

# phyloP conservation scores (browser-data-refactor.md Section 4)
PHYLOP_COLUMNS = [
    'phylop_scores_447way',
    'phylop_scores_100way',
]

# Constraint prediction column (from AOU_RGC_All_preds.ht)
# Format: Array of structs, one per variant at locus (REFACTORED from tuples)
# Each struct: {alt: str, pred: float32, n_pred: int32}
CONSTRAINT_PREDS_COLUMN = ['Constraint', 'Core', 'Complete']

# Stacked scores columns (variant-level data with struct format)
STACKED_SCORE_COLUMNS = [
    'AlphaMissense_stacked',  # array<struct{alt: str, score: float64, percentile: float64}>
    'ESM1b_stacked',          # array<struct{alt: str, score: float64, percentile: float64}>
]

# Variant consequences column
VARIANT_CONSEQUENCES_COLUMN = 'variant_consequences'  # array<struct{alt: str, csq: str}>

# Domains column (now array format)
DOMAIN_COLUMN = 'domains'

# Chromosome sort order for proper genomic ordering
CHROM_ORDER = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']


def detect_columns(df):
    """Detect standard column names in the dataframe."""
    chrom_col = 'chrom' if 'chrom' in df.columns else 'locus.contig'
    pos_col = 'pos' if 'pos' in df.columns else 'locus.position'
    
    gene_col = None
    for col in df.columns:
        if 'HGNC' in col:
            gene_col = col
            break
    
    return chrom_col, pos_col, gene_col


def get_columns_to_keep(df, chrom_col, pos_col, gene_col):
    """Determine which columns to include in the output."""
    columns_to_keep = ['filtered_idx', chrom_col, pos_col]
    
    # Add gene column if found
    if gene_col:
        columns_to_keep.append(gene_col)
    
    # Add transcript_id if present
    if 'transcript_id' in df.columns:
        columns_to_keep.append('transcript_id')
    
    # Add aa_pos if present
    if 'aa_pos' in df.columns:
        columns_to_keep.append('aa_pos')
    
    # Add ALL rgc_ columns (raw constraint metrics)
    rgc_columns = [col for col in df.columns if col.startswith('rgc_')]
    columns_to_keep.extend(rgc_columns)
    
    # Add ClinVar variants column (new format)
    if CLINVAR_VARIANTS_COLUMN in df.columns:
        columns_to_keep.append(CLINVAR_VARIANTS_COLUMN)
        clinvar_count = 1
    else:
        # Fall back to legacy format
        clinvar_cols = [col for col in CLINVAR_LEGACY_COLUMNS if col in df.columns]
        columns_to_keep.extend(clinvar_cols)
        clinvar_count = len(clinvar_cols)
    
    # Add Training columns
    training_cols = [col for col in TRAINING_COLUMNS if col in df.columns]
    columns_to_keep.extend(training_cols)
    
    # Add dbNSFP score columns
    dbnsfp_score_cols = [col for col in DBNSFP_SCORE_COLUMNS if col in df.columns]
    columns_to_keep.extend(dbnsfp_score_cols)
    
    # Add gnomAD coverage columns (12 total)
    gnomad_cov_cols = [col for col in GNOMAD_COVERAGE_COLUMNS if col in df.columns]
    columns_to_keep.extend(gnomad_cov_cols)

    # Add variant frequency array columns (3 total)
    variant_freq_cols = [col for col in VARIANT_FREQUENCY_COLUMNS if col in df.columns]
    columns_to_keep.extend(variant_freq_cols)

    # Add phyloP conservation scores (Section 4)
    phylop_cols = [col for col in PHYLOP_COLUMNS if col in df.columns]
    columns_to_keep.extend(phylop_cols)

    # Auto-detect gnomAD constraint columns (Task 5) - pattern: gnomad_*_oe, gnomad_*vir*
    # These columns have gnomad_ prefix and contain oe or vir metrics
    gnomad_constraint_cols = [col for col in df.columns
                              if col.startswith('gnomad_') and ('_oe' in col or 'vir' in col)]
    columns_to_keep.extend(gnomad_constraint_cols)

    # Add constraint predictions column (array of structs, REFACTORED from tuples)
    prediction_cols = [col for col in CONSTRAINT_PREDS_COLUMN if col in df.columns]
    columns_to_keep.extend(prediction_cols)

    # Add stacked score columns (variant-level with struct format)
    stacked_cols = [col for col in STACKED_SCORE_COLUMNS if col in df.columns]
    columns_to_keep.extend(stacked_cols)

    # Add variant consequences column
    if VARIANT_CONSEQUENCES_COLUMN in df.columns:
        columns_to_keep.append(VARIANT_CONSEQUENCES_COLUMN)

    # Auto-detect and add percentile columns (exome and cross-norm)
    perc_columns = [c for c in df.columns if c.endswith('_exome_perc') or c.endswith('_cross_norm_perc')]
    columns_to_keep.extend(perc_columns)

    # Add domains column
    if DOMAIN_COLUMN in df.columns:
        columns_to_keep.append(DOMAIN_COLUMN)

    # Filter to columns that exist and remove duplicates
    columns_to_keep = list(dict.fromkeys([
        col for col in columns_to_keep if col in df.columns
    ]))

    return columns_to_keep, {
        'rgc': len(rgc_columns),
        'clinvar': clinvar_count,
        'training': len(training_cols),
        'dbnsfp': len(dbnsfp_score_cols),
        'gnomad_coverage': len(gnomad_cov_cols),
        'variant_freq': len(variant_freq_cols),
        'phylop': len(phylop_cols),
        'gnomad_constraint': len(gnomad_constraint_cols),
        'constraint_preds': len(prediction_cols),
        'stacked_scores': len(stacked_cols),
        'variant_consequences': 1 if VARIANT_CONSEQUENCES_COLUMN in df.columns else 0,
        'percentiles': len(perc_columns),
        'domains': 1 if DOMAIN_COLUMN in df.columns else 0
    }


def process_filter(df, filter_id, filter_config, chrom_col, pos_col, gene_col, chrom_order_map, output_dir: Path):
    """Process a single filter and generate axis table + gene index."""
    print(f"\n{'='*60}")
    print(f"Processing filter: {filter_id}")
    print(f"  {filter_config['description']}")
    print(f"{'='*60}")
    
    # Apply filter
    print(f"\nApplying filter...")
    filtered_df = df.filter(filter_config['condition'])
    print(f"✓ Kept {len(filtered_df):,} positions")
    
    # Sort by chromosome (in genomic order) and position
    print(f"Sorting by chromosome and position...")
    sorted_df = filtered_df.with_columns([
        pl.col(chrom_col).replace_strict(chrom_order_map, default=99).alias('_chrom_order')
    ]).sort(['_chrom_order', pos_col]).drop('_chrom_order')
    
    # Add compressed coordinate index
    print(f"Generating compressed coordinates...")
    sorted_df = sorted_df.with_row_index('filtered_idx')
    
    # Get columns to keep
    columns_to_keep, col_counts = get_columns_to_keep(
        sorted_df, chrom_col, pos_col, gene_col
    )
    
    print(f"\nColumn breakdown:")
    print(f"  - Core (idx, chrom, pos, gene): 4")
    print(f"  - RGC raw metrics: {col_counts['rgc']}")
    print(f"  - ClinVar: {col_counts['clinvar']}")
    print(f"  - Training labels: {col_counts['training']}")
    print(f"  - dbNSFP scores: {col_counts['dbnsfp']}")
    print(f"  - gnomAD coverage: {col_counts['gnomad_coverage']}")
    print(f"  - Variant frequencies: {col_counts['variant_freq']}")
    print(f"  - phyloP: {col_counts['phylop']}")
    print(f"  - gnomAD constraint: {col_counts['gnomad_constraint']}")
    print(f"  - Constraint predictions: {col_counts['constraint_preds']}")
    print(f"  - Stacked scores: {col_counts['stacked_scores']}")
    print(f"  - Variant consequences: {col_counts['variant_consequences']}")
    print(f"  - Percentiles: {col_counts['percentiles']}")
    print(f"  - Domains: {col_counts['domains']}")
    print(f"  Total columns: {len(columns_to_keep)}")
    
    # Select and rename columns
    axis_table = sorted_df.select(columns_to_keep)
    
    rename_map = {
        chrom_col: 'chrom',
        pos_col: 'pos'
    }
    if gene_col:
        rename_map[gene_col] = 'gene_symbol'
    
    axis_table = axis_table.rename(rename_map)
    
    # Save axis table
    output_file = output_dir / f'{filter_id}.parquet'
    print(f"\nSaving axis table...")
    axis_table.write_parquet(output_file)
    print(f"✓ Saved: {output_file}")
    print(f"  Size: {output_file.stat().st_size / 1024 / 1024:.2f} MB")
    
    # Generate gene index
    print(f"\nGenerating gene index...")
    if 'gene_symbol' in axis_table.columns:
        gene_index = axis_table.group_by('gene_symbol').agg([
            pl.col('chrom').first(),
            pl.col('pos').min().alias('pos_start'),
            pl.col('pos').max().alias('pos_end'),
            pl.col('filtered_idx').min().alias('filtered_idx_start'),
            pl.col('filtered_idx').max().alias('filtered_idx_end'),
            pl.len().alias('num_positions')
        ]).sort('filtered_idx_start')
        
        gene_output = output_dir / f'gene_index_{filter_id}.parquet'
        gene_index.write_parquet(gene_output)
        print(f"✓ Saved: {gene_output}")
        print(f"  Genes: {len(gene_index):,}")
    else:
        print("  ⚠ No gene column found, skipping gene index")
        gene_index = None
    
    # Per-chromosome stats
    chrom_stats = axis_table.group_by('chrom').agg([
        pl.len().alias('count'),
        pl.col('filtered_idx').min().alias('start_idx'),
        pl.col('filtered_idx').max().alias('end_idx')
    ])
    chrom_stats = chrom_stats.with_columns([
        pl.col('chrom').replace_strict(chrom_order_map, default=99).alias('_order')
    ]).sort('_order').drop('_order')
    
    print(f"\nPer-chromosome breakdown:")
    for row in chrom_stats.iter_rows(named=True):
        print(f"  {row['chrom']:6s}: {row['count']:>10,} positions")
    
    return {
        'filter_id': filter_id,
        'positions': len(axis_table),
        'columns': len(axis_table.columns),
        'genes': len(gene_index) if gene_index is not None else 0
    }


def main(input_path: str = None, output_dir: str = None):
    print(f"{'='*60}")
    print(f"VARPRED BROWSER - DATA PREPROCESSING")
    print(f"Generating axis tables for all filter modes")
    print(f"{'='*60}\n")

    # Use provided paths or defaults
    input_parquet = Path(input_path) if input_path else Path(DEFAULT_INPUT)
    output_path = Path(output_dir) if output_dir else DEFAULT_OUTPUT

    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)

    # Load parquet file
    print(f"Loading data from {input_parquet}...")
    try:
        df = pl.read_parquet(input_parquet)
        print(f"Loaded {len(df):,} total positions")
        print(f"  Columns: {len(df.columns)}")
    except Exception as e:
        print(f"Error loading parquet: {e}")
        sys.exit(1)
    
    # Detect column names
    chrom_col, pos_col, gene_col = detect_columns(df)
    
    print(f"\nDetected columns:")
    print(f"  Chromosome: {chrom_col}")
    print(f"  Position: {pos_col}")
    print(f"  Gene: {gene_col or 'NOT FOUND'}")
    
    # Show chromosome distribution
    print(f"\nChromosome distribution (top 10):")
    chrom_counts = df.group_by(chrom_col).agg(
        pl.len().alias('count')
    ).sort('count', descending=True)
    for row in chrom_counts.head(10).iter_rows(named=True):
        print(f"  {row[chrom_col]}: {row['count']:,}")
    
    # Create chromosome order mapping
    chrom_order_map = {chrom: i for i, chrom in enumerate(CHROM_ORDER)}
    
    # Process each filter
    results = []
    for filter_id, filter_config in FILTERS.items():
        result = process_filter(
            df, filter_id, filter_config,
            chrom_col, pos_col, gene_col, chrom_order_map, output_path
        )
        results.append(result)
    
    # Final summary
    print(f"\n{'='*60}")
    print(f"FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"\nGenerated {len(results)} axis tables:\n")
    
    for result in results:
        print(f"  {result['filter_id']}:")
        print(f"    - Positions: {result['positions']:,}")
        print(f"    - Columns: {result['columns']}")
        print(f"    - Genes: {result['genes']:,}")
        print()
    
    print(f"Output directory: {output_path.absolute()}")
    print(f"{'='*60}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Preprocess browser data into axis tables')
    parser.add_argument('--input', '-i', default=None,
                        help='Input parquet file path (default: from config or ../rgc_browser_data_merged.parquet)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output directory (default: from config or data/)')

    args = parser.parse_args()
    main(input_path=args.input, output_dir=args.output)
