#!/usr/bin/env python3
"""
Calculate percentiles for browser scores using Polars.

Uses rank(method='average') for proper tie handling.

Two types of percentiles:
1. Exome-wide: Percentile across all non-null values for each score
2. Cross-normalized: Percentile only where ALL 6 scores are defined

Usage:
    python scripts/calculate_percentiles.py --input merged.parquet --output merged_with_perc.parquet
"""

import argparse
import polars as pl
from pathlib import Path


# Score columns to percentile
SCORE_COLUMNS = {
    # dbNSFP scores (flattened column names)
    'max_AlphaMissense_am_pathogenicity': 'AlphaMissense_am_pathogenicity',
    'max_ESM1b_score': 'ESM1b_score',
    'max_RGC_MTR_MTR': 'RGC_MTR_MTR',
}

# Stacked array columns - we'll extract max prediction per locus
STACKED_PRED_COLUMNS = ['Constraint', 'Core', 'Complete']


def get_oe_columns() -> dict:
    """Generate O/E column names for percentile calculation."""
    windows = ["3bp", "9bp", "21bp", "45bp", "93bp"]
    consequences = ["mis", "syn", "any"]
    af_suffix = "af0epos00"  # Only AF≥0 for O/E

    columns = {}
    for window in windows:
        for cons in consequences:
            col = f"rgc_{cons}_exomes_XX_XY_{window}_oe_{af_suffix}"
            columns[col] = col  # Use same name for output
    return columns


def get_vir_columns() -> dict:
    """Generate VIR column names for percentile calculation (no depth)."""
    consequences = ["mis", "syn", "any"]
    af_suffixes = ["af0epos00", "af1eneg04", "af1eneg06"]
    # Length, Expected μ, Mean Expected - NOT depth
    metrics = ["vir_length", "vir_mu_exp", "mean_vir_exp"]

    columns = {}
    for metric in metrics:
        for cons in consequences:
            for af in af_suffixes:
                col = f"rgc_{cons}_exomes_XX_XY_{metric}_{af}"
                columns[col] = col

    # AA Level (missense only) - only length, not depth
    for af in af_suffixes:
        col = f"rgc_mis_exomes_XX_XY_aa_vir_length_{af}"
        columns[col] = col

    return columns


def build_percentile_exprs(df: pl.DataFrame, columns: dict[str, str]) -> list[pl.Expr]:
    """
    Build percentile expressions for all columns at once.

    Returns list of expressions that can be passed to with_columns().
    """
    exprs = []
    for col, name in columns.items():
        if col in df.columns:
            # Count non-null values for this column
            non_null_count = df.select(pl.col(col).is_not_null().sum()).item()
            if non_null_count > 0:
                exprs.append(
                    (pl.col(col).rank(method='average') / non_null_count * 100)
                    .alias(f'{name}_exome_perc')
                )
    return exprs


def build_max_pred_exprs(df: pl.DataFrame, columns: list[str]) -> list[pl.Expr]:
    """
    Build expressions to extract max prediction from stacked array columns.

    Array format: [{alt: str, pred: float, n_pred: int}, ...]
    """
    exprs = []
    for col in columns:
        if col in df.columns:
            exprs.append(
                pl.col(col).list.eval(pl.element().struct.field('pred')).list.max()
                .alias(f'{col}_max_pred')
            )
    return exprs


def add_per_variant_percentiles(df: pl.DataFrame, stacked_col: str) -> pl.DataFrame:
    """
    Add percentile field to each variant in a stacked array column.

    Input format: [{alt: str, pred: float, n_pred: int}, ...]
    Output format: [{alt: str, pred: float, n_pred: int, percentile: float}, ...]

    Percentiles are calculated across ALL variants at ALL positions.
    """
    if stacked_col not in df.columns:
        return df

    print(f"    Processing {stacked_col}...")

    # Add row index for grouping back
    df = df.with_row_index('_row_idx')

    # Explode the array to get individual variants
    exploded = df.select(
        '_row_idx',
        pl.col(stacked_col).explode().alias('_variant')
    ).filter(
        pl.col('_variant').is_not_null()
    )

    if exploded.height == 0:
        print(f"      No variants found in {stacked_col}")
        df = df.drop('_row_idx')
        return df

    # Extract pred values and calculate percentiles
    exploded = exploded.with_columns([
        pl.col('_variant').struct.field('alt').alias('_alt'),
        pl.col('_variant').struct.field('pred').alias('_pred'),
        pl.col('_variant').struct.field('n_pred').alias('_n_pred'),
    ])

    # Count non-null predictions
    total_variants = exploded.filter(pl.col('_pred').is_not_null()).height
    print(f"      Total variants: {total_variants:,}")

    if total_variants > 0:
        # Calculate percentile rank
        exploded = exploded.with_columns(
            (pl.col('_pred').rank(method='average') / total_variants * 100)
            .alias('_percentile')
        )

        # Rebuild struct with percentile field
        exploded = exploded.with_columns(
            pl.struct([
                pl.col('_alt').alias('alt'),
                pl.col('_pred').alias('pred'),
                pl.col('_n_pred').alias('n_pred'),
                pl.col('_percentile').alias('percentile'),
            ]).alias('_variant_with_perc')
        )

        # Group back by row index to recreate arrays
        grouped = exploded.group_by('_row_idx').agg(
            pl.col('_variant_with_perc').alias(f'{stacked_col}_with_perc')
        )

        # Join back to original dataframe
        df = df.join(grouped, on='_row_idx', how='left')

        # Replace original column with new one containing percentiles
        df = df.drop(stacked_col)
        df = df.rename({f'{stacked_col}_with_perc': stacked_col})

    df = df.drop('_row_idx')
    return df


# Stacked columns for cross-normalization (must have per-variant scores)
CROSS_NORM_STACKED_COLUMNS = {
    # (column_name, score_field_name, output_field_name)
    'Constraint': ('pred', 'Constraint'),
    'Core': ('pred', 'Core'),
    'Complete': ('pred', 'Complete'),
    'AlphaMissense_stacked': ('score', 'AlphaMissense'),
    'ESM1b_stacked': ('score', 'ESM1b'),
}


def add_cross_norm_percentiles(df: pl.DataFrame) -> pl.DataFrame:
    """
    Add cross-normalized percentile to variants where ALL 5 scores are defined.

    Cross-normalization means: for each variant (position + alt allele),
    calculate percentile only among variants that have all 5 scores.

    This modifies the stacked arrays to add a 'cross_norm_perc' field.
    """
    print("\n=== Calculating Per-Variant Cross-Normalized Percentiles ===")

    # Check which stacked columns exist
    available_cols = {k: v for k, v in CROSS_NORM_STACKED_COLUMNS.items() if k in df.columns}
    print(f"  Available stacked columns: {list(available_cols.keys())}")

    if len(available_cols) < 2:
        print("  Not enough stacked columns for cross-normalization, skipping")
        return df

    # Add row index for tracking
    df = df.with_row_index('_row_idx')

    # Explode each stacked column and extract scores
    exploded_dfs = {}
    for col, (score_field, name) in available_cols.items():
        print(f"  Exploding {col}...")
        exploded = df.select(
            '_row_idx',
            pl.col(col).explode().alias('_variant')
        ).filter(
            pl.col('_variant').is_not_null()
        )

        if exploded.height == 0:
            print(f"    No variants in {col}")
            continue

        # Extract alt and score
        exploded = exploded.with_columns([
            pl.col('_variant').struct.field('alt').alias('alt'),
            pl.col('_variant').struct.field(score_field).alias(f'{name}_score'),
        ])

        exploded_dfs[name] = exploded.select('_row_idx', 'alt', f'{name}_score')

    if len(exploded_dfs) < 2:
        print("  Not enough valid stacked columns, skipping cross-norm")
        df = df.drop('_row_idx')
        return df

    # Join all exploded dataframes on (row_idx, alt)
    print(f"  Joining {len(exploded_dfs)} score tables by (row_idx, alt)...")
    names = list(exploded_dfs.keys())
    joined = exploded_dfs[names[0]]
    for name in names[1:]:
        joined = joined.join(
            exploded_dfs[name],
            on=['_row_idx', 'alt'],
            how='inner'  # Only keep variants with ALL scores
        )

    print(f"  Variants with all {len(names)} scores: {joined.height:,}")

    if joined.height == 0:
        print("  No variants have all scores defined, skipping cross-norm")
        df = df.drop('_row_idx')
        return df

    # Calculate cross-norm percentile for each score
    cross_norm_count = joined.height
    for name in names:
        score_col = f'{name}_score'
        perc_col = f'{name}_cross_norm_perc'
        joined = joined.with_columns(
            (pl.col(score_col).rank(method='average') / cross_norm_count * 100)
            .alias(perc_col)
        )
    print(f"  Calculated cross-norm percentiles for {len(names)} scores")

    # Now we need to add cross_norm_perc back to each stacked array
    # Create a lookup table: (row_idx, alt) -> cross_norm_perc for each score
    cross_norm_lookup = joined.select(
        '_row_idx', 'alt',
        *[f'{name}_cross_norm_perc' for name in names]
    )

    # For each stacked column, update the structs to include cross_norm_perc
    for col, (score_field, name) in available_cols.items():
        perc_col = f'{name}_cross_norm_perc'
        if perc_col not in cross_norm_lookup.columns:
            continue

        print(f"  Adding cross_norm_perc to {col}...")

        # Explode the array again to add cross_norm_perc
        exploded = df.select(
            '_row_idx',
            pl.col(col).explode().alias('_variant')
        )

        # Add sequence number within each group to maintain order
        exploded = exploded.with_columns(
            pl.col('_variant').struct.field('alt').alias('_alt')
        )

        # Join with cross_norm_lookup
        exploded = exploded.join(
            cross_norm_lookup.select('_row_idx', 'alt', perc_col),
            left_on=['_row_idx', '_alt'],
            right_on=['_row_idx', 'alt'],
            how='left'
        )

        # Rebuild struct with cross_norm_perc field
        # Get existing fields from the struct
        if col in ['Constraint', 'Core', 'Complete']:
            # Prediction format: {alt, pred, n_pred, percentile}
            exploded = exploded.with_columns(
                pl.struct([
                    pl.col('_variant').struct.field('alt').alias('alt'),
                    pl.col('_variant').struct.field('pred').alias('pred'),
                    pl.col('_variant').struct.field('n_pred').alias('n_pred'),
                    pl.col('_variant').struct.field('percentile').alias('percentile'),
                    pl.col(perc_col).alias('cross_norm_perc'),
                ]).alias('_new_variant')
            )
        else:
            # dbNSFP format: {alt, score, percentile}
            exploded = exploded.with_columns(
                pl.struct([
                    pl.col('_variant').struct.field('alt').alias('alt'),
                    pl.col('_variant').struct.field('score').alias('score'),
                    pl.col('_variant').struct.field('percentile').alias('percentile'),
                    pl.col(perc_col).alias('cross_norm_perc'),
                ]).alias('_new_variant')
            )

        # Filter out nulls and group back
        grouped = exploded.filter(
            pl.col('_new_variant').is_not_null()
        ).group_by('_row_idx').agg(
            pl.col('_new_variant').alias(f'{col}_with_xnorm')
        )

        # Join back and replace column
        df = df.join(grouped, on='_row_idx', how='left')
        df = df.drop(col)
        df = df.rename({f'{col}_with_xnorm': col})

    df = df.drop('_row_idx')
    print(f"  Cross-normalization complete")
    return df


def main(input_path: str, output_path: str):
    print(f"Loading data from {input_path}...")
    df = pl.read_parquet(input_path)
    print(f"Loaded {len(df):,} rows, {len(df.columns)} columns")

    # Collect all column mappings
    all_columns = {}
    all_columns.update(SCORE_COLUMNS)
    all_columns.update(get_oe_columns())
    all_columns.update(get_vir_columns())

    # Report which columns are available
    available = [col for col in all_columns if col in df.columns]
    missing = [col for col in all_columns if col not in df.columns]
    print(f"\n=== Column Availability ===")
    print(f"  Available: {len(available)}/{len(all_columns)}")
    if missing:
        print(f"  Missing: {missing[:5]}{'...' if len(missing) > 5 else ''}")

    # --- EXTRACT MAX PREDICTIONS FROM STACKED ARRAYS ---
    print("\n=== Extracting Max Predictions from Stacked Arrays ===")
    available_stacked = [col for col in STACKED_PRED_COLUMNS if col in df.columns]
    print(f"  Stacked columns: {available_stacked}")

    max_pred_exprs = build_max_pred_exprs(df, STACKED_PRED_COLUMNS)
    if max_pred_exprs:
        df = df.with_columns(max_pred_exprs)
        print(f"  Extracted {len(max_pred_exprs)} max_pred columns in parallel")

    # Add max_pred columns to percentile calculation
    for col in STACKED_PRED_COLUMNS:
        max_col = f'{col}_max_pred'
        if max_col in df.columns:
            all_columns[max_col] = col

    # --- PER-VARIANT PERCENTILES FOR STACKED ARRAYS ---
    print("\n=== Calculating Per-Variant Percentiles for Stacked Arrays ===")
    for stacked_col in STACKED_PRED_COLUMNS:
        df = add_per_variant_percentiles(df, stacked_col)

    # --- PER-VARIANT CROSS-NORMALIZED PERCENTILES ---
    # This adds cross_norm_perc to variants that have ALL 5 stacked scores
    df = add_cross_norm_percentiles(df)

    # --- EXOME-WIDE PERCENTILES (ALL AT ONCE) ---
    print("\n=== Calculating Exome-Wide Percentiles (Parallel) ===")
    percentile_exprs = build_percentile_exprs(df, all_columns)
    print(f"  Building {len(percentile_exprs)} percentile expressions...")

    if percentile_exprs:
        df = df.with_columns(percentile_exprs)
        print(f"  Calculated {len(percentile_exprs)} percentiles in parallel")

    # --- CROSS-NORMALIZED PERCENTILES ---
    print("\n=== Calculating Cross-Normalized Percentiles ===")

    # Define the 6 scores for cross-normalization
    cross_norm_scores = [
        ('max_AlphaMissense_am_pathogenicity', 'AlphaMissense_am_pathogenicity'),
        ('max_ESM1b_score', 'ESM1b_score'),
        ('max_RGC_MTR_MTR', 'RGC_MTR_MTR'),
        ('Constraint_max_pred', 'Constraint'),
        ('Core_max_pred', 'Core'),
        ('Complete_max_pred', 'Complete'),
    ]

    # Check which columns exist
    available_scores = [(col, name) for col, name in cross_norm_scores if col in df.columns]
    print(f"  Available scores for cross-norm: {len(available_scores)}/6")

    if len(available_scores) == 6:
        # Create filter for rows where ALL scores are defined
        cross_norm_filter = pl.lit(True)
        for col, _ in available_scores:
            cross_norm_filter = cross_norm_filter & pl.col(col).is_not_null()

        # Count cross-norm positions
        cross_norm_count = df.filter(cross_norm_filter).height
        print(f"  Cross-norm positions: {cross_norm_count:,}")

        # Build all cross-norm percentile expressions at once
        cross_norm_exprs = []
        for col, name in available_scores:
            perc_col = f'{name}_cross_norm_perc'
            cross_norm_exprs.append(
                pl.when(cross_norm_filter)
                .then(
                    pl.col(col).rank(method='average').over(cross_norm_filter)
                    / cross_norm_count * 100
                )
                .otherwise(None)
                .alias(perc_col)
            )

        df = df.with_columns(cross_norm_exprs)
        print(f"  Calculated {len(cross_norm_exprs)} cross-norm percentiles in parallel")
    else:
        print("  WARNING: Not all 6 scores available, skipping cross-norm percentiles")
        missing = [name for col, name in cross_norm_scores if col not in df.columns]
        print(f"  Missing: {missing}")

    # --- SAVE OUTPUT ---
    print(f"\n=== Saving Output ===")
    print(f"  Output: {output_path}")
    print(f"  Columns: {len(df.columns)}")

    # Show new percentile columns
    perc_cols = [c for c in df.columns if '_perc' in c or '_max_pred' in c]
    print(f"  New computed columns: {len(perc_cols)}")
    for c in sorted(perc_cols):
        print(f"    - {c}")

    df.write_parquet(output_path)
    print(f"\nDone! Output size: {Path(output_path).stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate percentiles for browser scores')
    parser.add_argument('--input', '-i', required=True, help='Input parquet file')
    parser.add_argument('--output', '-o', required=True, help='Output parquet file')

    args = parser.parse_args()
    main(args.input, args.output)
