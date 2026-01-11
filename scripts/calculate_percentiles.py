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
    # Constraint predictions (need to extract from stacked arrays - use max per locus)
    # These are handled separately since they're arrays
}

# Stacked array columns - we'll extract max prediction per locus
STACKED_PRED_COLUMNS = ['Constraint', 'Core', 'Complete']


def calculate_exome_percentile(df: pl.DataFrame, col: str, output_name: str) -> pl.DataFrame:
    """
    Calculate exome-wide percentile for a column using average rank method.

    Percentile = (rank / count_non_null) * 100
    """
    non_null_count = df.filter(pl.col(col).is_not_null()).height

    if non_null_count == 0:
        return df.with_columns(pl.lit(None).alias(f'{output_name}_exome_perc'))

    # Use rank with average method for ties
    return df.with_columns(
        (pl.col(col).rank(method='average') / non_null_count * 100)
        .alias(f'{output_name}_exome_perc')
    )


def extract_max_pred(df: pl.DataFrame, col: str) -> pl.DataFrame:
    """
    Extract max prediction value from stacked array column.

    Array format: [{alt: str, pred: float, n_pred: int}, ...]
    """
    # Check if column exists and has data
    if col not in df.columns:
        return df.with_columns(pl.lit(None).cast(pl.Float64).alias(f'{col}_max_pred'))

    # Extract max pred value from array of structs
    # Struct fields are named: alt, pred, n_pred (from hl.struct() in Hail)
    return df.with_columns(
        pl.col(col).list.eval(pl.element().struct.field('pred')).list.max()
        .alias(f'{col}_max_pred')
    )


def main(input_path: str, output_path: str):
    print(f"Loading data from {input_path}...")
    df = pl.read_parquet(input_path)
    print(f"Loaded {len(df):,} rows, {len(df.columns)} columns")

    # --- EXOME-WIDE PERCENTILES ---
    print("\n=== Calculating Exome-Wide Percentiles ===")

    for col, name in SCORE_COLUMNS.items():
        if col in df.columns:
            print(f"  {col} -> {name}_exome_perc")
            df = calculate_exome_percentile(df, col, name)
        else:
            print(f"  {col} - NOT FOUND, skipping")

    # Extract max predictions from stacked arrays
    print("\n=== Extracting Max Predictions from Stacked Arrays ===")
    for col in STACKED_PRED_COLUMNS:
        if col in df.columns:
            print(f"  {col} -> {col}_max_pred")
            df = extract_max_pred(df, col)
        else:
            print(f"  {col} - NOT FOUND, skipping")

    # Calculate exome percentiles for extracted predictions
    for col in STACKED_PRED_COLUMNS:
        max_col = f'{col}_max_pred'
        if max_col in df.columns:
            print(f"  {max_col} -> {col}_exome_perc")
            df = calculate_exome_percentile(df, max_col, col)

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

        # Calculate cross-norm percentiles
        for col, name in available_scores:
            perc_col = f'{name}_cross_norm_perc'
            print(f"  {col} -> {perc_col}")

            # Rank only among cross-norm positions, null otherwise
            df = df.with_columns(
                pl.when(cross_norm_filter)
                .then(
                    pl.col(col).rank(method='average').over(cross_norm_filter)
                    / cross_norm_count * 100
                )
                .otherwise(None)
                .alias(perc_col)
            )
    else:
        print("  WARNING: Not all 6 scores available, skipping cross-norm percentiles")
        missing = [name for col, name in cross_norm_scores if col not in df.columns]
        print(f"  Missing: {missing}")

    # --- SAVE OUTPUT ---
    print(f"\n=== Saving Output ===")
    print(f"  Output: {output_path}")
    print(f"  Columns: {len(df.columns)}")

    # Show new percentile columns
    perc_cols = [c for c in df.columns if '_perc' in c]
    print(f"  New percentile columns: {len(perc_cols)}")
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
