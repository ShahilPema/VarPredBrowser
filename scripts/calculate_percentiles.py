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
