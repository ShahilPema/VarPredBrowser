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

# Stacked array columns - we'll extract max prediction per locus and add per-variant percentiles
# Format: column_name -> (score_field, extra_fields_to_preserve)
# score_field: the field containing the score to percentile ('pred' or 'score')
# extra_fields: other fields to preserve when rebuilding the struct
STACKED_PRED_COLUMNS = {
    'Constraint': ('pred', ['n_pred']),
    'Core': ('pred', ['n_pred']),
    'Complete': ('pred', ['n_pred']),
    'AlphaMissense_stacked': ('score', []),
    'ESM1b_stacked': ('score', []),
}


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


def build_max_pred_exprs(df: pl.DataFrame, columns: dict[str, tuple]) -> list[pl.Expr]:
    """
    Build expressions to extract max score from stacked array columns.

    Args:
        df: Input DataFrame
        columns: Dict of column_name -> (score_field, extra_fields)

    Array format: [{alt: str, <score_field>: float, ...}, ...]
    """
    exprs = []
    for col, (score_field, _) in columns.items():
        if col in df.columns:
            exprs.append(
                pl.col(col).list.eval(pl.element().struct.field(score_field)).list.max()
                .alias(f'{col}_max_pred')
            )
    return exprs


def add_per_variant_percentiles(
    df: pl.DataFrame,
    stacked_col: str,
    score_field: str,
    extra_fields: list[str]
) -> pl.DataFrame:
    """
    Add percentile field to each variant in a stacked array column.

    Args:
        df: Input DataFrame
        stacked_col: Name of the stacked array column
        score_field: Name of the score field to percentile ('pred' or 'score')
        extra_fields: List of additional fields to preserve (e.g., ['n_pred'])

    Input format: [{alt: str, <score_field>: float, ...}, ...]
    Output format: [{alt: str, <score_field>: float, ..., percentile: float}, ...]

    Percentiles are calculated across ALL variants at ALL positions.
    """
    if stacked_col not in df.columns:
        return df

    print(f"    Processing {stacked_col} (score_field={score_field})...")

    # Add row index for grouping back
    df = df.with_row_index('_row_idx')

    # Explode the array to get individual variants
    # Note: explode() must be called on DataFrame, not inside select()
    exploded = df.select(
        '_row_idx',
        pl.col(stacked_col).alias('_variant')
    ).explode('_variant').filter(
        pl.col('_variant').is_not_null()
    )

    if exploded.height == 0:
        print(f"      No variants found in {stacked_col}")
        df = df.drop('_row_idx')
        return df

    # Extract score values - always extract 'alt' and the score field
    extract_exprs = [
        pl.col('_variant').struct.field('alt').alias('_alt'),
        pl.col('_variant').struct.field(score_field).alias('_score'),
    ]
    # Also extract any extra fields
    for field in extra_fields:
        extract_exprs.append(
            pl.col('_variant').struct.field(field).alias(f'_{field}')
        )
    exploded = exploded.with_columns(extract_exprs)

    # Count non-null scores
    total_variants = exploded.filter(pl.col('_score').is_not_null()).height
    print(f"      Total variants: {total_variants:,}")

    if total_variants > 0:
        # Calculate percentile rank
        exploded = exploded.with_columns(
            (pl.col('_score').rank(method='average') / total_variants * 100)
            .alias('_percentile')
        )

        # Rebuild struct with percentile field
        # Build the struct fields dynamically
        struct_fields = [
            pl.col('_alt').alias('alt'),
            pl.col('_score').alias(score_field),
        ]
        for field in extra_fields:
            struct_fields.append(pl.col(f'_{field}').alias(field))
        struct_fields.append(pl.col('_percentile').alias('percentile'))

        exploded = exploded.with_columns(
            pl.struct(struct_fields).alias('_variant_with_perc')
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

# MTR column used as gating filter for cross-normalization
# Sites without MTR are excluded from ALL cross-normalized percentiles
MTR_GATING_COLUMN = 'max_RGC_MTR_MTR'


def add_cross_norm_percentiles(df: pl.DataFrame) -> pl.DataFrame:
    """
    Add cross-normalized percentile to variants where ALL 5 scores are defined
    AND MTR is available at that position.

    Cross-normalization means: for each variant (position + alt allele),
    calculate percentile only among variants that have all 5 scores.

    MTR gating: Sites without MTR are excluded from cross-normalization entirely.
    This ensures consistency - if a site isn't in MTR, it won't have cross-norm
    percentiles for any of the other metrics either.

    This modifies the stacked arrays to add a 'cross_norm_perc' field.
    """
    print("\n=== Calculating Per-Variant Cross-Normalized Percentiles ===")

    # Check MTR gating column
    if MTR_GATING_COLUMN not in df.columns:
        print(f"  WARNING: MTR gating column '{MTR_GATING_COLUMN}' not found, skipping cross-norm")
        return df

    mtr_defined_count = df.filter(pl.col(MTR_GATING_COLUMN).is_not_null()).height
    print(f"  MTR gating: {mtr_defined_count:,} / {len(df):,} positions have MTR defined")

    # Check which stacked columns exist
    available_cols = {k: v for k, v in CROSS_NORM_STACKED_COLUMNS.items() if k in df.columns}
    print(f"  Available stacked columns: {list(available_cols.keys())}")

    if len(available_cols) < 2:
        print("  Not enough stacked columns for cross-normalization, skipping")
        return df

    # Add row index for tracking
    df = df.with_row_index('_row_idx')

    # Get row indices where MTR is defined (for filtering)
    mtr_valid_rows = df.filter(
        pl.col(MTR_GATING_COLUMN).is_not_null()
    ).select('_row_idx')

    # Explode each stacked column and extract scores (only for MTR-valid positions)
    exploded_dfs = {}
    for col, (score_field, name) in available_cols.items():
        print(f"  Exploding {col} (MTR-gated)...")
        exploded = df.select(
            '_row_idx',
            pl.col(col).alias('_variant')
        ).explode('_variant').filter(
            pl.col('_variant').is_not_null()
        ).join(
            mtr_valid_rows, on='_row_idx', how='inner'  # Only keep MTR-valid positions
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

    print(f"  Variants with all {len(names)} scores (MTR-gated): {joined.height:,}")

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
            pl.col(col).alias('_variant')
        ).explode('_variant')

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
        # Use STACKED_PRED_COLUMNS to get the correct field structure
        pred_score_field, extra_fields = STACKED_PRED_COLUMNS.get(col, (score_field, []))

        struct_fields = [
            pl.col('_variant').struct.field('alt').alias('alt'),
            pl.col('_variant').struct.field(pred_score_field).alias(pred_score_field),
        ]
        for field in extra_fields:
            struct_fields.append(pl.col('_variant').struct.field(field).alias(field))
        struct_fields.append(pl.col('_variant').struct.field('percentile').alias('percentile'))
        struct_fields.append(pl.col(perc_col).alias('cross_norm_perc'))

        exploded = exploded.with_columns(
            pl.struct(struct_fields).alias('_new_variant')
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

    # Calculate MTR cross-norm percentile (position-level)
    # MTR percentile is calculated only among positions with all 5 predictor scores
    print(f"  Calculating MTR cross-norm percentile...")
    cross_norm_positions = joined.select('_row_idx').unique()
    cross_norm_position_count = cross_norm_positions.height
    print(f"    Cross-norm positions for MTR: {cross_norm_position_count:,}")

    # Calculate MTR percentile among these positions
    mtr_at_cross_norm = df.join(
        cross_norm_positions, on='_row_idx', how='inner'
    ).select('_row_idx', MTR_GATING_COLUMN)

    mtr_with_perc = mtr_at_cross_norm.with_columns(
        (pl.col(MTR_GATING_COLUMN).rank(method='average') / cross_norm_position_count * 100)
        .alias('_mtr_cross_norm_perc')
    ).select('_row_idx', '_mtr_cross_norm_perc')

    # Join back to main dataframe
    df = df.join(mtr_with_perc, on='_row_idx', how='left')
    df = df.rename({'_mtr_cross_norm_perc': 'RGC_MTR_MTR_cross_norm_perc'})

    df = df.drop('_row_idx')
    print(f"  Cross-normalization complete (including MTR)")
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
    for stacked_col, (score_field, extra_fields) in STACKED_PRED_COLUMNS.items():
        df = add_per_variant_percentiles(df, stacked_col, score_field, extra_fields)

    # --- PER-VARIANT CROSS-NORMALIZED PERCENTILES ---
    # This adds cross_norm_perc to variants where:
    # 1. MTR is defined at the position (gating filter)
    # 2. ALL 5 stacked predictor scores exist (Constraint, Core, Complete, AlphaMissense, ESM1b)
    # Also calculates MTR cross-norm percentile among positions meeting these criteria
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
