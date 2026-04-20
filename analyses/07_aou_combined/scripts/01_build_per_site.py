#!/usr/bin/env python3
"""Workbench-side: build per_site_{gnomad_only,aou_only,combined}.parquet.

Pure Polars (no Hail / no dataproc). Reads:
  - base_table.parquet (shipped from BCM) — gnomAD AC/PASS, scores, roulette
    MR, ENST/HGNC, gnomAD-derived f_lowcov / f_ab_low.
  - AoU observed_snps parquet at $AOU_OBSERVED_PARQUET — locus.contig,
    locus.position, alleles[]. Row presence = variant observed AND passing
    in All of Us (no counts / no pass flag needed).

Joins AoU observed flag onto the base table by reconstructing (locus_str,
alleles_str) keys, then for each dataset writes per-site
(locus_str, alleles_str, transcript, gene_symbol, observed, expected) where:
  - **gnomAD QC filters apply to ALL three datasets** (a position dropped
    by f_lowcov or f_ab_low is dropped from every leg).
  - observed:
      gnomad_only : 1 iff gnomad_ac > 0
      aou_only    : 1 iff aou_observed
      combined    : 1 iff either above
  - expected: 1 - exp(-k * roulette_AR_MR), with k fit per-(dataset, region)
    on synonymous scaling sites.

Also writes per_site/synonymous_augmented.parquet for `02_diagnostic_plots.py`.

Runtime expectation on a 96-vCPU workbench Jupyter VM: ~5 min total.
"""
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import (
    BASE_TABLE_PARQUET, PER_SITE_DIR, INPUTS_DIR,
    DATASETS, AOU_CDR_VERSION, per_site_path,
)


def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)


# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
# AoU observed_snps parquet (workspace bucket). Row present => variant observed
# AND passing in AoU; no AC / pass columns needed. Schema:
#   locus.contig (str), locus.position (int32), alleles (list[str])
AOU_OBSERVED_PARQUET = os.environ.get(
    'AOU_OBSERVED_PARQUET',
    'gs://fc-secure-57a37491-77aa-498c-b0f0-2163b5feb12e/data/observed_snps.parquet/*.parquet',
)

# Local cached AoU+gnomAD-augmented base table; built once and reused.
AUGMENTED_BASE = INPUTS_DIR / 'base_with_aou.parquet'
SYN_OUT = PER_SITE_DIR / 'synonymous_augmented.parquet'


# -----------------------------------------------------------------------------
# Stage 1: join AoU VAT onto base table
# -----------------------------------------------------------------------------

def build_augmented_base():
    if AUGMENTED_BASE.exists():
        log(f'AUGMENTED base already at {AUGMENTED_BASE}; skipping join')
        return

    log(f'Loading AoU observed_snps: {AOU_OBSERVED_PARQUET}')
    aou = (pl.scan_parquet(AOU_OBSERVED_PARQUET)
           .with_columns([
               (pl.col('locus.contig') + pl.lit(':') +
                pl.col('locus.position').cast(pl.Utf8)).alias('locus_str'),
               (pl.lit('["') + pl.col('alleles').list.get(0) + pl.lit('","') +
                pl.col('alleles').list.get(1) + pl.lit('"]')).alias('alleles_str'),
           ])
           .select(['locus_str', 'alleles_str'])
           # Row presence => observed AND passing in AoU.
           .with_columns(pl.lit(True).alias('aou_observed')))

    log(f'Loading base table: {BASE_TABLE_PARQUET}')
    base_path = (f'{BASE_TABLE_PARQUET}/*.parquet'
                 if BASE_TABLE_PARQUET.is_dir() else str(BASE_TABLE_PARQUET))
    base = pl.scan_parquet(base_path)

    log('Joining base + AoU observed (left, on locus_str + alleles_str)...')
    augmented = (base.join(aou, on=['locus_str', 'alleles_str'], how='left')
                     .with_columns(pl.col('aou_observed').fill_null(False)))

    log(f'Sinking augmented base to {AUGMENTED_BASE}')
    augmented.sink_parquet(str(AUGMENTED_BASE), compression='zstd', compression_level=3)

    n_rows = pl.scan_parquet(str(AUGMENTED_BASE)).select(pl.len()).collect()['len'][0]
    log(f'  augmented rows: {n_rows:,}')


# -----------------------------------------------------------------------------
# Stage 2: synonymous slice for diagnostics
# -----------------------------------------------------------------------------

def write_synonymous_slice():
    if SYN_OUT.exists():
        log(f'Synonymous slice already at {SYN_OUT}; skipping')
        return
    log('Building synonymous_augmented.parquet for diagnostics...')
    syn_cols = [
        'locus_str', 'transcript', 'roulette_AR_MR',
        'is_chrX_nonPAR', 'roulette_syn_scaling_site',
        'gnomad_in_table', 'gnomad_ac', 'gnomad_is_pass',
        'aou_observed',
        'f_lowcov', 'f_ab_low',
    ]
    (pl.scan_parquet(str(AUGMENTED_BASE))
       .filter(pl.col('most_deleterious_consequence_cds') == 'synonymous_variant')
       .filter(pl.col('roulette_AR_MR') > 0)
       .select(syn_cols)
       .sink_parquet(str(SYN_OUT), compression='zstd', compression_level=3))
    n = pl.scan_parquet(str(SYN_OUT)).select(pl.len()).collect()['len'][0]
    log(f'  syn rows: {n:,}  -> {SYN_OUT}')


# -----------------------------------------------------------------------------
# Stage 3: per-dataset Poisson scaling factor on synonymous sites
# -----------------------------------------------------------------------------

def fit_k(syn_df: pl.DataFrame, observed_expr: pl.Expr, region: str) -> float:
    """Region: 'autosomes_par' or 'chrX_nonpar'."""
    if region == 'chrX_nonpar':
        sub = syn_df.filter(pl.col('is_chrX_nonPAR'))
    else:
        sub = syn_df.filter(~pl.col('is_chrX_nonPAR') &
                            pl.col('roulette_syn_scaling_site'))
    if sub.height == 0:
        return float('nan')
    sub = sub.with_columns(observed_expr.alias('_obs'))
    agg = (sub.group_by('roulette_AR_MR')
              .agg(pl.col('_obs').sum().alias('polymorphic'),
                   pl.len().alias('total')))
    rates = agg['roulette_AR_MR'].to_numpy().astype(np.float64)
    totals = agg['total'].to_numpy().astype(np.float64)
    observed_count = float(agg['polymorphic'].sum())

    def objective(x):
        expected = np.sum((1 - np.exp(-x[0] * rates)) * totals)
        return abs(expected - observed_count)

    res = minimize(objective, [1.0], method='Nelder-Mead',
                   options={'xatol': 1e-3, 'maxiter': 10000})
    return float(res.x[0])


# -----------------------------------------------------------------------------
# Stage 4: per-dataset extraction
# -----------------------------------------------------------------------------

# gnomAD's f_lowcov + f_ab_low apply to ALL three datasets (per user spec —
# if a position is filtered out of gnomAD it's filtered out of every leg).
COMMON_FILTER = ~pl.col('f_lowcov') & ~pl.col('f_ab_low')

OBSERVED_EXPR = {
    'gnomad_only': (pl.col('gnomad_ac') > 0).cast(pl.Int32),
    'aou_only':    pl.col('aou_observed').cast(pl.Int32),
    'combined':    ((pl.col('gnomad_ac') > 0) | pl.col('aou_observed')).cast(pl.Int32),
}


def extract_one(dataset: str, syn_df: pl.DataFrame):
    out_path = per_site_path(dataset)
    if out_path.exists():
        log(f'SKIP {dataset}: {out_path} exists')
        return

    log(f'=== {dataset} ===')
    obs_expr = OBSERVED_EXPR[dataset]

    log(f'  fitting k_auto / k_chrX on synonymous sites...')
    k_auto = fit_k(syn_df, obs_expr, 'autosomes_par')
    k_chrX = fit_k(syn_df, obs_expr, 'chrX_nonpar')
    log(f'  k_auto={k_auto:.6f}  k_chrX={k_chrX:.6f}')

    expected_expr = (
        pl.when(pl.col('is_chrX_nonPAR'))
          .then(1.0 - (-pl.col('roulette_AR_MR') * k_chrX).exp())
          .otherwise(1.0 - (-pl.col('roulette_AR_MR') * k_auto).exp())
    )

    out = (pl.scan_parquet(str(AUGMENTED_BASE))
           .filter(COMMON_FILTER)
           .filter(pl.col('most_deleterious_consequence_cds') == 'missense_variant')
           .filter(pl.col('roulette_AR_MR') > 0)
           .with_columns([
               obs_expr.alias('observed'),
               expected_expr.alias('expected'),
           ])
           .filter(pl.col('expected') > 0)
           .select(['locus_str', 'alleles_str', 'transcript', 'gene_symbol',
                    'observed', 'expected']))

    out.sink_parquet(str(out_path), compression='zstd', compression_level=3)
    n = pl.scan_parquet(str(out_path)).select(pl.len()).collect()['len'][0]
    log(f'  wrote {out_path}  ({n:,} rows)')


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--datasets', default=','.join(DATASETS),
                    help=f'comma-separated subset of {DATASETS}')
    args = ap.parse_args()

    requested = args.datasets.split(',')
    for d in requested:
        if d not in DATASETS:
            raise ValueError(f'unknown dataset {d}; expected {DATASETS}')

    log(f'AoU CDR pin: {AOU_CDR_VERSION}')

    build_augmented_base()
    write_synonymous_slice()

    log(f'Loading synonymous slice for k fitting: {SYN_OUT}')
    syn_df = pl.read_parquet(f'{SYN_OUT}')
    # apply common filter to syn slice as well so k is fit on the same site set
    syn_df = syn_df.filter((~pl.col('f_lowcov')) & (~pl.col('f_ab_low')))
    log(f'  syn rows after gnomAD QC filter: {syn_df.height:,}')

    for d in requested:
        extract_one(d, syn_df)

    log(f'Done. Per-site parquets under {PER_SITE_DIR}')


if __name__ == '__main__':
    main()
