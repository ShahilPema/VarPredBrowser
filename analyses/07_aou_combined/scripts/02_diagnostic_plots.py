#!/usr/bin/env python3
"""Diagnostic plots that validate the per-dataset mutation rate scaling.

Three outputs (all written under output/figures/):

  k_summary.csv
    table of (dataset, region, k) — Poisson scaling factors fit on synonymous
    sites. Sanity: k_combined should sit between max(k_gnomad, k_aou) and
    k_gnomad + k_aou (additive in the limit of no overlap, capped at union
    saturation).

  combined_k_sanity.png
    Per-dataset binned (mean(observed) vs mean(expected)) scatter on
    synonymous sites by roulette_AR_MR decile. Each dataset should hug y=x;
    if combined drifts off-diagonal, the union construction has a problem.

  mean_expected_by_percentile.png
    Mirror of analyses/03_mutation_rate_comparison/scripts/03_mean_expected_
    by_percentile.py. 8-panel grid (one per score tag), each panel overlays
    the three datasets' bin_mean(expected) / global_mean(expected) vs
    pathogenicity percentile.

Pure Polars/numpy/matplotlib — no Hail.
"""
import os
os.environ.setdefault('POLARS_MAX_THREADS', '32')

import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import (
    BASE_TABLE_PARQUET, FIG_DIR, OUTPUT_DIR, PER_SITE_DIR, DATASETS, ALL_TAGS,
    TAG_TO_SCORE_COL, per_site_path,
)

SYN_SLICE = PER_SITE_DIR / 'synonymous_augmented.parquet'


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


COLORS = {
    'gnomad_only': '#1f77b4',
    'aou_only':    '#2ca02c',
    'combined':    '#d62728',
}


# -----------------------------------------------------------------------------
# Re-fit k from per-site parquet (cheap; just optimize on aggregated synonymous)
# -----------------------------------------------------------------------------

def fit_k_from_per_site(per_site_df, region):
    """Fit Poisson k on synonymous scaling sites for a region.

    Note: the per-site parquet only contains MISSENSE rows (filtered by
    01_build_per_site.py). For diagnostics we re-load the synonymous slice
    from the base table — see compute_k_for_dataset below.
    """
    raise NotImplementedError('Use compute_k_for_dataset which loads syn from base.')


def load_synonymous_slice():
    """Load synonymous slice produced by 01_build_per_site.py with both gnomAD
    and AoU annotations attached (cannot be derived from the bare base_table
    since AoU columns live only on-workbench)."""
    if not SYN_SLICE.exists():
        raise FileNotFoundError(
            f'{SYN_SLICE} missing. Run scripts/01_build_per_site.py first '
            f'(it writes the synonymous_augmented slice as a side output).'
        )
    log(f'Loading synonymous slice: {SYN_SLICE}')
    df = pl.read_parquet(str(SYN_SLICE))
    log(f'  syn rows: {df.height:,}')
    return df


def compute_observed_for_dataset(syn_df, dataset):
    """Replicate the per-dataset observed rule for synonymous sites.

    Mirrors 01_build_per_site.py:OBSERVED_EXPR. gnomAD's f_lowcov / f_ab_low
    apply to all three legs and are assumed to have been pre-filtered out of
    syn_df by main() before we get here.
    """
    if dataset == 'gnomad_only':
        return (syn_df['gnomad_ac'] > 0).cast(pl.Int32)
    elif dataset == 'aou_only':
        return syn_df['aou_observed'].cast(pl.Int32)
    elif dataset == 'combined':
        return ((syn_df['gnomad_ac'] > 0) | syn_df['aou_observed']).cast(pl.Int32)
    raise ValueError(dataset)


def fit_k_poisson(rates, totals, observed_count):
    rates = np.asarray(rates, dtype=np.float64)
    totals = np.asarray(totals, dtype=np.float64)

    def objective(x):
        expected = np.sum((1 - np.exp(-x[0] * rates)) * totals)
        return abs(expected - observed_count)

    res = minimize(objective, [1.0], method='Nelder-Mead',
                   options={'xatol': 1e-3, 'maxiter': 10000})
    return float(res.x[0])


def compute_k_for_dataset(syn_df, dataset, region):
    """region: 'autosomes_par' or 'chrX_nonpar'."""
    if region == 'chrX_nonpar':
        sub = syn_df.filter(pl.col('is_chrX_nonPAR') == True)
    else:
        sub = syn_df.filter(pl.col('is_chrX_nonPAR') == False).filter(
            pl.col('roulette_syn_scaling_site') == True)
    if sub.height == 0:
        return float('nan')
    obs = compute_observed_for_dataset(sub, dataset)
    sub = sub.with_columns(obs.alias('_obs'))
    agg = (sub.group_by('roulette_AR_MR')
              .agg(pl.col('_obs').sum().alias('polymorphic'),
                   pl.len().alias('total')))
    return fit_k_poisson(
        agg['roulette_AR_MR'].to_numpy(),
        agg['total'].to_numpy(),
        float(agg['polymorphic'].sum()),
    )


# -----------------------------------------------------------------------------
# Plot 1: combined_k_sanity — observed vs expected on synonymous deciles
# -----------------------------------------------------------------------------

def plot_combined_k_sanity(syn_df, k_table, out_path):
    log('Plot: combined_k_sanity ...')
    syn_auto = syn_df.filter(pl.col('is_chrX_nonPAR') == False).filter(
        pl.col('roulette_syn_scaling_site') == True)
    deciles = np.quantile(syn_auto['roulette_AR_MR'].to_numpy(),
                          np.linspace(0, 1, 21))

    fig, ax = plt.subplots(figsize=(6, 6))
    for ds in DATASETS:
        obs_series = compute_observed_for_dataset(syn_auto, ds)
        df = syn_auto.with_columns(
            obs_series.alias('_obs'),
            pl.lit(0).alias('_dummy'),
        )
        rates = df['roulette_AR_MR'].to_numpy()
        # k_auto for this dataset
        k = float(k_table.filter(
            (pl.col('dataset') == ds) & (pl.col('region') == 'autosomes_par')
        )['k'][0])
        exp = 1.0 - np.exp(-k * rates)
        obs = df['_obs'].to_numpy()

        bin_idx = np.digitize(rates, deciles[1:-1])
        bin_obs, bin_exp = [], []
        for b in range(len(deciles) - 1):
            mask = bin_idx == b
            if mask.sum() < 50:
                continue
            bin_obs.append(obs[mask].mean())
            bin_exp.append(exp[mask].mean())
        ax.scatter(bin_exp, bin_obs, color=COLORS[ds], s=30, alpha=0.85,
                   label=f'{ds} (k={k:.3f})')

    lims = [0, max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', lw=1, label='y = x')
    ax.set_xlabel('mean(expected) per Roulette MR decile (synonymous, autosomes+PAR)')
    ax.set_ylabel('mean(observed) per Roulette MR decile')
    ax.set_title('Combined-cohort k sanity: synonymous deciles, observed vs expected')
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log(f'  wrote {out_path}')


# -----------------------------------------------------------------------------
# Plot 2: mean_expected_by_percentile — mirrors 03_mutation_rate_comparison
# -----------------------------------------------------------------------------

def binned_mean_expected(score, expected, n_bins=100):
    order = np.argsort(-score)
    e_ord = expected[order]
    n = len(e_ord)
    if n == 0:
        return np.array([]), np.array([])
    global_mean = e_ord.mean()
    edges = np.linspace(0, n, n_bins + 1, dtype=int)
    xs, ys = [], []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        if hi <= lo:
            continue
        bin_mean = e_ord[lo:hi].mean()
        pct_center = 100.0 * (1.0 - (lo + hi) / (2.0 * n))
        xs.append(pct_center)
        ys.append(float(bin_mean / global_mean))
    return np.array(xs), np.array(ys)


def plot_mean_expected_by_percentile(score_df_by_tag, per_site_by_dataset, out_path):
    log('Plot: mean_expected_by_percentile (8 panels x 3 datasets) ...')
    fig, axes = plt.subplots(2, 4, figsize=(20, 9), sharex=True, sharey=True)
    axes = axes.ravel()

    for k, tag in enumerate(ALL_TAGS):
        ax = axes[k]
        scol = TAG_TO_SCORE_COL[tag]
        for ds in DATASETS:
            ps = per_site_by_dataset.get(ds)
            if ps is None or ps.height == 0:
                continue
            joined = ps.join(
                score_df_by_tag.select(['locus_str', 'alleles_str', 'transcript', scol])
                              .filter(pl.col(scol).is_not_null()),
                on=['locus_str', 'alleles_str', 'transcript'], how='inner',
            )
            if joined.height == 0:
                continue
            score = joined[scol].to_numpy()
            expected = joined['expected'].to_numpy().astype(np.float64)
            x, y = binned_mean_expected(score, expected)
            ax.plot(x, y, color=COLORS[ds], lw=1.4, label=ds)
        ax.axhline(1.0, color='grey', lw=0.5, ls='--')
        ax.set_title(tag, fontsize=10)
        ax.grid(True, alpha=0.3)
        if k % 4 == 0:
            ax.set_ylabel('bin mean(E) / global mean(E)')
        if k >= 4:
            ax.set_xlabel('pathogenicity percentile (high = more pathogenic)')
    axes[-1].legend(fontsize=8, loc='best')
    fig.suptitle('Mean Roulette `expected` per score percentile — three datasets',
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    log(f'  wrote {out_path}')


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    syn_df = load_synonymous_slice()

    # gnomAD's f_lowcov + f_ab_low apply to all legs (per 01_build_per_site.py).
    # Apply the common QC filter so k-fitting and the diagnostic plots use the
    # same site set that the per-site parquets were built from.
    syn_df = syn_df.filter((~pl.col('f_lowcov')) & (~pl.col('f_ab_low')))
    log(f'  syn rows after gnomAD QC filter: {syn_df.height:,}')

    # k table
    log('Fitting k for each (dataset, region)...')
    rows = []
    for ds in DATASETS:
        for region in ('autosomes_par', 'chrX_nonpar'):
            try:
                k = compute_k_for_dataset(syn_df, ds, region)
            except Exception as e:
                log(f'  {ds}/{region}: ERROR {e}')
                k = float('nan')
            rows.append({'dataset': ds, 'region': region, 'k': k})
            log(f'  {ds:14s}  {region:14s}  k={k:.6f}')
    k_df = pl.DataFrame(rows)
    out_csv = FIG_DIR / 'k_summary.csv'
    k_df.write_csv(out_csv)
    log(f'wrote {out_csv}')

    # k sanity scatter
    plot_combined_k_sanity(syn_df, k_df, FIG_DIR / 'combined_k_sanity.png')

    # mean_expected_by_percentile across 3 datasets x 8 tags
    log('Loading per-site parquets for all datasets...')
    per_site_by_dataset = {}
    for ds in DATASETS:
        psp = per_site_path(ds)
        if not psp.exists():
            log(f'  SKIP {ds}: per-site missing')
            continue
        per_site_by_dataset[ds] = pl.read_parquet(str(psp))
        log(f'  {ds}: {per_site_by_dataset[ds].height:,} rows')

    log('Loading scores from base table...')
    score_cols = list(TAG_TO_SCORE_COL.values())
    base_path = (f'{BASE_TABLE_PARQUET}/*.parquet'
                 if BASE_TABLE_PARQUET.is_dir() else str(BASE_TABLE_PARQUET))
    score_df = pl.read_parquet(
        base_path,
        columns=['locus_str', 'alleles_str', 'transcript', *score_cols],
    )
    log(f'  scores: {score_df.height:,} rows')

    plot_mean_expected_by_percentile(
        score_df, per_site_by_dataset,
        FIG_DIR / 'mean_expected_by_percentile.png',
    )

    log('Done.')
    log(f'Figures + table under {FIG_DIR}')


if __name__ == '__main__':
    main()
