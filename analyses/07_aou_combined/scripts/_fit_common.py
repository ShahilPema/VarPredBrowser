"""Shared fit/apply skeleton used by 03..06 (IRLS-QP family).

popEVE GP (07) has a different per-tx training loop and lives in its own script.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('NUMBA_NUM_THREADS', '1')

import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fitters as F  # noqa
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import (  # noqa
    PER_SITE_DIR, CURVES_DIR, OUTPUT_DIR, BASE_TABLE_PARQUET,
    DATASETS, INHOUSE_TAGS, EXTERNAL_TAGS, ALL_TAGS, TAG_TO_SCORE_COL,
    per_site_path,
)


MIN_VARIANTS = 20
N_KNOTS = 24

# Quadprog needs a windowed O/E pre-step (Poisson-style aggregation along the
# sorted-by-score axis). Defaults match 02_poisson_quadprog_4models/scripts/
# 03_fit_and_apply.py exactly.
QUADPROG_MAX_WINDOW = 150
QUADPROG_MIN_EXPECTED = 20.0

# Hyperparameter defaults — match the canonical analysis dir for each method.
DEFAULT_LAMBDA = {
    'cloglog':   200.0,   # 04_bernoulli_fitter / 07_cloglog_sweep.py
    'quadprog':  0.5,     # 02_poisson_quadprog_4models grid
    'poisson':   500.0,   # 02_poisson_quadprog_4models grid
    'plm_alone': 200.0,   # 04_bernoulli_fitter / 04_lambda_sweep_plm_alone.py
}


def compute_fixed_window_oe(obs, exp, max_window=QUADPROG_MAX_WINDOW,
                            min_expected=QUADPROG_MIN_EXPECTED):
    """numpy prefix-sum windowed O/E on a sorted-by-score array (quadprog input)."""
    n = len(obs)
    cum_obs = np.concatenate([[0.0], np.cumsum(obs)])
    cum_exp = np.concatenate([[0.0], np.cumsum(exp)])
    idx = np.arange(n)
    lo = np.maximum(idx - max_window, 0)
    hi = np.minimum(idx + max_window + 1, n)
    w_obs = cum_obs[hi] - cum_obs[lo]
    w_exp = cum_exp[hi] - cum_exp[lo]
    oe = np.where(w_exp >= min_expected, w_obs / w_exp, np.nan)
    return oe, w_exp

# Which tags each method runs on
TAGS_PER_METHOD = {
    'cloglog':   ALL_TAGS,
    'quadprog':  INHOUSE_TAGS,
    'poisson':   INHOUSE_TAGS,
    'plm_alone': ALL_TAGS,
}


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


# -----------------------------------------------------------------------------
# Universe (per-tx sorted score arrays for percentile rank at predict time)
# -----------------------------------------------------------------------------

def load_scores_and_universe(score_col):
    """Load (locus_str, alleles_str, transcript, score_col) from base_table and
    build {tx -> sorted_score_array} for percentile axis at predict time."""
    log(f'Loading scores ({score_col}) from base_table...')
    df = (pl.scan_parquet(f'{BASE_TABLE_PARQUET}/*.parquet')
            .select(['locus_str', 'alleles_str', 'transcript', score_col])
            .filter(pl.col(score_col).is_not_null())
            .collect())
    log(f'  scored rows: {df.height:,}')
    universe = F.compute_sorted_universe(df, score_col)
    return df, universe


# -----------------------------------------------------------------------------
# Build per-tx records
# -----------------------------------------------------------------------------

def build_records(per_site_df, scores_df, score_col, universe):
    df = per_site_df.join(
        scores_df.select(['locus_str', 'alleles_str', score_col]),
        on=['locus_str', 'alleles_str'], how='inner',
    ).filter(pl.col(score_col).is_not_null())
    groups = df.partition_by('transcript', maintain_order=False)
    records = []
    for g in groups:
        if g.height < MIN_VARIANTS:
            continue
        g2 = g.sort(score_col)
        tx = g2['transcript'][0]
        sorted_full = universe.get(tx)
        if sorted_full is None or len(sorted_full) < MIN_VARIANTS:
            continue
        s = g2[score_col].to_numpy().astype(np.float64)
        obs = g2['observed'].fill_null(0).cast(pl.Float64).to_numpy()
        exp = g2['expected'].fill_null(0.0).cast(pl.Float64).to_numpy()
        records.append((tx, s, obs, exp, sorted_full))
    return records


# -----------------------------------------------------------------------------
# Fit chunk dispatch
# -----------------------------------------------------------------------------

def _fit_chunk(method, records, lam):
    out = {}
    for tx, s, obs, exp, sorted_full in records:
        try:
            if method == 'cloglog':
                keep = exp > 0
                if keep.sum() < MIN_VARIANTS:
                    continue
                m = F.fit_bernoulli_cloglog_irls_qp(
                    s[keep], obs[keep], exp[keep],
                    lam=lam, n_knots=N_KNOTS, sorted_universe=sorted_full)
            elif method == 'poisson':
                keep = exp > 0
                if keep.sum() < MIN_VARIANTS:
                    continue
                m = F.fit_poisson_irls_qp(
                    s[keep], obs[keep], exp[keep],
                    lam=lam, n_knots=N_KNOTS, sorted_universe=sorted_full)
            elif method == 'quadprog':
                # Records are already sorted by score in build_records; window
                # uses cumulative sums along that axis.
                oe, w_exp = compute_fixed_window_oe(obs, exp)
                valid = ~np.isnan(oe) & (w_exp > 0)
                if valid.sum() < MIN_VARIANTS:
                    continue
                m = F.fit_qp_quadprog(
                    s[valid], oe[valid], w_exp[valid],
                    lam=lam, n_knots=N_KNOTS, sorted_universe=sorted_full)
            elif method == 'plm_alone':
                # No offset — fits P(observed | score). Use all sites with score.
                m = F.fit_bernoulli_no_offset_irls_qp(
                    s, obs,
                    lam=lam, n_knots=N_KNOTS, sorted_universe=sorted_full)
            else:
                raise ValueError(method)
        except Exception:
            m = None
        if m is not None:
            out[tx] = m
    return out


def fit_all(method, records, lam, n_workers):
    chunks = [records[i::n_workers] for i in range(n_workers)]
    t0 = time.perf_counter()
    results = Parallel(n_jobs=n_workers, prefer='processes')(
        delayed(_fit_chunk)(method, c, lam) for c in chunks
    )
    elapsed = time.perf_counter() - t0
    models = {}
    for r in results:
        models.update(r)
    log(f'  fit {len(models)} / {len(records)} transcripts in {elapsed:.1f}s')
    return models


# -----------------------------------------------------------------------------
# Output paths
# -----------------------------------------------------------------------------

def curves_path(method, tag, dataset):
    return CURVES_DIR / f'curves_{method}_{tag}_{dataset}.pkl'


# -----------------------------------------------------------------------------
# Run loop
# -----------------------------------------------------------------------------

def run_method(method, datasets, tags, lam, n_workers, force=False):
    valid_tags = TAGS_PER_METHOD[method]
    requested_tags = [t for t in tags if t in valid_tags]
    skipped_tags = [t for t in tags if t not in valid_tags]
    if skipped_tags:
        log(f'note: {method} does not run on tags {skipped_tags}; skipping those')

    for tag in requested_tags:
        score_col = TAG_TO_SCORE_COL[tag]
        scores_df = None
        universe = None
        for dataset in datasets:
            cp = curves_path(method, tag, dataset)
            if cp.exists() and not force:
                log(f'SKIP {method}/{tag}/{dataset}: {cp.name} exists')
                continue

            psp = per_site_path(dataset)
            if not (psp / '_SUCCESS').exists() and not psp.exists():
                log(f'SKIP {method}/{tag}/{dataset}: per-site missing ({psp})')
                continue

            if scores_df is None:
                scores_df, universe = load_scores_and_universe(score_col)

            log(f'=== {method} / {tag} / {dataset} (lam={lam}) ===')
            per_site_df = pl.read_parquet(f'{psp}/*.parquet')
            log(f'  per_site: {per_site_df.height:,} rows')

            records = build_records(per_site_df, scores_df, score_col, universe)
            log(f'  {len(records)} transcripts with >= {MIN_VARIANTS} scored sites')
            del per_site_df
            if not records:
                continue

            models = fit_all(method, records, lam, n_workers)
            with open(cp, 'wb') as f:
                pickle.dump(models, f)
            log(f'  wrote {cp}')
            del records, models


def standard_argparser(method):
    ap = argparse.ArgumentParser()
    ap.add_argument('--datasets', default=','.join(DATASETS))
    ap.add_argument('--tags', default=','.join(TAGS_PER_METHOD[method]))
    ap.add_argument('--lam', type=float, default=DEFAULT_LAMBDA[method])
    ap.add_argument('--workers', type=int, default=80,
                    help='joblib n_jobs. AoU workbench cap is 96 vCPU; '
                         'default 80 leaves headroom for OS / notebook process.')
    ap.add_argument('--force', action='store_true')
    return ap


def standard_main(method):
    ap = standard_argparser(method)
    args = ap.parse_args()
    datasets = args.datasets.split(',')
    tags = args.tags.split(',')
    log(f'method={method}  datasets={datasets}  tags={tags}  '
        f'lam={args.lam}  workers={args.workers}')
    run_method(method, datasets, tags, args.lam, args.workers, args.force)
    log(f'Done. Curves under {CURVES_DIR}')
