#!/usr/bin/env python3
"""BCM-side: apply pickled curves to the local score table.

Pulls all curves_*.pkl from output/curves/ and produces per-variant
fitted parquets next to them. Run this at BCM after extracting the
tarball produced by 08_pack_outputs.py at the workbench.

Two pickle families:
  - IRLS-QP family ({cloglog, quadprog}): pickled numpy arrays + scalars.
    Apply uses fitters.predict() -> writes fitted_oe (rate-scale O/E).
  - popEVE GP: pickled torch state_dicts. Apply uses _popeve_apply ->
    writes (percentile, gp_mean, mean_prob). Requires torch + gpytorch
    + the vendored popEVE module (all already installed at BCM).

Run: python3 -u 09_apply_curves_local.py [--workers 64]
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')

import argparse
import pickle
import re
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
    BASE_TABLE_PARQUET, CURVES_DIR, TAG_TO_SCORE_COL,
)


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


CURVES_PATTERN = re.compile(
    r'^curves_(?P<method>cloglog|quadprog|popeve_gp)_'
    r'(?P<tag>[^_]+(?:_[^_]+)*?)_'
    r'(?P<dataset>gnomad_only|aou_only|combined)\.pkl$'
)


def parse_pkl(name):
    m = CURVES_PATTERN.match(name)
    if not m:
        return None
    return m.group('method'), m.group('tag'), m.group('dataset')


def _apply_irls_qp(pkl_path, score_col, scores_df, fitted_path):
    """cloglog / quadprog: predict() returns fitted_oe (rate-scale O/E)."""
    log(f'Loading {pkl_path}')
    with open(pkl_path, 'rb') as f:
        models = pickle.load(f)
    if not models:
        log('  empty pickle; skipping')
        return

    df = (scores_df
          .filter(pl.col(score_col).is_not_null())
          .filter(pl.col('transcript').is_in(list(models.keys()))))
    if df.height == 0:
        log('  no rows to score')
        return

    groups = df.partition_by('transcript', maintain_order=False)
    out_frames = []
    for g in groups:
        tx = g['transcript'][0]
        m = models.get(tx)
        if m is None:
            continue
        scores = g[score_col].to_numpy().astype(np.float64)
        try:
            oe = F.predict(m, scores)
        except Exception as e:
            log(f'  predict failed for tx={tx}: {e}')
            continue
        out_frames.append(
            g.select(['locus_str', 'alleles_str', 'transcript'])
             .with_columns(pl.Series('fitted_oe', oe.astype(np.float32)))
        )

    if not out_frames:
        log('  no per-tx outputs; skipping write')
        return

    out_df = pl.concat(out_frames)
    out_df = out_df.with_columns(
        pl.when(pl.col('fitted_oe').is_nan()).then(None)
          .otherwise(pl.col('fitted_oe')).alias('fitted_oe')
    )
    out_df.write_parquet(str(fitted_path))
    log(f'  wrote {fitted_path} ({out_df.height:,} rows)')


def _apply_popeve_gp(pkl_path, score_col, scores_df, fitted_path, n_workers):
    """popeve_gp: torch state_dicts -> (percentile, gp_mean, mean_prob)."""
    from _popeve_apply import apply_curves as gp_apply_curves
    log(f'Loading {pkl_path}')
    with open(pkl_path, 'rb') as f:
        fits = pickle.load(f)
    if not fits:
        log('  empty pickle; skipping')
        return
    label = pkl_path.stem.replace('curves_', '')
    gp_apply_curves(scores_df, score_col, fits, str(fitted_path), label, n_workers)


def apply_one(pkl_path, score_col, scores_df, gp_workers=8):
    fitted_path = pkl_path.with_name(
        pkl_path.name.replace('curves_', 'fitted_').replace('.pkl', '.parquet')
    )
    if fitted_path.exists():
        log(f'SKIP {pkl_path.name} (fitted exists)')
        return

    method, _, _ = parse_pkl(pkl_path.name)
    if method == 'popeve_gp':
        _apply_popeve_gp(pkl_path, score_col, scores_df, fitted_path, gp_workers)
    else:
        _apply_irls_qp(pkl_path, score_col, scores_df, fitted_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=8,
                    help='Parallel IRLS-QP pickles run in parallel (default 8). Each '
                         'worker reads a copy of scores_df, so memory is workers x scores_df_size.')
    ap.add_argument('--gp-workers', type=int, default=28,
                    help='Workers per popeve_gp apply (run serially across pickles to '
                         'avoid nested process pools). Default 28 to match shared-node norms.')
    args = ap.parse_args()

    pkls = sorted(CURVES_DIR.glob('curves_*.pkl'))
    pkls = [p for p in pkls if parse_pkl(p.name) is not None]
    log(f'Found {len(pkls)} pickles total')

    score_cols_needed = set()
    parsed = []
    for p in pkls:
        meta = parse_pkl(p.name)
        method, tag, dataset = meta
        score_cols_needed.add(TAG_TO_SCORE_COL[tag])
        parsed.append((p, meta))

    log(f'Loading score table from {BASE_TABLE_PARQUET}')
    base_path = (f'{BASE_TABLE_PARQUET}/*.parquet'
                 if BASE_TABLE_PARQUET.is_dir() else str(BASE_TABLE_PARQUET))
    scores_df = pl.read_parquet(
        base_path,
        columns=['locus_str', 'alleles_str', 'transcript', *sorted(score_cols_needed)],
    )
    log(f'  scores: {scores_df.height:,} rows')

    irls = [(p, m) for p, m in parsed if m[0] != 'popeve_gp']
    gp   = [(p, m) for p, m in parsed if m[0] == 'popeve_gp']
    log(f'  IRLS-QP pickles: {len(irls)}   popeve_gp pickles: {len(gp)}')

    # IRLS-QP: serial inside, parallel across pickles
    if irls:
        log(f'Applying IRLS-QP pickles ({args.workers}-way across pickles)...')
        if args.workers <= 1:
            for p, (method, tag, dataset) in irls:
                apply_one(p, TAG_TO_SCORE_COL[tag], scores_df)
        else:
            Parallel(n_jobs=args.workers, prefer='processes')(
                delayed(apply_one)(p, TAG_TO_SCORE_COL[tag], scores_df)
                for p, (method, tag, dataset) in irls
            )

    # popeve_gp: parallel inside (joblib in _popeve_apply), serial across pickles
    if gp:
        log(f'Applying popeve_gp pickles serially (each uses {args.gp_workers} workers internally)...')
        for p, (method, tag, dataset) in gp:
            apply_one(p, TAG_TO_SCORE_COL[tag], scores_df, gp_workers=args.gp_workers)

    log(f'Done. Per-variant fitted parquets under {CURVES_DIR}')


if __name__ == '__main__':
    main()
