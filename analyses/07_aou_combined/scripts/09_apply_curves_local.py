#!/usr/bin/env python3
"""BCM-side: apply pickled curves to the local score table.

Pulls all curves_*.pkl from output/curves/ and produces per-variant
fitted_oe parquets next to them. Run this at BCM after extracting the
tarball produced by 08_pack_outputs.py at the workbench.

For IRLS-QP family pickles ({cloglog, quadprog, poisson, plm_alone}):
  fitted_oe = predict(model, score) -> exp(B @ coeffs), the rate-scale
  O/E modifier (matches the prior Poisson application path).

For popEVE GP pickles, this script defers to a separate apply path
since GP prediction needs gpytorch and chunked dispatch — see
07_fit_popeve_gp.py:apply_curves which produces the parquet inline at
fit time. If you only have the pickle and need to re-apply, run
07_fit_popeve_gp.py with --force on a stripped-down score-only input.

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
    r'^curves_(?P<method>cloglog|quadprog|poisson|plm_alone)_'
    r'(?P<tag>[^_]+(?:_[^_]+)*?)_'
    r'(?P<dataset>gnomad_only|aou_only|combined)\.pkl$'
)


def parse_pkl(name):
    m = CURVES_PATTERN.match(name)
    if not m:
        return None
    return m.group('method'), m.group('tag'), m.group('dataset')


def apply_one(pkl_path, score_col, scores_df):
    fitted_path = pkl_path.with_name(
        pkl_path.name.replace('curves_', 'fitted_').replace('.pkl', '.parquet')
    )
    if fitted_path.exists():
        log(f'SKIP {pkl_path.name} (fitted exists)')
        return

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=8,
                    help='Parallel pickles processed in parallel (default 8). Each worker '
                         'reads a copy of scores_df, so memory is workers x scores_df_size.')
    args = ap.parse_args()

    pkls = sorted(CURVES_DIR.glob('curves_*.pkl'))
    pkls = [p for p in pkls if parse_pkl(p.name) is not None]
    log(f'Found {len(pkls)} IRLS-QP family pickles to apply')

    score_cols_needed = set()
    parsed = []
    for p in pkls:
        meta = parse_pkl(p.name)
        method, tag, dataset = meta
        score_cols_needed.add(TAG_TO_SCORE_COL[tag])
        parsed.append((p, meta))

    log(f'Loading score table from {BASE_TABLE_PARQUET}')
    scores_df = pl.read_parquet(
        f'{BASE_TABLE_PARQUET}/*.parquet',
        columns=['locus_str', 'alleles_str', 'transcript', *sorted(score_cols_needed)],
    )
    log(f'  scores: {scores_df.height:,} rows')

    if args.workers <= 1:
        for p, (method, tag, dataset) in parsed:
            apply_one(p, TAG_TO_SCORE_COL[tag], scores_df)
    else:
        Parallel(n_jobs=args.workers, prefer='processes')(
            delayed(apply_one)(p, TAG_TO_SCORE_COL[tag], scores_df)
            for p, (method, tag, dataset) in parsed
        )

    log(f'Done. Per-variant fitted_oe parquets under {CURVES_DIR}')


if __name__ == '__main__':
    main()
