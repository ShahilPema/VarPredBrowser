#!/usr/bin/env python3
"""Fit popEVE-style GP calibration curves for all 8 score tags on the three
07_aou_combined datasets: {gnomad_only, aou_only, combined}.

Method: popEVE's ApproximateGP with Polya-Gamma Bernoulli likelihood
(GPModel + PGLikelihood vendored under scripts/popEVE/), fit per-transcript.
Input transform is gene-level PERCENTILE rank (not popEVE's min-max rescale)
so the comparison against the Bernoulli-cloglog P-spline fitter isolates
GP-vs-spline rather than confounding it with the input transform. No
mutation-rate offset (native popEVE: answers P(observed | score)).

Percentile axis definition: each gene's percentile ladder is built from
ALL variants in scores_df with a non-null score for that transcript --
not just the filtered training slice. This ensures the calibrated-score
mapping is defined over the full score range the gene ever sees (so
filtered-out variants get a true rank rather than being pinned to the
training endpoints). The GP is still trained only on the filtered
(observed, expected) subset; what changes is where `p` is defined.

Per-transcript training (verbatim from popEVE's train_popEVE.py):
  - 20 inducing points linearly spaced on [0, 1], locations learned
  - ScaleKernel(RBFKernel), lengthscale init 0.2
  - ZeroMean, PG Bernoulli likelihood
  - NGD lr=0.1 for variational params, Adam lr=0.05 for hyperparameters
  - 6000 epochs (with patience-based early stop on ELBO)

Outputs per (dataset, score_tag):
  output/curves/curves_popeve_gp_{tag}_{dataset}.pkl      # per-tx fits

Per-variant predictions are NOT generated here -- they're produced at BCM
by 09_apply_curves_local.py, which imports the GP apply path from
_popeve_apply.py. Keeps the workbench egress small (just pickles).

Runtime expectation on 80 vCPU (workbench standard): ~30 min per
(score, dataset), ~10 h total for 9 fits. Run per-dataset overnight
with `--datasets gnomad_only` etc.
"""
# --- THREAD GUARDRAILS (must run before any torch / numpy import) -----------
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMBA_NUM_THREADS'] = '1'
# -----------------------------------------------------------------------------

import argparse
import io
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from joblib import Parallel, delayed


# =============================================================================
# Paths
# =============================================================================

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config.paths import (  # noqa
    BASE_TABLE_PARQUET, PER_SITE_DIR, CURVES_DIR, DATASETS,
    INHOUSE_TAGS, EXTERNAL_TAGS, ALL_TAGS, TAG_TO_SCORE_COL, per_site_path,
)

# popEVE module is vendored next to this script (scripts/popEVE/)
POPEVE_REPO = Path(__file__).resolve().parent / 'popEVE'


# =============================================================================
# Config
# =============================================================================

MIN_VARIANTS = 20
M_INDUCING = 20
NUM_EPOCHS = 6000
LR_NGD = 0.1
LR_HYPER = 0.05
LENGTHSCALE_INIT = 0.2
CONV_TOL = 1e-5
CONV_PATIENCE = 100       # epochs with no ELBO improvement -> stop
MIN_EPOCHS = 500          # always run at least this many before early-stop
# Predict returns the GP posterior MEAN at each point (a single calibrated
# logit per variant). No sampling -- the mean is cheap (O(N*M)) while
# sampling requires O(N^3) Cholesky of the posterior covariance which was
# the apply-step bottleneck in earlier runs.

ALL_MODELS = [(tag, TAG_TO_SCORE_COL[tag]) for tag in ALL_TAGS]


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


# =============================================================================
# Worker-side fit. Helpers (_init_worker, _build_model, scores_to_pct) are
# imported from _popeve_apply so the apply path at BCM (09_apply_curves_local.py)
# uses the IDENTICAL model constructor — state_dicts load cleanly across both.
# =============================================================================

from _popeve_apply import (  # noqa
    _init_worker, _build_model, scores_to_pct, M_INDUCING as _M_INDUCING_CHECK,
)
assert _M_INDUCING_CHECK == M_INDUCING, 'M_INDUCING mismatch with _popeve_apply'


def _fit_one(tx, s, obs, sorted_full, seed=42):
    """Fit popEVE GP on one transcript.

    `sorted_full` is this gene's full sorted score universe (all scored
    variants with a non-null score in scores_df for this transcript),
    used as the percentile axis. Training uses (s, obs) from the filtered
    per-site set; `p` is the rank of each filtered-set score within the
    full gene universe.
    """
    _init_worker()
    import torch
    import gpytorch

    torch.manual_seed(seed)

    # Gene-level percentile against the FULL per-gene score universe
    p = scores_to_pct(s, sorted_full)
    train_x = torch.tensor(p.reshape(-1, 1), dtype=torch.float32)
    train_y = torch.tensor(obs, dtype=torch.float32)

    inducing = torch.linspace(0.0, 1.0, M_INDUCING, dtype=train_x.dtype).unsqueeze(-1)
    model, likelihood = _build_model(inducing)

    variational_opt = gpytorch.optim.NGD(
        model.variational_parameters(), num_data=train_y.size(0), lr=LR_NGD)
    hyper_opt = torch.optim.Adam(
        [{'params': model.hyperparameters()}, {'params': likelihood.parameters()}],
        lr=LR_HYPER)

    model.train()
    likelihood.train()
    mll = gpytorch.mlls.VariationalELBO(likelihood, model, num_data=train_y.size(0))

    best_loss = float('inf')
    stale = 0
    final_loss = float('nan')
    n_iter = 0
    for i in range(NUM_EPOCHS):
        variational_opt.zero_grad()
        hyper_opt.zero_grad()
        out = model(train_x)
        loss = -mll(out, train_y)
        loss.backward()
        variational_opt.step()
        hyper_opt.step()
        loss_val = float(loss.item())
        final_loss = loss_val
        n_iter = i + 1
        if i < MIN_EPOCHS:
            best_loss = min(best_loss, loss_val)
            continue
        if loss_val + CONV_TOL < best_loss:
            best_loss = loss_val
            stale = 0
        else:
            stale += 1
        if stale >= CONV_PATIENCE:
            break

    # serialize state_dict to bytes so it pickles cleanly across workers
    buf = io.BytesIO()
    torch.save({k: v.detach().cpu() for k, v in model.state_dict().items()}, buf)
    state_bytes = buf.getvalue()

    return {
        'transcript': tx,
        'state_bytes': state_bytes,
        'sorted_train_scores': np.asarray(sorted_full, dtype=np.float64),
        'final_loss': final_loss,
        'n_iter': n_iter,
        'n_train': int(len(s)),
        'n_observed': int(obs.sum()),
        'n_universe': int(len(sorted_full)),
        'lengthscale': float(model.covar_module.base_kernel.lengthscale.item()),
        'outputscale': float(model.covar_module.outputscale.item()),
    }


# =============================================================================
# Record building
# =============================================================================

def compute_sorted_universe(scores_df, score_col):
    """Per-transcript sorted array of all non-null scores in scores_df.

    This is the percentile axis the GP will be defined on. It is the full
    per-gene score universe -- every variant scores_df has a non-null score
    for, regardless of filter/consequence.
    """
    t0 = time.perf_counter()
    scored = (scores_df
              .select(['transcript', score_col])
              .filter(pl.col(score_col).is_not_null()))
    # Sorted per-transcript in one shot via Polars groupby-agg
    agg = (scored.group_by('transcript')
                 .agg(pl.col(score_col).sort().alias('sorted_scores')))
    out = {}
    for row in agg.iter_rows(named=True):
        arr = np.asarray(row['sorted_scores'], dtype=np.float64)
        if len(arr) >= 2:
            out[row['transcript']] = arr
    log(f'    universe sorted: {len(out):,} txs in {time.perf_counter()-t0:.1f}s')
    return out


def build_records(per_site_df, scores_df, score_col, sorted_by_tx, max_tx=None):
    """Filtered-set (score, observed) pairs plus each gene's full sorted universe."""
    df = per_site_df.join(
        scores_df.select(['locus_str', 'alleles_str', 'transcript', score_col]),
        on=['locus_str', 'alleles_str', 'transcript'], how='inner',
    ).filter(pl.col(score_col).is_not_null())
    groups = df.partition_by('transcript', maintain_order=False)
    records = []
    for g in groups:
        if g.height < MIN_VARIANTS:
            continue
        g2 = g.sort(score_col)
        tx = g2['transcript'][0]
        sorted_full = sorted_by_tx.get(tx)
        if sorted_full is None or len(sorted_full) < MIN_VARIANTS:
            continue
        s = g2[score_col].to_numpy().astype(np.float64)
        obs = g2['observed'].fill_null(0).cast(pl.Float64).to_numpy()
        records.append((tx, s, obs, sorted_full))
    if max_tx is not None:
        records = records[:max_tx]
    return records


# =============================================================================
# Fit for a single (tag, dataset). The apply step lives in 09_apply_curves_local.py
# and reuses _popeve_apply.
# =============================================================================

def fit_all(records, n_workers):
    t0 = time.perf_counter()
    results = Parallel(n_jobs=n_workers, prefer='processes', verbose=0)(
        delayed(_fit_one)(tx, s, obs, sorted_full)
        for tx, s, obs, sorted_full in records
    )
    elapsed = time.perf_counter() - t0
    log(f'    fit {len(results)} / {len(records)} transcripts in {elapsed:.1f}s '
        f'({elapsed / max(1, len(results)):.2f}s / tx)')
    return {r['transcript']: r for r in results if r is not None}


def run_one(dataset, tag, score_col, scores_df, per_site_df, n_workers,
            sorted_by_tx, force=False, max_tx=None, curves_suffix=''):
    curves_path = CURVES_DIR / f'curves_popeve_gp_{tag}_{dataset}{curves_suffix}.pkl'
    if curves_path.exists() and not force:
        log(f'  SKIP {tag} / {dataset} (curves exist)')
        return

    log(f'=== popeve_gp / {tag} / {dataset} (score={score_col}) ===')
    records = build_records(per_site_df, scores_df, score_col, sorted_by_tx,
                            max_tx=max_tx)
    log(f'  {len(records)} transcripts with >= {MIN_VARIANTS} scored sites'
        f'{f" (CAPPED to {max_tx})" if max_tx else ""}')
    if not records:
        return

    fits = fit_all(records, n_workers)
    with open(curves_path, 'wb') as f:
        pickle.dump(fits, f, protocol=pickle.HIGHEST_PROTOCOL)
    log(f'  wrote {curves_path}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=80,
                    help='joblib n_jobs. AoU workbench cap is 96 vCPU; '
                         'default 80 leaves headroom for OS / notebook process.')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--tags', type=str, default=None,
                    help='Comma-separated subset of score tags to run')
    ap.add_argument('--datasets', type=str, default=None,
                    help='Comma-separated subset of datasets to run')
    ap.add_argument('--max-tx', type=int, default=None,
                    help='SMOKE TEST: cap at this many transcripts per (tag, dataset)')
    ap.add_argument('--smoke-suffix', type=str, default='',
                    help='Suffix appended to output filenames (e.g. "_smoke")')
    args = ap.parse_args()

    datasets = args.datasets.split(',') if args.datasets else list(DATASETS)
    wanted_tags = set(args.tags.split(',')) if args.tags else None

    log(f'popEVE GP fit: workers={args.workers}, epochs={NUM_EPOCHS} '
        f'(patience={CONV_PATIENCE}), datasets={datasets}')

    # Load all 8 score columns from the consolidated base table once.
    log(f'Loading scores from {BASE_TABLE_PARQUET}')
    score_cols = [TAG_TO_SCORE_COL[t] for t, _ in ALL_MODELS]
    scores_df = pl.read_parquet(
        f'{BASE_TABLE_PARQUET}/*.parquet',
        columns=['locus_str', 'alleles_str', 'transcript', *score_cols],
    )
    log(f'  scores: {scores_df.height:,} rows')

    # Pre-compute per-tx sorted score universes once per tag (dataset-independent).
    log('Pre-computing per-gene full sorted score universes...')
    universe_by_tag = {}
    for tag, col in ALL_MODELS:
        if wanted_tags and tag not in wanted_tags:
            continue
        log(f'  universe for {tag} ({col})')
        universe_by_tag[tag] = compute_sorted_universe(scores_df, col)

    for dataset in datasets:
        psp = per_site_path(dataset)
        if not psp.exists():
            log(f'SKIP {dataset}: per-site parquet missing ({psp})')
            continue
        log(f'Loading per_site: {psp}')
        per_site_df = pl.read_parquet(f'{psp}/*.parquet')
        log(f'  per_site: {per_site_df.height:,} rows  '
            f'(O rate: {per_site_df["observed"].mean():.4f})')

        for tag, col in ALL_MODELS:
            if wanted_tags and tag not in wanted_tags:
                continue
            run_one(dataset, tag, col, scores_df, per_site_df,
                    args.workers, universe_by_tag[tag], force=args.force,
                    max_tx=args.max_tx, curves_suffix=args.smoke_suffix)

        del per_site_df

    log('All done.')


if __name__ == '__main__':
    main()
