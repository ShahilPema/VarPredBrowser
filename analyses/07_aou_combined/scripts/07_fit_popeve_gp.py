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
  output/curves/fitted_popeve_gp_{tag}_{dataset}.parquet  # per-variant preds

Runtime expectation on 60 vCPU (workbench standard): ~30 min per
(score, dataset), ~12 h total for 24 fits. Run per-dataset overnight
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
# Per-gene percentile transform (matches fitters.py:scores_to_percentiles)
# =============================================================================

def scores_to_pct(x, sorted_train_scores):
    return np.searchsorted(sorted_train_scores, x, side='right') / len(sorted_train_scores)


# =============================================================================
# Worker-side fit / predict. Imports torch / gpytorch LAZILY inside the worker
# so the thread-guardrail env vars are applied before torch initializes its
# threadpool.
# =============================================================================

def _init_worker():
    """Set single-threaded torch in every worker process."""
    sys.path.insert(0, str(POPEVE_REPO))
    import torch
    torch.set_num_threads(1)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        # set_num_interop_threads can only be called once per process
        pass


def _build_model(inducing_points):
    """Constructor used both at fit time (scratch) and predict time (from state).
    Kept identical so state_dict loads cleanly."""
    from popEVE.popEVE import GPModel, PGLikelihood  # noqa
    model = GPModel(inducing_points=inducing_points)
    model.covar_module.base_kernel.initialize(lengthscale=LENGTHSCALE_INIT)
    likelihood = PGLikelihood()
    return model, likelihood


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


def _predict_chunk(chunk, seed=42):
    """Score a chunk of (fit_dict, scores_array) pairs. Returns a list of
    per-tx prediction dicts in the same order as the input chunk.

    Outputs the GP posterior MEAN at each point (a single calibrated logit
    per variant) plus its sigmoid. No posterior sampling -- we only need
    one value per input score.

    Chunked dispatch (vs per-tx) reduces joblib IPC overhead: with 18K
    per-tx tasks the main process bottlenecks on pickling ~225KB fit
    payloads into the dispatch queue. Chunks of ~50 tx reduce that to
    ~360 dispatches of ~11 MB each.
    """
    _init_worker()
    import torch

    torch.manual_seed(seed)
    inducing = torch.linspace(0.0, 1.0, M_INDUCING, dtype=torch.float32).unsqueeze(-1)
    model, likelihood = _build_model(inducing)
    model.eval()
    likelihood.eval()

    results = []
    for fit, scores in chunk:
        sorted_s = fit['sorted_train_scores']
        p = scores_to_pct(scores, sorted_s)
        p = np.clip(p, 0.0, 1.0)
        x = torch.tensor(p.reshape(-1, 1), dtype=torch.float32)

        state = torch.load(io.BytesIO(fit['state_bytes']), weights_only=True)
        model.load_state_dict(state)

        with torch.no_grad():
            out = model(x)
            # out.mean is the posterior marginal mean at each input.
            # Computing it is O(N*M) -- kernel block + one matmul -- no
            # O(N^3) Cholesky of the posterior covariance required.
            gp_mean = out.mean.numpy().astype(np.float32)
            mean_prob = torch.sigmoid(out.mean).numpy().astype(np.float32)

        results.append({
            'transcript': fit['transcript'],
            'percentile': p.astype(np.float32),
            'gp_mean': gp_mean,
            'mean_prob': mean_prob,
        })
    return results


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
# Fit + apply for a single (tag, dataset)
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


def apply_curves(scores_df, score_col, fits, out_path, label, n_workers):
    df = (scores_df
          .filter(pl.col(score_col).is_not_null())
          .filter(pl.col('transcript').is_in(list(fits.keys()))))
    if df.height == 0:
        log(f'    [{label}] no rows to score')
        return

    t_part = time.perf_counter()
    groups = df.partition_by('transcript', maintain_order=False)
    log(f'    [{label}] partitioned {df.height:,} rows into '
        f'{len(groups)} tx groups in {time.perf_counter()-t_part:.1f}s')

    # Build per-tx records: (fit_dict, scores_array, keys_df)
    records = []
    for g in groups:
        tx = g['transcript'][0]
        fit = fits.get(tx)
        if fit is None:
            continue
        records.append((fit, g[score_col].to_numpy().astype(np.float64),
                        g.select(['locus_str', 'alleles_str', 'transcript'])))
    log(f'    [{label}] built {len(records)} predict records')

    # Chunk records. With 18K tx and ~180 workers we want ~2x workers in
    # chunks so fast workers can pick up extra chunks when slower ones
    # finish. Each chunk then carries a slab of fits + arrays (~10-20 MB)
    # which is cheap to pickle.
    n_chunks = min(n_workers * 2, max(1, len(records)))
    idxs = np.arange(len(records))
    chunk_indices = [idxs[i::n_chunks] for i in range(n_chunks)]
    chunk_payloads = [
        [(records[i][0], records[i][1]) for i in ixs] for ixs in chunk_indices
    ]

    t0 = time.perf_counter()
    chunk_preds = Parallel(n_jobs=n_workers, prefer='processes', verbose=0)(
        delayed(_predict_chunk)(payload) for payload in chunk_payloads
    )
    log(f'    [{label}] scored {sum(len(cp) for cp in chunk_preds)} txs '
        f'in {time.perf_counter()-t0:.1f}s ({n_chunks} chunks)')

    # Stitch: preds come back in chunk order; each chunk is in chunk_indices order.
    out_frames = [None] * len(records)
    for ixs, preds in zip(chunk_indices, chunk_preds):
        for i, pred in zip(ixs, preds):
            _fit, _sc, keys = records[i]
            out_frames[i] = keys.with_columns([
                pl.Series('percentile', pred['percentile']),
                pl.Series('gp_mean',    pred['gp_mean']),
                pl.Series('mean_prob',  pred['mean_prob']),
            ])
    out_df = pl.concat(out_frames)
    out_df.write_parquet(out_path)
    log(f'    [{label}] wrote {out_path} ({out_df.height:,} rows)')


def run_one(dataset, tag, score_col, scores_df, per_site_df, n_workers,
            sorted_by_tx, force=False, max_tx=None, curves_suffix=''):
    curves_path = CURVES_DIR / f'curves_popeve_gp_{tag}_{dataset}{curves_suffix}.pkl'
    fitted_path = CURVES_DIR / f'fitted_popeve_gp_{tag}_{dataset}{curves_suffix}.parquet'
    if curves_path.exists() and fitted_path.exists() and not force:
        log(f'  SKIP {tag} / {dataset} (outputs present)')
        return

    log(f'=== popeve_gp / {tag} / {dataset} (score={score_col}) ===')
    records = build_records(per_site_df, scores_df, score_col, sorted_by_tx,
                            max_tx=max_tx)
    log(f'  {len(records)} transcripts with >= {MIN_VARIANTS} scored sites'
        f'{f" (CAPPED to {max_tx})" if max_tx else ""}')
    if not records:
        return

    if not curves_path.exists() or force:
        fits = fit_all(records, n_workers)
        with open(curves_path, 'wb') as f:
            pickle.dump(fits, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(curves_path, 'rb') as f:
            fits = pickle.load(f)
        log(f'  loaded cached curves ({len(fits)} txs)')

    if not fitted_path.exists() or force:
        apply_curves(scores_df, score_col, fits, str(fitted_path),
                     f'popeve_gp/{tag}/{dataset}', n_workers)


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
