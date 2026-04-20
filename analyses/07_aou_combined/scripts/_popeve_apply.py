"""popEVE GP apply path — used by both the workbench fit script (07) and the
BCM-side per-variant scorer (09). Importing this module touches torch /
gpytorch / popEVE only when its functions are actually called (lazy imports
inside _init_worker / _build_model), so it's safe to import on systems
without those dependencies installed (you just can't call the functions).
"""
import io
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from joblib import Parallel, delayed

# Hyperparameters that must match _fit_one in 07_fit_popeve_gp.py so loaded
# state_dicts reconstruct cleanly.
M_INDUCING = 20
LENGTHSCALE_INIT = 0.2

POPEVE_REPO = Path(__file__).resolve().parent / 'popEVE'


def log(m):
    print(f'[{time.strftime("%H:%M:%S")}] {m}', flush=True)


def scores_to_pct(x, sorted_train_scores):
    return np.searchsorted(sorted_train_scores, x, side='right') / len(sorted_train_scores)


def _init_worker():
    """Set single-threaded torch in every worker process."""
    sys.path.insert(0, str(POPEVE_REPO))
    import torch
    torch.set_num_threads(1)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass


def _build_model(inducing_points):
    """Constructor used both at fit time (scratch) and predict time (from state).
    Kept identical so state_dict loads cleanly."""
    from popEVE.popEVE import GPModel, PGLikelihood  # noqa
    model = GPModel(inducing_points=inducing_points)
    model.covar_module.base_kernel.initialize(lengthscale=LENGTHSCALE_INIT)
    likelihood = PGLikelihood()
    return model, likelihood


def _predict_chunk(chunk, seed=42):
    """Score a chunk of (fit_dict, scores_array) pairs. Returns a list of
    per-tx prediction dicts in chunk order."""
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
            gp_mean = out.mean.numpy().astype(np.float32)
            mean_prob = torch.sigmoid(out.mean).numpy().astype(np.float32)

        results.append({
            'transcript': fit['transcript'],
            'percentile': p.astype(np.float32),
            'gp_mean': gp_mean,
            'mean_prob': mean_prob,
        })
    return results


def apply_curves(scores_df, score_col, fits, out_path, label, n_workers):
    """Apply pickled GP fits to scores_df, write per-variant parquet."""
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

    records = []
    for g in groups:
        tx = g['transcript'][0]
        fit = fits.get(tx)
        if fit is None:
            continue
        records.append((fit, g[score_col].to_numpy().astype(np.float64),
                        g.select(['locus_str', 'alleles_str', 'transcript'])))
    log(f'    [{label}] built {len(records)} predict records')

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
