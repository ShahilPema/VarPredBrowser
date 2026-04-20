"""Monotone P-spline fitters — Bernoulli replaces Poisson.

The per-site "expected" in this project is P(observe >= 1 variant | cohort),
not a Poisson rate, and observed is strictly {0,1}. The correct likelihood is
Bernoulli:

    O_i ~ Bernoulli(p_i)
    logit(p_i) = logit(E_i) + g(score_i)
    g(score) = B(score) @ theta

with theta penalized by lam * ||D2 theta||^2 and constrained monotone
non-increasing (Delta theta <= 0) via quadprog. logit(E_i) enters as a
fixed per-site offset and is NOT absorbed into theta.

Shared machinery (basis, knots, penalty, monotonicity, quadprog solver) is
identical to analyses/poisson_quadprog_4models/scripts/fitters.py; only the
IRLS inner-step working response/weights and the deviance change.

fit_poisson_irls_qp and fit_qp_quadprog are retained unchanged so this file
can be used for side-by-side comparisons.
"""
import time
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import minimize as sp_minimize, LinearConstraint
from scipy.special import expit  # logistic

try:
    import quadprog  # noqa
    HAVE_QUADPROG = True
except ImportError:
    HAVE_QUADPROG = False


N_KNOTS = 24      # grid_search_phase2 best row
LAMBDA = 500.0    # Poisson default (Bernoulli default TBD by sweep)
ORDER = 3

P_CLIP = 1e-10    # clip for p and 1-p before forming weights / working response


# =============================================================================
# Shared basis construction
# =============================================================================

def scores_to_percentiles(x, sorted_train_scores):
    return np.searchsorted(sorted_train_scores, x, side='right') / len(sorted_train_scores)


def compute_sorted_universe(scores_df, score_col, transcript_col='transcript'):
    """Per-transcript sorted array of all non-null scores in scores_df.

    This is the percentile axis used both at fit time and at predict time:
    every variant that has a score for this transcript gets its true rank,
    independent of whether it survived the per-site training filter.

    Returns {transcript: np.ndarray} (float64, sorted ascending). Transcripts
    with fewer than 2 unique scores are omitted (degenerate basis).
    """
    import polars as pl  # local import to keep fitters.py import-light
    scored = (scores_df
              .select([transcript_col, score_col])
              .filter(pl.col(score_col).is_not_null()))
    agg = (scored.group_by(transcript_col)
                 .agg(pl.col(score_col).sort().alias('sorted_scores')))
    out = {}
    for row in agg.iter_rows(named=True):
        arr = np.asarray(row['sorted_scores'], dtype=np.float64)
        if len(arr) >= 2:
            out[row[transcript_col]] = arr
    return out


def build_basis(x, n_knots=N_KNOTS, order=ORDER, sorted_universe=None):
    """Build monotone P-spline basis in gene-level percentile space.

    sorted_universe: if provided, a sorted 1D array of all scores in the
    gene's full universe (not just training rows). Percentiles are then ranked
    against it, so out-of-training scores get their true full-universe rank at
    apply time. Interior knots are placed at quantiles of the training-set
    percentile values (not evenly spaced on [p_min, p_max]) so filter-induced
    training-mass skew does not leave knot intervals starved.

    When sorted_universe is None, falls back to legacy behavior: percentiles
    are ranked against x itself (training-set only) with evenly spaced knots.
    This path is retained for back-compat with saved curves but should not be
    used for new production fits.
    """
    if sorted_universe is not None:
        sorted_scores = np.asarray(sorted_universe, dtype=np.float64)
        p = scores_to_percentiles(x, sorted_scores)
        p_min, p_max = float(p.min()), float(p.max())
        if p_max - p_min < 1e-10:
            return None
        inner = np.quantile(p, np.linspace(0.0, 1.0, n_knots + 2)[1:-1])
        inner = np.unique(inner)
        if len(inner) < 3:
            # pathological ties — fall back to evenly spaced in observed range
            inner = np.linspace(p_min, p_max, n_knots + 2)[1:-1]
    else:
        sorted_scores = np.sort(x)
        p = scores_to_percentiles(x, sorted_scores)
        p_min, p_max = float(p.min()), float(p.max())
        if p_max - p_min < 1e-10:
            return None
        inner = np.linspace(p_min, p_max, n_knots + 2)[1:-1]

    knots = np.concatenate([
        np.repeat(p_min, order + 1), inner, np.repeat(p_max, order + 1),
    ])
    n_bases = len(knots) - order - 1

    B = np.zeros((len(p), n_bases))
    for i in range(n_bases):
        c = np.zeros(n_bases)
        c[i] = 1.0
        vals = BSpline(knots, c, order, extrapolate=False)(p)
        vals[np.isnan(vals)] = 0.0
        B[:, i] = vals
    return {
        'B': B, 'knots': knots, 'order': order,
        'n_bases': n_bases, 'p': p, 'p_min': p_min, 'p_max': p_max,
        'sorted_train_scores': sorted_scores,
    }


def eval_basis(model, x):
    p = scores_to_percentiles(x, model['sorted_train_scores'])
    p = np.clip(p, model['p_min'], model['p_max'])
    n_bases = model['n_bases']
    B = np.zeros((len(p), n_bases))
    for i in range(n_bases):
        c = np.zeros(n_bases)
        c[i] = 1.0
        vals = BSpline(model['knots'], c, model['order'], extrapolate=False)(p)
        vals[np.isnan(vals)] = 0.0
        B[:, i] = vals
    return B


def mono_constraint_matrix(n_bases):
    A = np.zeros((n_bases - 1, n_bases))
    for i in range(n_bases - 1):
        A[i, i] = 1.0
        A[i, i + 1] = -1.0
    return A


def d2_matrix(n_bases):
    return np.diff(np.eye(n_bases), n=2, axis=0)


def d_matrix(n_bases, order):
    return np.diff(np.eye(n_bases), n=order, axis=0)


# =============================================================================
# quadprog inner solver
# =============================================================================

def _solve_qp_quadprog(G, a, A_mono):
    G = 0.5 * (G + G.T)
    jitter = 1e-10 * np.trace(G) / G.shape[0]
    for _ in range(8):
        try:
            sol, *_ = quadprog.solve_qp(G, a, A_mono.T, np.zeros(A_mono.shape[0]), 0)
            return sol
        except ValueError:
            G = G + jitter * np.eye(G.shape[0])
            jitter *= 10
    raise RuntimeError('quadprog failed after jittering')


# =============================================================================
# Fitter: Gaussian QP on windowed O/E (unchanged, kept for reference)
# =============================================================================

def fit_qp_quadprog(scores, y, weights, lam=LAMBDA, n_knots=N_KNOTS,
                    sorted_universe=None):
    if not HAVE_QUADPROG:
        raise RuntimeError('quadprog not installed')

    basis = build_basis(scores, n_knots=n_knots, sorted_universe=sorted_universe)
    if basis is None:
        return None
    B = basis['B']
    n_bases = basis['n_bases']
    w = weights / weights.sum()
    D2 = d2_matrix(n_bases)
    A_mono = mono_constraint_matrix(n_bases)

    WB = B * w[:, None]
    G = 2.0 * (B.T @ WB + lam * D2.T @ D2)
    a = 2.0 * (B.T @ (w * y))

    t0 = time.perf_counter()
    theta = _solve_qp_quadprog(G, a, A_mono)
    elapsed = time.perf_counter() - t0

    return {
        'method': 'qp_quadprog',
        'coeffs': theta,
        'knots': basis['knots'], 'order': basis['order'],
        'n_bases': n_bases, 'p_min': basis['p_min'], 'p_max': basis['p_max'],
        'sorted_train_scores': basis['sorted_train_scores'],
        'output_scale': 'oe',
        'fit_time': elapsed,
        'n_iter': 1,
        'converged': True,
    }


# =============================================================================
# Fitter: Poisson IRLS-QP (unchanged — retained for side-by-side comparison)
# =============================================================================

def fit_poisson_irls_qp(scores, observed, expected, lam=LAMBDA, n_knots=N_KNOTS,
                        max_iter=50, tol=1e-7, monotone=True, penalty_order=2,
                        sorted_universe=None):
    """log(mu_i) = log(E_i) + B theta, PIRLS with quadprog inner step."""
    if not HAVE_QUADPROG:
        raise RuntimeError('quadprog not installed')

    basis = build_basis(scores, n_knots=n_knots, sorted_universe=sorted_universe)
    if basis is None:
        return None
    B = basis['B']
    n_bases = basis['n_bases']
    D2 = d_matrix(n_bases, penalty_order)
    A_mono = mono_constraint_matrix(n_bases)

    log_E = np.log(np.maximum(expected, 1e-12))
    O = observed.astype(np.float64)

    global_oe = max(O.sum() / max(expected.sum(), 1e-12), 1e-6)
    theta = np.full(n_bases, np.log(global_oe))

    t0 = time.perf_counter()
    prev_dev = np.inf
    n_iter = 0
    converged = False
    for it in range(max_iter):
        eta = B @ theta
        mu = np.exp(eta + log_E)
        mu = np.maximum(mu, 1e-12)

        z = eta + (O - mu) / mu
        W = mu

        WB = B * W[:, None]
        G = 2.0 * (B.T @ WB + lam * D2.T @ D2)
        a = 2.0 * (B.T @ (W * z))

        if monotone:
            theta_new = _solve_qp_quadprog(G, a, A_mono)
        else:
            try:
                theta_new = np.linalg.solve(0.5 * G, 0.5 * a)
            except np.linalg.LinAlgError:
                theta_new = theta

        mu_new = np.exp(B @ theta_new + log_E)
        mu_new = np.maximum(mu_new, 1e-12)
        with np.errstate(divide='ignore', invalid='ignore'):
            term = np.where(O > 0, O * np.log(O / mu_new), 0.0)
        dev = 2.0 * np.sum(term - (O - mu_new))

        theta = theta_new
        n_iter = it + 1
        if np.isfinite(prev_dev) and abs(prev_dev - dev) < tol * max(1.0, abs(prev_dev)):
            converged = True
            break
        prev_dev = dev
    elapsed = time.perf_counter() - t0

    return {
        'method': 'poisson_irls_qp' if monotone else 'poisson_irls_unconstrained',
        'coeffs': theta,
        'knots': basis['knots'], 'order': basis['order'],
        'n_bases': n_bases, 'p_min': basis['p_min'], 'p_max': basis['p_max'],
        'sorted_train_scores': basis['sorted_train_scores'],
        'output_scale': 'log_oe',   # predict returns exp(B theta)
        'fit_time': elapsed,
        'n_iter': n_iter,
        'converged': converged,
        'final_deviance': float(prev_dev),
    }


# =============================================================================
# Fitter: Bernoulli IRLS-QP with logit(E) offset  (NEW — the point of this dir)
# =============================================================================

def fit_bernoulli_irls_qp(scores, observed, expected, lam=LAMBDA, n_knots=N_KNOTS,
                          max_iter=50, tol=1e-7, monotone=True, penalty_order=2,
                          sorted_universe=None):
    """logit(p_i) = logit(E_i) + B theta.

    PIRLS with Bernoulli working response, quadprog inner step.
    Basis, knots, D2 penalty, monotonicity, quadprog call — all identical to
    fit_poisson_irls_qp. Only the working weights, working response, and
    deviance change.

    Returns a model dict with output_scale='bernoulli_logit_oe', meaning the
    stored coeffs parameterize g(score) = B theta on the log-odds scale (a
    deviation from neutral). To recover per-site O/E = p/E, callers must
    supply `expected` at predict time via predict_bernoulli_oe().
    """
    if not HAVE_QUADPROG:
        raise RuntimeError('quadprog not installed')

    basis = build_basis(scores, n_knots=n_knots, sorted_universe=sorted_universe)
    if basis is None:
        return None
    B = basis['B']
    n_bases = basis['n_bases']
    D2 = d_matrix(n_bases, penalty_order)
    A_mono = mono_constraint_matrix(n_bases)

    E = np.clip(expected.astype(np.float64), P_CLIP, 1.0 - P_CLIP)
    O = observed.astype(np.float64)
    logit_E = np.log(E) - np.log1p(-E)   # log(E/(1-E))

    theta = np.zeros(n_bases)   # init at neutral: p = E

    t0 = time.perf_counter()
    prev_dev = np.inf
    n_iter = 0
    converged = False
    for it in range(max_iter):
        g = B @ theta                # deviation on logit scale
        eta = logit_E + g            # full linear predictor
        p = expit(eta)
        p = np.clip(p, P_CLIP, 1.0 - P_CLIP)

        W = p * (1.0 - p)             # working weights
        # working response, full scale, then subtract offset so QP solves for g
        z_full = eta + (O - p) / W
        z = z_full - logit_E          # target for B theta

        WB = B * W[:, None]
        G = 2.0 * (B.T @ WB + lam * D2.T @ D2)
        a = 2.0 * (B.T @ (W * z))

        if monotone:
            theta_new = _solve_qp_quadprog(G, a, A_mono)
        else:
            try:
                theta_new = np.linalg.solve(0.5 * G, 0.5 * a)
            except np.linalg.LinAlgError:
                theta_new = theta

        # Bernoulli deviance at theta_new
        p_new = expit(logit_E + B @ theta_new)
        p_new = np.clip(p_new, P_CLIP, 1.0 - P_CLIP)
        with np.errstate(divide='ignore', invalid='ignore'):
            term1 = np.where(O > 0, O * np.log(O / p_new), 0.0)
            term0 = np.where(O < 1, (1.0 - O) * np.log((1.0 - O) / (1.0 - p_new)), 0.0)
        dev = 2.0 * float(np.sum(term1 + term0))

        theta = theta_new
        n_iter = it + 1
        if np.isfinite(prev_dev) and abs(prev_dev - dev) < tol * max(1.0, abs(prev_dev)):
            converged = True
            break
        prev_dev = dev
    elapsed = time.perf_counter() - t0

    # Diagnostic: did any p hit the clip boundary at the final theta?
    p_final = expit(logit_E + B @ theta)
    p_hit_hi = int(np.sum(p_final >= 1.0 - 1e-6))
    p_hit_lo = int(np.sum(p_final <= 1e-6))

    return {
        'method': 'bernoulli_irls_qp' if monotone else 'bernoulli_irls_unconstrained',
        'coeffs': theta,
        'knots': basis['knots'], 'order': basis['order'],
        'n_bases': n_bases, 'p_min': basis['p_min'], 'p_max': basis['p_max'],
        'sorted_train_scores': basis['sorted_train_scores'],
        'output_scale': 'bernoulli_logit_oe',
        'fit_time': elapsed,
        'n_iter': n_iter,
        'converged': converged,
        'final_deviance': float(prev_dev),
        'p_hit_hi_boundary': p_hit_hi,
        'p_hit_lo_boundary': p_hit_lo,
        'n_sites': int(len(O)),
    }


# =============================================================================
# Fitter: Bernoulli IRLS-QP with cloglog link and cloglog(E) offset  (NEW)
# =============================================================================
#
# Motivation: the data-generating process IS Poisson-thresholded-to-binary:
#   E_i = 1 - exp(-lambda_neutral_i)
# cloglog is the canonical link for this process:
#   cloglog(p_i) = log(-log(1 - p_i)) = cloglog(E_i) + g(score_i)
# Equivalently on the rate scale:
#   lambda_fitted_i = lambda_neutral_i * exp(g(score_i))
# So rate-scale O/E = exp(g), same ranker as the Poisson fitter — but the
# probability stays bounded in [0, 1] and the offset is the correct neutral
# rate -log(1-E) rather than the probability E itself.
#
# Parameterized via lambda = exp(eta) throughout so all derivative
# expressions are unambiguously positive and avoid log(1-p) cancellation.

def fit_bernoulli_cloglog_irls_qp(scores, observed, expected, lam=LAMBDA,
                                  n_knots=N_KNOTS, max_iter=50, tol=1e-7,
                                  monotone=True, penalty_order=2,
                                  eta_clip=30.0, sorted_universe=None):
    """Monotone Bernoulli P-spline with cloglog link and cloglog(E) offset.

    IRLS internals:
        eta     = cloglog(E) + B theta = log(-log(1-E)) + g
        lambda  = exp(eta)                  # Poisson rate, > 0
        p       = 1 - exp(-lambda)          # cloglog inverse
        dp_deta = (1-p) * lambda            # always > 0
        W       = dp_deta^2 / (p*(1-p)) = exp(-lambda)*lambda^2 / p
        z_g     = g + (O - p) / dp_deta     # offset-removed working response
        deviance: 2 sum[ O log(O/p) + (1-O) log((1-O)/(1-p)) ]

    output_scale='log_oe' so predict() returns exp(g), identical semantics to
    the Poisson fitter's ranker.
    """
    if not HAVE_QUADPROG:
        raise RuntimeError('quadprog not installed')

    basis = build_basis(scores, n_knots=n_knots, sorted_universe=sorted_universe)
    if basis is None:
        return None
    B = basis['B']
    n_bases = basis['n_bases']
    D2 = d_matrix(n_bases, penalty_order)
    A_mono = mono_constraint_matrix(n_bases)

    E = np.clip(expected.astype(np.float64), P_CLIP, 1.0 - P_CLIP)
    O = observed.astype(np.float64)
    lam_neutral = -np.log1p(-E)           # > 0, = -log(1-E)
    cloglog_E = np.log(lam_neutral)       # fixed offset

    theta = np.zeros(n_bases)             # g=0 -> p=E (neutral init)

    def _state(theta_):
        g_ = B @ theta_
        eta_ = np.clip(cloglog_E + g_, -eta_clip, eta_clip)
        lam_fit = np.exp(eta_)
        one_minus_p = np.clip(np.exp(-lam_fit), P_CLIP, 1.0)
        p_ = np.clip(1.0 - one_minus_p, P_CLIP, 1.0 - P_CLIP)
        return g_, lam_fit, p_, one_minus_p

    t0 = time.perf_counter()
    prev_dev = np.inf
    n_iter = 0
    converged = False
    for it in range(max_iter):
        g_cur, lam_fit, p, one_minus_p = _state(theta)

        dp_deta = one_minus_p * lam_fit                    # >0
        W = (dp_deta * dp_deta) / (p * one_minus_p)         # >0
        z = g_cur + (O - p) / dp_deta                       # g-scale target

        WB = B * W[:, None]
        G = 2.0 * (B.T @ WB + lam * D2.T @ D2)
        a = 2.0 * (B.T @ (W * z))

        if monotone:
            theta_new = _solve_qp_quadprog(G, a, A_mono)
        else:
            try:
                theta_new = np.linalg.solve(0.5 * G, 0.5 * a)
            except np.linalg.LinAlgError:
                theta_new = theta

        _, _, p_new, _ = _state(theta_new)
        with np.errstate(divide='ignore', invalid='ignore'):
            term1 = np.where(O > 0, O * np.log(O / p_new), 0.0)
            term0 = np.where(O < 1, (1.0 - O) * np.log((1.0 - O) / (1.0 - p_new)), 0.0)
        dev = 2.0 * float(np.sum(term1 + term0))

        theta = theta_new
        n_iter = it + 1
        if np.isfinite(prev_dev) and abs(prev_dev - dev) < tol * max(1.0, abs(prev_dev)):
            converged = True
            break
        prev_dev = dev
    elapsed = time.perf_counter() - t0

    _, _, p_final, _ = _state(theta)
    p_hit_hi = int(np.sum(p_final >= 1.0 - 1e-6))
    p_hit_lo = int(np.sum(p_final <= 1e-6))

    return {
        'method': 'bernoulli_cloglog_irls_qp' if monotone else 'bernoulli_cloglog_unconstrained',
        'coeffs': theta,
        'knots': basis['knots'], 'order': basis['order'],
        'n_bases': n_bases, 'p_min': basis['p_min'], 'p_max': basis['p_max'],
        'sorted_train_scores': basis['sorted_train_scores'],
        'output_scale': 'log_oe',   # predict() -> exp(g) — matches Poisson
        'fit_time': elapsed,
        'n_iter': n_iter,
        'converged': converged,
        'final_deviance': float(prev_dev),
        'p_hit_hi_boundary': p_hit_hi,
        'p_hit_lo_boundary': p_hit_lo,
        'n_sites': int(len(O)),
    }


def predict_cloglog_p(model, x, expected):
    """Exact cloglog per-site probability."""
    g = _basis_eval_fast(model, x)
    E = np.clip(np.asarray(expected, dtype=np.float64), P_CLIP, 1.0 - P_CLIP)
    lam_neutral = -np.log1p(-E)
    lam_fit = lam_neutral * np.exp(g)
    return -np.expm1(-lam_fit)


# =============================================================================
# Post-hoc: recover converged IRLS working weights for diagnostic plots
# =============================================================================

def weights_at_convergence(method, model, scores, observed, expected):
    """Return the IRLS working weight vector W_i at the fitted theta.

    method in {'poisson', 'bernoulli_logit', 'bernoulli_cloglog'}.
    Recomputes W using the exact formulas inside each fitter's iteration.
    """
    s = np.asarray(scores, dtype=np.float64)
    E = np.asarray(expected, dtype=np.float64)

    B = eval_basis(model, s)
    g = B @ model['coeffs']

    if method == 'poisson':
        log_E = np.log(np.maximum(E, 1e-12))
        mu = np.exp(log_E + g)
        return np.maximum(mu, 1e-12)

    if method == 'bernoulli_logit':
        Ec = np.clip(E, P_CLIP, 1.0 - P_CLIP)
        logit_E = np.log(Ec) - np.log1p(-Ec)
        p = expit(logit_E + g)
        p = np.clip(p, P_CLIP, 1.0 - P_CLIP)
        return p * (1.0 - p)

    if method == 'bernoulli_cloglog':
        Ec = np.clip(E, P_CLIP, 1.0 - P_CLIP)
        cloglog_E = np.log(-np.log1p(-Ec))
        eta = np.clip(cloglog_E + g, -30.0, 30.0)
        lam_fit = np.exp(eta)
        one_minus_p = np.clip(np.exp(-lam_fit), P_CLIP, 1.0)
        p = np.clip(1.0 - one_minus_p, P_CLIP, 1.0 - P_CLIP)
        dp_deta = one_minus_p * lam_fit
        return (dp_deta * dp_deta) / (p * one_minus_p)

    raise ValueError(f'unknown method: {method}')


# =============================================================================
# Prediction
# =============================================================================

def _basis_eval_fast(model, x):
    """BSpline-based evaluation of g(score) = B theta (single spline call)."""
    x = np.asarray(x, dtype=np.float64)
    p = scores_to_percentiles(x, model['sorted_train_scores'])
    p = np.clip(p, model['p_min'], model['p_max'])
    spline = BSpline(model['knots'], model['coeffs'], model['order'],
                     extrapolate=False)
    y = spline(p)
    np.nan_to_num(y, copy=False, nan=0.0)
    return y


def predict(model, x):
    """Prediction. Semantics vary by output_scale:

    - 'log_oe' or 'bernoulli_logit_oe': returns exp(B theta)  (rate O/E)
    - 'logit_p': returns sigmoid(B theta)  (fitted probability p, for the
      no-offset PLM-alone ablation). Note this is NOT an O/E — it is the
      raw fitted P(observed | score). Rank ascending for pathogenicity.
    - 'oe' (Gaussian): returns B theta clipped >=0.
    """
    y = _basis_eval_fast(model, x)
    scale = model.get('output_scale')
    if scale in ('log_oe', 'bernoulli_logit_oe'):
        return np.exp(np.clip(y, -30.0, 30.0))
    if scale == 'logit_p':
        return expit(np.clip(y, -30.0, 30.0))
    return np.maximum(y, 0.0)


# =============================================================================
# Fitter: Bernoulli-logit with NO offset (PLM-alone ablation, Model A)
# =============================================================================
#
# Purpose: fit P(observed | score) using only the PLM score — no mutation
# rate offset, no E. g(score) absorbs everything that correlates with both
# the score and variant observation (selection + mutation rate + sequence
# context). Used to quantify the marginal value of the mutation rate offset
# when compared against fit_bernoulli_cloglog_irls_qp.
#
#     eta_i = B theta                    (no offset)
#     p_i   = sigmoid(eta_i)
#     W_i   = p_i (1 - p_i)
#     z_i   = eta_i + (O_i - p_i) / W_i
#
# Shared machinery (basis, knots, D2 penalty, monotonicity, quadprog QP) is
# unchanged from all other fitters in this file. Init theta = 0 puts p = 0.5
# everywhere; first iteration moves to the global prevalence level.

def fit_bernoulli_no_offset_irls_qp(scores, observed, lam=LAMBDA,
                                    n_knots=N_KNOTS, max_iter=60, tol=1e-7,
                                    monotone=True, penalty_order=2,
                                    sorted_universe=None):
    """Monotone Bernoulli-logit P-spline without offset."""
    if not HAVE_QUADPROG:
        raise RuntimeError('quadprog not installed')

    basis = build_basis(scores, n_knots=n_knots, sorted_universe=sorted_universe)
    if basis is None:
        return None
    B = basis['B']
    n_bases = basis['n_bases']
    D2 = d_matrix(n_bases, penalty_order)
    A_mono = mono_constraint_matrix(n_bases)

    O = observed.astype(np.float64)
    theta = np.zeros(n_bases)     # p = 0.5 at init

    t0 = time.perf_counter()
    prev_dev = np.inf
    n_iter = 0
    converged = False
    for it in range(max_iter):
        eta = B @ theta
        p = np.clip(expit(eta), P_CLIP, 1.0 - P_CLIP)
        W = p * (1.0 - p)
        z = eta + (O - p) / W

        WB = B * W[:, None]
        G = 2.0 * (B.T @ WB + lam * D2.T @ D2)
        a = 2.0 * (B.T @ (W * z))

        if monotone:
            theta_new = _solve_qp_quadprog(G, a, A_mono)
        else:
            try:
                theta_new = np.linalg.solve(0.5 * G, 0.5 * a)
            except np.linalg.LinAlgError:
                theta_new = theta

        p_new = np.clip(expit(B @ theta_new), P_CLIP, 1.0 - P_CLIP)
        with np.errstate(divide='ignore', invalid='ignore'):
            term1 = np.where(O > 0, O * np.log(O / p_new), 0.0)
            term0 = np.where(O < 1, (1.0 - O) * np.log((1.0 - O) / (1.0 - p_new)), 0.0)
        dev = 2.0 * float(np.sum(term1 + term0))

        theta = theta_new
        n_iter = it + 1
        if np.isfinite(prev_dev) and abs(prev_dev - dev) < tol * max(1.0, abs(prev_dev)):
            converged = True
            break
        prev_dev = dev
    elapsed = time.perf_counter() - t0

    p_final = np.clip(expit(B @ theta), P_CLIP, 1.0 - P_CLIP)
    p_hit_hi = int(np.sum(p_final >= 1.0 - 1e-6))
    p_hit_lo = int(np.sum(p_final <= 1e-6))

    return {
        'method': 'bernoulli_no_offset_irls_qp' if monotone else 'bernoulli_no_offset_unconstrained',
        'coeffs': theta,
        'knots': basis['knots'], 'order': basis['order'],
        'n_bases': n_bases, 'p_min': basis['p_min'], 'p_max': basis['p_max'],
        'sorted_train_scores': basis['sorted_train_scores'],
        'output_scale': 'logit_p',   # predict() -> sigmoid(B theta)
        'fit_time': elapsed,
        'n_iter': n_iter,
        'converged': converged,
        'final_deviance': float(prev_dev),
        'p_hit_hi_boundary': p_hit_hi,
        'p_hit_lo_boundary': p_hit_lo,
        'n_sites': int(len(O)),
    }


def predict_bernoulli_oe(model, x, expected):
    """Exact Bernoulli per-site O/E:

        p = sigmoid(logit(E) + B theta)
        O/E = p / E

    This is the function to use when applying a Bernoulli curve to new sites.
    """
    g = _basis_eval_fast(model, x)
    E = np.clip(np.asarray(expected, dtype=np.float64), P_CLIP, 1.0 - P_CLIP)
    logit_E = np.log(E) - np.log1p(-E)
    p = expit(logit_E + g)
    return p / E


def predict_bernoulli_p(model, x, expected):
    g = _basis_eval_fast(model, x)
    E = np.clip(np.asarray(expected, dtype=np.float64), P_CLIP, 1.0 - P_CLIP)
    logit_E = np.log(E) - np.log1p(-E)
    return expit(logit_E + g)
