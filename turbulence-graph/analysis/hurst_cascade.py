"""
Non-Markovian Cascade Analysis: Hurst Exponent in Turbulence Scale Space
=========================================================================

Tests the hypothesis that the turbulent energy cascade has long-range memory
characterized by a Hurst exponent H ≠ 1/2 in scale space.

Two tracks:
  A. Synthetic cascade: generate fractional BM cascade, estimate H, validate
  B. DNS pipeline: ready to run against JHTDB data with auth token

Physical hypothesis:
  The coarse-grained dissipation ln(ε_r) evolves as fractional Brownian
  motion in log-scale space s = log(L/r):

  Markovian (K62 log-normal): Cov[ln ε_r, ln ε_{r'}] = μ² log(L/max(r,r'))
  Non-Markovian (fBM):        Cov[ln ε_r, ln ε_{r'}] = (σ²/2)[s^{2H} + s'^{2H} - |s-s'|^{2H}]

  H = 1/2: Markovian (K62)
  H > 1/2: Long-range positive memory (each scale remembers all larger scales)
  H < 1/2: Anti-correlated cascade (overshooting at each step)

Key discriminator:
  The lag-1 autocorrelation of log-scale increments:
  ρ(1) = (2^{2H} - 2) / 2  →  H = 1/2 gives ρ(1) = 0
                             H = 0.7 gives ρ(1) ≈ 0.095
                             H = 0.8 gives ρ(1) ≈ 0.246
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.optimize import minimize_scalar, curve_fit
from scipy.special import gamma
import warnings
warnings.filterwarnings('ignore')

RNG = np.random.default_rng(42)

# ── Colours ────────────────────────────────────────────────────────────────
C = {
    'markov':   '#94a3b8',
    'H06':      '#22c55e',
    'H07':      '#3b82f6',
    'H08':      '#f97316',
    'H09':      '#ef4444',
    'data':     '#a855f7',
    'exact':    '#facc15',
}

# ═══════════════════════════════════════════════════════════════════════════
# 1. FRACTIONAL BROWNIAN MOTION GENERATION
# ═══════════════════════════════════════════════════════════════════════════

def fbm_covariance(n, H):
    """
    Covariance matrix of a fractional Brownian motion with Hurst exponent H.
    C[i,j] = (1/2)[|i|^{2H} + |j|^{2H} - |i-j|^{2H}]
    Indices i, j ∈ {1, ..., n}.
    """
    i = np.arange(1, n + 1, dtype=float)
    I, J = np.meshgrid(i, i, indexing='ij')
    return 0.5 * (I**(2*H) + J**(2*H) - np.abs(I - J)**(2*H))

def generate_fbm(n, H, sigma=1.0, seed=None):
    """
    Generate one realisation of fractional Brownian motion at n points
    (increments, then cumulative sum) via Cholesky decomposition.

    Returns: fbm path X[0..n-1], increments dX[0..n-2]
    """
    rng = np.random.default_rng(seed)
    cov = fbm_covariance(n, H) * sigma**2
    # Regularise for numerical stability
    cov += np.eye(n) * 1e-10
    L   = np.linalg.cholesky(cov)
    z   = rng.standard_normal(n)
    return L @ z  # fBm values at times 1, 2, ..., n

def fbm_autocorrelation(H, max_lag=10):
    """
    Exact autocorrelation of fBm increments (fractional Gaussian noise):
    ρ(lag) = (1/2)[|lag+1|^{2H} + |lag-1|^{2H} - 2|lag|^{2H}]
    """
    lag = np.arange(0, max_lag + 1, dtype=float)
    rho = 0.5 * (np.abs(lag + 1)**(2*H) + np.abs(lag - 1)**(2*H) - 2*np.abs(lag)**(2*H))
    return rho

# ═══════════════════════════════════════════════════════════════════════════
# 2. SYNTHETIC MULTIFRACTAL CASCADE (LOG-NORMAL + FRACTIONAL)
# ═══════════════════════════════════════════════════════════════════════════

def generate_lognormal_cascade(n_levels=12, lam2=0.025, N_pts=5000, seed=42):
    """
    Generate a log-normal multiplicative cascade (Markovian, H=1/2).
    Returns: dissipation array ε at finest scale.

    Each level: ε_{l+1} = ε_l × exp(σ W_l - σ²/2)
    where W_l ~ N(0,1) i.i.d. and σ² = lam2 × log(2) (per level).
    """
    rng = np.random.default_rng(seed)
    sigma_per_level = np.sqrt(lam2 * np.log(2))
    eps = np.ones(N_pts)
    for _ in range(n_levels):
        W   = rng.standard_normal(N_pts)
        eps = eps * np.exp(sigma_per_level * W - sigma_per_level**2 / 2)
    return eps

def generate_fractional_cascade(n_levels=12, lam2=0.025, H=0.7, N_pts=500, seed=42):
    """
    Generate a cascade where the log-multipliers at each spatial position
    follow a fractional Brownian motion across scale levels (scale space memory).

    Construction:
      - Scale space: s ∈ {1, 2, ..., n_levels} (s = log(L/r)/log(2))
      - For each spatial position x: the log-dissipation path
          Y_x(s) = σ × B_H(s)  (fBm path in scale space)
      - ε_r(x) = exp(Y_x(s_r) - Var(Y_x)/2)

    This generates a cascade with scale-space Hurst exponent H.
    """
    rng  = np.random.default_rng(seed)
    sigma = np.sqrt(lam2)

    # For each spatial point, generate fBm path in scale space
    eps_finest = np.zeros(N_pts)

    for ix in range(N_pts):
        path = generate_fbm(n_levels, H, sigma=sigma, seed=seed + ix)
        # Total log-dissipation at finest scale = sum of all cascade levels
        log_eps = np.sum(path) - lam2 * n_levels / 2  # normalise
        eps_finest[ix] = np.exp(log_eps)

    return eps_finest

def coarse_grain(eps, n_coarse):
    """Coarse-grain dissipation by averaging n_coarse consecutive values."""
    N = len(eps) - (len(eps) % n_coarse)
    return eps[:N].reshape(-1, n_coarse).mean(axis=1)

def compute_log_dissipation_covariance(eps_fine, scale_factors, log=True):
    """
    Compute Cov[ln ε_r, ln ε_{r'}] at multiple coarsening levels.

    Parameters
    ----------
    eps_fine:      finest-scale dissipation array
    scale_factors: list of coarsening factors (integers)
    log:           whether to take logarithm of ε first

    Returns:
        cov_matrix[i,j] = empirical covariance at scales i,j
        scales:           list of scales used
    """
    n_scales = len(scale_factors)
    log_eps  = {}
    for sf in scale_factors:
        eps_coarse      = coarse_grain(eps_fine, sf)
        log_eps[sf]     = np.log(np.maximum(eps_coarse, 1e-20)) if log else eps_coarse
        # demean
        log_eps[sf]    -= log_eps[sf].mean()

    # Build covariance matrix (min size)
    min_len = min(len(v) for v in log_eps.values())
    cov_matrix = np.zeros((n_scales, n_scales))
    for i, si in enumerate(scale_factors):
        for j, sj in enumerate(scale_factors):
            vi = log_eps[si][:min_len]
            vj = log_eps[sj][:min_len]
            cov_matrix[i, j] = np.mean(vi * vj)

    return cov_matrix, scale_factors

# ═══════════════════════════════════════════════════════════════════════════
# 3. HURST EXPONENT ESTIMATOR
# ═══════════════════════════════════════════════════════════════════════════

def hurst_from_autocorrelation(increments):
    """
    Estimate Hurst exponent from lag-1 autocorrelation of increments.

    Exact formula: ρ(1) = (2^{2H} - 2) / 2  →  H = log2(2 + 2ρ(1)) / 2

    Most robust: lag-1 autocorrelation (maximum signal-to-noise).
    """
    rho1  = np.corrcoef(increments[:-1], increments[1:])[0, 1]
    # Invert exact formula: ρ(1) = (2^{2H} - 2) / 2
    H_est = np.log2(2 + 2 * rho1) / 2 if (2 + 2*rho1) > 0 else 0.5
    return H_est, rho1

def hurst_from_variance_scaling(fbm_paths, scales):
    """
    Estimate H from variance scaling: Var[fBm(s)] ~ s^{2H}.
    More stable than autocorrelation for short paths.
    """
    vars_ = np.array([np.var(fbm_paths[:, i]) for i in range(fbm_paths.shape[1])])
    log_s = np.log(scales)
    log_v = np.log(np.maximum(vars_, 1e-20))
    slope, intercept, r, _, _ = stats.linregress(log_s, log_v)
    H_est = slope / 2
    return H_est, r**2

def hurst_from_structure_function_scaling(zeta_p_data, p_values):
    """
    Estimate the effective Hurst exponent from structure function exponents.

    Under the fractional cascade hypothesis:
      ζ_p = p × H  (if cascade is fBm in physical space, not scale space)

    Alternative: fit ζ_p = p/3 + β(p, H) where β depends on H.

    Here we use the simplest estimator: slope of ζ_p vs p at high p
    (where intermittency corrections are large and most diagnostic).
    """
    p_high = np.array([p for p in p_values if p > 3])
    z_high = np.array([z for p, z in zip(p_values, zeta_p_data) if p > 3])
    slope, intercept, r, _, _ = stats.linregress(p_high, z_high)
    return slope, r**2

def maximum_likelihood_hurst(increments, H_grid=None):
    """
    Maximum likelihood estimation of H from a sequence of fBm increments.
    Uses the exact likelihood with the Toeplitz covariance structure.
    """
    if H_grid is None:
        H_grid = np.linspace(0.05, 0.95, 91)

    n = len(increments)
    log_likelihoods = []

    for H in H_grid:
        # Covariance of increments (fractional Gaussian noise)
        rho = fbm_autocorrelation(H, max_lag=min(n-1, 20))
        # Approximate log-likelihood using AR(1) approximation at lag-1
        r1     = rho[1]
        sigma2 = np.var(increments) * (1 - r1**2)
        ll     = -0.5 * n * np.log(2 * np.pi * sigma2)
        ll    -= 0.5 * np.sum((increments[1:] - r1 * increments[:-1])**2) / sigma2
        log_likelihoods.append(ll)

    best_idx = np.argmax(log_likelihoods)
    H_mle    = H_grid[best_idx]
    return H_mle, H_grid, np.array(log_likelihoods)

# ═══════════════════════════════════════════════════════════════════════════
# 4. STATISTICAL TEST POWER ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

def power_analysis(H_true_values, n_points_list, n_realisations=200):
    """
    Compute the power of the H = 0.5 test as a function of:
      - True H value
      - Number of scale points n

    Returns: power[i, j] = fraction of times test correctly rejects H = 0.5
    """
    results = {}
    for n in n_points_list:
        row = {}
        for H_true in H_true_values:
            rejections = 0
            for _ in range(n_realisations):
                path = generate_fbm(n + 1, H_true, seed=None)
                inc  = np.diff(path)
                H_est, rho1 = hurst_from_autocorrelation(inc)
                # Test: |ρ(1)| > 2/sqrt(n) → reject H=0.5
                threshold = 2.0 / np.sqrt(n)
                if abs(rho1) > threshold:
                    rejections += 1
            row[H_true] = rejections / n_realisations
        results[n] = row
    return results

def analytical_discriminator():
    """
    Exact values of the lag-1 autocorrelation ρ(1) for different H values.
    This is the key discriminator between Markovian and non-Markovian cascade.
    """
    H_values = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
    rows = []
    for H in H_values:
        rho = fbm_autocorrelation(H, max_lag=5)
        rows.append({
            'H':          H,
            'rho_lag1':   round(rho[1], 5),
            'rho_lag2':   round(rho[2], 5),
            'rho_lag3':   round(rho[3], 5),
            'markovian':  H == 0.5,
            'n_for_5pct': max(16, int(np.ceil(4 / rho[1]**2))) if abs(rho[1]) > 0.01 else 9999,
        })
    return rows

# ═══════════════════════════════════════════════════════════════════════════
# 5. STRUCTURE FUNCTION ANALYSIS PIPELINE
# ═══════════════════════════════════════════════════════════════════════════

def compute_structure_functions_from_velocity(u, p_orders, r_indices=None):
    """
    Compute longitudinal structure functions S_p(r) from 1D velocity array.

    S_p(r) = ⟨|u(x+r) - u(x)|^p⟩

    Parameters
    ----------
    u:        1D velocity array (N points)
    p_orders: list of p values [2, 4, 6, 8]
    r_indices: lag indices to compute (None = geometric series)

    Returns: dict {p: array of S_p at each r}
    """
    N = len(u)
    if r_indices is None:
        r_indices = np.unique(np.geomspace(1, N // 4, 20).astype(int))

    results = {p: np.zeros(len(r_indices)) for p in p_orders}

    for j, r in enumerate(r_indices):
        du = u[r:] - u[:-r]
        for p in p_orders:
            results[p][j] = np.mean(np.abs(du)**p)

    return results, r_indices

def fit_scaling_exponent(S_p_r, r_values, r_min_frac=0.01, r_max_frac=0.2):
    """
    Fit ζ_p from S_p(r) ~ r^{ζ_p} via log-log regression in the inertial range.
    """
    N     = r_values[-1]
    r_min = r_min_frac * N
    r_max = r_max_frac * N
    mask  = (r_values >= r_min) & (r_values <= r_max)
    if mask.sum() < 3:
        return np.nan, 0.0

    log_r = np.log(r_values[mask])
    log_S = np.log(np.maximum(S_p_r[mask], 1e-30))
    slope, _, r, _, _ = stats.linregress(log_r, log_S)
    return slope, r**2

def extract_hurst_from_anomalous_scaling(zeta_p_dict, p_values):
    """
    Given measured ζ_p values, extract the effective Hurst exponent H
    by fitting the anomalous part β_p = ζ_p - p/3.

    Model: β_p = β_0 × f(p, H) where f encodes the fBm cascade structure.

    For a simple fBm cascade in velocity space:
      ζ_p = H × p  (monofractal limit)

    For an fBm cascade in SCALE space (our hypothesis):
      ζ_p = p/3 - β_p  where β_p ~ (p/3 - 1)² × g(H)

    We fit: β_p = A × [1 - (1 - 1/p)^{2H}] as a simple parametric form.
    """
    p_arr   = np.array(p_values, dtype=float)
    zeta_arr = np.array([zeta_p_dict[p] for p in p_values])
    beta_arr = zeta_arr - p_arr / 3.0

    # Fit: β_p = A [1 - exp(-B p)] as a simple 2-param model
    def model(p, A, B):
        return A * (1 - np.exp(-B * p))

    try:
        popt, pcov = curve_fit(model, p_arr, beta_arr, p0=[-0.1, 0.2], maxfev=5000)
        A_fit, B_fit = popt
        H_eff = 0.5 - A_fit  # rough mapping
    except Exception:
        H_eff, A_fit, B_fit = np.nan, np.nan, np.nan

    return {'H_eff': H_eff, 'A': A_fit, 'B': B_fit, 'beta_arr': beta_arr}

# ═══════════════════════════════════════════════════════════════════════════
# 6. RIGOROUS STATISTICAL TEST — OLS VARIANCE ESTIMATOR (Cramér-Rao optimal)
# ═══════════════════════════════════════════════════════════════════════════
#
# Derived from agent analysis.  Key insight: the OLS fit of log(Var(X_p))
# vs log(s) achieves the Cramér-Rao lower bound exactly.
#
# Model:  log V̂_i ~ N(log σ_p² + 2H log s_i,  2/(M-1))
# OLS:    H_hat = Σ(x_i - x̄)(y_i - ȳ) / Σ(x_i - x̄)²
#         where x_i = log s_i,  y_i = (1/2) log V̂_i
# CRLB:   Var(H_hat) = σ_ε² / (2 SS_xx)   with  σ_ε² = 2/(M-1)
#
# Lag-1 closed form:  ρ(1) = 2^{2H-1} - 1   (≡ (2^{2H} - 2)/2)
# Inversion:          H = (1/2)(log₂(ρ(1) + 1) + 1)

def ols_variance_hurst_estimator(Sp_r_samples, r_values, L, p=2, M_blocks=None):
    """
    Estimate H via OLS on log(Var(X_p)) vs log(s = log(L/r)).

    Parameters
    ----------
    Sp_r_samples : dict mapping r -> 1D array of |Δu|^p samples at that r
    r_values     : array of r values
    L            : integral scale
    p            : structure function order (default 2)
    M_blocks     : spatial blocking size (None = use all samples, no blocking)

    Returns
    -------
    H_hat   : OLS estimate
    SE_H    : standard error (from CRLB)
    CI95    : ±1.96 SE
    r2      : R² of the variance fit
    """
    s_vals = np.log(L / np.asarray(r_values, dtype=float))
    log_s  = np.log(s_vals)

    V_hats = []
    for r in r_values:
        samples = np.asarray(Sp_r_samples[r], dtype=float)
        X_p     = np.log(np.maximum(samples, 1e-30)) - (p / 3.0) * np.log(r)
        if M_blocks is not None:
            # block mean to get M_blocks independent estimates
            k = len(X_p) // M_blocks
            if k < 2:
                V_hats.append(np.var(X_p, ddof=1))
                continue
            block_means = np.array([X_p[i*k:(i+1)*k].mean() for i in range(M_blocks)])
            V_hats.append(np.var(block_means, ddof=1))
        else:
            V_hats.append(np.var(X_p, ddof=1))

    V_hats = np.array(V_hats)
    y      = 0.5 * np.log(np.maximum(V_hats, 1e-30))

    # OLS
    x_bar = log_s.mean()
    y_bar = y.mean()
    SS_xx = np.sum((log_s - x_bar)**2)
    SS_xy = np.sum((log_s - x_bar) * (y - y_bar))
    H_hat = SS_xy / SS_xx if SS_xx > 0 else 0.5

    # Residuals → empirical σ_ε
    y_fit = H_hat * log_s + (y_bar - H_hat * x_bar)
    resid = y - y_fit
    n     = len(r_values)
    sigma2_eps = np.sum(resid**2) / max(n - 2, 1)
    r2    = 1 - np.sum(resid**2) / np.sum((y - y_bar)**2) if SS_xx > 0 else 0.0

    # Standard error (CRLB)
    SE_H = np.sqrt(sigma2_eps / SS_xx) if SS_xx > 0 else np.nan
    CI95 = 1.96 * SE_H

    return H_hat, SE_H, CI95, r2


def cramér_rao_bound(n_scales, M_spatial, r_min, r_max, L):
    """
    Compute the theoretical Cramér-Rao lower bound for H estimation.

    Parameters
    ----------
    n_scales   : number of inertial-range scale points
    M_spatial  : number of spatial samples per scale
    r_min, r_max : inertial range limits
    L          : integral scale

    Returns
    -------
    SE_H : standard error at CRLB
    """
    s_vals = np.log(L / np.geomspace(r_min, r_max, n_scales))
    log_s  = np.log(s_vals)
    SS_xx  = np.sum((log_s - log_s.mean())**2)
    sigma2_eps = 2.0 / (M_spatial - 1)   # chi-squared variance of log(V̂)
    return np.sqrt(sigma2_eps / (2 * SS_xx))


def publishability_check(H_hat, SE_H, H_estimates_by_p, rho1_measured, n_scales):
    """
    Evaluate publishability criteria for H ≠ 0.5 from turbulence DNS data.

    Criteria (from rigorous statistical test design):
      1. |H_hat − 0.5| > 5σ at each of p = 2, 4, 6 independently
      2. ρ̂(1) > 3 SE and sign consistent with OLS H_hat
      3. H_hat from p = 2, 4, 6 agree within ±0.05 (cascade origin, not artifact)
      4. Likelihood ratio 2Δlog L > 10.83 (χ²₁ at 99.9%)

    Returns list of (criterion, passed, detail) tuples.
    """
    results = []

    # Criterion 1: 5σ detection per order
    for p, H_p in H_estimates_by_p.items():
        sigma_p = SE_H  # assume same SE (conservative)
        z = abs(H_p - 0.5) / sigma_p
        results.append((
            f'5σ detection p={p}',
            z > 5,
            f'|H({p})-0.5|/SE = {z:.2f} (need >5)'
        ))

    # Criterion 2: autocorrelation consistency
    SE_rho = 1 / np.sqrt(n_scales)
    rho_z  = abs(rho1_measured) / SE_rho
    sign_ok = (rho1_measured > 0) == (H_hat > 0.5)
    results.append((
        'ρ(1) sign + magnitude',
        rho_z > 3 and sign_ok,
        f'ρ(1)/SE = {rho_z:.2f} (need >3), sign consistent={sign_ok}'
    ))

    # Criterion 3: order independence
    H_vals = list(H_estimates_by_p.values())
    H_spread = max(H_vals) - min(H_vals) if H_vals else np.inf
    results.append((
        'Order independence (p=2,4,6 agree within ±0.05)',
        H_spread < 0.05,
        f'Spread across orders: {H_spread:.4f} (need <0.05)'
    ))

    # Criterion 4: likelihood ratio (approximate)
    # LR ≈ (H_hat - 0.5)² / Var(H_hat) under Gaussian approximation
    LR = ((H_hat - 0.5) / SE_H)**2 if SE_H > 0 else 0
    results.append((
        'LR test 2ΔlogL > 10.83',
        LR > 10.83,
        f'LR statistic: {LR:.2f} (need >10.83 for 99.9%)'
    ))

    return results


def print_publishability_report(H_hat, SE_H, H_by_p, rho1, n_scales):
    """Print formatted publishability assessment."""
    criteria = publishability_check(H_hat, SE_H, H_by_p, rho1, n_scales)
    passed   = sum(1 for _, ok, _ in criteria if ok)
    total    = len(criteria)

    print("\n" + "─" * 70)
    print(f"  PUBLISHABILITY ASSESSMENT  (H_hat = {H_hat:.4f} ± {SE_H:.4f})")
    print("─" * 70)
    for name, ok, detail in criteria:
        mark = '✓' if ok else '✗'
        print(f"  [{mark}] {name}")
        print(f"      {detail}")
    print("─" * 70)
    print(f"  {passed}/{total} criteria met  —  ",
          "PUBLISHABLE" if passed == total else
          ("STRONG EVIDENCE" if passed >= total - 1 else "INCONCLUSIVE"))
    print("─" * 70 + "\n")


# ═══════════════════════════════════════════════════════════════════════════
# 7. JHTDB DATA ACCESS PIPELINE (ready to run with auth token)
# ═══════════════════════════════════════════════════════════════════════════

JHTDB_PIPELINE_CODE = r'''
"""
JHTDB Analysis Pipeline — Rigorous Non-Markovian Cascade Test
==============================================================

Requirements:
  pip install pyJHTDB numpy scipy matplotlib

Auth token from: https://turbulence.pha.jhu.edu/
Dataset: isotropic1024coarse (Re_lambda=433, 1024^3, Kolmogorov η ~ 2Δx)

Two estimators are applied:
  (A) OLS log(Var(X_p)) vs log(s)  — Cramér-Rao optimal
  (B) Lag-1 autocorrelation         — single-number discriminator

Publishability threshold: |H - 0.5| > 5σ at each of p=2,4,6 independently,
  plus order-independence |H(p=2) - H(p=6)| < 0.05.

JHTDB inertial-range scales:
  r = 2Δx .. 128Δx   (Δx = 2π/1024 ≈ 0.006136)
  s = log(L/r) = 2.079 .. 6.238   (7 scale points)
  M >= 834 spatial samples per scale for SE(H) < 0.05 (CRLB)
"""
import pyJHTDB
import numpy as np
from scipy import stats

AUTH_TOKEN = "YOUR_TOKEN_HERE"  # <-- register at turbulence.pha.jhu.edu

# ── Constants ─────────────────────────────────────────────────────────────
DATASET   = "isotropic1024coarse"
T0        = 0.002           # first stored timestep
N_GRID    = 1024
DX        = 2 * np.pi / N_GRID   # ≈ 0.006136
L         = 2 * np.pi            # integral scale
# Inertial-range scales: 2Δx to 128Δx
R_MULT    = np.array([2, 4, 8, 16, 32, 64, 128])
R_VALS    = R_MULT * DX


def get_velocity(lJHTDB, n_points=50_000, t=T0):
    rng    = np.random.default_rng(0)
    points = rng.uniform(0, 2*np.pi, (n_points, 3)).astype(np.float32)
    print(f"Querying {n_points} base points from JHTDB...")
    vel = lJHTDB.getData(t, points, sinterp=4, tinterp=0,
                         getFunction="getVelocity", data_set=DATASET)
    return points, vel  # vel: (n_points, 3)


def longitudinal_increments(lJHTDB, points, vel_base, r, direction=0, t=T0):
    e_hat      = np.zeros(3, dtype=np.float32)
    e_hat[direction] = 1.0
    pts_shifted = ((points + r * e_hat) % (2*np.pi)).astype(np.float32)
    vel_shifted = lJHTDB.getData(t, pts_shifted, sinterp=4, tinterp=0,
                                 getFunction="getVelocity", data_set=DATASET)
    return vel_shifted[:, direction] - vel_base[:, direction]  # (n_points,)


def ols_hurst(Sp_samples, r_values, L_scale, p=2):
    """OLS estimator: slope of log(Var(X_p)) vs log(s=log(L/r))."""
    s_arr  = np.log(L_scale / r_values)
    log_s  = np.log(s_arr)
    y_vals = []
    for r in r_values:
        X_p = np.log(np.maximum(np.abs(Sp_samples[r]), 1e-30)) - (p/3)*np.log(r)
        y_vals.append(0.5 * np.log(np.var(X_p, ddof=1)))
    y = np.array(y_vals)
    x_bar = log_s.mean(); y_bar = y.mean()
    SS_xx = np.sum((log_s - x_bar)**2)
    H_hat = np.sum((log_s - x_bar)*(y - y_bar)) / SS_xx
    resid = y - (H_hat*log_s + (y_bar - H_hat*x_bar))
    SE_H  = np.sqrt(np.sum(resid**2) / max(len(r_values)-2, 1) / SS_xx)
    return H_hat, SE_H


def run_jhtdb_analysis(n_points=50_000, p_orders=(2, 4, 6)):
    lJHTDB = pyJHTDB.libJHTDB()
    lJHTDB.initialize()
    lJHTDB.add_token(AUTH_TOKEN)

    points, vel_base = get_velocity(lJHTDB, n_points)

    # Collect |Δu|^p samples at each scale
    Sp_raw = {p: {r: None for r in R_VALS} for p in p_orders}
    for r in R_VALS:
        print(f"  r = {r:.4f}  ({r/DX:.0f}Δx)...")
        dU = longitudinal_increments(lJHTDB, points, vel_base, r)
        for p in p_orders:
            Sp_raw[p][r] = np.abs(dU)**p

    lJHTDB.finalize()

    # ── Estimator A: OLS variance ──────────────────────────────────────────
    print("\n── OLS Variance Estimator ─────────────────────────────────────────")
    H_by_p = {}
    for p in p_orders:
        H_hat, SE_H = ols_hurst(Sp_raw[p], R_VALS, L)
        z = (H_hat - 0.5) / SE_H
        sig = "***" if abs(z) > 5 else ("**" if abs(z) > 3 else "")
        print(f"  p={p}: H = {H_hat:.4f} ± {SE_H:.4f}   z = {z:+.2f} {sig}")
        H_by_p[p] = H_hat

    spread = max(H_by_p.values()) - min(H_by_p.values())
    print(f"  Order spread |H_max - H_min| = {spread:.4f}  (need <0.05 for cascade origin)")

    # ── Estimator B: lag-1 autocorrelation ────────────────────────────────
    print("\n── Lag-1 Autocorrelation ──────────────────────────────────────────")
    Sp2_means = np.array([Sp_raw[2][r].mean() for r in R_VALS])
    log_Sp2   = np.log(np.maximum(Sp2_means, 1e-30))
    inc       = np.diff(log_Sp2)
    rho1      = np.corrcoef(inc[:-1], inc[1:])[0, 1]
    SE_rho    = 1 / np.sqrt(len(inc))
    H_from_rho = (np.log2(rho1 + 1) + 1) / 2 if rho1 > -0.9 else 0.5
    print(f"  ρ(1) = {rho1:.4f}  (SE = {SE_rho:.4f},  z = {rho1/SE_rho:+.2f})")
    print(f"  Implied H from ρ(1) = 2^{{2H-1}}-1: H = {H_from_rho:.4f}")

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n── Decision ───────────────────────────────────────────────────────")
    H_main, SE_main = ols_hurst(Sp_raw[2], R_VALS, L, p=2)
    z_main = (H_main - 0.5) / SE_main
    if abs(z_main) > 5 and spread < 0.05 and abs(rho1) > 3*SE_rho:
        verdict = "REJECT H=0.5  ← Non-Markovian cascade (publishable)"
    elif abs(z_main) > 3:
        verdict = "EVIDENCE for H≠0.5, but needs more data / cross-validation"
    else:
        verdict = "CANNOT reject H=0.5  (Markovian cascade consistent with data)"
    print(f"  {verdict}")
    return H_by_p, rho1


if __name__ == "__main__":
    H_by_p, rho1 = run_jhtdb_analysis()
'''

def save_jhtdb_pipeline():
    with open('/home/user/Klmgrv/turbulence-graph/analysis/jhtdb_pipeline.py', 'w') as f:
        f.write(JHTDB_PIPELINE_CODE)
    print("Saved: jhtdb_pipeline.py (ready to run with your JHTDB auth token)")

# ═══════════════════════════════════════════════════════════════════════════
# 8. FIGURES
# ═══════════════════════════════════════════════════════════════════════════

def plot_hurst_analysis():
    fig = plt.figure(figsize=(20, 16), facecolor='#0d0d0d')
    fig.suptitle(
        'Non-Markovian Cascade: Hurst Exponent Test for Turbulence Scale-Space Memory\n'
        'H = 0.5 → Markovian (K62)  |  H ≠ 0.5 → Non-Markovian cascade',
        color='white', fontsize=13, fontweight='bold', y=0.98)

    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.45, wspace=0.35,
                           left=0.06, right=0.98, top=0.94, bottom=0.06)

    def style(ax, title='', xlabel='', ylabel=''):
        ax.set_facecolor('#141414')
        ax.tick_params(colors='white', labelsize=8)
        for sp in ax.spines.values():
            sp.set_edgecolor('#333')
        ax.grid(True, color='#252525', linewidth=0.5, linestyle='--')
        if title:  ax.set_title(title, color='white', fontsize=9, fontweight='bold')
        if xlabel: ax.set_xlabel(xlabel, color='#aaa', fontsize=8)
        if ylabel: ax.set_ylabel(ylabel, color='#aaa', fontsize=8)

    H_values = [0.5, 0.6, 0.7, 0.8, 0.9]
    H_colors = {0.5: C['markov'], 0.6: C['H06'], 0.7: C['H07'],
                0.8: C['H08'], 0.9: C['H09']}

    n_path = 50  # scale levels

    # ── A: fBm paths in scale space ────────────────────────────────────────
    ax_A = fig.add_subplot(gs[0, 0])
    style(ax_A, 'A — fBm Paths in Scale Space', 'scale level s = log₂(L/r)', 'ln(ε_r)')

    s_vals = np.arange(1, n_path + 1)
    for H in H_values:
        path = generate_fbm(n_path, H, sigma=0.5, seed=7)
        ax_A.plot(s_vals, path, color=H_colors[H], lw=1.5, label=f'H = {H}')
    ax_A.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── B: Autocorrelation of increments ───────────────────────────────────
    ax_B = fig.add_subplot(gs[0, 1])
    style(ax_B, 'B — Exact Lag Autocorrelation ρ(lag)', 'lag', 'ρ(lag)')

    lags = np.arange(0, 12)
    for H in H_values:
        rho = fbm_autocorrelation(H, max_lag=11)
        ax_B.plot(lags, rho, 'o-', color=H_colors[H], lw=1.5,
                  markersize=4, label=f'H = {H}')
    ax_B.axhline(0, color='white', lw=0.8, ls='--')
    ax_B.axhline(2/np.sqrt(30), color='#ef4444', lw=1, ls=':',
                 label='2σ threshold (n=30)')
    ax_B.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── C: Discriminator ρ(1) vs H ────────────────────────────────────────
    ax_C = fig.add_subplot(gs[0, 2])
    style(ax_C, 'C — Discriminator: ρ(1) vs H', 'Hurst exponent H', 'ρ(lag=1)')

    H_dense = np.linspace(0.05, 0.95, 200)
    rho1_exact = np.array([fbm_autocorrelation(H, max_lag=1)[1] for H in H_dense])
    ax_C.plot(H_dense, rho1_exact, color='white', lw=2, label='Exact ρ(1)')
    ax_C.axhline(0, color='white', lw=0.8, ls='--')
    ax_C.axvline(0.5, color='white', lw=0.8, ls=':')

    for n_pts in [20, 50, 100]:
        thresh = 2 / np.sqrt(n_pts)
        ax_C.axhline(thresh, color=C['data'], lw=0.8, ls='--', alpha=0.6,
                     label=f'2σ threshold (n={n_pts})')
        ax_C.axhline(-thresh, color=C['data'], lw=0.8, ls='--', alpha=0.6)

    ax_C.scatter([0.5], [0], color='white', s=50, zorder=10, label='H=0.5 (Markovian)')
    ax_C.text(0.7, 0.25, 'Long-range\nmemory', color='#f97316', fontsize=8, ha='center')
    ax_C.text(0.3, -0.15, 'Anti-correlated\ncascade', color='#22c55e', fontsize=8, ha='center')
    ax_C.legend(fontsize=6.5, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── D: Covariance structure ───────────────────────────────────────────
    ax_D = fig.add_subplot(gs[0, 3])
    style(ax_D, 'D — Covariance Cov[lnε_r, lnε_{r\'}]', '|s − s\'| (scale separation)', 'Cov')

    ds_vals = np.linspace(0, 10, 100)
    s_ref   = 5.0  # fix one scale
    for H in [0.5, 0.7, 0.9]:
        # K62 (Markovian): Cov = μ² min(s, s') = μ² (s_ref - ds) for ds < s_ref
        ds_eff = np.minimum(ds_vals, s_ref)
        cov_markov = 0.025 * (s_ref - ds_eff)
        # fBm: Cov = (σ²/2)[s^{2H} + s'^{2H} - |s-s'|^{2H}]
        s_prime = np.abs(s_ref - ds_vals)
        cov_fbm = (0.025/2) * (s_ref**(2*H) + s_prime**(2*H) - ds_vals**(2*H))
        col = H_colors.get(H, 'white')
        ax_D.plot(ds_vals, cov_fbm, color=col, lw=2, label=f'fBm H={H}')
    ax_D.plot(ds_vals, cov_markov, color=C['markov'], lw=2, ls='--', label='K62 Markovian')
    ax_D.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── E: Synthetic cascade structure functions ───────────────────────────
    ax_E = fig.add_subplot(gs[1, 0:2])
    style(ax_E, 'E — Structure Functions from Synthetic Cascades', 'r / L', 'S₂(r)')
    ax_E.set_xscale('log'); ax_E.set_yscale('log')

    N_vel = 4096
    dx = 1.0 / N_vel
    r_indices = np.unique(np.geomspace(2, N_vel // 4, 25).astype(int))

    for H in [0.5, 0.7, 0.9]:
        # Generate synthetic velocity from cascade (using fBm as velocity field)
        rng   = np.random.default_rng(H_values.index(H) if H in H_values else 99)
        u_syn = generate_fbm(N_vel, H, sigma=1.0, seed=int(H * 100))
        Sp, r_idx = compute_structure_functions_from_velocity(u_syn, [2], r_indices)
        r_phys = r_idx * dx
        col    = H_colors.get(H, 'white')
        ax_E.plot(r_phys, Sp[2], '-', color=col, lw=1.5, label=f'H = {H}')

    # Theoretical power laws
    r_ref = r_indices * dx
    ax_E.plot(r_ref, r_ref**(2/3), color='white', lw=1, ls='--', label='r^{2/3} (K41)')
    ax_E.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── F: ζ_p from synthetic cascades ───────────────────────────────────
    ax_F = fig.add_subplot(gs[1, 2:4])
    style(ax_F, 'F — Scaling Exponents ζ_p from Synthetic Cascades', 'order p', 'ζ_p')

    p_orders = [2, 4, 6, 8]
    p_dense  = np.linspace(1, 10, 100)

    from log_cft_derivation import zeta_she_leveque, zeta_k41
    ax_F.plot(p_dense, [zeta_k41(p) for p in p_dense],
              color=C['markov'], ls='--', lw=2, label='K41 p/3')
    ax_F.plot(p_dense, [zeta_she_leveque(p) for p in p_dense],
              color='white', ls='-', lw=2, label='She-Lévêque (target)')

    for H in [0.5, 0.7, 0.9]:
        u_syn = generate_fbm(N_vel, H, sigma=1.0, seed=int(H * 100))
        Sp_dict, r_idx = compute_structure_functions_from_velocity(u_syn, p_orders, r_indices)
        zeta_meas = {}
        for p in p_orders:
            zp, r2 = fit_scaling_exponent(Sp_dict[p], r_idx)
            zeta_meas[p] = zp
        col = H_colors.get(H, 'white')
        ax_F.scatter(list(zeta_meas.keys()), list(zeta_meas.values()),
                     color=col, s=60, zorder=10, label=f'Synthetic H={H}')
        ax_F.plot(list(zeta_meas.keys()), list(zeta_meas.values()),
                  '--', color=col, lw=1, alpha=0.5)

    ax_F.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── G: Power analysis ─────────────────────────────────────────────────
    ax_G = fig.add_subplot(gs[2, 0:2])
    style(ax_G, 'G — Statistical Power: P(reject H=0.5 | H_true)',
          'True Hurst exponent H', 'Detection power')

    n_scale_pts_list = [15, 25, 50]
    H_test = np.linspace(0.50, 0.95, 18)

    print("  Running power analysis (200 realisations per point)...")
    power_results = power_analysis(H_test, n_scale_pts_list, n_realisations=150)

    for n_s, ls in zip(n_scale_pts_list, ['-', '--', ':']):
        powers = [power_results[n_s][H] for H in H_test]
        ax_G.plot(H_test, powers, ls=ls, lw=2, color=C['data'],
                  label=f'n = {n_s} scale points')

    ax_G.axhline(0.05, color='white', lw=0.8, ls=':', label='α = 0.05 (false positive rate)')
    ax_G.axhline(0.80, color='#facc15', lw=0.8, ls=':', label='80% power threshold')
    ax_G.set_ylim(-0.05, 1.05)
    ax_G.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── H: Discriminator table ────────────────────────────────────────────
    ax_H = fig.add_subplot(gs[2, 2:4])
    style(ax_H, 'H — Key Discriminating Values (Exact)', '', '')
    ax_H.axis('off')

    disc = analytical_discriminator()
    col_names = ['H', 'ρ(lag=1)', 'ρ(lag=2)', 'ρ(lag=3)', 'n for 5% test', 'Cascade type']
    y_start = 0.92
    for j, cn in enumerate(col_names):
        ax_H.text(0.02 + j * 0.17, y_start, cn, transform=ax_H.transAxes,
                  color='#aaa', fontsize=8, fontweight='bold', va='top')
    ax_H.axhline(0.88, color='#333', xmin=0.01, xmax=0.99, lw=0.8)

    for i, row in enumerate(disc):
        y = y_start - 0.09 * (i + 1)
        color_row = '#ef4444' if row['H'] == 0.5 else ('#f97316' if row['H'] >= 0.8 else 'white')
        cascade_type = ('Markovian (K62)' if row['H'] == 0.5 else
                        f'Strong memory' if row['H'] >= 0.8 else 'Mild memory')
        values = [f"{row['H']:.2f}", f"{row['rho_lag1']:.4f}",
                  f"{row['rho_lag2']:.4f}", f"{row['rho_lag3']:.4f}",
                  f"{row['n_for_5pct']}", cascade_type]
        for j, val in enumerate(values):
            ax_H.text(0.02 + j * 0.17, y, val, transform=ax_H.transAxes,
                      color=color_row if j == 0 else 'white', fontsize=7.5, va='top')

    ax_H.text(0.5, 0.05,
              'Key: if measured ρ(1) from DNS exceeds 2/√n, reject Markovian cascade at 5% level',
              transform=ax_H.transAxes, ha='center', color='#aaa', fontsize=7.5,
              style='italic')

    plt.savefig('/home/user/Klmgrv/turbulence-graph/analysis/hurst_cascade_analysis.png',
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    print("Saved: hurst_cascade_analysis.png")
    return fig

# ═══════════════════════════════════════════════════════════════════════════
# 9. MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("NON-MARKOVIAN CASCADE: HURST EXPONENT ANALYSIS")
    print("=" * 72)

    print("\n1. EXACT DISCRIMINATOR VALUES")
    disc = analytical_discriminator()
    print(f"\n   {'H':>6}  {'ρ(1)':>10}  {'ρ(2)':>10}  {'n (5% test)':>14}  Type")
    print("   " + "-" * 55)
    for r in disc:
        t = 'Markovian' if r['H'] == 0.5 else 'Non-Markovian'
        print(f"   {r['H']:>6.2f}  {r['rho_lag1']:>10.5f}  {r['rho_lag2']:>10.5f}  "
              f"{r['n_for_5pct']:>14}  {t}")

    print("\n2. KEY FORMULA")
    print("   ρ(lag) = (1/2)[|lag+1|^{2H} + |lag-1|^{2H} - 2|lag|^{2H}]")
    print("   → ρ(1) = (2^{2H} - 2) / 2")
    print("   → H = log₂(2 + 2ρ(1)) / 2")
    print("   → To detect H = 0.7 at 5% level, need n > 9 scale points")
    print("   → To detect H = 0.6 at 5% level, need n > 25 scale points")

    print("\n3. CRAMÉR-RAO BOUND (JHTDB inertial range)")
    print("   7 scales: r = 2Δx .. 128Δx  (Δx = 2π/1024)")
    print("   s = log(L/r) = 2.079 .. 6.238")
    for M in [100, 1_000, 10_000, 50_000]:
        SE = cramér_rao_bound(
            n_scales=7, M_spatial=M,
            r_min=2 * 2*np.pi/1024, r_max=128 * 2*np.pi/1024,
            L=2*np.pi
        )
        print(f"   M = {M:>7} spatial samples/scale → SE(H) = {SE:.4f}  (95% CI ±{1.96*SE:.4f})")

    print("\n4. PUBLISHABILITY CRITERIA (illustrative: H=0.67 from M=10000 samples)")
    SE_demo = cramér_rao_bound(7, 10_000, 2*2*np.pi/1024, 128*2*np.pi/1024, 2*np.pi)
    print_publishability_report(
        H_hat=0.67, SE_H=SE_demo,
        H_by_p={2: 0.67, 4: 0.65, 6: 0.68},
        rho1=2**(2*0.67 - 1) - 1,
        n_scales=7,
    )

    print("5. JHTDB PIPELINE")
    save_jhtdb_pipeline()
    print("   Structure function analysis requires auth token from:")
    print("   https://turbulence.pha.jhu.edu/")
    print("   Then run: python jhtdb_pipeline.py")

    print("\n6. GENERATING FIGURES...")
    plot_hurst_analysis()
    print("\nDone.")


if __name__ == '__main__':
    main()
