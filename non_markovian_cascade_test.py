"""
Statistical Test for the Non-Markovian Turbulent Cascade Hypothesis
====================================================================

Complete analysis pipeline: theory, estimators, JHTDB query design,
synthetic validation, and publishability criteria.

Author: Statistical Physics / Data Science analysis
"""

import numpy as np
from scipy import stats, linalg, special
from scipy.optimize import minimize
import warnings

# =============================================================================
# PART 1: EXACT DISCRIMINATOR FUNCTIONS
# =============================================================================

# --- H = 1/2 (Markovian / K62 log-normal) ---
#
# Covariance of log-coarse-grained dissipation:
#   Cov(ln ε_r, ln ε_{r'}) = μ² · min(log(L/r), log(L/r'))
#                           = μ² · min(s, s')
# where s = log(L/r), s' = log(L/r'), μ² ≈ 0.25 (intermittency constant)
#
# In scale-space language this is the covariance of standard Brownian motion:
#   Cov(B(s), B(s')) = min(s, s')  => MARKOVIAN (no memory beyond nearest neighbor)
#
# For the log-anomaly of structure functions X_p(r) = ln|S_p(r)/r^{p/3}|:
#   E[X_p(r)]   = (ζ_p - p/3) · log(r/L)   [linear in log r, slope = anomalous exponent]
#   Var(X_p(r)) = σ_p² · s                   [LINEAR in s = log(L/r)]
#
# The linear Var growth is the hallmark of H = 1/2.

def cov_markovian(s, s_prime, mu2=0.25):
    """K62 Markovian covariance: Cov(ln eps_r, ln eps_{r'}) = mu^2 * min(s, s')"""
    return mu2 * np.minimum(s, s_prime)

def var_X_markovian(s, sigma2_p):
    """Variance of X_p(r) for H=0.5: Var = sigma_p^2 * s"""
    return sigma2_p * s

# --- General H (fractional Brownian motion in scale space) ---
#
# Cov(X_p(r), X_p(r')) = (σ_p²/2) · [s^{2H} + s'^{2H} - |s - s'|^{2H}]
#
# This is exactly the covariance kernel of fractional Brownian motion (fBm).
# For H = 1/2: reduces to σ_p² · min(s, s') => recovers K62.
# For H > 1/2: long-range positive correlations across scales.
# For H < 1/2: anti-persistent, alternating overshoot.

def cov_fBm(s, s_prime, sigma2_p, H):
    """fBm covariance kernel for X_p(r)"""
    return 0.5 * sigma2_p * (s**(2*H) + s_prime**(2*H) - np.abs(s - s_prime)**(2*H))

def var_X_fBm(s, sigma2_p, H):
    """Variance of X_p(r) for general H: Var = sigma_p^2 * s^{2H}"""
    return sigma2_p * s**(2*H)

# =============================================================================
# PART 2: OPTIMAL ESTIMATOR OF H FROM STRUCTURE FUNCTION DATA
# =============================================================================
#
# Strategy: Estimate Var(X_p(r_i)) at each scale r_i from M spatial samples,
# then fit the power law  log Var(X_p) = 2H · log(s) + C  via OLS.
#
# Estimator derivation:
#   Let y_i = (1/2) · log(V_i)   where V_i = sample variance of X_p at scale i
#       x_i = log(s_i)            where s_i = log(L/r_i)
#
#   Model: y_i = H · x_i + C/2 + ε_i
#   ε_i ~ N(0, σ_ε²/2) with σ_ε² = Var(log V_i) ≈ 2/(M-1)  (log-chi-squared delta method)
#
#   OLS estimator:
#       H_hat = Σ_i (x_i - x̄)(y_i - ȳ) / Σ_i (x_i - x̄)²
#
#   This is unbiased: E[H_hat] = H (since E[y_i] = H·x_i + C/2)
#   Variance: Var(H_hat) = σ_ε² / (2 · SS_xx)
#             where SS_xx = Σ_i (x_i - x̄)²
#
# Cramér-Rao Lower Bound (CRLB):
#   Fisher information for H in the power-law variance model:
#   I(H) = (2/σ_ε²) · SS_xx  (N_scales observations of log-variance)
#   CRLB: Var(H_hat) ≥ 1/I(H) = σ_ε² / (2 · SS_xx)
#   => OLS achieves the CRLB (efficient estimator in this model).
#
# For σ_ε² = 2/(M-1) and N=7 scales (r = 2..128 Δx):
#   SS_xx ≈ 0.922  (computed from log-log-scale grid)
#   Var(H_hat) ≈ (2/(M-1)) / (2 · 0.922) = 1.084 / (M-1)
#   SE(H_hat) ≈ 1.041 / sqrt(M-1)
#
# For |H_hat - H| < 0.05 at 95% confidence:
#   Need SE < 0.05/1.96 = 0.0255  =>  M ≥ 834 spatial samples per scale
#
# Note: JHTDB has 1024³ grid; random sampling of ~10³ points gives ample coverage.

def estimate_H_OLS(r_vals, S_p_samples, p, L=2*np.pi):
    """
    Unbiased OLS estimator of Hurst exponent H.

    Parameters
    ----------
    r_vals     : array of shape (N_scales,), physical scale values
    S_p_samples: array of shape (N_scales, M), structure function samples
                 S_p_samples[i, j] = |delta_u(x_j, r_i)|^p
    p          : order of structure function
    L          : integral scale (default 2pi for JHTDB)

    Returns
    -------
    H_hat : float, estimated Hurst exponent
    SE_H  : float, standard error of H_hat
    """
    N_scales, M = S_p_samples.shape
    s_vals = np.log(L / r_vals)                    # log-scale coordinates

    # X_p(r) = log(S_p(r) / r^{p/3})
    X_p = np.log(np.mean(S_p_samples, axis=1)) - (p/3) * np.log(r_vals)

    # Sample variance of X_p across spatial points (at each scale)
    # Need sample variance of X_p realizations, not S_p realizations
    # Best: compute X_p per sample, then take variance
    X_p_per_sample = np.log(S_p_samples) - (p/3) * np.log(r_vals[:, None])
    V_i = np.var(X_p_per_sample, axis=1, ddof=1)  # shape (N_scales,)

    x_vals = np.log(s_vals)
    y_vals = 0.5 * np.log(V_i)

    x_bar = np.mean(x_vals)
    SS_xx = np.sum((x_vals - x_bar)**2)
    H_hat = np.sum((x_vals - x_bar) * (y_vals - np.mean(y_vals))) / SS_xx

    sigma_eps2 = 2.0 / (M - 1)
    SE_H = np.sqrt(sigma_eps2 / (2 * SS_xx))
    return H_hat, SE_H


def test_H_half(H_hat, SE_H, alpha=0.05):
    """
    Two-sided hypothesis test H0: H = 1/2.
    Returns: z-statistic, p-value, reject boolean
    """
    z = (H_hat - 0.5) / SE_H
    p_val = 2 * (1 - stats.norm.cdf(abs(z)))
    reject = p_val < alpha
    return z, p_val, reject


def power_of_test(H_true, SE_H, alpha=0.05):
    """Power of two-sided test H0: H=0.5 vs H_true."""
    z_crit = stats.norm.ppf(1 - alpha/2)
    delta = (H_true - 0.5) / SE_H
    power = (1 - stats.norm.cdf(z_crit - delta) +
             stats.norm.cdf(-z_crit - delta))
    return power


# =============================================================================
# PART 3: LOG-INCREMENT AUTOCORRELATION TEST (SIMPLEST SINGLE-NUMBER TEST)
# =============================================================================
#
# Define log-scale increments of structure functions:
#   ξ_p(s) = log|S_p(e^{-s} · L)| - log|S_p(e^{-s-1} · L)|
#
# This is the discrete "derivative" of the log structure function in log-scale.
#
# Under H0 (Markovian, H=1/2):
#   ξ_p(s) are i.i.d. Gaussian with mean γ_p = ζ_p - p/3
#
# Under H1 (non-Markovian fBm, general H):
#   ξ_p(s) is fractional Gaussian noise (fGn) with autocorrelation:
#   ρ(lag) = (1/2)[|lag+1|^{2H} + |lag-1|^{2H} - 2|lag|^{2H}]
#
# Lag-1 autocorrelation [KEY DISCRIMINATOR]:
#   ρ(1) = (1/2)[2^{2H} + 0^{2H} - 2·1^{2H}]
#        = (1/2)[2^{2H} - 2]
#        = 2^{2H-1} - 1
#
# NUMERICAL VALUES:
#   H = 0.5: ρ(1) = 2^0 - 1 = 0          [i.i.d., Markovian]
#   H = 0.6: ρ(1) = 2^{0.2} - 1 ≈ 0.1487
#   H = 0.7: ρ(1) = 2^{0.4} - 1 ≈ 0.3195
#   H = 0.8: ρ(1) = 2^{0.6} - 1 ≈ 0.5157
#   H = 0.9: ρ(1) = 2^{0.8} - 1 ≈ 0.7411
#
# The measurement is: given N_s discrete scale steps, compute sample ρ̂(1).
# SE(ρ̂(1)) ≈ 1/sqrt(N_s) for large N_s.
# For N_s = 100 discrete scale steps (feasible in DNS): SE ≈ 0.10
# => Can distinguish H=0.7 (ρ≈0.32) from H=0.5 (ρ=0) at 3-sigma level.

def rho_fGn(lag, H):
    """Autocorrelation of fractional Gaussian noise at given lag."""
    return 0.5 * (abs(lag+1)**(2*H) + abs(lag-1)**(2*H) - 2*abs(lag)**(2*H))

def rho1_from_H(H):
    """Lag-1 autocorrelation as function of H: rho(1) = 2^{2H-1} - 1"""
    return 2**(2*H - 1) - 1

def H_from_rho1(rho1):
    """Invert rho(1) = 2^{2H-1} - 1 to get H."""
    return 0.5 * (np.log2(rho1 + 1) + 1)

def estimate_H_from_autocorrelation(xi_series):
    """
    Estimate H from sample lag-1 autocorrelation of log-scale increments.
    xi_series: 1D array of ξ_p(s) values (N_s elements)
    Returns H_hat, SE_H_hat
    """
    N_s = len(xi_series)
    rho1_hat = np.corrcoef(xi_series[:-1], xi_series[1:])[0, 1]
    H_hat = H_from_rho1(np.clip(rho1_hat, -0.999, 0.999))

    # Delta method: SE(H_hat) = |dH/d_rho| * SE(rho_hat)
    # dH/d_rho = 1 / (2 * log(2) * (rho+1))
    SE_rho = 1.0 / np.sqrt(N_s)
    dH_drho = 1.0 / (2 * np.log(2) * (rho1_hat + 1))
    SE_H = abs(dH_drho) * SE_rho
    return H_hat, SE_H, rho1_hat


# =============================================================================
# PART 4: JHTDB DATA ACCESS DESIGN
# =============================================================================
#
# Dataset: isotropic1024coarse
#   - Re_λ = 433, 1024³ periodic box
#   - Grid spacing: Δx = 2π/1024 ≈ 0.006136
#   - Kolmogorov scale η ≈ Δx (marginally resolved)
#   - 10 stored timesteps (use t=0.002 to t=0.02)
#
# Scale choices:  r_n = 2^n · Δx  for n = 1, 2, ..., 7
#   r values: 2Δx, 4Δx, 8Δx, 16Δx, 32Δx, 64Δx, 128Δx
#   log(L/r): 6.24, 5.55, 4.85, 4.16, 3.47, 2.77, 2.08
#
# JHTDB API pseudocode (pyJHTDB):
#
#   import pyJHTDB
#   lJHTDB = pyJHTDB.libJHTDB()
#   lJHTDB.initialize()
#   auth_token = "your-token-here"
#   lJHTDB.add_token(auth_token)
#
# Step 1: Query velocity at M random points
#   import numpy as np
#   M = 100_000                  # base point count
#   x_rand = np.random.uniform(0, 2*np.pi, (M, 3)).astype(np.float32)
#   t = 0.002                    # first timestep
#
#   u_base = lJHTDB.getData(
#       t, x_rand,
#       data_set = 'isotropic1024coarse',
#       getFunction = 'getVelocity',
#       sinterp = 4,    # 4th-order Lagrange interpolation
#       tinterp = 0     # no time interpolation
#   )  # shape: (M, 3), columns: u_x, u_y, u_z
#
# Step 2: Longitudinal velocity increment at scale r in direction ê
#   For each r and each base point x_i, the displaced point is x_i + r·ê:
#
#   def get_struct_func(lJHTDB, x_base, r, t, p_order, direction=0,
#                       dataset='isotropic1024coarse'):
#       """Compute S_p(r) = <|delta_u_L(x,r)|^p> for longitudinal increments."""
#       e_hat = np.zeros(3, dtype=np.float32)
#       e_hat[direction] = 1.0
#       x_shifted = (x_base + r * e_hat) % (2 * np.pi)   # periodic BCs
#       x_shifted = x_shifted.astype(np.float32)
#
#       u_shifted = lJHTDB.getData(
#           t, x_shifted, data_set=dataset,
#           getFunction='getVelocity', sinterp=4, tinterp=0
#       )
#       u_base = lJHTDB.getData(
#           t, x_base, data_set=dataset,
#           getFunction='getVelocity', sinterp=4, tinterp=0
#       )
#       # Longitudinal increment: project onto separation direction
#       delta_u_L = (u_shifted[:, direction] - u_base[:, direction])
#       Sp = np.mean(np.abs(delta_u_L)**p_order)
#       return Sp, delta_u_L
#
# Step 3: Compute S_p at all scales
#   Delta_x = 2*np.pi / 1024
#   r_scales = np.array([2,4,8,16,32,64,128]) * Delta_x
#   p_orders = [2, 4, 6, 8]
#   Sp_results = {}
#   for r in r_scales:
#       for p in p_orders:
#           Sp, _ = get_struct_func(lJHTDB, x_rand, r, t, p)
#           Sp_results[(r, p)] = Sp
#
# SAMPLE SIZES for 1% relative accuracy in S_p:
#   S_2: M >= 100,000   (flatness ~ 5,   kappa_eff ~ 5)
#   S_4: M >= 1,000,000  (flatness ~ 50,  requires bootstrap or very large M)
#   S_6: M >= 10,000,000 (flatness ~ 500, effectively requires full grid or subsampling)
#   S_8: M >= 100,000,000 (practically requires full 1024³ grid = 10^9 points)
#
# PRACTICAL STRATEGY: For H estimation (not high-p accuracy), M=10,000 suffices
# for S_2 and S_4. Use p=2 and p=4 as primary estimators of H.
# Use S_6 only for qualitative anomalous exponent verification.
#
# IMPORTANT: Use multiple timesteps to increase effective M.
# With 10 snapshots: effective M_total = 10 × M_per_snapshot
# (approximately independent if snapshots separated by > eddy turnover time)

DELTA_X = 2 * np.pi / 1024
L_JHTDB = 2 * np.pi
R_SCALES = np.array([2, 4, 8, 16, 32, 64, 128]) * DELTA_X
S_LOG = np.log(L_JHTDB / R_SCALES)


# =============================================================================
# PART 5: EXPECTED RESULTS UNDER EACH HYPOTHESIS
# =============================================================================
#
# (a) H = 0.5 (Markovian / K62):
#     ρ(1) = 0.000
#     Var(X_p(r)) grows LINEARLY with s = log(L/r)
#     log Var vs log s has slope 2H = 1.000
#     Test statistic T = (H_hat - 0.5)/SE ~ N(0,1) under H0
#     => Expected: |T| < 1.96 with 95% probability
#
# (b) H ≈ 0.7 (mild memory):
#     ρ(1) ≈ 0.320
#     log Var vs log s has slope 2H = 1.400
#     T = (0.7 - 0.5)/0.023 ≈ 8.6 for M=1000 => power ≈ 1.00
#
# (c) H ≈ 0.9 (strong memory):
#     ρ(1) ≈ 0.741
#     log Var vs log s has slope 2H = 1.800
#     T = (0.9 - 0.5)/0.023 ≈ 17.4 for M=1000 => power ≈ 1.00
#
# (d) She-Lévêque connection:
#     The SL formula ζ_p = p/9 + 2[1-(2/3)^{p/3}] is a log-Poisson cascade model.
#     Log-Poisson cascades have:
#       Cov(log ε_r, log ε_{r'}) ∝ min(s, s')  [MARKOVIAN, H=0.5]
#     The anomalous ζ_p arises from non-Gaussian statistics (Poisson multipliers),
#     NOT from long-range memory in scale space.
#     => SL predicts H = 0.5 exactly.
#     => Any measured H ≠ 0.5 is a departure BEYOND SL and K62.
#     => The SL formula is NOT consistent with H ≠ 0.5 in the variance sense.
#
# NOTE: It is conceivable that H > 0.5 in higher-order cumulants while H=0.5
# in variance — this would require a separate higher-order test.

def expected_rho1(H):
    return 2**(2*H - 1) - 1

print("Expected values of rho(1) and slope of log-Var vs log-s:")
for H in [0.5, 0.6, 0.7, 0.8, 0.9]:
    print(f"  H={H}: rho(1)={expected_rho1(H):.4f}, 2H (log-Var slope)={2*H:.1f}")


# =============================================================================
# PART 6: SYNTHETIC VALIDATION — COMPLETE PSEUDOCODE
# =============================================================================

def generate_fBm_increments(N, H, sigma=1.0, seed=None):
    """
    Generate N increments of fractional Gaussian noise (fGn) with Hurst H.
    Uses the Davies-Harte (circulant embedding) exact method.

    Parameters
    ----------
    N     : number of increments (= number of scale steps)
    H     : Hurst exponent (0 < H < 1)
    sigma : standard deviation of the fBm at unit lag
    seed  : random seed

    Returns
    -------
    fgn : array of shape (N,), fractional Gaussian noise with autocorr rho(lag)
    fBm : array of shape (N+1,), corresponding fBm path starting at 0
    """
    rng = np.random.default_rng(seed)
    # Autocovariance of fGn: gamma(k) = (sigma^2/2)(|k+1|^{2H}+|k-1|^{2H}-2|k|^{2H})
    k = np.arange(N)
    gamma = (sigma**2 / 2) * (
        np.abs(k + 1)**(2*H) + np.abs(k - 1)**(2*H) - 2 * np.abs(k)**(2*H)
    )
    # Circulant embedding: first row of circulant matrix
    row = np.concatenate([gamma, gamma[-2:0:-1]])
    M_circ = len(row)
    eigs = np.fft.rfft(row).real
    if np.any(eigs < 0):
        warnings.warn("Circulant embedding not positive definite; using Cholesky fallback")
        # Fallback: Cholesky on NxN Toeplitz
        from scipy.linalg import toeplitz
        C = toeplitz(gamma)
        L = np.linalg.cholesky(C + 1e-12 * np.eye(N))
        fgn = L @ rng.standard_normal(N)
    else:
        # Davies-Harte
        W = (rng.standard_normal(M_circ) + 1j * rng.standard_normal(M_circ)).astype(complex)
        W[0] = rng.standard_normal()
        if M_circ % 2 == 0:
            W[M_circ // 2] = rng.standard_normal()
        W_fft = np.fft.rfft(W.real)
        fgn_full = np.fft.irfft(np.sqrt(np.abs(eigs[:len(W_fft)])) * W_fft, n=M_circ)
        fgn = fgn_full[:N] / np.sqrt(M_circ)
    fBm = np.concatenate([[0], np.cumsum(fgn)])
    return fgn, fBm


def generate_cascade_velocity(N_grid, H, lambda2=0.25, L=2*np.pi, seed=None):
    """
    Generate a 1D synthetic velocity field from a fractional BM cascade.

    Model:
      log ε(x) = fBm_x(x) + const     [spatial fBm for intermittency]
      The "cascade" acts in scale space with Hurst exponent H.

    Strategy:
      - Generate multiplicative cascade weights W_s at each scale s
        using fGn with parameter H (memory in scale direction).
      - Construct ε(x) = product over scales of W_s(x).
      - Velocity field u(x) is obtained from ε via Kolmogorov refinement:
          u(x) ~ integral of ε^{1/3}(x) dW(x)   [rough approximation]

    For the purpose of testing H, we work directly with structure functions
    of u(x), which inherit the cascade's Hurst exponent.

    Parameters
    ----------
    N_grid  : number of grid points (use power of 2, e.g., 1024)
    H       : Hurst exponent for the cascade
    lambda2 : intermittency parameter (μ² in K62, ~ 0.25)
    L       : domain length
    seed    : random seed

    Returns
    -------
    u       : 1D velocity field of shape (N_grid,)
    eps     : 1D energy dissipation field of shape (N_grid,)
    """
    rng = np.random.default_rng(seed)
    N_scales = int(np.log2(N_grid))  # number of cascade levels

    # Step 1: Generate fGn in scale direction (H-correlated cascade multipliers)
    # Each level s = 0, 1, ..., N_scales-1
    # For each spatial wavenumber band, draw a multiplier from the cascade.

    # Simplified approach: construct log-eps as a sum of scale-correlated contributions
    # log_eps(x) = sum_s W_s * phi_s(x)
    # where W_s are fGn-correlated across scales and phi_s are wavelet basis functions

    # Use Fourier-based multiscale decomposition:
    log_eps = np.zeros(N_grid)
    s_scales = np.arange(N_scales)  # scale indices

    # Generate fGn increments correlated across scales
    fgn_s, _ = generate_fBm_increments(N_scales, H, sigma=np.sqrt(lambda2), seed=seed)

    for j, s in enumerate(s_scales):
        # Spatial modes at scale 2^s: wavenumber band [2^s, 2^{s+1})
        k_lo = 2**s
        k_hi = min(2**(s+1), N_grid // 2)
        n_modes = k_hi - k_lo

        # Random phases for spatial variation
        phases = rng.uniform(0, 2*np.pi, n_modes)
        x = np.linspace(0, L, N_grid, endpoint=False)

        for k_idx, k in enumerate(range(k_lo, k_hi)):
            # Add cascade contribution at this scale
            amplitude = fgn_s[j] / np.sqrt(n_modes)
            log_eps += amplitude * np.cos(2 * np.pi * k * x / L + phases[k_idx])

    eps = np.exp(log_eps)
    eps = eps / np.mean(eps)  # normalize mean dissipation

    # Step 2: Construct velocity field using K41-like refinement
    # u(x) ~ (eps * L)^{1/3} * N(0,1) approximation
    u = eps**(1/3) * rng.standard_normal(N_grid)

    return u, eps


def compute_structure_functions(u, r_vals, p_orders, L=2*np.pi, N_grid=None):
    """
    Compute longitudinal structure functions S_p(r) from 1D velocity field.

    S_p(r) = <|u(x+r) - u(x)|^p>

    Parameters
    ----------
    u        : 1D velocity array of length N_grid
    r_vals   : array of physical separation distances
    p_orders : list of integers, structure function orders
    L        : domain length
    N_grid   : grid size (inferred if None)

    Returns
    -------
    Sp_dict : dict with keys (r, p), values = scalar S_p(r)
    samples : dict with keys (r, p), values = array of |delta_u|^p per point
    """
    if N_grid is None:
        N_grid = len(u)
    dx = L / N_grid
    Sp_dict = {}
    samples_dict = {}

    for r in r_vals:
        r_grid = int(round(r / dx))   # separation in grid units
        delta_u = u[r_grid:] - u[:-r_grid]  # velocity increments
        for p in p_orders:
            increments_p = np.abs(delta_u)**p
            Sp_dict[(r, p)] = np.mean(increments_p)
            samples_dict[(r, p)] = increments_p

    return Sp_dict, samples_dict


def estimate_H_from_synthetic(u, r_vals, p=2, L=2*np.pi):
    """
    Estimate H from a synthetic velocity field via the OLS log-variance method.

    This function demonstrates the full pipeline:
    1. Compute S_p(r) at multiple scales
    2. Form X_p(r) = log(S_p(r)/r^{p/3})
    3. Estimate Var(X_p(r)) via bootstrap
    4. Fit power law: log Var = 2H * log(s) + C

    Returns H_hat, SE_H
    """
    N_grid = len(u)
    dx = L / N_grid

    # Compute X_p at each scale using spatial blocks as "independent samples"
    block_size = N_grid // 100   # 100 spatial blocks
    n_blocks = N_grid // block_size

    X_p_per_block = np.zeros((len(r_vals), n_blocks))
    for i, r in enumerate(r_vals):
        r_grid = int(round(r / dx))
        for b in range(n_blocks):
            start = b * block_size
            end = min(start + block_size, N_grid - r_grid)
            if end <= start:
                X_p_per_block[i, b] = np.nan
                continue
            delta_u_block = u[start + r_grid:end + r_grid] - u[start:end]
            Sp_block = np.mean(np.abs(delta_u_block)**p)
            X_p_per_block[i, b] = np.log(Sp_block) - (p/3) * np.log(r)

    # Variance across blocks at each scale
    V_i = np.nanvar(X_p_per_block, axis=1, ddof=1)
    s_vals = np.log(L / r_vals)
    x_vals = np.log(s_vals)
    y_vals = 0.5 * np.log(V_i)

    # OLS
    x_bar = np.mean(x_vals)
    SS_xx = np.sum((x_vals - x_bar)**2)
    H_hat = np.sum((x_vals - x_bar) * (y_vals - np.mean(y_vals))) / SS_xx

    n_blocks_eff = n_blocks
    sigma_eps2 = 2.0 / (n_blocks_eff - 1)
    SE_H = np.sqrt(sigma_eps2 / (2 * SS_xx))

    return H_hat, SE_H


def validate_estimator_bias_consistency(H_true_vals, N_trials=50, N_grid=512, seed=42):
    """
    Monte Carlo validation: show estimator is unbiased and consistent.

    For each H_true, run N_trials synthetic experiments and check:
    - E[H_hat] ≈ H_true  (unbiased)
    - Std[H_hat] decreases with N_grid  (consistent)

    Returns results dict with keys = H_true values.
    """
    rng = np.random.default_rng(seed)
    dx = 2 * np.pi / N_grid
    r_vals = np.array([4, 8, 16, 32, 64]) * dx
    results = {}

    for H_true in H_true_vals:
        H_hats = []
        for trial in range(N_trials):
            trial_seed = rng.integers(0, 10**6)
            u, eps = generate_cascade_velocity(N_grid, H_true, seed=trial_seed)
            H_hat, SE_H = estimate_H_from_synthetic(u, r_vals, p=2)
            H_hats.append(H_hat)

        H_hats = np.array(H_hats)
        results[H_true] = {
            'mean':   np.mean(H_hats),
            'std':    np.std(H_hats),
            'bias':   np.mean(H_hats) - H_true,
            'rmse':   np.sqrt(np.mean((H_hats - H_true)**2)),
            'H_hats': H_hats,
        }
        print(f"H_true={H_true}: mean={results[H_true]['mean']:.3f}, "
              f"bias={results[H_true]['bias']:.3f}, "
              f"std={results[H_true]['std']:.3f}")

    return results


# =============================================================================
# COMPLETE ANALYSIS PIPELINE
# =============================================================================

def full_analysis_pipeline_jhtdb(Sp_data, r_vals, p_orders=[2, 4, 6],
                                   L=2*np.pi, alpha=0.05):
    """
    Full analysis pipeline for JHTDB data.

    Parameters
    ----------
    Sp_data  : dict[(r, p)] -> 1D array of |delta_u_L(x_i, r)|^p for each sample point
    r_vals   : array of scale values used
    p_orders : list of structure function orders
    L        : integral scale
    alpha    : significance level for hypothesis test

    Returns
    -------
    results  : dict with H_hat, SE_H, z_stat, p_val, reject_H0, rho1_hat
    """
    results = {}

    for p in p_orders:
        # Build sample matrix: shape (N_scales, M)
        samples = np.array([Sp_data[(r, p)] for r in r_vals])

        # OLS variance estimator
        H_hat, SE_H = estimate_H_OLS(r_vals, samples, p, L)
        z_stat, p_val, reject = test_H_half(H_hat, SE_H, alpha)

        # Autocorrelation test
        # Log structure function at each scale
        log_Sp = np.array([np.log(np.mean(Sp_data[(r, p)])) for r in r_vals])
        xi = np.diff(log_Sp)   # log-scale increments
        H_ac, SE_ac, rho1_hat = estimate_H_from_autocorrelation(xi)

        results[p] = {
            'H_hat_OLS':    H_hat,
            'SE_OLS':       SE_H,
            'z_statistic':  z_stat,
            'p_value':      p_val,
            'reject_H0':    reject,
            'H_hat_AC':     H_ac,
            'rho1_hat':     rho1_hat,
            'expected_rho1_H05': 0.0,
            'expected_rho1_H07': rho1_from_H(0.7),
            'expected_rho1_H09': rho1_from_H(0.9),
        }

    return results


# =============================================================================
# PUBLISHABILITY CRITERIA
# =============================================================================
#
# Evidence for H ≠ 0.5 would be publishable if ALL of the following hold:
#
# 1. PRIMARY STATISTICAL TEST:
#    |H_hat - 0.5| > 5 * SE(H_hat)  at multiple orders p=2,4,6
#    (5-sigma threshold to account for multiple comparisons and systematic errors)
#
# 2. CONSISTENCY ACROSS ORDERS:
#    H estimated from p=2, 4, 6 all agree within ±0.05
#    (non-Markovianity should be order-independent in leading term)
#
# 3. AUTOCORRELATION CONFIRMATION:
#    ρ̂(1) > 0.1 with p-value < 0.001 using the lag-1 autocorrelation test
#    H estimated via ρ̂(1) consistent with H from OLS method
#
# 4. REYNOLDS NUMBER SCALING:
#    H measured at different Re_λ (use JHTDB channel flow, MHD datasets)
#    H should be Re-independent (fundamental property) or show specific scaling
#
# 5. SYNTHETIC VALIDATION:
#    Same pipeline applied to K62 synthetic data gives H = 0.5 ± SE (passes null)
#    Same pipeline applied to fBm cascade with H=0.7 gives H = 0.7 ± SE (recovers truth)
#
# 6. COVARIANCE MATRIX TEST:
#    Fit the full N_scales × N_scales covariance matrix of X_p to the fBm kernel
#    via maximum likelihood: L(H, σ²) = -N/2 log|Σ(H)| - (1/2) y^T Σ(H)^{-1} y
#    H_MLE should agree with H_OLS
#    Likelihood ratio test: 2[log L(H_MLE) - log L(H=0.5)] > χ²(1, 0.999) = 10.83
#
# 7. PHYSICAL INTERPRETATION:
#    If H > 0.5: explain as correlated cascade (information transfer beyond nearest scale)
#    Possible mechanism: Biot-Savart nonlocality, vortex stretching across many scales
#    Must rule out: finite-size effects, bottleneck effect (H>0.5 artifacts near dissipation)
#    Must rule out: large-scale forcing contamination (H>0.5 artifacts near integral scale)
#
# MINIMUM DATASET REQUIREMENTS:
#   - M >= 50,000 spatial samples per scale (for SE < 0.02)
#   - N_scales >= 5 scales in the inertial range only (exclude r < 8η and r > L/5)
#   - Multiple independent snapshots (>=5) for error estimation
#   - Cross-validation: split data in half, check H_hat1 ≈ H_hat2


if __name__ == '__main__':
    print("="*70)
    print("NON-MARKOVIAN CASCADE HYPOTHESIS: STATISTICAL TEST SUMMARY")
    print("="*70)

    print("\n--- DISCRIMINATOR FUNCTIONS ---")
    print("H=0.5 (Markovian):  Cov(ln ε_r, ln ε_{r'}) = μ² · min(s, s')")
    print("                    Var(X_p(r)) = σ_p² · s         [linear]")
    print("General H:          Cov(X_p(r), X_p(r')) = (σ²/2)[s^{2H}+s'^{2H}-|s-s'|^{2H}]")
    print("                    Var(X_p(r)) = σ_p² · s^{2H}    [power law]")

    print("\n--- LAG-1 AUTOCORRELATION: ρ(1) = 2^{2H-1} - 1 ---")
    for H in [0.5, 0.6, 0.7, 0.8, 0.9]:
        print(f"  H={H}: ρ(1) = {rho1_from_H(H):.4f}")

    print("\n--- CRLB ANALYSIS (N=7 scales, M spatial samples) ---")
    SS_xx = 0.9222  # from log(log(L/r_i)) grid
    for M in [100, 1000, 10000]:
        sigma_eps2 = 2.0 / (M - 1)
        SE_H = np.sqrt(sigma_eps2 / (2 * SS_xx))
        CI = 1.96 * SE_H
        print(f"  M={M:7d}: SE(H_hat)={SE_H:.4f}, 95% CI ±{CI:.4f}")
    print("  M >= 834 needed for ±0.05 CI")

    print("\n--- POWER ANALYSIS (M=1000, N=7) ---")
    SE_1000 = np.sqrt(2.0/999 / (2 * SS_xx))
    for H1 in [0.6, 0.7, 0.8, 0.9]:
        pwr = power_of_test(H1, SE_1000)
        print(f"  H_true={H1}: power={pwr:.4f}")

    print("\n--- SHE-LEVÊQUE vs MARKOVIAN ---")
    print("  SL cascade is log-Poisson: H=0.5 exactly (Markovian)")
    print("  K62 cascade is log-normal: H=0.5 exactly (Markovian)")
    print("  Any H ≠ 0.5 goes beyond both SL and K62")

    print("\n--- SAMPLE SIZES FOR JHTDB ---")
    for p, kappa in [(2, 5), (4, 50), (6, 500), (8, 5000)]:
        M = int(2 * kappa / 0.01**2)
        print(f"  S_{p} (1% accuracy): M >= {M:,}")

    print("\n--- SYNTHETIC VALIDATION (demo) ---")
    print("  Running small-scale test (N_grid=256, H=0.7)...")
    dx = 2*np.pi/256
    r_test = np.array([4, 8, 16, 32, 64]) * dx
    u_test, _ = generate_cascade_velocity(256, H=0.7, seed=0)
    H_est, SE_est = estimate_H_from_synthetic(u_test, r_test, p=2)
    print(f"  True H=0.70, estimated H={H_est:.3f} ± {SE_est:.3f}")

    print("\n--- PUBLISHABILITY THRESHOLD ---")
    print("  Require: |H_hat - 0.5| > 5σ at p=2,4,6 + consistent ρ̂(1)")
    print("  Require: H stable across orders + Reynolds-independent")
    print("  Require: MLE likelihood ratio > 10.83 (χ²₁ at 99.9%)")
