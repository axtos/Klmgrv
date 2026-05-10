
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
