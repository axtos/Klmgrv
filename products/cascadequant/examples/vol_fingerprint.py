"""
vol_fingerprint.py
==================
End-to-end demo: fit a log-Poisson cascade to synthetic return data,
run the Hankel fingerprint test, price options, compute VaR.

Run with: python examples/vol_fingerprint.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cascadequant import (
    she_leveque,
    CascadeFingerprintTest,
    CascadeVolatility,
    fit_to_returns,
    cascade_var,
    EXPERIMENTAL_ZETA,
)

rng = np.random.default_rng(42)

# ── 1. SL model at a glance ────────────────────────────────────────────────
print("=" * 60)
print("She-Lévêque cascade model")
print("=" * 60)
mdl = she_leveque()
print(f"  q  = {mdl.q:.6f}   (cascade ratio)")
print(f"  μ  = {mdl.mu:.6f}   (KS entropy)")
print(f"  λ  = {mdl.lam:.1f}     (Poisson rate)")
print(f"  γ∞ = {mdl._h_inf:.6f}  (1/9, vortex filaments)")
print(f"  ζ₃ = {mdl.zeta(3.0):.8f}  (should be 1.0)")
print(f"  Hankel ratio = {mdl.hankel_ratio():.8f}  (should be 1.0)")

print("\nStructure-function exponents vs experiment:")
print(f"  {'p':>3}  {'SL':>8}  {'Exp':>8}  {'Err':>8}")
for p, (z_exp, z_err) in EXPERIMENTAL_ZETA.items():
    z_sl = float(mdl.zeta(float(p)))
    print(f"  {p:>3}  {z_sl:>8.4f}  {z_exp:>8.4f}  {abs(z_sl-z_exp):>8.4f}")

# ── 2. Fingerprint test: log-Poisson vs Gaussian ──────────────────────────
print("\n" + "=" * 60)
print("Cascade fingerprint test (n=20,000 each)")
print("=" * 60)

s_cascade = 1.0
# Log-Poisson increments (the "fingerprint" signal lives in log|increment|)
N_lp = rng.poisson(mdl.lam * s_cascade, 20_000)
log_inc_lp = -s_cascade * mdl._h_inf + N_lp * (-mdl.mu)

# Gaussian with same mean and variance (lognormal baseline — should FAIL)
log_inc_gauss = rng.normal(log_inc_lp.mean(), log_inc_lp.std(), 20_000)

result_lp    = CascadeFingerprintTest(log_inc_lp,    bootstrap_B=300).run()
result_gauss = CascadeFingerprintTest(log_inc_gauss, bootstrap_B=300).run()

print("\nLog-Poisson (SL cascade):")
print(result_lp.summary())

print("\nGaussian / lognormal (K41):")
print(result_gauss.summary())

# ── 3. Options pricing ─────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("Options pricing: cascade vs Black-Scholes")
print("  (s=1 cascade step, 20% base vol)")
print("=" * 60)

cv = CascadeVolatility(s=0.25, sigma=0.10)
# sigma_spread = sigma_annual * sqrt(T) = 0.20 * 0.50 = 0.10 (total vol, not annualized)
S, T, r, sigma_spread = 100.0, 0.25, 0.0, 0.10

print(f"\n{'K/S':>6}  {'Cascade IV':>12}  {'BS flat':>10}  {'Wing edge':>12}")
for mn in [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]:
    K = S * mn
    otype = "call" if mn >= 1 else "put"
    iv_c = cv.implied_vol(S, K, T, r, otype, sigma_spread=sigma_spread)
    bs_iv = 0.20
    edge = (iv_c - bs_iv) * 100
    print(f"  {mn:.2f}  {iv_c*100:>10.2f}%  {bs_iv*100:>8.1f}%  {edge:>+10.2f}%")

# ── 4. Risk metrics ────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("VaR and Expected Shortfall (cascade vs historical)")
print("  (100k simulated daily log-returns, t(4) noise)")
print("=" * 60)

# Realistic daily returns: t(4) fat tails, scale ~1%
daily_returns = rng.standard_t(df=4, size=100_000) * 0.01
risk = cascade_var(daily_returns, confidence=0.99)
print(f"  Cascade VaR 99%      = {risk['cascade_var_99']:.4f}")
print(f"  Historical VaR 99%   = {risk['historical_var_99']:.4f}")
print(f"  Cascade ES 99%       = {risk['cascade_es_99']:.4f}")
print(f"  Historical ES 99%    = {risk['historical_es_99']:.4f}")
print(f"  Implied μ            = {risk['implied_mu']:.5f}")
print(f"  SL μ                 = {risk['sl_mu']:.5f}")

# ── 5. Parameter fitting from returns ─────────────────────────────────────
print("\n" + "=" * 60)
print("Parameter fitting from financial returns")
print("=" * 60)
params = fit_to_returns(daily_returns)
print(f"  Fitted μ        = {params['mu']:.5f}")
print(f"  Fitted s        = {params['s']:.4f}")
print(f"  Hankel ratio    = {params['hankel_ratio']:.4f}  (1.0 → log-Poisson)")
print(f"  μ / μ_SL        = {params['mu'] / risk['sl_mu']:.3f}")

print("\nAll done. Figures can be generated with:")
print("  from cascadequant.plots import generate_all_figures")
print("  generate_all_figures('figures/')")
