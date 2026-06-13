"""
cascadequant.statistics
=======================
Cascade fingerprint tests on empirical time-series data.

Primary test: Hankel identity  κ₃² = κ₂·κ₄  for log-increments.

This is a *parameter-free* discriminator between cascade models:
  - Lognormal (K62):     κ₃=κ₄=0                    → ratio = NaN
  - Log-Poisson (SL):   κ₃²  = κ₂·κ₄  exactly      → ratio = 1.0
  - Log-stable:          divergent / different ratio

Usage
-----
    from cascadequant.statistics import CascadeFingerprintTest
    result = CascadeFingerprintTest(log_increments).run()
    print(result.summary())
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Optional


# ── helpers ────────────────────────────────────────────────────────────────

def _log_cumulants(x: np.ndarray, n_max: int = 4) -> dict[int, float]:
    """
    Estimate log-cumulants κ₁..κ_{n_max} from data x via moment-cumulant relations.
    Uses unbiased central moments.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    m1 = x.mean()
    xc = x - m1
    moments = {1: m1}
    for k in range(2, n_max + 1):
        moments[k] = np.mean(xc ** k)

    # Moment → cumulant conversion (standard)
    kappas = {1: m1}
    kappas[2] = moments[2]
    kappas[3] = moments[3]
    kappas[4] = moments[4] - 3 * moments[2] ** 2
    if n_max >= 5:
        kappas[5] = moments[5] - 10 * moments[3] * moments[2]
    if n_max >= 6:
        kappas[6] = (moments[6] - 15 * moments[4] * moments[2]
                     - 10 * moments[3] ** 2 + 30 * moments[2] ** 3)
    return kappas


def _bootstrap_hankel(x: np.ndarray, B: int = 500) -> tuple[float, float]:
    """Bootstrap mean and std of the Hankel ratio κ₃²/(κ₂κ₄)."""
    n = len(x)
    ratios = []
    for _ in range(B):
        sample = x[np.random.randint(0, n, n)]
        k = _log_cumulants(sample, 4)
        denom = k[2] * k[4]
        if abs(denom) > 1e-20:
            ratios.append(k[3] ** 2 / denom)
    arr = np.array(ratios)
    return float(arr.mean()), float(arr.std())


# ── main result container ───────────────────────────────────────────────────

@dataclass
class FingerprintResult:
    n_samples: int
    kappa: dict[int, float]           # estimated cumulants
    hankel_ratio: float               # κ₃²/(κ₂κ₄)
    hankel_std: float                 # bootstrap std
    mu_implied: Optional[float]       # |κ₃/κ₂| = μ if log-Poisson
    lam_implied: Optional[float]      # κ₂/μ² = λ·s if log-Poisson
    is_log_poisson: bool              # ratio within 2σ of 1.0
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"Cascade Fingerprint Test  (n={self.n_samples:,})",
            "─" * 50,
            f"  κ₂ (log-variance)   = {self.kappa[2]:.6f}",
            f"  κ₃ (log-skewness)   = {self.kappa[3]:.6f}",
            f"  κ₄ (log-kurtosis)   = {self.kappa[4]:.6f}",
            "",
            f"  Hankel ratio κ₃²/(κ₂κ₄) = {self.hankel_ratio:.4f}  ±  {self.hankel_std:.4f}",
            f"  Log-Poisson prediction:      1.0000",
            f"  Is log-Poisson?              {'YES ✓' if self.is_log_poisson else 'NO  ✗'}",
        ]
        if self.mu_implied is not None:
            from .core import MU_SL
            lines += [
                "",
                f"  Implied μ  = |κ₃/κ₂| = {self.mu_implied:.6f}",
                f"  SL μ (theory)          = {MU_SL:.6f}",
                f"  Relative deviation     = {abs(self.mu_implied/MU_SL - 1)*100:.1f}%",
                f"  Implied λ·s = κ₂/μ²  = {self.lam_implied:.4f}",
            ]
        for note in self.notes:
            lines.append(f"  [note] {note}")
        return "\n".join(lines)


# ── primary test class ──────────────────────────────────────────────────────

class CascadeFingerprintTest:
    """
    Test whether a log-increment series is consistent with a
    log-Poisson (She-Lévêque) cascade.

    Parameters
    ----------
    log_increments : 1-D array
        Log of absolute increments:  log|X_{t+lag} - X_t|
        (log-velocity increments for turbulence, log-return for finance).
    bootstrap_B : int
        Bootstrap replications for uncertainty in Hankel ratio.
    """

    def __init__(self, log_increments: np.ndarray, bootstrap_B: int = 500):
        self.x = np.asarray(log_increments, dtype=float)
        self.x = self.x[np.isfinite(self.x)]
        self.B = bootstrap_B

    def run(self) -> FingerprintResult:
        kappas = _log_cumulants(self.x, 4)
        k2, k3, k4 = kappas[2], kappas[3], kappas[4]

        denom = k2 * k4
        if abs(denom) < 1e-20:
            hankel = float("nan")
            hankel_std = float("nan")
        else:
            hankel = k3 ** 2 / denom
            _, hankel_std = _bootstrap_hankel(self.x, self.B)

        mu_implied = abs(k3 / k2) if abs(k2) > 1e-20 else None
        lam_implied = k2 / mu_implied ** 2 if mu_implied else None
        is_lp = (not np.isnan(hankel)) and (abs(hankel - 1.0) < 2 * hankel_std)

        notes = []
        if len(self.x) < 500:
            notes.append("Small sample (<500) — bootstrap uncertainty may be large")
        if abs(k3) < 1e-10:
            notes.append("Near-zero skewness — consistent with lognormal, not log-Poisson")

        return FingerprintResult(
            n_samples=len(self.x),
            kappa=kappas,
            hankel_ratio=hankel,
            hankel_std=hankel_std,
            mu_implied=mu_implied,
            lam_implied=lam_implied,
            is_log_poisson=is_lp,
            notes=notes,
        )


# ── structure-function scaling estimator ───────────────────────────────────

def estimate_scaling_exponents(
    series: np.ndarray,
    orders: list[int],
    lags: list[int],
    fit_range: Optional[tuple[int, int]] = None,
) -> dict[int, float]:
    """
    Estimate ζ_p from  log S_p(r) ~ ζ_p · log r  over a range of lags.

    Parameters
    ----------
    series  : time series (velocity, price, etc.)
    orders  : moment orders p
    lags    : list of lag values in samples
    fit_range : (min_lag, max_lag) to use for the power-law fit

    Returns
    -------
    dict p -> ζ_p
    """
    lags = np.array(sorted(lags))
    if fit_range is not None:
        mask = (lags >= fit_range[0]) & (lags <= fit_range[1])
    else:
        mask = np.ones(len(lags), bool)

    log_r = np.log(lags[mask])
    exponents = {}
    for p in orders:
        sf = np.array([np.mean(np.abs(np.diff(series, n=1))[:-lag + 1 if lag > 1 else None] ** p)
                       for lag in lags])
        log_sf = np.log(sf[mask] + 1e-300)
        slope, _ = np.polyfit(log_r, log_sf, 1)
        exponents[p] = slope
    return exponents


# ── multiscale Hankel test (across lag/scale) ──────────────────────────────

def multiscale_fingerprint(
    series: np.ndarray,
    lags: list[int],
    bootstrap_B: int = 200,
) -> dict[int, FingerprintResult]:
    """
    Run the fingerprint test at multiple lags.
    Returns dict: lag -> FingerprintResult.
    Expected: Hankel ratio ≈ 1.0 at all lags for a log-Poisson cascade.
    """
    results = {}
    for lag in lags:
        increments = np.diff(series, n=1)
        increments = increments[::lag]       # sub-sampled at this lag
        log_inc = np.log(np.abs(increments[increments != 0]))
        test = CascadeFingerprintTest(log_inc, bootstrap_B)
        results[lag] = test.run()
    return results
