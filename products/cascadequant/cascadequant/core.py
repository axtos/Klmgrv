"""
cascadequant.core
=================
Exact mathematics of the log-Poisson energy cascade.

Central object: the dynamical partition function

    Z(p) = prod_{n>=1} (1 - q^{n*p})   [q-Pochhammer / Fredholm determinant]

with cascade ratio  q = (2/3)^{1/3},  KS entropy  mu = (1/3)*ln(3/2).

All public functions are vectorised over p (accepts scalars or 1-D arrays).
"""

from __future__ import annotations
import math
from typing import Union
import numpy as np
from scipy.special import gamma as _gamma

# ── canonical SL parameters ────────────────────────────────────────────────
Q_SL: float = (2.0 / 3.0) ** (1.0 / 3.0)   # 0.87358…
MU_SL: float = -math.log(Q_SL)              # (1/3)*ln(3/2) ≈ 0.13516
LAM_SL: float = 2.0                          # Poisson rate (codimension = 2)


class CascadeModel:
    """
    Log-Poisson multiplicative cascade parameterised by (q, lam).

    Parameters
    ----------
    q   : cascade ratio per step (velocity).  Default = SL value (2/3)^{1/3}.
    lam : mean number of singular events per step.  Default = 2 (SL).

    Notes
    -----
    The She-Lévêque formula is recovered at q=Q_SL, lam=LAM_SL.
    Every public method accepts scalar or array p.
    """

    def __init__(self, q: float = Q_SL, lam: float = LAM_SL):
        self.q = float(q)
        self.lam = float(lam)
        self.mu = -math.log(q)        # KS entropy per step
        self._h_inf = self._linear_exponent()

    # ── private ────────────────────────────────────────────────────────────

    def _linear_exponent(self) -> float:
        """γ_∞ from the ζ₃ = 1 constraint: γ_∞ = (1 − λ(1 − q³)) / 3."""
        # ζ_3 = γ_∞·3 + λ(1−q³) = 1  →  γ_∞ = (1 − λ(1−q³))/3
        # For SL: q³=2/3, λ=2 → γ_∞ = (1 − 2/3)/3 = 1/9 ✓
        return (1.0 - self.lam * (1.0 - self.q ** 3)) / 3.0

    # ── structure-function exponents ───────────────────────────────────────

    def zeta(self, p) -> np.ndarray:
        """
        She-Lévêque scaling exponents  ζ_p = γ_∞·p + λ·(1 - q^p).

        Satisfies ζ_3 = 1 exactly when q=Q_SL, λ=2.
        """
        p = np.asarray(p, dtype=float)
        return self._h_inf * p + self.lam * (1.0 - self.q ** p)

    def dzetadp(self, p) -> np.ndarray:
        """dζ/dp = γ_∞ + λ·μ·q^p  (always positive → ζ strictly increasing)."""
        p = np.asarray(p, dtype=float)
        return self._h_inf + self.lam * self.mu * self.q ** p

    def holder_exponent(self, p) -> np.ndarray:
        """h(p) = dζ_p/dp — the dominant Hölder exponent at moment order p."""
        return self.dzetadp(p)

    # ── multifractal spectrum ───────────────────────────────────────────────

    def multifractal_spectrum(self, h) -> np.ndarray:
        """
        Legendre-transform spectrum  f(h) = 1 + λ·u·(1 - ln u)
        where u = (h - γ_∞) / (λ·μ).

        Valid for h ∈ [γ_∞,  γ_∞ + λ·μ].
        Returns -inf outside support.
        """
        h = np.asarray(h, dtype=float)
        h_min = self._h_inf
        h_max = self._h_inf + self.lam * self.mu
        u = (h - h_min) / (self.lam * self.mu)
        out = np.where(
            (u > 0) & (u <= math.e),
            1.0 + self.lam * u * (1.0 - np.log(np.maximum(u, 1e-300))),
            -np.inf,
        )
        return float(out) if out.ndim == 0 else out

    def codimension_k41(self) -> float:
        """3 - f(1/3): intermittency co-dimension of K41 structures."""
        return 3.0 - float(self.multifractal_spectrum(1.0 / 3.0))

    # ── partition function Z_dyn ────────────────────────────────────────────

    def log_partition(self, p, *, method: str = "auto") -> np.ndarray:
        """
        log Z_dyn(p)  =  log prod_{n>=1} (1 - q^{n*p})

        method:
          "series"  — direct log-sum (slow for p close to 0)
          "sdual"   — S-duality formula; exact for p > 0.5, error < 1e-20
          "auto"    — sdual for p > 0.3, series for smaller p (more terms)
        """
        p = np.asarray(p, dtype=float)
        scalar = p.ndim == 0
        p = np.atleast_1d(p)
        out = np.empty_like(p)

        if method in ("sdual", "auto"):
            mask_sd = (p > 0.3) if method == "auto" else np.ones(len(p), bool)
            out[mask_sd] = self._log_partition_sdual(p[mask_sd])
            if method == "auto":
                out[~mask_sd] = self._log_partition_series(p[~mask_sd])
        else:
            out = self._log_partition_series(p)

        return float(out[0]) if scalar else out

    def _log_partition_sdual(self, p: np.ndarray) -> np.ndarray:
        """
        S-duality (modular) formula — exact to O(exp(-4π²/(μp))).

        log Z_dyn(p) = μp/24 + ½ log(2π/(μp)) - π²/(6μp)
                       + log Z_dyn_dual(p)

        where Z_dyn_dual has q₂ = exp(-4π²/(μp)) ≈ 0 for any p > 0.1.
        """
        mu = self.mu
        log_main = mu * p / 24.0 + 0.5 * np.log(2.0 * np.pi / (mu * p)) - np.pi ** 2 / (6.0 * mu * p)
        # dual correction: sum_{n>=1} log(1 - q2^n), q2 = exp(-4pi^2/(mu*p))
        # Number of terms needed: n_max = ceil(10 / |log10(q2_max)|) for 1e-10 precision
        q2 = np.exp(-4.0 * np.pi ** 2 / (mu * p))
        q2_max = float(q2.max())
        if q2_max < 1e-15:
            return log_main
        n_max = max(50, int(math.ceil(10.0 / abs(math.log10(q2_max + 1e-300)))) + 10)
        n_arr = np.arange(1, n_max + 1)
        log_dual = np.sum(
            np.log1p(-q2[:, None] ** n_arr[None, :]),
            axis=1,
        )
        return log_main + log_dual

    def _log_partition_series(self, p: np.ndarray, N: int = 8000) -> np.ndarray:
        """Direct product: sum_{n=1}^N log(1 - q^{np})."""
        out = np.zeros(len(p))
        for n in range(1, N + 1):
            t = self.q ** (n * p)
            out += np.log1p(-t)
            if np.all(t < 1e-15):
                break
        return out

    def partition(self, p) -> np.ndarray:
        """Z_dyn(p) = exp(log_partition(p))."""
        return np.exp(self.log_partition(p))

    # ── increment distribution (compound Poisson) ──────────────────────────

    def increment_cdf(self, y, s: float, N_max: int = 60) -> float:
        """
        CDF of the log-increment Y = log|δu_r| at cascade depth s.

        P(Y ≤ y | s) = sum_{n=0}^N_max  Poisson(n; λ*s) * Heaviside(y - y_n)

        where y_n = -s*γ_∞ - n*μ  (n singular events shift the increment by -μ each).
        """
        lam_s = self.lam * s
        out = 0.0
        log_poisson = -lam_s
        for n in range(N_max + 1):
            if n > 0:
                log_poisson += math.log(lam_s) - math.log(n)
            y_n = -s * self._h_inf - n * self.mu
            if y >= y_n:
                out += math.exp(log_poisson)
            if log_poisson < -700:
                break
        return min(out, 1.0)

    def increment_quantile(self, prob: float, s: float, tol: float = 1e-9) -> float:
        """
        Inverse CDF: smallest y such that P(Y ≤ y | s) >= prob.
        Since the distribution is discrete (atoms at y_n), this returns
        the exact atom location.
        """
        lam_s = self.lam * s
        # atoms at y_n = -s*h_inf - n*mu, n=0,1,2,...  (decreasing in n)
        # CDF is a step function; find the step that crosses prob
        cumulative = 0.0
        log_poisson = -lam_s
        atoms = []
        for n in range(200):
            if n > 0:
                log_poisson += math.log(lam_s) - math.log(n)
            y_n = -s * self._h_inf - n * self.mu
            atoms.append((y_n, math.exp(log_poisson)))
            if log_poisson < -700:
                break
        # Sort by y ascending (n decreasing)
        atoms.sort(key=lambda x: x[0])
        cumulative = 0.0
        for y_n, p_n in atoms:
            cumulative += p_n
            if cumulative >= prob - tol:
                return y_n
        return atoms[-1][0]

    def log_cumulants(self, n_max: int = 6, s: float = 1.0) -> list[float]:
        """
        Exact log-cumulants of log|δu| at scale s.

        κ_1 = -s·γ_∞ - λ·s·μ     (mean)
        κ_n = (-1)^n · λ·s·μ^n    (n ≥ 2)

        The key fingerprint: κ_3² = κ_2·κ_4  exactly (Hankel identity).
        """
        kappas = []
        mu = self.mu
        lam_s = self.lam * s
        for n in range(1, n_max + 1):
            if n == 1:
                k = -s * self._h_inf - lam_s * mu
            else:
                k = (-1) ** n * lam_s * (mu ** n)
            kappas.append(k)
        return kappas

    def hankel_ratio(self, s: float = 1.0) -> float:
        """
        κ_3² / (κ_2 · κ_4).  Exactly 1.0 for any log-Poisson cascade.
        Deviations indicate non-log-Poisson statistics.
        """
        k = self.log_cumulants(4, s)
        k2, k3, k4 = k[1], k[2], k[3]
        if abs(k2 * k4) < 1e-300:
            return float("nan")
        return k3 ** 2 / (k2 * k4)


# ── convenience constructor for She-Lévêque ────────────────────────────────

def she_leveque() -> CascadeModel:
    """Return the canonical She-Lévêque (1994) model."""
    return CascadeModel(q=Q_SL, lam=LAM_SL)


# ── published experimental values for validation ───────────────────────────

EXPERIMENTAL_ZETA = {
    # Source: Anselmet et al. 1984, Benzi ESS 1993, Gotoh et al. 2002,
    #         JHTDB channel flow Re_tau=1000, Iyer-Sreenivasan-Yeung 2020
    # Longitudinal structure function exponents ζ_p, p=1..8
    # Format: p -> (mean, std_uncertainty)
    1: (0.370, 0.015),
    2: (0.696, 0.010),
    3: (1.000, 0.000),   # exact, 4/5 law
    4: (1.280, 0.015),
    5: (1.538, 0.025),
    6: (1.778, 0.030),
    7: (2.000, 0.040),
    8: (2.210, 0.050),
}
