"""
cascadequant.finance
====================
Cascade-based financial volatility analytics.

The key insight: if log-volatility follows a log-Poisson cascade, then:

  1. The Hankel identity  κ₃² = κ₂·κ₄  holds in the realized-vol time series.
  2. The cascade PDF gives better tail estimates than Black-Scholes or t-distributions.
  3. Option prices can be computed analytically from the compound Poisson PDF.

Entry points
------------
  CascadeVolatility  — fit cascade model to returns and price options
  cascade_var        — compute Value-at-Risk from cascade tail
  fit_to_returns     — estimate cascade parameters from return data
"""

from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm
from typing import Optional

from .core import CascadeModel, she_leveque, Q_SL, MU_SL, LAM_SL


# ── parameter fitting ───────────────────────────────────────────────────────

def fit_to_returns(
    log_returns: np.ndarray,
    timescale_days: float = 1.0,
    cascade_days: float = 260.0,
) -> dict:
    """
    Fit a log-Poisson cascade to financial log-returns.

    The cascade model over s cascade steps (s = log_horizon / log(1/q)):
        log|r| ~ N(μ_drift, σ_base²) * cascade_multiplier

    We fit:
      - σ_base  from variance after cascade correction
      - cascade timescale s  (single free parameter)
      - q, λ  fixed at SL values (can relax)

    Returns dict with: q, lam, mu, s, sigma_base, fit_quality
    """
    r = np.asarray(log_returns, dtype=float)
    r = r[np.isfinite(r)]
    # Log-vol cumulants
    log_abs = np.log(np.abs(r[r != 0]))
    var_log = float(np.var(log_abs))
    skew_log = float(np.mean((log_abs - log_abs.mean()) ** 3))
    kurt_log = float(np.mean((log_abs - log_abs.mean()) ** 4)) - 3 * var_log ** 2

    # From  κ₂ = λ·s·μ²  and  κ₃ = -λ·s·μ³  → μ = |κ₃/κ₂|, λ·s = κ₂/μ²
    mu_fit = abs(skew_log / var_log) if abs(var_log) > 1e-12 else MU_SL
    lam_s = var_log / mu_fit ** 2 if mu_fit > 1e-10 else LAM_SL
    q_fit = math.exp(-mu_fit)

    # Timescale in cascade units
    s_fit = lam_s / LAM_SL  # with lam fixed at 2, s = lam_s/2

    # Sigma base from residual variance
    total_var = float(np.var(r))
    sigma_base = math.sqrt(max(total_var - var_log, 1e-10))

    # Hankel quality check
    from .statistics import _log_cumulants
    k = _log_cumulants(log_abs, 4)
    denom = k[2] * k[4]
    hankel = k[3] ** 2 / denom if abs(denom) > 1e-20 else float("nan")

    return {
        "q": q_fit,
        "lam": LAM_SL,
        "mu": mu_fit,
        "s": s_fit,
        "sigma_base": sigma_base,
        "hankel_ratio": hankel,
        "mu_sl_deviation_pct": abs(mu_fit / MU_SL - 1) * 100,
    }


# ── main class ──────────────────────────────────────────────────────────────

class CascadeVolatility:
    """
    Options pricing and volatility analytics using the log-Poisson cascade.

    Parameters
    ----------
    q     : cascade ratio (default = SL)
    lam   : Poisson rate  (default = SL)
    s     : cascade depth in steps for the target horizon
    sigma : base volatility per cascade step
    """

    def __init__(
        self,
        q: float = Q_SL,
        lam: float = LAM_SL,
        s: float = 1.0,
        sigma: float = 0.20,
    ):
        self.model = CascadeModel(q=q, lam=lam)
        self.s = float(s)
        self.sigma = float(sigma)

    @classmethod
    def from_returns(cls, log_returns: np.ndarray, **kwargs) -> "CascadeVolatility":
        """Construct by fitting cascade parameters to empirical returns."""
        params = fit_to_returns(log_returns)
        return cls(
            q=params["q"],
            lam=params["lam"],
            s=params["s"],
            sigma=params["sigma_base"],
            **kwargs,
        )

    # ── return distribution ─────────────────────────────────────────────────

    def return_pdf_grid(
        self,
        n_atoms: int = 80,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Return (y_grid, pdf_weights) for the compound Poisson return distribution.

        The distribution of log-return Y = log|δS/S| over cascade depth s is
        a discrete mixture of Gaussians:

          P(Y) = sum_{n=0}^∞  Poisson(n; λs) · N(Y; y_n, σ²)

        where  y_n = -s·γ_∞ - n·μ  (mean shifts down by μ per singular event).
        """
        lam_s = self.model.lam * self.s
        mu = self.model.mu
        h_inf = self.model._h_inf
        sigma = self.sigma

        y_vals = []
        weights = []
        log_pois = -lam_s
        for n in range(n_atoms):
            if n > 0:
                log_pois += math.log(lam_s) - math.log(n)
            y_n = -self.s * h_inf - n * mu
            y_vals.append(y_n)
            weights.append(math.exp(log_pois))
            if log_pois < -700:
                break

        return np.array(y_vals), np.array(weights)

    def return_cdf(self, y_query: np.ndarray, sigma_spread: float = 0.0) -> np.ndarray:
        """
        CDF of log-return.  If sigma_spread > 0, atoms are Gaussian-broadened.
        """
        y_query = np.asarray(y_query, dtype=float)
        y_vals, weights = self.return_pdf_grid()
        if sigma_spread > 0:
            # Gaussian CDF at each atom
            result = np.zeros_like(y_query)
            for y_n, w in zip(y_vals, weights):
                result += w * norm.cdf(y_query, loc=y_n, scale=sigma_spread)
            return result
        else:
            # Step function
            result = np.zeros_like(y_query)
            for y_n, w in zip(y_vals, weights):
                result += w * (y_query >= y_n).astype(float)
            return result

    # ── option pricing ──────────────────────────────────────────────────────

    def option_price(
        self,
        S: float,
        K: float,
        T: float,
        r: float = 0.0,
        option_type: str = "call",
        sigma_spread: float = 0.0,
    ) -> float:
        """
        European option price via cascade PDF.

        Uses numerical integration over the compound Poisson distribution.
        For sigma_spread > 0, each cascade atom is broadened by a Gaussian
        (accounts for residual diffusion within a step).

        Parameters
        ----------
        S    : spot price
        K    : strike
        T    : time to expiry (years)
        r    : risk-free rate
        option_type : 'call' or 'put'
        sigma_spread : per-step diffusive volatility broadening
        """
        df = math.exp(-r * T)
        y_vals, weights = self.return_pdf_grid()

        # Risk-neutral drift correction: shift atoms so that E[S*e^Y] = S*e^{rT}.
        # Without this the compound-Poisson forward price is below S (cascade drift).
        if sigma_spread > 0:
            log_fwd_model = math.log(
                float(np.dot(weights, np.exp(y_vals + 0.5 * sigma_spread ** 2)))
            )
        else:
            log_fwd_model = math.log(float(np.dot(weights, np.exp(y_vals))) + 1e-300)
        rn_shift = r * T - log_fwd_model   # adjust to risk-neutral forward
        y_vals = y_vals + rn_shift

        if sigma_spread > 0:
            # Convolve each atom with Gaussian
            if option_type == "call":
                price = 0.0
                for y_n, w in zip(y_vals, weights):
                    # E[max(S*e^Y - K, 0)] where Y ~ N(y_n, sigma_spread^2)
                    fwd = S * math.exp(y_n + 0.5 * sigma_spread ** 2)
                    d1 = (y_n + sigma_spread ** 2 - math.log(K / S)) / sigma_spread
                    d2 = d1 - sigma_spread
                    atom_val = fwd * norm.cdf(d1) - K * norm.cdf(d2)
                    price += w * atom_val
            else:
                # Put via put-call parity (valid now that forward is risk-neutral)
                call = self.option_price(S, K, T, r, "call", sigma_spread)
                price = call - S + K * df
        else:
            # Pure discrete: each atom contributes max(S*e^y_n - K, 0)
            payoffs = np.maximum(S * np.exp(y_vals) - K, 0.0)
            if option_type == "put":
                payoffs = np.maximum(K - S * np.exp(y_vals), 0.0)
            price = float(np.dot(weights, payoffs))

        return df * price

    def implied_vol(
        self,
        S: float,
        K: float,
        T: float,
        r: float = 0.0,
        option_type: str = "call",
        sigma_spread: float = 0.02,
    ) -> float:
        """
        Black-Scholes implied volatility of the cascade option price.
        (Converts the cascade price into the 'BS-equivalent' vol.)
        """
        cascade_px = self.option_price(S, K, T, r, option_type, sigma_spread)
        return _bs_implied_vol(cascade_px, S, K, T, r, option_type)

    def vol_surface(
        self,
        S: float,
        moneyness: np.ndarray,   # K/S values
        maturities: np.ndarray,  # T values in years
        r: float = 0.0,
        sigma_spread: float = 0.02,
    ) -> np.ndarray:
        """
        Return implied-vol surface  shape (len(moneyness), len(maturities)).
        """
        out = np.full((len(moneyness), len(maturities)), float("nan"))
        for j, T in enumerate(maturities):
            for i, mn in enumerate(moneyness):
                K = S * mn
                try:
                    out[i, j] = self.implied_vol(S, K, T, r,
                                                  "call" if mn >= 1 else "put",
                                                  sigma_spread)
                except Exception:
                    pass
        return out

    # ── risk metrics ────────────────────────────────────────────────────────

    def var(self, confidence: float = 0.99, sigma_spread: float = 0.0) -> float:
        """
        Value at Risk at given confidence level.
        Returns the (1-confidence) quantile of log-return (negative = loss).
        """
        target = 1.0 - confidence
        y_vals, weights = self.return_pdf_grid()
        # Sort ascending
        idx = np.argsort(y_vals)
        y_sorted = y_vals[idx]
        w_sorted = weights[idx]
        cumulative = np.cumsum(w_sorted)
        # Find crossing
        crossing = np.searchsorted(cumulative, target)
        return float(y_sorted[min(crossing, len(y_sorted) - 1)])

    def expected_shortfall(
        self, confidence: float = 0.99, sigma_spread: float = 0.0
    ) -> float:
        """CVaR / Expected Shortfall at given confidence."""
        var_val = self.var(confidence)
        y_vals, weights = self.return_pdf_grid()
        tail_mask = y_vals <= var_val
        if not np.any(tail_mask):
            return var_val
        tail_prob = float(weights[tail_mask].sum())
        if tail_prob < 1e-15:
            return var_val
        return float(np.dot(weights[tail_mask], y_vals[tail_mask])) / tail_prob


# ── Black-Scholes helpers ───────────────────────────────────────────────────

def _bs_price(S, K, T, r, sigma, option_type="call"):
    if T <= 0 or sigma <= 0:
        return max(S - K, 0.0) if option_type == "call" else max(K - S, 0.0)
    d1 = (math.log(S / K) + (r + 0.5 * sigma ** 2) * T) / (sigma * math.sqrt(T))
    d2 = d1 - sigma * math.sqrt(T)
    if option_type == "call":
        return S * norm.cdf(d1) - K * math.exp(-r * T) * norm.cdf(d2)
    else:
        return K * math.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)


def _bs_implied_vol(
    market_price: float,
    S: float, K: float, T: float, r: float = 0.0,
    option_type: str = "call",
    tol: float = 1e-6,
) -> float:
    """Brent-method BS implied vol."""
    df = math.exp(-r * T)
    if option_type == "call":
        intrinsic = max(S - K * df, 0.0)
    else:
        intrinsic = max(K * df - S, 0.0)
    if market_price <= intrinsic + 1e-10:
        return 0.0
    try:
        iv = brentq(
            lambda s: _bs_price(S, K, T, r, s, option_type) - market_price,
            1e-6, 20.0, xtol=tol,
        )
    except ValueError:
        iv = float("nan")
    return iv


# ── standalone VaR function ─────────────────────────────────────────────────

def cascade_var(
    returns: np.ndarray,
    horizon_days: int = 1,
    confidence: float = 0.99,
    fit_cascade: bool = True,
) -> dict:
    """
    Compute cascade-based VaR and compare with historical simulation.

    Returns dict with keys: cascade_var, historical_var, cascade_es, hist_es.
    """
    r = np.asarray(returns, dtype=float)
    r = r[np.isfinite(r)]

    # Historical
    hist_var = float(np.percentile(r, (1 - confidence) * 100))
    hist_es = float(r[r <= hist_var].mean()) if (r <= hist_var).any() else hist_var

    # Cascade
    if fit_cascade:
        cv = CascadeVolatility.from_returns(r)
    else:
        cv = CascadeVolatility()
    casc_var = cv.var(confidence)
    casc_es = cv.expected_shortfall(confidence)

    return {
        "cascade_var_99": casc_var,
        "historical_var_99": hist_var,
        "cascade_es_99": casc_es,
        "historical_es_99": hist_es,
        "implied_mu": cv.model.mu,
        "sl_mu": MU_SL,
    }
