"""
cascadequant.plots
==================
Publication-quality figures for cascade analysis.
All functions return (fig, axes) — caller controls saving.
"""

from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from typing import Optional

from .core import CascadeModel, she_leveque, EXPERIMENTAL_ZETA


PALETTE = {
    "cascade": "#1f77b4",
    "k41":     "#d62728",
    "lognorm": "#ff7f0e",
    "data":    "#2ca02c",
    "neutral": "#7f7f7f",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


# ── Figure 1: Structure-function exponents vs experiment ───────────────────

def fig_zeta_comparison(
    models: Optional[dict[str, CascadeModel]] = None,
    save_path: Optional[str] = None,
) -> tuple:
    """
    ζ_p from SL cascade + K41 vs published experimental data.
    """
    if models is None:
        models = {"She-Lévêque (exact)": she_leveque()}

    p_arr = np.linspace(0, 10, 300)
    fig, ax = plt.subplots(figsize=(7, 5))

    for label, mdl in models.items():
        ax.plot(p_arr, mdl.zeta(p_arr), lw=2.5, label=label, color=PALETTE["cascade"])

    # K41
    ax.plot(p_arr, p_arr / 3, "--", lw=1.5, color=PALETTE["k41"], label="K41  (p/3)")

    # Experimental data
    ps_exp = sorted(EXPERIMENTAL_ZETA.keys())
    z_exp = [EXPERIMENTAL_ZETA[p][0] for p in ps_exp]
    z_err = [EXPERIMENTAL_ZETA[p][1] for p in ps_exp]
    ax.errorbar(
        ps_exp, z_exp, yerr=z_err,
        fmt="o", color=PALETTE["data"], ms=7, capsize=4,
        label="Experiment (Anselmet 1984, Benzi 1993,\n  Gotoh 2002, Iyer–Sreenivasan–Yeung 2020)",
        zorder=5,
    )

    ax.set_xlabel("Moment order $p$", fontsize=13)
    ax.set_ylabel(r"$\zeta_p$", fontsize=13)
    ax.set_title("Velocity structure-function scaling exponents", fontsize=12)
    ax.legend(fontsize=9, framealpha=0.9)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3.5)
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig, ax


# ── Figure 2: Multifractal spectrum f(h) ───────────────────────────────────

def fig_multifractal_spectrum(save_path: Optional[str] = None) -> tuple:
    mdl = she_leveque()
    h_arr = np.linspace(mdl._h_inf + 1e-6, mdl._h_inf + mdl.lam * mdl.mu - 1e-6, 400)
    f_arr = mdl.multifractal_spectrum(h_arr)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(h_arr, f_arr, lw=2.5, color=PALETTE["cascade"], label="SL spectrum")
    ax.axvline(1.0 / 9.0, color=PALETTE["neutral"], ls=":", lw=1.2,
               label="$h_{min}=1/9$ (vortex filaments, $f=1$)")
    ax.axvline(1.0 / 3.0, color=PALETTE["k41"], ls="--", lw=1.5,
               label=f"K41: $h=1/3$,  $f=$ {mdl.multifractal_spectrum(1/3):.3f}")
    ax.axhline(1.0, color=PALETTE["neutral"], ls=":", lw=0.8)
    ax.axhline(3.0, color=PALETTE["neutral"], ls=":", lw=0.8)
    ax.fill_between(h_arr, f_arr, alpha=0.12, color=PALETTE["cascade"])
    ax.set_xlabel("Hölder exponent $h$", fontsize=13)
    ax.set_ylabel("Fractal dimension $f(h)$", fontsize=13)
    ax.set_title("Multifractal spectrum of She-Lévêque turbulence", fontsize=12)
    ax.legend(fontsize=9)
    ax.set_ylim(0.5, 3.2)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig, ax


# ── Figure 3: Compound Poisson increment PDF ───────────────────────────────

def fig_increment_pdf(
    scales: list[float] = (0.1, 0.5, 1.0, 2.0),
    y_range: tuple = (-3.0, 1.0),
    save_path: Optional[str] = None,
) -> tuple:
    mdl = she_leveque()
    y_fine = np.linspace(*y_range, 600)
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(scales)))

    for s, col in zip(scales, colors):
        cdf = np.array([mdl.increment_cdf(y, s) for y in y_fine])
        pdf_approx = np.gradient(cdf, y_fine)
        ax.plot(y_fine, pdf_approx, lw=2, color=col, label=f"$s={s}$")

    ax.set_xlabel(r"$\log|\delta u_r|$", fontsize=13)
    ax.set_ylabel("Density (approximate)", fontsize=13)
    ax.set_title("Compound Poisson increment distribution at cascade depths $s$", fontsize=11)
    ax.legend(fontsize=10, title="Cascade depth $s$")
    ax.set_xlim(*y_range)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig, ax


# ── Figure 4: Hankel fingerprint test on simulated cascade vs lognormal ────

def fig_hankel_test(
    n_samples: int = 5000,
    n_bootstrap: int = 200,
    save_path: Optional[str] = None,
) -> tuple:
    from .statistics import CascadeFingerprintTest
    rng = np.random.default_rng(42)
    mdl = she_leveque()
    mu = mdl.mu
    lam = mdl.lam
    s_steps = [0.5, 1.0, 2.0, 4.0, 8.0]

    # Simulate log-Poisson increments at each scale
    ratios_lp, ratios_ln = [], []
    for s in s_steps:
        # Log-Poisson: Y = -s*h_inf + sum_{i=1}^N (-mu) where N~Poisson(lam*s)
        N = rng.poisson(lam * s, n_samples)
        Y_lp = -s * mdl._h_inf + N * (-mu)
        # Lognormal: same mean and variance
        mean_lp = float(np.mean(Y_lp))
        std_lp = float(np.std(Y_lp))
        Y_ln = rng.normal(mean_lp, std_lp, n_samples)

        r_lp = CascadeFingerprintTest(Y_lp, n_bootstrap).run().hankel_ratio
        r_ln = CascadeFingerprintTest(Y_ln, n_bootstrap).run().hankel_ratio
        ratios_lp.append(r_lp)
        ratios_ln.append(r_ln)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(s_steps, ratios_lp, "o-", lw=2, ms=8, color=PALETTE["cascade"],
            label="Log-Poisson cascade (SL)")
    ax.plot(s_steps, ratios_ln, "s--", lw=1.5, ms=7, color=PALETTE["lognorm"],
            label="Lognormal (K62)")
    ax.axhline(1.0, color="black", lw=1, ls="--", alpha=0.5,
               label="Exact prediction: ratio = 1")
    ax.set_xlabel("Cascade depth $s$", fontsize=13)
    ax.set_ylabel(r"Hankel ratio  $\kappa_3^2 / (\kappa_2\,\kappa_4)$", fontsize=12)
    ax.set_title("Cascade fingerprint test: log-Poisson vs lognormal", fontsize=11)
    ax.set_ylim(-0.5, 2.5)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig, ax


# ── Figure 5: Implied volatility smile — cascade vs Black-Scholes ──────────

def fig_vol_smile(
    S: float = 100.0,
    T: float = 0.25,
    r: float = 0.0,
    sigma_bs: float = 0.20,
    s_cascade: float = 0.25,
    sigma_spread: float = 0.10,
    save_path: Optional[str] = None,
) -> tuple:
    """
    sigma_spread is the *total* (not annual) log-return standard deviation.
    sigma_spread = sigma_annual * sqrt(T).  At T=0.25, sigma_spread=0.10 ↔ 20% annual.
    """
    from .finance import CascadeVolatility, _bs_price

    moneyness = np.linspace(0.70, 1.30, 61)
    Ks = S * moneyness
    cv = CascadeVolatility(s=s_cascade, sigma=sigma_spread)

    ivs_cascade = []
    ivs_bs = []
    for K in Ks:
        otype = "call" if K >= S else "put"
        iv_c = cv.implied_vol(S, K, T, r, otype, sigma_spread)
        ivs_cascade.append(iv_c * 100)
        ivs_bs.append(sigma_bs * 100)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(moneyness, ivs_cascade, lw=2.5, color=PALETTE["cascade"],
            label="Cascade model (log-Poisson)")
    ax.plot(moneyness, ivs_bs, "--", lw=1.5, color=PALETTE["k41"],
            label="Black-Scholes (flat smile)")
    ax.fill_between(moneyness, ivs_cascade, ivs_bs,
                    alpha=0.15, color=PALETTE["cascade"],
                    label="Edge region (wing mispricing)")
    ax.set_xlabel("Moneyness  $K/S$", fontsize=13)
    ax.set_ylabel("Implied volatility (%)", fontsize=12)
    ax.set_title("Cascade vol smile vs Black-Scholes", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig, ax


# ── Figure 6: Partition function and S-duality ─────────────────────────────

def fig_partition_sdual(save_path: Optional[str] = None) -> tuple:
    mdl = she_leveque()
    p_arr = np.linspace(0.05, 6, 200)
    log_z = mdl.log_partition(p_arr, method="sdual")
    mu = mdl.mu
    main_3term = mu * p_arr / 24 + 0.5 * np.log(2 * np.pi / (mu * p_arr)) - np.pi ** 2 / (6 * mu * p_arr)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    ax = axes[0]
    ax.plot(p_arr, log_z, lw=2.5, color=PALETTE["cascade"], label=r"$\log Z_{\rm dyn}(p)$")
    ax.plot(p_arr, main_3term, "--", lw=1.5, color=PALETTE["k41"],
            label=r"S-dual: $\mu p/24 + \frac{1}{2}\log(2\pi/\mu p) - \pi^2/(6\mu p)$")
    ax.set_xlabel("Moment order $p$", fontsize=12)
    ax.set_ylabel(r"$\log Z_{\rm dyn}(p)$", fontsize=12)
    ax.set_title("Partition function + S-dual formula", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    err = np.abs(log_z - main_3term)
    ax.semilogy(p_arr, err + 1e-25, lw=2, color=PALETTE["cascade"])
    ax.set_xlabel("$p$", fontsize=12)
    ax.set_ylabel("Absolute error", fontsize=12)
    ax.set_title(r"Error: 3-term S-dual formula vs exact", fontsize=11)
    ax.grid(True, alpha=0.25)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    return fig, axes


def generate_all_figures(out_dir: str) -> list[str]:
    """Generate all article figures and save to out_dir. Returns list of paths."""
    import os
    os.makedirs(out_dir, exist_ok=True)
    paths = []

    figs = [
        ("fig1_zeta.png",        lambda p: fig_zeta_comparison(save_path=p)),
        ("fig2_spectrum.png",     lambda p: fig_multifractal_spectrum(save_path=p)),
        ("fig3_pdf.png",          lambda p: fig_increment_pdf(save_path=p)),
        ("fig4_hankel.png",       lambda p: fig_hankel_test(save_path=p)),
        ("fig5_smile.png",        lambda p: fig_vol_smile(save_path=p)),
        ("fig6_partition.png",    lambda p: fig_partition_sdual(save_path=p)),
    ]
    for fname, fn in figs:
        path = os.path.join(out_dir, fname)
        fn(path)
        plt.close("all")
        paths.append(path)
        print(f"  saved {path}")
    return paths
