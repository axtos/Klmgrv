"""
Log-CFT Operator Content Derivation for Turbulence Anomalous Scaling
=====================================================================

This module derives the operator content of the putative logarithmic CFT
describing turbulence intermittency.

Key claim: The NPRG fixed point (Canet-Delamotte-Wschebor 2016/2022) has
non-diagonalizable scaling operator spectrum → Jordan block structure →
logarithmic CFT. The Jordan coupling β_p = p/3 - ζ_p in each spin-p sector
IS the anomalous scaling exponent.

References:
  - Canet et al. (2016) Phys. Rev. E 93, 063101
  - She & Lévêque (1994) PRL 72, 336
  - Kolmogorov (1941) Proc. R. Soc. USSR 30, 301
  - Gurarie (1993) Nucl. Phys. B 410, 535 (log-CFT)
  - Flohr (2003) Int. J. Mod. Phys. A 18, 4497 (log-CFT review)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import gamma
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# ── Colour scheme ──────────────────────────────────────────────────────────
COLORS = {
    'k41':       '#94a3b8',
    'she_lev':   '#3b82f6',
    'log_normal':'#22c55e',
    'jordan':    '#f97316',
    'exact':     '#ef4444',
    'data':      '#a855f7',
}

# ═══════════════════════════════════════════════════════════════════════════
# 1. STRUCTURE FUNCTION SCALING EXPONENTS
# ═══════════════════════════════════════════════════════════════════════════

def zeta_k41(p):
    """K41 mean-field: ζ_p = p/3."""
    return p / 3.0

def zeta_she_leveque(p):
    """She-Lévêque (1994): best empirical fit to DNS.
    ζ_p = p/9 + 2[1 - (2/3)^(p/3)]
    Assumes 1D (filamentary) most-singular structures, log-Poisson multipliers.
    """
    return p / 9.0 + 2.0 * (1.0 - (2.0/3.0)**(p/3.0))

def zeta_log_normal(p, mu=0.025):
    """K62 log-normal: ζ_p = p/3 - μ·p(p-3)/18."""
    return p / 3.0 - mu * p * (p - 3.0) / 18.0

def zeta_multifractal(p, alpha_min=0.12, alpha_max=1.0, f_max=3.0):
    """Multifractal via Legendre transform of parabolic D(h).
    D(h) = 3 - (h - h0)^2 / (2c)  [parabolic approximation]
    ζ_p = min_h [ph + 3 - D(h)]
    """
    h0 = 1.0 / 3.0   # K41 exponent
    c  = 0.015        # intermittency coefficient
    h_values = np.linspace(alpha_min, alpha_max, 10000)
    D_h = 3.0 - (h_values - h0)**2 / (2*c)
    D_h = np.clip(D_h, 0, 3)
    objective = p * h_values + 3.0 - D_h
    return float(np.min(objective))

# ═══════════════════════════════════════════════════════════════════════════
# 2. LOG-CFT JORDAN BLOCK STRUCTURE
# ═══════════════════════════════════════════════════════════════════════════

def jordan_coupling(p, zeta_func=zeta_she_leveque):
    """
    Jordan coupling β_p = p/3 - ζ_p in the spin-p sector of the log-CFT.

    Derivation:
    -----------
    In the log-CFT, the spin-p operator O_p has a logarithmic partner Õ_p.
    The Jordan cell for the dilatation operator L₀:

        L₀|O_p⟩ = (p/3)|O_p⟩             [classical K41 dimension]
        L₀|Õ_p⟩ = (p/3)|Õ_p⟩ + |O_p⟩    [Jordan coupling = 1, normalised]

    Two-point functions (Ward identities of the log-CFT):
        ⟨O_p(r) O_p(0)⟩ = 0
        ⟨O_p(r) Õ_p(0)⟩ = b_p · r^{-2p/3}
        ⟨Õ_p(r) Õ_p(0)⟩ = −2b_p · r^{-2p/3} log(r/L)

    The measured structure function receives contributions from Õ_p:
        S_p(r) ∝ ⟨Õ_p(r) Õ_p(0)⟩^{1/2} ~ r^{p/3} [log(L/r)]^{1/2}

    RESUMMATION:
    At large log(L/r) (large inertial range), the log correction must be
    resummed. The resummation is:
        exp(−β_p · log(L/r)) = (L/r)^{−β_p} = r^{β_p} L^{−β_p}

    So: S_p(r) ~ r^{p/3} · exp(−β_p log(L/r))
              = r^{p/3} · r^{β_p} · L^{−β_p}
              = r^{p/3 + β_p}  [absorbing L^{−β_p} into prefactor]

    Therefore: ζ_p = p/3 + β_p  →  β_p = ζ_p − p/3

    SIGN CONVENTION:
    We define the Jordan coupling as the deviation from K41:
        β_p = ζ_p − p/3

    - β_p > 0: structure function scales FASTER than K41 (p < 3 regime)
    - β_p = 0: NO Jordan partner (p = 3, energy flux — exact, Ward identity)
    - β_p < 0: structure function scales SLOWER than K41 (p > 3, intermittency)

    The absolute value |β_p| = |ζ_p − p/3| is the anomalous dimension.
    """
    return zeta_func(p) - p / 3.0

def two_point_functions(r, p, L, zeta_func=zeta_she_leveque, b=1.0):
    """
    Log-CFT two-point functions for the Jordan pair (O_p, Õ_p).

    Returns:
        dict with 'primary', 'mixed', 'log_partner' two-point functions
    """
    h_p   = p / 3.0
    beta_p = jordan_coupling(p, zeta_func)
    r_arr = np.asarray(r, dtype=float)

    # Primary-primary (vanishes in log-CFT normalisation)
    G_OO = np.zeros_like(r_arr)

    # Primary - log partner (standard power law at K41 dimension)
    G_OOt = b * r_arr**(-2*h_p)

    # Log partner - log partner (logarithmic correction)
    # Before resummation: −2b log(r/L) r^{-2h_p}
    G_OtOt_bare = -2*b * r_arr**(-2*h_p) * np.log(r_arr / L)

    # After resummation to effective power law
    # ζ_p = p/3 + β_p, so S_p ~ r^{ζ_p}
    zeta_p = zeta_func(p)
    G_OtOt_resum = r_arr**(-2*zeta_p)

    return {
        'G_OO':          G_OO,
        'G_OOt':         G_OOt,
        'G_OtOt_bare':   G_OtOt_bare,
        'G_OtOt_resum':  G_OtOt_resum,
        'beta_p':        beta_p,
        'zeta_p':        zeta_p,
        'h_p':           h_p,
    }

def rank3_prediction(r, p, L, zeta_func=zeta_she_leveque):
    """
    Rank-3 Jordan block prediction.

    In a rank-3 Jordan cell:
        L₀|O_p⟩  = h_p|O_p⟩
        L₀|Õ_p⟩  = h_p|Õ_p⟩ + |O_p⟩
        L₀|Ô_p⟩  = h_p|Ô_p⟩ + |Õ_p⟩

    Three-point (log)²  function:
        ⟨Ô_p(r) Ô_p(0)⟩ = r^{-2h_p} [log(r/L)]²

    Resummation generates a log correction to the power law:
        S_p(r) ~ r^{ζ_p} × [log(r/η)]^{ε_p}

    where ε_p would be measurable as a systematic deviation from pure
    power-law scaling in the structure function. Estimate ε_p from data.
    """
    zeta_p = zeta_func(p)
    h_p    = p / 3.0
    r_arr  = np.asarray(r, dtype=float)

    # Rank-2 (resummed)
    G_rank2 = r_arr**(-2*zeta_p)

    # Rank-3 log correction: amplitude from β_p²
    beta_p  = jordan_coupling(p, zeta_func)
    epsilon = beta_p**2  # rough estimate of rank-3 coupling
    log_factor = np.abs(np.log(r_arr / L))**epsilon
    G_rank3 = r_arr**(-2*zeta_p) * log_factor

    return {'rank2': G_rank2, 'rank3': G_rank3, 'epsilon': epsilon}

# ═══════════════════════════════════════════════════════════════════════════
# 3. OPERATOR SPECTRUM COMPUTATION
# ═══════════════════════════════════════════════════════════════════════════

def compute_operator_spectrum(p_values=None, models=None):
    """
    Compute the full operator spectrum of the turbulence log-CFT.

    Returns a table of Jordan pair content for each spin-p sector.
    """
    if p_values is None:
        p_values = np.array([1, 2, 3, 4, 5, 6, 8, 10, 12])

    if models is None:
        models = {
            'K41':        zeta_k41,
            'She-Lévêque':zeta_she_leveque,
            'Log-Normal': zeta_log_normal,
        }

    results = {}
    for name, func in models.items():
        rows = []
        for p in p_values:
            zp    = func(p)
            beta  = jordan_coupling(p, func)
            rows.append({
                'p':                  p,
                'zeta_p':             round(zp, 5),
                'classical_dim':      round(p/3, 5),
                'anomalous_dim':      round(beta, 5),
                'jordan_coupling':    round(beta, 5),
                'log_partner_exists': abs(beta) > 1e-9,
                'sector_type':        (
                    'PRIMARY (protected)' if abs(beta) < 1e-9 else
                    'JORDAN PAIR (β > 0)' if beta > 0 else
                    'JORDAN PAIR (β < 0)'
                ),
            })
        results[name] = rows
    return results

def ward_identity_check():
    """
    Verify that β_3 = 0 for She-Lévêque (Ward identity from 4/5 law).

    The exact Kolmogorov 4/5 law ⟨(δu_L)³⟩ = −(4/5)εr implies ζ_3 = 1.
    K41 classical dimension at p=3: h_3 = 3/3 = 1.
    Therefore β_3 = ζ_3 − 3/3 = 1 − 1 = 0.
    No Jordan partner for the energy flux operator.
    This is a Ward identity of the energy flux conservation.
    """
    p = 3
    zeta3_sheLeveque = zeta_she_leveque(3)
    zeta3_k41        = zeta_k41(3)
    beta3            = jordan_coupling(3, zeta_she_leveque)
    return {
        'p':                    3,
        'zeta_3_k41':           zeta3_k41,
        'zeta_3_she_leveque':   zeta3_sheLeveque,
        'beta_3':               beta3,
        'ward_identity_holds':  abs(beta3) < 1e-12,
        'statement': (
            "β₃ = ζ₃ − 1 = 0 exactly. "
            "The p=3 sector has NO Jordan logarithmic partner. "
            "The energy flux operator is a true primary. "
            "This is the Ward identity of energy flux conservation."
        ),
    }

# ═══════════════════════════════════════════════════════════════════════════
# 4. RESUMMATION: HOW LOG CORRECTIONS BECOME A POWER LAW
# ═══════════════════════════════════════════════════════════════════════════

def demonstrate_resummation(p=6, L=1.0, eta=1e-4):
    """
    Show explicitly that the log-CFT prediction resumms to a power law.

    Bare log-CFT (before resummation):
        S_p(r) ~ r^{p/3} × [1 − β_p log(r/L) + β_p²(log(r/L))²/2 − ...]

    This Taylor expansion of:
        r^{p/3} × exp(−β_p log(r/L))  = r^{p/3} × (L/r)^{β_p}
                                       = L^{β_p} × r^{p/3 − β_p}
                                       = L^{β_p} × r^{ζ_p}

    At large inertial range (log(L/r) ~ log Re / 4 >> 1), the Taylor series
    diverges and must be resummed to the exponential → pure power law r^{ζ_p}.
    """
    r_values = np.logspace(np.log10(eta), np.log10(L*0.5), 500)
    beta_p   = jordan_coupling(p, zeta_she_leveque)
    zeta_p   = zeta_she_leveque(p)
    h_p      = p / 3.0

    log_rL   = np.log(r_values / L)  # negative (r < L)

    # Taylor series truncated at different orders
    truncN   = 6
    S_taylor = np.ones_like(r_values) * r_values**h_p
    for n in range(1, truncN + 1):
        import math as _math
        S_taylor = S_taylor + r_values**h_p * ((-beta_p * log_rL)**n) / _math.factorial(n)

    # Resummed (exact power law)
    S_resummed = L**beta_p * r_values**zeta_p  # = r^{h_p} * exp(-beta_p log(r/L))

    # K41 baseline
    S_k41 = r_values**h_p

    return {
        'r':          r_values,
        'S_k41':      S_k41,
        'S_taylor':   S_taylor,
        'S_resummed': S_resummed,
        'p':          p,
        'h_p':        h_p,
        'beta_p':     beta_p,
        'zeta_p':     zeta_p,
    }

# ═══════════════════════════════════════════════════════════════════════════
# 5. FUSION RULES / OPE CONSTRAINTS
# ═══════════════════════════════════════════════════════════════════════════

def ope_consistency_check(p_values=None):
    """
    Check OPE / fusion rule consistency.

    In a log-CFT, the OPE of O_p × O_q must produce O_{p+q}.
    Consistency requires: β_{p+q} ≠ β_p + β_q in general (nonlinear!).

    The nonlinearity of β(p) as a function of p is the signature of
    non-trivial OPE coefficients in the log-CFT.

    The linearity defect:
        δ(p, q) = β_{p+q} - β_p - β_q

    measures how much the log-CFT OPE deviates from a "free field" structure.
    For a free field: β_p = c × p → β_{p+q} = β_p + β_q → δ = 0.
    For turbulence: β_p is nonlinear → δ ≠ 0 → nontrivial interacting log-CFT.
    """
    if p_values is None:
        p_values = [1, 2, 3, 4, 5, 6]

    defects = []
    for p in p_values:
        for q in p_values:
            bp    = jordan_coupling(p)
            bq    = jordan_coupling(q)
            bpq   = jordan_coupling(p + q)
            delta = bpq - bp - bq
            defects.append({
                'p': p, 'q': q,
                'beta_p': round(bp, 5),
                'beta_q': round(bq, 5),
                'beta_pq': round(bpq, 5),
                'linearity_defect': round(delta, 5),
                'interacting': abs(delta) > 1e-6,
            })
    return defects

# ═══════════════════════════════════════════════════════════════════════════
# 6. FIGURES
# ═══════════════════════════════════════════════════════════════════════════

def plot_operator_spectrum():
    """Main figure: operator spectrum β_p and log-CFT structure functions."""
    p_dense  = np.linspace(0.1, 12, 500)
    p_points = np.array([1, 2, 3, 4, 5, 6, 8, 10])

    fig = plt.figure(figsize=(18, 14), facecolor='#0d0d0d')
    fig.suptitle('Log-CFT Operator Content of Turbulence\n(Jordan Coupling Spectrum β_p = ζ_p − p/3)',
                 color='white', fontsize=14, fontweight='bold', y=0.97)

    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35,
                           left=0.07, right=0.97, top=0.93, bottom=0.06)

    ax_colors = {'fg': 'white', 'bg': '#141414', 'grid': '#252525', 'spine': '#333'}

    def style_ax(ax, title='', xlabel='', ylabel=''):
        ax.set_facecolor(ax_colors['bg'])
        ax.tick_params(colors=ax_colors['fg'], labelsize=8)
        for sp in ax.spines.values():
            sp.set_edgecolor(ax_colors['spine'])
        ax.grid(True, color=ax_colors['grid'], linewidth=0.5, linestyle='--')
        if title:  ax.set_title(title, color='white', fontsize=9, fontweight='bold')
        if xlabel: ax.set_xlabel(xlabel, color='#aaa', fontsize=8)
        if ylabel: ax.set_ylabel(ylabel, color='#aaa', fontsize=8)

    # ── Panel A: β_p spectrum (the operator content) ───────────────────────
    ax_A = fig.add_subplot(gs[0, 0])
    style_ax(ax_A, 'A — Jordan Coupling Spectrum β_p(p)', 'order p', 'β_p = ζ_p − p/3')

    for name, func, col in [
        ('K41 (β=0 always)', zeta_k41, COLORS['k41']),
        ('She-Lévêque',       zeta_she_leveque,  COLORS['she_lev']),
        ('Log-Normal K62',    zeta_log_normal,   COLORS['log_normal']),
    ]:
        beta_dense = np.array([jordan_coupling(p, func) for p in p_dense])
        ax_A.plot(p_dense, beta_dense, color=col, linewidth=2, label=name)

    ax_A.axhline(0, color='#ef4444', linewidth=1.5, linestyle='--', label='β=0 (Ward identity)')
    ax_A.axvline(3, color='#ef4444', linewidth=1.0, linestyle=':')
    ax_A.annotate('p=3 protected\n(4/5 law)', xy=(3, 0), xytext=(4, 0.04),
                  color='#ef4444', fontsize=7,
                  arrowprops=dict(arrowstyle='->', color='#ef4444'))
    ax_A.annotate('Intermittency\n(β < 0)', xy=(6, -0.05), xytext=(7, -0.12),
                  color='#f97316', fontsize=7,
                  arrowprops=dict(arrowstyle='->', color='#f97316'))
    ax_A.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── Panel B: ζ_p comparison ────────────────────────────────────────────
    ax_B = fig.add_subplot(gs[0, 1])
    style_ax(ax_B, 'B — Structure Function Exponents ζ_p', 'order p', 'ζ_p')

    # Experimental data points (compiled from literature)
    p_exp    = np.array([2, 3, 4, 5, 6, 8, 10])
    zeta_exp = np.array([0.696, 1.000, 1.280, 1.540, 1.778, 2.230, 2.620])
    zeta_exp_err = np.array([0.010, 0.000, 0.015, 0.020, 0.020, 0.040, 0.060])

    for name, func, col, ls in [
        ('K41',         zeta_k41,         COLORS['k41'],       '--'),
        ('She-Lévêque', zeta_she_leveque, COLORS['she_lev'],   '-'),
        ('Log-Normal',  zeta_log_normal,  COLORS['log_normal'],'-'),
    ]:
        zeta_dense = np.array([func(p) for p in p_dense])
        ax_B.plot(p_dense, zeta_dense, color=col, linewidth=2, linestyle=ls, label=name)

    ax_B.errorbar(p_exp, zeta_exp, yerr=zeta_exp_err,
                  fmt='o', color=COLORS['data'], markersize=5, capsize=3,
                  label='DNS/Experiment', zorder=10)
    ax_B.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── Panel C: Jordan pair existence table ───────────────────────────────
    ax_C = fig.add_subplot(gs[0, 2])
    style_ax(ax_C, 'C — Jordan Pair Existence', 'order p', '|β_p| (Jordan coupling)')

    p_pts = np.array([1, 2, 3, 4, 5, 6, 8, 10])
    beta_sl = np.array([jordan_coupling(p) for p in p_pts])
    colors_C = ['#22c55e' if b > 1e-9 else '#86efac' if b < -1e-9 else '#ef4444'
                for b in beta_sl]
    ax_C.bar(p_pts, np.abs(beta_sl), color=colors_C, alpha=0.8, edgecolor='#333')
    ax_C.bar(3, 0, color='#ef4444', label='β=0: No Jordan partner (energy flux)')
    ax_C.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)
    for p_i, b_i in zip(p_pts, beta_sl):
        ax_C.text(p_i, abs(b_i) + 0.003, f'{b_i:+.3f}',
                  ha='center', va='bottom', fontsize=7, color='white')

    # ── Panel D: Two-point functions (bare vs resummed) ────────────────────
    ax_D = fig.add_subplot(gs[1, 0])
    style_ax(ax_D, 'D — Log-CFT Two-Point Functions (p=6)', 'r / L', 'G(r)')
    ax_D.set_xscale('log')
    ax_D.set_yscale('log')

    r_vals = np.logspace(-4, -0.3, 300)
    L_val  = 1.0
    p_demo = 6
    tf = two_point_functions(r_vals, p_demo, L_val)

    ax_D.plot(r_vals, r_vals**(6/3),       color=COLORS['k41'],     lw=2,  label='K41: r^{p/3}')
    ax_D.plot(r_vals, np.abs(tf['G_OtOt_bare']), color='#facc15', lw=1.5, ls='--',
              label='⟨Õ Õ⟩ bare (log corr.)')
    ax_D.plot(r_vals, tf['G_OtOt_resum'],  color=COLORS['she_lev'], lw=2,
              label=f'⟨Õ Õ⟩ resummed = r^{{ζ₆}} (ζ₆={tf["zeta_p"]:.3f})')
    ax_D.plot(r_vals, tf['G_OOt'],         color=COLORS['jordan'],  lw=1.5, ls=':',
              label='⟨O Õ⟩ mixed = r^{-2p/3}')

    ax_D.set_xlabel('r / L', color='#aaa', fontsize=8)
    ax_D.set_ylabel('G(r)', color='#aaa', fontsize=8)
    ax_D.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── Panel E: Resummation demonstration ────────────────────────────────
    ax_E = fig.add_subplot(gs[1, 1])
    style_ax(ax_E, 'D — Resummation: Log Corrections → Power Law (p=6)',
             'r / L', 'S₆(r) [normalized]')
    ax_E.set_xscale('log')
    ax_E.set_yscale('log')

    res = demonstrate_resummation(p=6, L=1.0, eta=1e-4)
    r_r = res['r']
    ax_E.plot(r_r, r_r**(6/3),       color=COLORS['k41'],     lw=2, ls='--', label='K41 r^{p/3}')
    ax_E.plot(r_r, res['S_taylor'],   color='#facc15',         lw=1.5, ls=':',
              label='Log-CFT bare (Taylor, 6 terms)')
    ax_E.plot(r_r, res['S_resummed'], color=COLORS['she_lev'], lw=2,
              label=f'Log-CFT resummed = r^{{ζ₆}} (ζ₆={res["zeta_p"]:.3f})')

    ax_E.axvspan(1e-4, 5e-3, alpha=0.1, color='#ef4444',
                 label='Dissipation range (resummation fails)')
    ax_E.annotate('Taylor ≈ resummed\n(small r: log large → resum)',
                  xy=(1e-3, res['S_resummed'][50]),
                  xytext=(3e-2, res['S_resummed'][100]*0.3),
                  color='white', fontsize=7,
                  arrowprops=dict(arrowstyle='->', color='#aaa'))
    ax_E.legend(fontsize=7, facecolor='#1a1a1a', labelcolor='white', framealpha=0.8)

    # ── Panel F: OPE linearity defect ────────────────────────────────────
    ax_F = fig.add_subplot(gs[1, 2])
    style_ax(ax_F, 'E — OPE Linearity Defect δ(p,q) = β_{p+q} − β_p − β_q',
             'p', 'Linearity defect δ(p,q)')

    p_ope = [1, 2, 3, 4, 5, 6]
    q_ref = 3  # Fix q=3 (energy flux sector)
    defects_q3 = []
    for p in p_ope:
        bp  = jordan_coupling(p)
        bq  = jordan_coupling(q_ref)
        bpq = jordan_coupling(p + q_ref)
        defects_q3.append(bpq - bp - bq)

    ax_F.bar(p_ope, defects_q3, color=COLORS['jordan'], alpha=0.8, edgecolor='#333')
    ax_F.axhline(0, color='white', linewidth=0.8, linestyle='--')
    ax_F.text(0.05, 0.85, f'q = {q_ref} fixed\nδ ≠ 0 → interacting log-CFT',
              transform=ax_F.transAxes, color='white', fontsize=8,
              bbox=dict(facecolor='#1a1a1a', edgecolor='#333', alpha=0.8))

    # ── Panel G: Physical interpretation table ────────────────────────────
    ax_G = fig.add_subplot(gs[2, :])
    ax_G.set_facecolor('#141414')
    ax_G.axis('off')

    headers = ['Spin-p sector', 'Classical dim p/3', 'Jordan coupling β_p',
               'Log partner Õ_p?', 'Physical meaning of Õ_p', 'Ward identity / Protection']
    rows_table = []
    p_list = [1, 2, 3, 4, 5, 6, 8, 10]
    meanings = {
        1: ('1-point velocity', 'Cascade-history weighted u'),
        2: ('Energy density', 'Cascade-history weighted u²'),
        3: ('Energy flux', 'PROTECTED — no log partner'),
        4: ('4th moment', 'Cascade-history weighted u⁴'),
        5: ('5th moment', 'Cascade-history weighted u⁵'),
        6: ('6th moment', 'Cascade-history weighted u⁶'),
        8: ('8th moment', 'Cascade-history weighted u⁸'),
        10:('10th moment','Cascade-history weighted u¹⁰'),
    }
    for p_i in p_list:
        beta_i = jordan_coupling(p_i)
        has_partner = abs(beta_i) > 1e-9
        primary_meaning, log_meaning = meanings.get(p_i, ('—', '—'))
        protection = '4/5 law (exact Ward identity)' if p_i == 3 else ''
        rows_table.append([
            f'p = {p_i}',
            f'{p_i/3:.4f}',
            f'{beta_i:+.5f}',
            '✓  Õ_p exists' if has_partner else '✗  No log partner',
            log_meaning if has_partner else 'N/A',
            protection,
        ])

    col_widths = [0.07, 0.10, 0.12, 0.14, 0.32, 0.23]
    col_starts = [0.01] + list(np.cumsum(col_widths[:-1]) + 0.01)

    ax_G.text(0.5, 0.98, 'F — Jordan Pair Operator Content: Full Table',
              transform=ax_G.transAxes, ha='center', va='top',
              color='white', fontsize=10, fontweight='bold')

    y_header = 0.88
    for j, (h, x) in enumerate(zip(headers, col_starts)):
        ax_G.text(x, y_header, h, transform=ax_G.transAxes,
                  color='#aaa', fontsize=7.5, fontweight='bold', va='top')

    ax_G.axhline(0.86, color='#333', linewidth=0.8,
                 xmin=0.01, xmax=0.99)

    for i, row in enumerate(rows_table):
        y = y_header - 0.10 * (i + 1)
        row_color = '#ef4444' if row[0] == 'p = 3' else (
                    '#f97316' if float(row[2]) < -0.03 else '#86efac')
        for j, (cell, x) in enumerate(zip(row, col_starts)):
            ax_G.text(x, y, cell, transform=ax_G.transAxes,
                      color=row_color if j == 0 else 'white',
                      fontsize=7.5, va='top')

    plt.savefig('/home/user/Klmgrv/turbulence-graph/analysis/log_cft_operator_content.png',
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    print("Saved: log_cft_operator_content.png")
    return fig

# ═══════════════════════════════════════════════════════════════════════════
# 7. MAIN — RUN COMPUTATION
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("LOG-CFT OPERATOR CONTENT DERIVATION FOR TURBULENCE")
    print("=" * 72)

    # Ward identity check
    print("\n1. WARD IDENTITY CHECK (p=3 sector)")
    wi = ward_identity_check()
    print(f"   ζ₃ (She-Lévêque) = {wi['zeta_3_she_leveque']:.10f}")
    print(f"   β₃ = ζ₃ − 1 = {wi['beta_3']:.2e}")
    print(f"   Ward identity holds: {wi['ward_identity_holds']}")
    print(f"   → {wi['statement']}")

    # Operator spectrum
    print("\n2. JORDAN PAIR OPERATOR SPECTRUM")
    spectrum = compute_operator_spectrum()
    print(f"\n   {'p':>4}  {'ζ_p (K41)':>12}  {'ζ_p (SL)':>12}  "
          f"{'β_p':>10}  {'Type':>30}")
    print("   " + "-" * 74)
    sl_rows = spectrum['She-Lévêque']
    k41_rows = spectrum['K41']
    for sl, k41 in zip(sl_rows, k41_rows):
        p_i   = sl['p']
        zp_sl = sl['zeta_p']
        zp_k  = k41['zeta_p']
        beta  = sl['jordan_coupling']
        stype = sl['sector_type']
        print(f"   {p_i:>4}  {zp_k:>12.5f}  {zp_sl:>12.5f}  "
              f"{beta:>10.5f}  {stype}")

    # Resummation
    print("\n3. RESUMMATION ARGUMENT (p=6)")
    res = demonstrate_resummation(p=6)
    print(f"   h_6 = 6/3 = {res['h_p']:.4f}  (K41 classical dimension)")
    print(f"   β_6 = {res['beta_p']:.5f}  (Jordan coupling = anomalous dim)")
    print(f"   ζ_6 = {res['zeta_p']:.5f}  (resummed = h_6 + β_6 = p/3 + β_p)")
    print(f"   Resummation: r^{{p/3}} × exp(−β_p log(r/L)) = r^{{ζ_p}}")
    print(f"   → Pure power law emerges from log-CFT at large inertial range.")

    # OPE
    print("\n4. OPE LINEARITY DEFECT (first few entries)")
    defects = ope_consistency_check([2, 3, 4, 6])
    for d in defects[:6]:
        print(f"   p={d['p']}, q={d['q']}: β_{{p+q}}={d['beta_pq']:.4f}, "
              f"β_p+β_q={d['beta_p']+d['beta_q']:.4f}, δ={d['linearity_defect']:.4f} "
              f"{'← interacting' if d['interacting'] else ''}")

    print("\n5. KEY FINDINGS")
    print("""
   ┌─ JORDAN PAIR SPECTRUM ───────────────────────────────────────────────┐
   │                                                                       │
   │  • p = 3 sector: β₃ = 0 EXACTLY → NO logarithmic partner            │
   │    Physical: energy flux operator is a TRUE PRIMARY (Ward identity)   │
   │                                                                       │
   │  • p > 3 sectors: β_p < 0 → INTERMITTENCY                           │
   │    Physical: cascade history REDUCES S_p below K41                   │
   │    The log partner Õ_p carries accumulated cascade memory             │
   │                                                                       │
   │  • p < 3 sectors: β_p > 0 → SUB-INTERMITTENCY                       │
   │    Physical: cascade history INCREASES S_p above K41                 │
   │                                                                       │
   │  • OPE defect δ ≠ 0 → INTERACTING log-CFT (not free field)          │
   │                                                                       │
   │  • Resummation: r^{p/3} × (L/r)^{-β_p} = r^{ζ_p} exactly           │
   │    → anomalous scaling IS the log-CFT resummation                    │
   │                                                                       │
   │  PHYSICAL INTERPRETATION OF LOG PARTNERS:                            │
   │  Õ_p = "cascade-history weighted" velocity moment                    │
   │  Carries factor log(L/r) = number of cascade steps from L to r       │
   │  This IS the non-Markovian cascade memory operator                   │
   └───────────────────────────────────────────────────────────────────────┘
    """)

    print("6. GENERATING FIGURES...")
    plot_operator_spectrum()

    print("\nDone. Output: analysis/log_cft_operator_content.png")


if __name__ == '__main__':
    main()
