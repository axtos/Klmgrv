"""
Turbulence L-function: Connection to Riemann Zeta and Polylogarithm

R20: TURBULENCE L-FUNCTION

    L(s) = sum_{p=1}^inf  zeta_p / p^s   (Dirichlet series of SL exponents)

DECOMPOSITION THEOREM:
    L(s) = (1/9)*zeta_R(s-1) + 2*zeta_R(s) - 2*Li_s(q)

where:
- zeta_R(s) = Riemann zeta function
- Li_s(q) = polylogarithm = sum_{n=1}^inf q^n/n^s  at q = (2/3)^{1/3}
- 1/9 = K41 exponent (from zeta_p = p/9 + 2(1-q^p))
- coefficient 2 = log-Poisson intensity parameter

POLES: simple poles at s=2 (Res=1/9) and s=1 (Res=2). Li_s(q) is entire.

R21: SPECIAL VALUES — TURBULENCE MEETS NUMBER THEORY

    L(3) = pi^2/54 + 2*zeta(3) - 2*Li_3(q)

where pi^2/6 = zeta_R(2) (Euler's Basel problem) and
      zeta(3) = Apery's constant = 1.20205690315959...
      Li_3(q) = 1.01191958... at q = (2/3)^{1/3}

Physical: K41 kinetic energy (pi^2/54) + bulk cascade (2*Apery) + deformation

R22: TURBULENCE MERTENS CONSTANT

    C_T = lim_{N->inf} [sum_{p=1}^N zeta_p/p - (19/9)*ln(N)]
        = (19/9)*gamma_E + 2*ln(1-q)

where gamma_E = Euler-Mascheroni constant = 0.5772...
      19/9 = 1/9 + 2 (combined pole residue at s=1)
      q = (2/3)^{1/3} (cascade factor)

Numerically: C_T = (19/9)*gamma_E + 2*ln(1-(2/3)^{1/3}) ~ -2.917

R23: SUMMATORY FUNCTION AND TURBULENT PRIME NUMBER THEOREM

The summatory L-value Sigma(N) = sum_{p=1}^N zeta_p/p satisfies:
    Sigma(N) ~ (19/9)*ln(N) + C_T  as N -> inf

This is the turbulence analogue of Mertens' prime theorem sum 1/p ~ ln(ln N).

R24: LERCH-HURWITZ REPRESENTATION

    L(s) = (1/9)*zeta_R(s-1) + 2*sum_{n=1}^inf (1-q^n)/n^s

Since 2*(1-q^n) = zeta_sl(n) - n/9 (the anomalous part without K41),
the second sum is the "anomalous Dirichlet series" A(s) = sum (zetasl(n)-n/9)/n^s.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import zeta as riemann_zeta

# ============================================================
# Constants
# ============================================================
q = (2/3)**(1/3)
mu = abs(np.log(q))
EULER_GAMMA = 0.5772156649015329

def zeta_sl(p):
    return p/9 + 2*(1 - q**p)

print("=" * 65)
print("TURBULENCE L-FUNCTION: CONNECTION TO NUMBER THEORY")
print("=" * 65)
print()
print(f"  q = (2/3)^{{1/3}} = {q:.10f}")
print(f"  mu = (1/3)*ln(3/2) = {mu:.10f}")
print(f"  gamma_E = {EULER_GAMMA:.10f}")
print()

# ============================================================
# R20: Turbulence L-function decomposition
# ============================================================
def Li_polylog(s, z, N=3000):
    """Li_s(z) = sum_{k=1}^inf z^k / k^s  (polylogarithm)"""
    result = 0.0
    zk = z
    for k in range(1, N+1):
        term = zk / k**s
        result += term
        zk *= z
        if abs(term) < 1e-15 * abs(result):
            break
    return result

def L_formula(s):
    """L(s) = (1/9)*zeta_R(s-1) + 2*zeta_R(s) - 2*Li_s(q)"""
    return (1/9)*riemann_zeta(s-1) + 2*riemann_zeta(s) - 2*Li_polylog(s, q)

def L_direct(s, N=2000):
    """Direct Dirichlet series"""
    return sum(zeta_sl(p) / p**s for p in range(1, N+1))

print("─" * 65)
print("R20: L(s) = (1/9)*zeta_R(s-1) + 2*zeta_R(s) - 2*Li_s(q)")
print("─" * 65)
print()
print("Verification:")
print(f"  {'s':>5}  {'L_direct':>13}  {'L_formula':>13}  {'error':>10}")
max_err = 0
for s in [3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0]:
    ld = L_direct(s)
    lf = L_formula(s)
    err = abs(ld - lf)
    max_err = max(max_err, err)
    print(f"  {s:>5.1f}  {ld:>13.8f}  {lf:>13.8f}  {err:>10.2e}")
print(f"  Max error: {max_err:.2e}")
print()

# Pole structure
print("Pole structure:")
print(f"  Pole at s=2: Res = 1/9 = {1/9:.8f}  (K41 exponent)")
print(f"  Pole at s=1: Res = 2         (log-Poisson intensity)")
print(f"  Li_s(q) entire (|q|<1 ensures convergence)")
print()

# ============================================================
# R21: Special values connecting turbulence to number theory
# ============================================================
print("─" * 65)
print("R21: SPECIAL VALUES — TURBULENCE MEETS NUMBER THEORY")
print("─" * 65)
print()

pi = np.pi
Apery = riemann_zeta(3)  # = 1.2020569...
pi2_over_54 = pi**2 / 54
Li3_q = Li_polylog(3, q)
Li2_q = Li_polylog(2, q)

print(f"L(3) = pi^2/54 + 2*zeta(3) - 2*Li_3(q)")
print(f"     = {pi2_over_54:.10f}  +  {2*Apery:.10f}  +  {-2*Li3_q:.10f}")
print(f"     = {pi2_over_54 + 2*Apery - 2*Li3_q:.10f}")
print()
print(f"  pi^2/54   = {pi2_over_54:.10f}  [Euler's pi^2/6 = zeta_R(2), K41 sector]")
print(f"  2*zeta(3) = {2*Apery:.10f}  [2 * Apery's constant, bulk cascade]")
print(f"  Li_3(q)   = {Li3_q:.10f}  [polylogarithm at q=(2/3)^{{1/3}}]")
print()
print(f"L(4) = (1/9)*zeta(3) + 2*(pi^4/90) - 2*Li_4(q)")
Li4_q = Li_polylog(4, q)
L4 = L_formula(4)
print(f"     = {Apery/9:.8f} + {2*pi**4/90:.8f} - {2*Li4_q:.8f} = {L4:.8f}")
print()
print(f"L(6) = (1/9)*zeta(5) + 2*(pi^6/945) - 2*Li_6(q)")
Li6_q = Li_polylog(6, q)
L6 = L_formula(6)
pi6 = 2*pi**6/945
print(f"     = {riemann_zeta(5)/9:.8f} + {pi6:.8f} - {2*Li6_q:.8f} = {L6:.8f}")
print()

# Ward identity connection: zeta_3 = 1 forces specific contribution to L(s)
print("Ward identity / 4/5 law: zeta_3 = 1 contributes exactly 1/3^s to L(s)")
for s_val in [3, 4, 6]:
    contrib = 1 / 3**s_val
    print(f"  L({s_val}) Ward contribution: 1/{3}^{s_val} = {contrib:.8f}")
print()

# ============================================================
# R22: Turbulence Mertens Constant
# ============================================================
print("─" * 65)
print("R22: TURBULENCE MERTENS CONSTANT")
print("─" * 65)
print()
# Correct asymptotic:
# Sigma(N) = sum zeta_p/p = N/9 + 2*H_N - 2*Li_1(q,N)
#          ~ N/9 + 2*ln(N) + 2*gamma_E - 2*Li_1(q)
#          = N/9 + 2*ln(N) + 2*gamma_E + 2*ln(1-q)

C_T_formula = 2 * EULER_GAMMA + 2 * np.log(1 - q)
print(f"Sigma(N) = sum_{{p=1}}^N zeta_p/p  ~  N/9 + 2*ln(N) + C_T")
print()
print(f"C_T = 2*gamma_E + 2*ln(1-q)")
print(f"    = 2*{EULER_GAMMA:.8f} + 2*ln(1-{q:.8f})")
print(f"    = {2*EULER_GAMMA:.8f} + {2*np.log(1-q):.8f}")
print(f"    = {C_T_formula:.8f}")
print()
print("Decomposition:")
print(f"  K41 linear growth: N/9  (from zeta_p = p/9 + ...)")
print(f"  Log growth: 2*ln(N)     (from 2/(s-1) pole of L(s))")
print(f"  Constant C_T = {C_T_formula:.6f}  (Euler-Mascheroni + q-log)")
print()
print("Verification: Sigma(N) - N/9 - 2*ln(N)  ->  C_T")
print()
print(f"  {'N':>8}  {'Sigma(N)':>12}  {'N/9+2lnN+C_T':>14}  {'residual':>10}")
for N in [100, 1000, 10000, 100000]:
    sigma_N = sum(zeta_sl(p)/p for p in range(1, N+1))
    approx = N/9 + 2*np.log(N) + C_T_formula
    diff = sigma_N - approx
    print(f"  {N:>8}  {sigma_N:>12.4f}  {approx:>14.4f}  {diff:>10.2e}")
print()

# ============================================================
# R23: Asymptotic of partial sums
# ============================================================
print("─" * 65)
print("R23: ANOMALOUS MERTENS  A(N) = Sigma(N)-N/9 ~ 2*ln(N) + C_T")
print("─" * 65)
print()
print(f"  2 = residue at s=1 [= log-Poisson intensity]")
print(f"  C_T = 2*gamma_E + 2*ln(1-q) = {C_T_formula:.8f}")
print()
print("  Decomposition of C_T:")
print(f"  2*gamma_E  = {2*EULER_GAMMA:.6f}  [Euler-Mascheroni constant * 2]")
print(f"  2*ln(1-q) = {2*np.log(1-q):.6f}  [cascade q-log: 2*ln(1-(2/3)^{{1/3}})]")
print(f"  Note: -ln(1-q) = Li_1(q) = {Li_polylog(1, q):.6f}")
print()

# ============================================================
# R24: Anomalous Dirichlet series
# ============================================================
print("─" * 65)
print("R24: ANOMALOUS DIRICHLET SERIES")
print("─" * 65)
print()
print("A(s) = sum_{n>=1} [zeta_n - n/9] / n^s = 2*zeta_R(s) - 2*Li_s(q)")
print("     = 2 * [zeta_R(s) - Li_s(q)]")
print("     = 2 * sum_{n>=1} (1-q^n) / n^s")
print()
print("This is the Dirichlet series of the anomalous exponents 2*(1-q^n).")
print("A(s) has a single pole at s=1 with residue 2.")
print()
for s in [3.0, 4.0, 6.0]:
    A_direct = sum(2*(1-q**n)/n**s for n in range(1, 3000))
    A_formula = 2*(riemann_zeta(s) - Li_polylog(s, q))
    print(f"  A({s:.0f}) = {A_direct:.8f}  (formula: {A_formula:.8f},  err: {abs(A_direct-A_formula):.2e})")

print()
print("=" * 65)
print("SUMMARY: Turbulence-Number Theory Bridge")
print("=" * 65)
print()
print("R20: L(s) = (1/9)*zeta_R(s-1) + 2*zeta_R(s) - 2*Li_s(q)")
print("     Poles at s=1 (Res=2) and s=2 (Res=1/9)")
print()
print("R21: L(3) = pi^2/54 + 2*Apery - 2*Li_3((2/3)^{1/3})")
print(f"          = {L_formula(3):.8f}")
print()
print(f"R22: C_T = 2*gamma_E + 2*ln(1-q) = {C_T_formula:.8f}")
print("     Turbulence Mertens: Sigma(N) = sum zeta_p/p ~ N/9 + 2*ln(N) + C_T")
print()
print("R24: Anomalous: A(s) = 2*(zeta_R(s) - Li_s(q)),  single pole Res=2")
print()
print("Physical meaning of number-theory connection:")
print("  pi^2/54:   K41 contribution to p=3 structure function spectrum")
print("  2*Apery:   log-Poisson bulk cascade (Apery = integral zeta)")
print("  Li_s(q):   Cascade deformation at q=(2/3)^{1/3}")
print("  C_T:       Turbulence spectral asymptotic (Euler-Mascheroni + log cascade)")
print()

# ============================================================
# Figures
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle(r'Turbulence L-function: $L(s) = \frac{1}{9}\zeta_R(s-1) + 2\zeta_R(s) - 2\mathrm{Li}_s(q)$',
             fontsize=13, fontweight='bold')

# Panel 1: L(s) decomposition
ax = axes[0, 0]
s_vals = np.linspace(2.1, 10, 200)
L_vals = [L_formula(s) for s in s_vals]
riemann_part = [2*riemann_zeta(s) for s in s_vals]
Li_part = [-2*Li_polylog(s, q) for s in s_vals]
K41_part = [(1/9)*riemann_zeta(s-1) for s in s_vals]

ax.plot(s_vals, L_vals, 'k-', linewidth=2.5, label='$L(s)$')
ax.plot(s_vals, K41_part, 'b--', alpha=0.7, linewidth=1.5, label=r'$\frac{1}{9}\zeta_R(s-1)$ (K41)')
ax.plot(s_vals, riemann_part, 'g:', alpha=0.7, linewidth=1.5, label=r'$2\zeta_R(s)$ (bulk)')
ax.plot(s_vals, Li_part, 'r-.', alpha=0.7, linewidth=1.5, label=r'$-2\mathrm{Li}_s(q)$ (deformation)')
ax.set_xlabel('$s$', fontsize=12)
ax.set_ylabel('Value', fontsize=12)
ax.set_title('Decomposition of turbulence L-function', fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim([2.1, 10])

# Panel 2: Residues and pole structure
ax = axes[0, 1]
s_near_2 = np.linspace(2.02, 3.5, 200)
L_near_2 = [L_formula(s) for s in s_near_2]
pole_approx_2 = [(1/9)/(s-2) + L_formula(2.5) - (1/9)/(2.5-2) for s in s_near_2]
ax.plot(s_near_2, L_near_2, 'b-', linewidth=2, label='$L(s)$')
ax.axvline(x=2, color='red', linestyle='--', alpha=0.5, label='pole s=2')
ax.scatter([3], [L_formula(3)], s=100, c='red', zorder=5,
           label=f'$L(3) = \\pi^2/54 + 2\\zeta(3) - 2\\mathrm{{Li}}_3(q)$\n$={L_formula(3):.4f}$')
ax.set_xlabel('$s$', fontsize=12)
ax.set_ylabel('$L(s)$', fontsize=12)
ax.set_title('Near s=2 pole (Res = 1/9)', fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim([0, 5])

# Panel 3: Turbulence Mertens constant
ax = axes[1, 0]
N_vals = np.array([int(10**x) for x in np.linspace(1, 5, 50)])
sigma_vals = []
approx_vals = []
for N in N_vals:
    sigma = sum(zeta_sl(p)/p for p in range(1, N+1))
    sigma_vals.append(sigma)
    approx_vals.append(N/9 + 2*np.log(N) + C_T_formula)

diff_vals = [s - a for s, a in zip(sigma_vals, approx_vals)]
ax.semilogx(N_vals, diff_vals, 'b-', linewidth=2)
ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('$N$', fontsize=12)
ax.set_ylabel(r'$\Sigma(N) - N/9 - 2\ln N - C_T$', fontsize=11)
ax.set_title(f'Turbulence Mertens: $C_T = {C_T_formula:.4f}$', fontsize=10)
ax.grid(True, alpha=0.3)

# Panel 4: L(s) vs Li_s(q) for context
ax = axes[1, 1]
s_range = np.linspace(3, 20, 300)
Li_vals = [Li_polylog(s, q) for s in s_range]
ax.plot(s_range, Li_vals, 'b-', linewidth=2, label=r'$\mathrm{Li}_s(q)$')
# Asymptotic Li_s(q) ~ q as s -> inf
ax.axhline(y=q, color='r', linestyle='--', alpha=0.7, label=f'$q = {q:.4f}$')
ax.axhline(y=Li_polylog(3, q), color='g', linestyle=':', alpha=0.7,
           label=f'$\\mathrm{{Li}}_3(q) = {Li_polylog(3,q):.4f}$')
# Mark pi^2/6 / 2
ax.axhline(y=pi**2/6, color='purple', linestyle='-.', alpha=0.7,
           label=r'$\zeta_R(2) = \pi^2/6$')
ax.set_xlabel('$s$', fontsize=12)
ax.set_ylabel(r'$\mathrm{Li}_s(q)$', fontsize=12)
ax.set_title(r'Polylogarithm $\mathrm{Li}_s(q)$ at $q=(2/3)^{1/3}$', fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('l_function.png', dpi=150, bbox_inches='tight')
print("  Saved: l_function.png")
