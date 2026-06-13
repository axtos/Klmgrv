# cascadequant

**Log-Poisson energy cascade analytics for turbulence and financial volatility.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## What it does

The [She-Lévêque (1994)](https://doi.org/10.1103/PhysRevLett.72.336) turbulence cascade assigns a *compound Poisson distribution* to velocity increments:

```
log|δu_r| = −s/9 − N(s)·μ,   N(s) ~ Poisson(2s),   μ = (1/3)ln(3/2)
```

This package implements:

| Module | What it provides |
|--------|------------------|
| `core` | CascadeModel, SL exponents ζ_p, multifractal spectrum f(h), partition function Z_dyn, S-duality |
| `statistics` | Hankel fingerprint test κ₃²/(κ₂κ₄) = 1, bootstrap CI, multiscale test |
| `finance` | CascadeVolatility, compound Poisson option pricing, risk-neutral calibration |
| `plots` | 6 publication-quality figures |

## Quick start

```bash
pip install -e ".[dev]"
```

```python
from cascadequant import she_leveque, CascadeFingerprintTest, CascadeVolatility

# SL model
mdl = she_leveque()
print(mdl.zeta([1,2,3,4,5,6]))     # [0.364, 0.696, 1.000, 1.280, 1.538, 1.778]
print(mdl.hankel_ratio())           # 1.0 (exact)

# Fingerprint test on your data
result = CascadeFingerprintTest(log_increments).run()
print(result.summary())
# Hankel ratio κ₃²/(κ₂κ₄) = 0.97 ± 0.12  → log-Poisson ✓

# Options pricing
cv = CascadeVolatility(s=0.25, sigma=0.10)
iv = cv.implied_vol(100, 95, 0.25, 0.0, 'put', sigma_spread=0.10)
print(f"Cascade implied vol: {iv*100:.1f}%")
```

## The Hankel identity

The parameter-free discriminator between cascade models:

| Model | κ₃²/(κ₂κ₄) |
|-------|-------------|
| Lognormal (K62) | 0 |
| **Log-Poisson (SL)** | **1.0 (exact)** |
| Log-stable | ∞ |

## Figures

```python
from cascadequant.plots import generate_all_figures
generate_all_figures("figures/")
# → fig1_zeta.png, fig2_spectrum.png, fig3_pdf.png,
#   fig4_hankel.png, fig5_smile.png, fig6_partition.png
```

## Tests

```bash
pytest tests/ -v
# 45 passed
```

## References

- She & Lévêque (1994). *Phys. Rev. Lett.* 72:336
- Iyer, Sreenivasan & Yeung (2020). *Phys. Rev. Fluids* 5:054605
- See `articles/arxiv_preprint.md` for full mathematical derivations

## License

MIT — see [LICENSE](LICENSE)
