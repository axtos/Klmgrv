"""
cascadequant
============
Log-Poisson energy cascade: analytics for turbulence and finance.

Quick start
-----------
    from cascadequant import she_leveque, CascadeFingerprintTest, CascadeVolatility

    mdl = she_leveque()
    print(mdl.zeta([1,2,3,4,5]))          # SL scaling exponents
    print(mdl.hankel_ratio())              # should be 1.0

    # Test real data
    result = CascadeFingerprintTest(log_increments).run()
    print(result.summary())

    # Options pricing
    cv = CascadeVolatility(s=1.0, sigma=0.20)
    iv = cv.implied_vol(100, 105, 0.25)
"""

from .core import (
    CascadeModel,
    she_leveque,
    Q_SL,
    MU_SL,
    LAM_SL,
    EXPERIMENTAL_ZETA,
)
from .statistics import (
    CascadeFingerprintTest,
    FingerprintResult,
    multiscale_fingerprint,
    estimate_scaling_exponents,
)
from .finance import (
    CascadeVolatility,
    fit_to_returns,
    cascade_var,
)

__version__ = "0.1.0"
__author__  = "AxiomBlack Research"

__all__ = [
    # core
    "CascadeModel", "she_leveque",
    "Q_SL", "MU_SL", "LAM_SL", "EXPERIMENTAL_ZETA",
    # statistics
    "CascadeFingerprintTest", "FingerprintResult",
    "multiscale_fingerprint", "estimate_scaling_exponents",
    # finance
    "CascadeVolatility", "fit_to_returns", "cascade_var",
    # version
    "__version__",
]
