"""Tests for cascadequant.statistics."""
import numpy as np
import pytest
from cascadequant.statistics import (
    _log_cumulants,
    CascadeFingerprintTest,
    FingerprintResult,
)
from cascadequant.core import she_leveque, MU_SL, LAM_SL


def _generate_log_poisson(n=10000, s=2.0, rng=None):
    """Simulate log-Poisson increments matching the SL cascade."""
    if rng is None:
        rng = np.random.default_rng(0)
    mdl = she_leveque()
    N = rng.poisson(mdl.lam * s, n)
    return -s * mdl._h_inf + N * (-mdl.mu)


class TestLogCumulants:
    def test_gaussian_k3_k4_near_zero(self):
        rng = np.random.default_rng(1)
        x = rng.normal(0, 1, 50000)
        k = _log_cumulants(x, 4)
        assert abs(k[3]) < 0.05   # skewness near 0
        assert abs(k[4]) < 0.15   # excess kurtosis near 0

    def test_cumulant_count(self):
        x = np.arange(100, dtype=float)
        k = _log_cumulants(x, 4)
        assert set(k.keys()) >= {1, 2, 3, 4}


class TestCascadeFingerprintTest:
    def test_log_poisson_hankel_near_one(self):
        y = _generate_log_poisson(n=20000, s=2.0)
        result = CascadeFingerprintTest(y, bootstrap_B=100).run()
        # Should be close to 1 within 3-sigma
        assert abs(result.hankel_ratio - 1.0) < 3 * result.hankel_std + 0.3

    def test_gaussian_hankel_nan_or_different(self):
        rng = np.random.default_rng(2)
        x = rng.normal(0, 1, 5000)
        result = CascadeFingerprintTest(x, bootstrap_B=50).run()
        # Gaussian has near-zero k4 → ratio may be nan or very large
        assert np.isnan(result.hankel_ratio) or abs(result.hankel_ratio - 1.0) > 0.3

    def test_result_fields(self):
        y = _generate_log_poisson(n=5000)
        result = CascadeFingerprintTest(y, bootstrap_B=50).run()
        assert isinstance(result, FingerprintResult)
        assert result.n_samples == 5000
        assert result.mu_implied is not None
        assert result.lam_implied is not None

    def test_mu_implied_close_to_sl(self):
        y = _generate_log_poisson(n=50000, s=3.0)
        result = CascadeFingerprintTest(y, bootstrap_B=50).run()
        # mu_implied = |k3/k2| should be close to MU_SL
        if result.mu_implied is not None:
            assert abs(result.mu_implied / MU_SL - 1.0) < 0.15

    def test_summary_string(self):
        y = _generate_log_poisson(n=2000)
        result = CascadeFingerprintTest(y, bootstrap_B=20).run()
        summary = result.summary()
        assert "Hankel ratio" in summary
        assert "κ₂" in summary

    def test_small_sample_note(self):
        y = _generate_log_poisson(n=200)
        result = CascadeFingerprintTest(y, bootstrap_B=20).run()
        assert any("Small sample" in note for note in result.notes)
