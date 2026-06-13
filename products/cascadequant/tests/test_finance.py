"""Tests for cascadequant.finance."""
import math
import numpy as np
import pytest
from cascadequant.finance import (
    CascadeVolatility,
    fit_to_returns,
    cascade_var,
    _bs_price,
    _bs_implied_vol,
)


class TestBlackScholes:
    def test_call_put_parity(self):
        S, K, T, r, sig = 100, 100, 0.5, 0.05, 0.20
        call = _bs_price(S, K, T, r, sig, "call")
        put  = _bs_price(S, K, T, r, sig, "put")
        # C - P = S - K*e^{-rT}
        parity = S - K * math.exp(-r * T)
        assert abs(call - put - parity) < 1e-10

    def test_implied_vol_roundtrip(self):
        S, K, T, r, sig = 100, 105, 0.25, 0.0, 0.18
        px = _bs_price(S, K, T, r, sig, "call")
        iv = _bs_implied_vol(px, S, K, T, r, "call")
        assert abs(iv - sig) < 1e-5

    def test_deep_itm_call(self):
        px = _bs_price(100, 10, 1.0, 0.0, 0.20, "call")
        assert abs(px - 90.0) < 1.0

    def test_zero_time(self):
        px = _bs_price(100, 100, 0.0, 0.0, 0.20, "call")
        assert px == 0.0


class TestCascadeVolatility:
    def setup_method(self):
        self.cv = CascadeVolatility(s=1.0, sigma=0.20)

    def test_return_pdf_grid_weights_sum_to_one(self):
        _, weights = self.cv.return_pdf_grid()
        assert abs(weights.sum() - 1.0) < 1e-10

    def test_return_cdf_monotone(self):
        y_arr = np.linspace(-5, 5, 50)
        cdf = self.cv.return_cdf(y_arr, sigma_spread=0.1)
        assert np.all(np.diff(cdf) >= 0)

    def test_return_cdf_bounds(self):
        cdf_lo = self.cv.return_cdf(np.array([-100.0]), sigma_spread=0.1)[0]
        cdf_hi = self.cv.return_cdf(np.array([100.0]),  sigma_spread=0.1)[0]
        assert cdf_lo < 0.01
        assert cdf_hi > 0.99

    def test_option_price_call_positive(self):
        px = self.cv.option_price(100, 100, 0.25, 0.0, "call", sigma_spread=0.15)
        assert px > 0

    def test_put_call_parity_approx(self):
        S, K, T, r = 100, 100, 0.25, 0.0
        sigma_spread = 0.15
        call = self.cv.option_price(S, K, T, r, "call", sigma_spread)
        put  = self.cv.option_price(S, K, T, r, "put",  sigma_spread)
        # Approximate parity: C - P ≈ S - K (r=0)
        assert abs(call - put - (S - K)) < 5.0

    def test_implied_vol_otm_call(self):
        iv = self.cv.implied_vol(100, 110, 0.25, 0.0, "call", sigma_spread=0.15)
        assert 0.01 < iv < 2.0

    def test_var_99_negative(self):
        var_val = self.cv.var(0.99)
        assert var_val < 0

    def test_es_less_than_var(self):
        var_val = self.cv.var(0.99)
        es_val  = self.cv.expected_shortfall(0.99)
        assert es_val <= var_val


class TestFitToReturns:
    def setup_method(self):
        rng = np.random.default_rng(42)
        from cascadequant.core import she_leveque
        mdl = she_leveque()
        s = 1.5
        N = rng.poisson(mdl.lam * s, 8000)
        self.returns = -s * mdl._h_inf + N * (-mdl.mu) + rng.normal(0, 0.01, 8000)

    def test_returns_dict_keys(self):
        result = fit_to_returns(self.returns)
        for key in ["q", "lam", "mu", "s", "sigma_base", "hankel_ratio"]:
            assert key in result

    def test_fit_mu_reasonable(self):
        result = fit_to_returns(self.returns)
        from cascadequant.core import MU_SL
        assert 0.01 < result["mu"] < 0.5

    def test_fit_q_in_range(self):
        result = fit_to_returns(self.returns)
        assert 0.5 < result["q"] < 1.0


class TestCascadeVar:
    def setup_method(self):
        rng = np.random.default_rng(7)
        self.returns = rng.standard_t(df=4, size=5000) * 0.01

    def test_cascade_var_dict(self):
        result = cascade_var(self.returns)
        for key in ["cascade_var_99", "historical_var_99", "cascade_es_99", "historical_es_99"]:
            assert key in result

    def test_var_less_than_es(self):
        result = cascade_var(self.returns)
        assert result["cascade_es_99"] <= result["cascade_var_99"]
        assert result["historical_es_99"] <= result["historical_var_99"]
