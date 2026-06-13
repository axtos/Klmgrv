"""Tests for cascadequant.core."""
import math
import numpy as np
import pytest
from cascadequant.core import CascadeModel, she_leveque, Q_SL, MU_SL, LAM_SL


class TestSLParameters:
    def test_q_value(self):
        assert abs(Q_SL - (2/3)**(1/3)) < 1e-15

    def test_mu_value(self):
        assert abs(MU_SL - (-math.log(Q_SL))) < 1e-15

    def test_lam_value(self):
        assert LAM_SL == 2.0


class TestCascadeModel:
    def setup_method(self):
        self.mdl = she_leveque()

    def test_zeta_3_equals_one(self):
        assert abs(self.mdl.zeta(3.0) - 1.0) < 1e-14

    def test_zeta_array(self):
        p = np.array([1.0, 2.0, 3.0])
        z = self.mdl.zeta(p)
        assert z.shape == (3,)
        assert abs(z[2] - 1.0) < 1e-14

    def test_zeta_increasing(self):
        p = np.linspace(0.1, 10, 100)
        z = self.mdl.zeta(p)
        assert np.all(np.diff(z) > 0)

    def test_zeta_zero(self):
        assert abs(self.mdl.zeta(0.0)) < 1e-15

    def test_multifractal_h_min(self):
        mdl = self.mdl
        # f(h) at h_min + epsilon should be close to 1 (vortex filaments)
        h = mdl._h_inf + 1e-6
        f = mdl.multifractal_spectrum(h)
        assert abs(f - 1.0) < 0.01

    def test_multifractal_h_max_approaches_d(self):
        # At h_max, f -> 3 (space-filling)
        mdl = self.mdl
        h = mdl._h_inf + mdl.lam * mdl.mu - 1e-6
        f = mdl.multifractal_spectrum(h)
        assert abs(f - 3.0) < 0.01

    def test_multifractal_outside_support(self):
        f = self.mdl.multifractal_spectrum(-10.0)
        assert f == -np.inf or np.isneginf(f)

    def test_log_partition_negative(self):
        # log Z_dyn < 0 for all p > 0
        p_arr = np.array([0.5, 1.0, 2.0, 5.0])
        lz = self.mdl.log_partition(p_arr)
        assert np.all(lz < 0)

    def test_log_partition_sdual_vs_series(self):
        # Both methods agree to ~10^{-10} for p >= 0.5
        p_arr = np.array([0.5, 1.0, 2.0])
        lz_sd = self.mdl.log_partition(p_arr, method="sdual")
        lz_ser = self.mdl.log_partition(p_arr, method="series")
        np.testing.assert_allclose(lz_sd, lz_ser, rtol=0, atol=1e-9)

    def test_increment_cdf_bounds(self):
        mdl = self.mdl
        assert mdl.increment_cdf(-100.0, 1.0) == 0.0
        assert abs(mdl.increment_cdf(100.0, 1.0) - 1.0) < 1e-10

    def test_increment_cdf_monotone(self):
        mdl = self.mdl
        ys = np.linspace(-5, 2, 20)
        cdfs = [mdl.increment_cdf(y, 1.0) for y in ys]
        assert all(cdfs[i] <= cdfs[i+1] for i in range(len(cdfs)-1))

    def test_log_cumulants_hankel(self):
        # κ₃² = κ₂·κ₄ exactly (up to float64 rounding ~1e-20)
        k = self.mdl.log_cumulants(4, s=1.0)
        k2, k3, k4 = k[1], k[2], k[3]
        assert abs(k3**2 - k2*k4) < 1e-19

    def test_hankel_ratio_exact_one(self):
        assert abs(self.mdl.hankel_ratio(1.0) - 1.0) < 1e-14
        assert abs(self.mdl.hankel_ratio(2.0) - 1.0) < 1e-14

    def test_h_inf_one_ninth(self):
        assert abs(self.mdl._h_inf - 1.0/9.0) < 1e-14

    def test_functional_equation(self):
        # Z_dyn(4π²/(μ²p)) = exp(π²/(6μp) - μp/24) * sqrt(μp/2π) * Z_dyn(p)
        mdl = self.mdl
        mu = mdl.mu
        p = 1.5
        p_dual = 4 * math.pi**2 / (mu**2 * p)

        lhs = mdl.log_partition(p_dual)
        rhs = (math.pi**2 / (6*mu*p) - mu*p/24
               + 0.5*math.log(mu*p / (2*math.pi))
               + mdl.log_partition(p))
        assert abs(lhs - rhs) < 1e-10


class TestExperimentalZeta:
    def test_zeta_3_exact(self):
        from cascadequant.core import EXPERIMENTAL_ZETA
        assert EXPERIMENTAL_ZETA[3] == (1.000, 0.000)

    def test_sl_within_experiment(self):
        from cascadequant.core import EXPERIMENTAL_ZETA
        mdl = she_leveque()
        for p, (z_exp, z_err) in EXPERIMENTAL_ZETA.items():
            z_sl = float(mdl.zeta(float(p)))
            # SL should be within 3-sigma of experiment for p=1..6
            if p <= 6:
                assert abs(z_sl - z_exp) < 3 * z_err + 0.01, (
                    f"p={p}: SL={z_sl:.4f}, exp={z_exp:.4f}±{z_err:.4f}"
                )
