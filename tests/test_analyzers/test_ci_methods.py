"""Tests for ci_methods.py -- percentile_ci, bca_ci, scipy_bca_ci.

Deterministic tests only: no bootstrap loops. Tolerance is 1e-12 for
exact computations, relaxed where scipy or floating-point order matters.
"""

from __future__ import annotations

import numpy as np
import pytest

from interrater.analyzers.ci_methods import (
    percentile_ci,
    bca_ci,
    scipy_bca_ci,
)


# ------------------------------------------------------------------ #
# Module-level _ERROR_PREFIX                                           #
# ------------------------------------------------------------------ #


class TestErrorPrefix:
    """_ERROR_PREFIX lives at module level, not inside functions."""

    def test_error_prefix_at_module_level(self):
        import interrater.analyzers.ci_methods as mod

        assert hasattr(mod, "_ERROR_PREFIX")
        assert mod._ERROR_PREFIX == "[ci_methods] Error"

    def test_error_prefix_appears_in_errors(self):
        with pytest.raises(ValueError, match=r"\[ci_methods\] Error"):
            percentile_ci(np.array([]), 0.05)


# ------------------------------------------------------------------ #
# percentile_ci                                                        #
# ------------------------------------------------------------------ #


class TestPercentileCi:
    """Tests for percentile_ci()."""

    def test_known_symmetric_distribution(self):
        """np.arange(1, 101): 2.5th percentile=3.475, 97.5th=97.525."""
        dist = np.arange(1, 101, dtype=float)
        lo, hi = percentile_ci(dist, 0.05)
        assert lo == pytest.approx(3.475, abs=1e-12)
        assert hi == pytest.approx(97.525, abs=1e-12)

    def test_all_same_values(self):
        """Constant distribution -> CI collapses to the single value."""
        dist = np.full(1000, 0.75)
        lo, hi = percentile_ci(dist, 0.05)
        assert lo == pytest.approx(0.75, abs=1e-12)
        assert hi == pytest.approx(0.75, abs=1e-12)

    def test_alpha_10_percent(self):
        """alpha=0.10 -> 5th and 95th percentile."""
        dist = np.arange(1, 101, dtype=float)
        lo, hi = percentile_ci(dist, 0.10)
        assert lo == pytest.approx(5.95, abs=1e-12)
        assert hi == pytest.approx(95.05, abs=1e-12)

    def test_small_alpha_gives_wider_ci(self):
        """Smaller alpha -> wider CI."""
        dist = np.arange(1, 1001, dtype=float)
        lo_05, hi_05 = percentile_ci(dist, 0.05)
        lo_01, hi_01 = percentile_ci(dist, 0.01)
        assert lo_01 < lo_05
        assert hi_01 > hi_05

    def test_single_element(self):
        """Single element: CI collapses."""
        dist = np.array([0.5])
        lo, hi = percentile_ci(dist, 0.05)
        assert lo == pytest.approx(0.5, abs=1e-12)
        assert hi == pytest.approx(0.5, abs=1e-12)

    def test_returns_tuple_of_floats(self):
        dist = np.arange(1, 101, dtype=float)
        result = percentile_ci(dist, 0.05)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)

    # --- Guards ---

    def test_empty_distribution_raises(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            percentile_ci(np.array([]), 0.05)

    def test_alpha_zero_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            percentile_ci(np.arange(100.0), 0.0)

    def test_alpha_one_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            percentile_ci(np.arange(100.0), 1.0)

    def test_alpha_negative_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            percentile_ci(np.arange(100.0), -0.1)

    def test_alpha_greater_than_one_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            percentile_ci(np.arange(100.0), 1.5)


# ------------------------------------------------------------------ #
# bca_ci                                                               #
# ------------------------------------------------------------------ #


class TestBcaCi:
    """Tests for bca_ci() -- custom BCa implementation."""

    def test_known_symmetric_close_to_percentile(self):
        """Symmetric dist -> BCa ≈ percentile (bias & skewness ≈ 0)."""
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.6, scale=0.1, size=5000)
        point = 0.6

        bca_lo, bca_hi = bca_ci(dist, point, 0.05)
        pct_lo, pct_hi = percentile_ci(dist, 0.05)

        # Should be close for a symmetric distribution
        assert bca_lo == pytest.approx(pct_lo, abs=0.02)
        assert bca_hi == pytest.approx(pct_hi, abs=0.02)

    def test_skewed_distribution_shifts_ci(self):
        """Skewed dist -> BCa adjusts bounds relative to percentile."""
        rng = np.random.default_rng(42)
        # Right-skewed distribution (exponential shifted)
        dist = rng.exponential(scale=0.1, size=5000) + 0.3
        point = 0.3 + 0.1  # Mean of exponential(0.1) is 0.1

        bca_lo, bca_hi = bca_ci(dist, point, 0.05)
        pct_lo, pct_hi = percentile_ci(dist, 0.05)

        # Both methods should give valid CIs
        assert bca_lo < bca_hi
        assert pct_lo < pct_hi
        # BCa should differ from percentile for skewed data
        # (not necessarily in a specific direction, just not identical)
        assert not (
            bca_lo == pytest.approx(pct_lo, abs=1e-12)
            and bca_hi == pytest.approx(pct_hi, abs=1e-12)
        )

    def test_constant_distribution_fallback(self):
        """All-same distribution -> acceleration=0, z0≈0, close to percentile."""
        dist = np.full(500, 0.75)
        point = 0.75
        lo, hi = bca_ci(dist, point, 0.05)
        # Should degenerate to the constant
        assert lo == pytest.approx(0.75, abs=1e-6)
        assert hi == pytest.approx(0.75, abs=1e-6)

    def test_returns_tuple_of_floats(self):
        rng = np.random.default_rng(42)
        dist = rng.normal(size=200)
        result = bca_ci(dist, 0.0, 0.05)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)

    def test_does_not_raise_with_extreme_inputs(self):
        """BCa should not raise even with extreme point/distribution mismatch.

        When the point is far from the distribution mass, BCa either
        falls back to percentile or produces extreme (but finite) bounds.
        Either way it should not crash.
        """
        dist = np.full(100, 0.99)
        dist[0] = 0.98
        point = 0.01
        lo, hi = bca_ci(dist, point, 0.05)
        # Should complete without error; bounds are numeric
        assert isinstance(lo, float)
        assert isinstance(hi, float)

    # --- Guards ---

    def test_empty_distribution_raises(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            bca_ci(np.array([]), 0.5, 0.05)

    def test_alpha_zero_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            bca_ci(np.arange(100.0), 0.5, 0.0)

    def test_alpha_one_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            bca_ci(np.arange(100.0), 0.5, 1.0)


# ------------------------------------------------------------------ #
# scipy_bca_ci                                                         #
# ------------------------------------------------------------------ #


class TestScipyBcaCi:
    """Tests for scipy_bca_ci() -- reference implementation."""

    def test_runs_without_error(self):
        """Smoke test: scipy BCa runs and returns a valid tuple."""
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=500)
        point = 0.5
        lo, hi = scipy_bca_ci(dist, point, 0.05)
        assert lo < hi
        assert isinstance(lo, float)
        assert isinstance(hi, float)

    def test_returns_finite_values(self):
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.7, scale=0.05, size=1000)
        lo, hi = scipy_bca_ci(dist, 0.7, 0.05)
        assert np.isfinite(lo)
        assert np.isfinite(hi)

    # --- Guards ---

    def test_empty_distribution_raises(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            scipy_bca_ci(np.array([]), 0.5, 0.05)

    def test_alpha_zero_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            scipy_bca_ci(np.arange(100.0), 0.5, 0.0)

    def test_alpha_one_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            scipy_bca_ci(np.arange(100.0), 0.5, 1.0)


# ------------------------------------------------------------------ #
# Cross-method consistency                                             #
# ------------------------------------------------------------------ #


class TestCrossMethodConsistency:
    """Sanity checks across the three CI methods."""

    def test_all_methods_return_lower_le_upper(self):
        """Every method: lower <= upper."""
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=1000)
        point = 0.5

        lo_p, hi_p = percentile_ci(dist, 0.05)
        lo_b, hi_b = bca_ci(dist, point, 0.05)
        lo_s, hi_s = scipy_bca_ci(dist, point, 0.05)

        assert lo_p <= hi_p
        assert lo_b <= hi_b
        assert lo_s <= hi_s

    def test_wider_alpha_gives_narrower_ci_percentile(self):
        """alpha=0.10 gives narrower CI than alpha=0.05."""
        dist = np.arange(1, 1001, dtype=float)
        lo_05, hi_05 = percentile_ci(dist, 0.05)
        lo_10, hi_10 = percentile_ci(dist, 0.10)
        assert (hi_05 - lo_05) > (hi_10 - lo_10)

    def test_wider_alpha_gives_narrower_ci_bca(self):
        """alpha=0.10 gives narrower CI than alpha=0.05 for BCa."""
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=2000)
        point = 0.5
        lo_05, hi_05 = bca_ci(dist, point, 0.05)
        lo_10, hi_10 = bca_ci(dist, point, 0.10)
        assert (hi_05 - lo_05) > (hi_10 - lo_10)
