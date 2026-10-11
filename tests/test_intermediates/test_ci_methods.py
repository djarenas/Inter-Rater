"""
Tests for ci_methods: percentile_ci, bca_ci, scipy_bca_ci.

These compute confidence intervals from bootstrap distributions. The
percentile method is the simplest (alpha/2 and 1-alpha/2 quantiles).
BCa (bias-corrected and accelerated) adjusts for bias and skewness.

Tests use small deterministic distributions whose percentiles can be
computed by hand.
"""

import numpy as np
import pytest

from interrater.intermediates.ci_methods import (
    percentile_ci,
    bca_ci,
    scipy_bca_ci,
)


# -- Test data --------------------------------------------------------------- #

# 10 values: [1, 2, ..., 10].
# alpha=0.05 -> lower_pct=2.5, upper_pct=97.5.
# numpy linear interpolation on 10 points:
#   2.5th percentile: index = 9 * 0.025 = 0.225 -> 1 + 0.225 = 1.225
#   97.5th percentile: index = 9 * 0.975 = 8.775 -> 9 + 0.775 = 9.775

DIST_10 = np.arange(1, 11, dtype=np.float64)
EXPECTED_LOWER_05 = 1.225
EXPECTED_UPPER_05 = 9.775


# ============================================================================ #
#                              percentile_ci                                   #
# ============================================================================ #


class TestPercentileCi:
    """Percentile confidence interval from bootstrap distribution."""

    def test_known_values(self):
        """Hand-computed percentiles for [1..10] at alpha=0.05."""
        lower, upper = percentile_ci(DIST_10, alpha=0.05)
        assert abs(lower - EXPECTED_LOWER_05) < 1e-12
        assert abs(upper - EXPECTED_UPPER_05) < 1e-12

    def test_narrower_at_higher_alpha(self):
        """alpha=0.10 gives a narrower interval than alpha=0.05."""
        lo_05, hi_05 = percentile_ci(DIST_10, alpha=0.05)
        lo_10, hi_10 = percentile_ci(DIST_10, alpha=0.10)
        assert lo_10 > lo_05
        assert hi_10 < hi_05

    def test_returns_floats(self):
        """Output is a tuple of Python floats, not np.float64."""
        lower, upper = percentile_ci(DIST_10, alpha=0.05)
        assert isinstance(lower, float)
        assert isinstance(upper, float)

    def test_single_value(self):
        """Bootstrap of one value: lower == upper."""
        lower, upper = percentile_ci(np.array([5.0]), alpha=0.05)
        assert lower == upper == 5.0

    def test_alpha_half_gives_50pct_ci(self):
        """alpha=0.5 -> 50% CI: 25th and 75th percentiles."""
        # [1..10], numpy linear interpolation:
        #   25th pct: index = 9 * 0.25 = 2.25 -> 3 + 0.25 = 3.25
        #   75th pct: index = 9 * 0.75 = 6.75 -> 7 + 0.75 = 7.75
        lower, upper = percentile_ci(DIST_10, alpha=0.5)
        assert abs(lower - 3.25) < 1e-12
        assert abs(upper - 7.75) < 1e-12

    def test_empty_raises(self):
        """Empty bootstrap distribution is rejected."""
        with pytest.raises(ValueError, match="cannot be empty"):
            percentile_ci(np.array([]), alpha=0.05)

    @pytest.mark.parametrize("bad_alpha", [0.0, 1.0, -0.1, 1.5])
    def test_alpha_out_of_range(self, bad_alpha):
        """Alpha must be in (0, 1), exclusive."""
        with pytest.raises(ValueError, match="alpha must be in"):
            percentile_ci(DIST_10, alpha=bad_alpha)

    def test_error_prefix(self):
        """Error message starts with [percentile_ci] Error."""
        with pytest.raises(ValueError, match=r"\[percentile_ci\] Error"):
            percentile_ci(np.array([]), alpha=0.05)


# ============================================================================ #
#                                 bca_ci                                       #
# ============================================================================ #


class TestBcaCi:
    """Bias-corrected and accelerated bootstrap CI."""

    def test_symmetric_close_to_percentile(self):
        """For a symmetric distribution centered on point, BCa ~ percentile."""
        # Symmetric around 50: [1, 2, ..., 99].
        dist = np.arange(1, 100, dtype=np.float64)
        point = 50.0
        bca_lo, bca_hi = bca_ci(dist, point, alpha=0.05)
        pct_lo, pct_hi = percentile_ci(dist, alpha=0.05)
        # With symmetric data and point at center, BCa should be close.
        assert abs(bca_lo - pct_lo) < 2.0
        assert abs(bca_hi - pct_hi) < 2.0

    def test_lower_leq_upper(self):
        """Lower bound is at most upper bound."""
        dist = np.random.default_rng(42).normal(0.5, 0.1, size=200)
        point = 0.5
        lower, upper = bca_ci(dist, point, alpha=0.05)
        assert lower <= upper

    def test_returns_floats(self):
        """Output is a tuple of Python floats."""
        lower, upper = bca_ci(DIST_10, point=5.5, alpha=0.05)
        assert isinstance(lower, float)
        assert isinstance(upper, float)

    def test_constant_distribution(self):
        """All identical values: acceleration denominator = 0, still valid."""
        dist = np.full(50, 0.7)
        lower, upper = bca_ci(dist, point=0.7, alpha=0.05)
        # All values identical -> CI collapses to that value
        assert abs(lower - 0.7) < 1e-12
        assert abs(upper - 0.7) < 1e-12

    def test_biased_point_shifts_interval(self):
        """Point below distribution center: BCa shifts interval downward."""
        # Distribution centered around 0.6
        rng = np.random.default_rng(99)
        dist = rng.normal(0.6, 0.05, size=500)
        # Point well below the center -> z0 < 0 -> lower tail gets more weight
        lo_biased, hi_biased = bca_ci(dist, point=0.45, alpha=0.05)
        lo_centered, hi_centered = bca_ci(dist, point=0.6, alpha=0.05)
        # Biased point should shift the interval downward relative to centered
        assert lo_biased < lo_centered

    @pytest.mark.parametrize("bad_alpha", [0.0, 1.0, -0.1, 1.5])
    def test_alpha_out_of_range(self, bad_alpha):
        """Alpha must be in (0, 1), exclusive."""
        with pytest.raises(ValueError, match="alpha must be in"):
            bca_ci(DIST_10, point=5.5, alpha=bad_alpha)

    def test_empty_raises(self):
        """Empty bootstrap distribution is rejected."""
        with pytest.raises(ValueError, match="cannot be empty"):
            bca_ci(np.array([]), point=0.5, alpha=0.05)

    def test_error_prefix(self):
        """Error message starts with [bca_ci] Error."""
        with pytest.raises(ValueError, match=r"\[bca_ci\] Error"):
            bca_ci(np.array([]), point=0.5, alpha=0.05)


# ============================================================================ #
#                              scipy_bca_ci                                    #
# ============================================================================ #


class TestScipyBcaCi:
    """Reference BCa implementation wrapping scipy.stats.bootstrap."""

    def test_returns_tuple_of_floats(self):
        """Happy path: returns (float, float)."""
        dist = np.random.default_rng(7).normal(0.5, 0.1, size=200)
        lower, upper = scipy_bca_ci(dist, point=0.5, alpha=0.05)
        assert isinstance(lower, float)
        assert isinstance(upper, float)

    def test_lower_leq_upper(self):
        """Lower bound is at most upper bound."""
        dist = np.random.default_rng(7).normal(0.5, 0.1, size=200)
        lower, upper = scipy_bca_ci(dist, point=0.5, alpha=0.05)
        assert lower <= upper

    def test_empty_raises(self):
        """Empty bootstrap distribution is rejected."""
        with pytest.raises(ValueError, match="cannot be empty"):
            scipy_bca_ci(np.array([]), point=0.5, alpha=0.05)

    def test_alpha_out_of_range(self):
        """Alpha must be in (0, 1)."""
        with pytest.raises(ValueError, match="alpha must be in"):
            scipy_bca_ci(DIST_10, point=5.5, alpha=0.0)

    def test_error_prefix(self):
        """Error message starts with [scipy_bca_ci] Error."""
        with pytest.raises(ValueError, match=r"\[scipy_bca_ci\] Error"):
            scipy_bca_ci(np.array([]), point=0.5, alpha=0.05)
