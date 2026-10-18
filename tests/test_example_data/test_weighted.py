"""
Integration tests for weighted analysis (linear + quadratic weights).

Tests that distance-based weighting changes the point estimate and exercises
both WeightableGroupConfig (Fleiss, B-P, Krippendorff) and
WeightableIndividualConfig (Cohen).

Weights only matter when K >= 3 -- with binary categories, linear and
quadratic collapse to identity (there is only one distance). All datasets
below have K=3.

Point estimates are asserted to 3 decimal places (abs=5e-4) and CI
structure (lower < point < upper) is checked as a sanity guard, using the
percentile method at B=250.

Containment is a real check under percentile and not an invariant under
BCa -- see the module docstring of test_reference_data.py, which also
holds the one end-to-end BCa test.

Every config here comes from a conftest fixture. These tests used to build
their own, which is how they kept BCa at B=50 after the shared config had
moved on.
"""

import pytest

from interrater import (
    FleissKappaAnalyzer,
    KrippendorffAlphaAnalyzer,
    CohenKappaAnalyzer,
    BrennanPredigerAnalyzer,
)
from interrater.analyzers.configs import (
    FleissKappaConfig,
    KrippendorffAlphaConfig,
    CohenKappaConfig,
    BrennanPredigerConfig,
)

# Tolerance: 3 decimal places.
TOL = 5e-4


def _assert_brackets_point(result) -> None:
    """
    Assert the point estimate falls strictly inside its own interval.

    An interval is a BootstrapCIResult, not a pair of floats, so this reads
    ``.lower`` and ``.upper`` rather than ``[0]`` and ``[1]``. It also checks
    the interval has bounds at all: an interval whose every draw failed
    carries NaN, and ``nan < point`` is False, so the comparison alone would
    fail with an unhelpful message instead of naming the real problem.
    """
    ci = result.confidence_interval
    assert ci is not None, "no interval was computed"
    assert ci.is_defined, f"interval has no bounds: {ci}"
    assert ci.lower < result.point < ci.upper, (
        f"point {result.point:.4f} outside its own interval {ci}"
    )


# ── SR1 prelim (ratings, 4 raters, K=3, balanced, ordinal -1/0/1) ────── #
#
# Two weight types on the same data, exercising WeightableGroupConfig
# with two different analyzers.  Categories are ordinal (-1, 0, 1), so
# -1/0 is closer than -1/1 under distance-based weights.
#
# Unweighted reference values (from test_research_data.py):
#   Fleiss kappa  = 0.493
#   Krippendorff  = 0.493  (approx; 0.493269 from ratings pathway)


class TestSR1PrelimWeighted:
    """Weighted analysis on SR1 prelim (199 subjects, 4 raters, K=3)."""

    def test_fleiss_kappa_quadratic(self, sr1_prelim_ratings, weighted_config, fast_fleiss_config):
        """Fleiss' kappa with quadratic weights = 0.593.

        Quadratic weights (w_ij = 1 - (i-j)^2 / (K-1)^2) upweight partial
        agreements, increasing the estimate from 0.493 to 0.593.
        Exercises WeightableGroupConfig.
        """
        config = weighted_config(FleissKappaConfig, "quadratic")
        analyzer = FleissKappaAnalyzer(sr1_prelim_ratings, config)
        result = analyzer.estimate()

        # Point estimate
        assert result.point == pytest.approx(0.593, abs=TOL)

        # CI structure
        _assert_brackets_point(result)

        # Structural: weights changed the estimate
        result_uw = FleissKappaAnalyzer(sr1_prelim_ratings, fast_fleiss_config).estimate()
        assert result.point != pytest.approx(result_uw.point, abs=1e-6)

        # Result reports the correct weights type
        assert result.weights_type == "quadratic"

    def test_krippendorff_alpha_linear(self, sr1_prelim_ratings, weighted_config, fast_krippendorff_config):
        """Krippendorff's alpha with linear weights = 0.555.

        Linear weights (w_ij = 1 - |i-j| / (K-1)) are less aggressive
        than quadratic, but still upweight partial agreements.
        Different weight type on same data as Fleiss quadratic above.
        Exercises WeightableGroupConfig.
        """
        config = weighted_config(KrippendorffAlphaConfig, "linear")
        analyzer = KrippendorffAlphaAnalyzer(sr1_prelim_ratings, config)
        result = analyzer.estimate()

        # Point estimate
        assert result.point == pytest.approx(0.555, abs=TOL)

        # CI structure
        _assert_brackets_point(result)

        # Structural: weights changed the estimate
        result_uw = KrippendorffAlphaAnalyzer(sr1_prelim_ratings, fast_krippendorff_config).estimate()
        assert result.point != pytest.approx(result_uw.point, abs=1e-6)

        # Result reports the correct weights type
        assert result.weights_type == "linear"


# ── Gwet cont3x3 (ratings, 2 raters, K=3, balanced) ─────────────────── #
#
# Exercises WeightableIndividualConfig (Cohen is an individual-marginal
# analyzer).  irrCAC may provide an external cross-check.
#
# Unweighted reference (from test_reference_data.py):
#   Cohen kappa = 0.796


class TestGwetCont3x3Weighted:
    """Weighted analysis on Gwet cont3x3 (100 subjects, 2 raters, K=3)."""

    def test_cohen_kappa_quadratic(self, gwet_ratings, weighted_config, fast_cohen_config):
        """Cohen's kappa with quadratic weights = 0.892.

        Quadratic weights increase the estimate from 0.796 to 0.892 because
        adjacent-category disagreements get partial credit.
        Exercises WeightableIndividualConfig.
        """
        config = weighted_config(CohenKappaConfig, "quadratic")
        analyzer = CohenKappaAnalyzer(gwet_ratings, config)
        result = analyzer.estimate()

        # Point estimate
        assert result.point == pytest.approx(0.892, abs=TOL)

        # CI structure
        _assert_brackets_point(result)

        # Structural: weights changed the estimate
        result_uw = CohenKappaAnalyzer(gwet_ratings, fast_cohen_config).estimate()
        assert result.point != pytest.approx(result_uw.point, abs=1e-6)

        # Result reports the correct weights type
        assert result.weights_type == "quadratic"


# ── Synthetic 2-rater (ratings, 2 raters, K=3, balanced, ordinal) ────── #
#
# Second WeightableGroupConfig analyzer with weights, different dataset.
#
# Unweighted reference (from test_reference_data.py):
#   Brennan-Prediger = 0.670  (synthetic_2rater B-P via Scott's equivalence)


class TestSynthetic2RaterWeighted:
    """Weighted analysis on Synthetic 2-rater (80 subjects, 2 raters, K=3)."""

    def test_brennan_prediger_linear(self, synthetic_ratings, weighted_config, fast_bp_config):
        """Brennan-Prediger with linear weights = 0.651.

        Linear weights lower the estimate from 0.670 to 0.651 here because
        the B-P chance correction (Pe = 1/K) interacts differently with
        partial credit than the observed-agreement correction in kappa.
        Exercises WeightableGroupConfig.
        """
        config = weighted_config(BrennanPredigerConfig, "linear")
        analyzer = BrennanPredigerAnalyzer(synthetic_ratings, config)
        result = analyzer.estimate()

        # Point estimate
        assert result.point == pytest.approx(0.651, abs=TOL)

        # CI structure
        _assert_brackets_point(result)

        # Structural: weights changed the estimate
        result_uw = BrennanPredigerAnalyzer(synthetic_ratings, fast_bp_config).estimate()
        assert result.point != pytest.approx(result_uw.point, abs=1e-6)

        # Result reports the correct weights type
        assert result.weights_type == "linear"
