"""
Integration tests for weighted analysis (linear + quadratic weights).

Tests that distance-based weighting changes the point estimate and exercises
both WeightableGroupConfig (Fleiss, B-P, Krippendorff) and
WeightableIndividualConfig (Cohen).

Weights only matter when K >= 3 -- with binary categories, linear and
quadratic collapse to identity (there is only one distance). All datasets
below have K=3.

Point estimates are asserted to 3 decimal places (abs=5e-4) and CI
structure (lower < point < upper) is checked as a sanity guard with
n_resamples=50.
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

    def test_fleiss_kappa_quadratic(self, sr1_prelim_ratings):
        """Fleiss' kappa with quadratic weights = 0.593.

        Quadratic weights (w_ij = 1 - (i-j)^2 / (K-1)^2) upweight partial
        agreements, increasing the estimate from 0.493 to 0.593.
        Exercises WeightableGroupConfig.
        """
        config = FleissKappaConfig(
            n_resamples=50, random_state=42, ci_method="bca",
            weights_type="quadratic",
        )
        analyzer = FleissKappaAnalyzer(sr1_prelim_ratings, config)
        result = analyzer.estimate()

        # Point estimate
        assert result.point == pytest.approx(0.593, abs=TOL)

        # CI structure
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

        # Structural: weights changed the estimate
        config_uw = FleissKappaConfig(
            n_resamples=50, random_state=42, ci_method="bca",
        )
        result_uw = FleissKappaAnalyzer(sr1_prelim_ratings, config_uw).estimate()
        assert result.point != pytest.approx(result_uw.point, abs=1e-6)

        # Result reports the correct weights type
        assert result.weights_type == "quadratic"

    def test_krippendorff_alpha_linear(self, sr1_prelim_ratings):
        """Krippendorff's alpha with linear weights = 0.555.

        Linear weights (w_ij = 1 - |i-j| / (K-1)) are less aggressive
        than quadratic, but still upweight partial agreements.
        Different weight type on same data as Fleiss quadratic above.
        Exercises WeightableGroupConfig.
        """
        config = KrippendorffAlphaConfig(
            n_resamples=50, random_state=42, ci_method="bca",
            weights_type="linear",
        )
        analyzer = KrippendorffAlphaAnalyzer(sr1_prelim_ratings, config)
        result = analyzer.estimate()

        # Point estimate
        assert result.point == pytest.approx(0.555, abs=TOL)

        # CI structure
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

        # Structural: weights changed the estimate
        config_uw = KrippendorffAlphaConfig(
            n_resamples=50, random_state=42, ci_method="bca",
        )
        result_uw = KrippendorffAlphaAnalyzer(sr1_prelim_ratings, config_uw).estimate()
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

    def test_cohen_kappa_quadratic(self, gwet_ratings):
        """Cohen's kappa with quadratic weights = 0.892.

        Quadratic weights increase the estimate from 0.796 to 0.892 because
        adjacent-category disagreements get partial credit.
        Exercises WeightableIndividualConfig.
        """
        config = CohenKappaConfig(
            n_resamples=50, random_state=42, ci_method="bca",
            weights_type="quadratic",
        )
        analyzer = CohenKappaAnalyzer(gwet_ratings, config)
        result = analyzer.estimate()

        # Point estimate
        assert result.point == pytest.approx(0.892, abs=TOL)

        # CI structure
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

        # Structural: weights changed the estimate
        config_uw = CohenKappaConfig(
            n_resamples=50, random_state=42, ci_method="bca",
        )
        result_uw = CohenKappaAnalyzer(gwet_ratings, config_uw).estimate()
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

    def test_brennan_prediger_linear(self, synthetic_ratings):
        """Brennan-Prediger with linear weights = 0.651.

        Linear weights lower the estimate from 0.670 to 0.651 here because
        the B-P chance correction (Pe = 1/K) interacts differently with
        partial credit than the observed-agreement correction in kappa.
        Exercises WeightableGroupConfig.
        """
        config = BrennanPredigerConfig(
            n_resamples=50, random_state=42, ci_method="bca",
            weights_type="linear",
        )
        analyzer = BrennanPredigerAnalyzer(synthetic_ratings, config)
        result = analyzer.estimate()

        # Point estimate
        assert result.point == pytest.approx(0.651, abs=TOL)

        # CI structure
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

        # Structural: weights changed the estimate
        config_uw = BrennanPredigerConfig(
            n_resamples=50, random_state=42, ci_method="bca",
        )
        result_uw = BrennanPredigerAnalyzer(synthetic_ratings, config_uw).estimate()
        assert result.point != pytest.approx(result_uw.point, abs=1e-6)

        # Result reports the correct weights type
        assert result.weights_type == "linear"
