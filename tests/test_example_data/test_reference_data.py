"""
Integration tests against datasets with externally published reference values.

Each test asserts the point estimate to 3 decimal places (abs=5e-4) and
checks CI structure (lower < point < upper) as a sanity guard, using the
percentile method at B=250.

The method matters to that guard. Under percentile the point estimate
always falls inside its own interval, so asserting it is a real check:
a failure means the bootstrap drew from a different set of subjects than
the point estimate was computed from. Under BCa it is not an invariant --
strong bias correction can place the interval entirely to one side of the
point, legitimately -- so the same assertion there would be asserting
something allowed to be false. BCa is covered separately below.

Datasets:
    Fleiss 1971     -- nij pathway, 6 raters, 5 categories, balanced
    Gwet cont3x3   -- ratings + nij pathways, 2 raters, 3 categories, balanced
    Synthetic 2-rater -- ratings pathway, 2 raters, 3 categories, balanced
"""

import pytest

from interrater import (
    FleissKappaAnalyzer,
    GwetAC1Analyzer,
    BrennanPredigerAnalyzer,
    KrippendorffAlphaAnalyzer,
    CohenKappaAnalyzer,
    ScottsAnalyzer,
)

# Tolerance: 3 decimal places.  Looser than the unit test tolerance (1e-12)
# because these are end-to-end checks, not formula-level assertions.
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


# ── Fleiss 1971 (nij pathway, 6 raters, K=5, balanced) ──────────────────── #
#
# Source: Fleiss, J. L. (1971). Measuring nominal scale agreement among
#         many raters. Psychological Bulletin, 76(5), 378–382.


class TestFleiss1971:
    """Fleiss 1971 psychiatric diagnosis dataset (30 subjects, 6 raters, K=5)."""

    def test_fleiss_kappa(self, fleiss1971_nij, fast_fleiss_config):
        """Fleiss' kappa = 0.430 (Fleiss 1971 original)."""
        analyzer = FleissKappaAnalyzer.build_from_nij_matrix(
            fleiss1971_nij, fast_fleiss_config,
        )
        result = analyzer.estimate_from_nij_matrix()

        assert result.point == pytest.approx(0.430, abs=TOL)
        _assert_brackets_point(result)

    def test_gwet_ac1(self, fleiss1971_nij, fast_gwet_config):
        """Gwet AC1 = 0.448 (V2 validated)."""
        analyzer = GwetAC1Analyzer.build_from_nij_matrix(
            fleiss1971_nij, fast_gwet_config,
        )
        result = analyzer.estimate_from_nij_matrix()

        assert result.point == pytest.approx(0.448, abs=TOL)
        _assert_brackets_point(result)

    def test_brennan_prediger(self, fleiss1971_nij, fast_bp_config):
        """Brennan-Prediger = 0.444 (V2 validated, Pe = 1/5)."""
        analyzer = BrennanPredigerAnalyzer.build_from_nij_matrix(
            fleiss1971_nij, fast_bp_config,
        )
        result = analyzer.estimate_from_nij_matrix()

        assert result.point == pytest.approx(0.444, abs=TOL)
        _assert_brackets_point(result)


# ── Gwet cont3x3 (2 raters, K=3, balanced) ──────────────────────────────── #
#
# Source: Gwet, K. L. (2014). Handbook of Inter-Rater Reliability (4th ed.).
#         Table 4.1 (cont3x3abstractors).
# Cross-check: irrCAC: Cohen κ = 0.79641, AC1 = 0.84933.


class TestGwetCont3x3:
    """Gwet cont3x3 dataset (100 subjects, 2 raters, K=3)."""

    def test_cohen_kappa_ratings(self, gwet_ratings, fast_cohen_config):
        """Cohen κ = 0.796 via ratings pathway (irrCAC: 0.79641)."""
        analyzer = CohenKappaAnalyzer(gwet_ratings, fast_cohen_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.796, abs=TOL)
        _assert_brackets_point(result)

    def test_gwet_ac1_ratings(self, gwet_ratings, fast_gwet_config):
        """Gwet AC1 = 0.849 via ratings pathway (irrCAC: 0.84933)."""
        analyzer = GwetAC1Analyzer(gwet_ratings, fast_gwet_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.849, abs=TOL)
        _assert_brackets_point(result)

    def test_krippendorff_alpha_nij(self, gwet_nij, fast_krippendorff_config):
        """Krippendorff α = 0.797 via nij pathway (V2 validated)."""
        analyzer = KrippendorffAlphaAnalyzer.build_from_nij_matrix(
            gwet_nij, fast_krippendorff_config,
        )
        result = analyzer.estimate_from_nij_matrix()

        assert result.point == pytest.approx(0.797, abs=TOL)
        _assert_brackets_point(result)


# ── Synthetic 2-rater (ratings pathway, 2 raters, K=3, balanced) ─────────── #
#
# Synthetic dataset.  Cohen κ = Fleiss κ = Scott's π to 3 decimal places
# on this balanced 2-rater data (individual marginals are close enough
# that Pe_Cohen ≈ Pe_Fleiss).


class TestSynthetic2Rater:
    """Synthetic 2-rater dataset (100 subjects, 2 raters, K=3)."""

    def test_cohen_kappa(self, synthetic_ratings, fast_cohen_config):
        """Cohen κ = 0.641 (V2 validated)."""
        analyzer = CohenKappaAnalyzer(synthetic_ratings, fast_cohen_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.641, abs=TOL)
        _assert_brackets_point(result)

    def test_scotts_pi(self, synthetic_ratings, fast_scotts_config):
        """Scott's π = 0.641 (V2 validated, = Fleiss at n=2)."""
        analyzer = ScottsAnalyzer(synthetic_ratings, fast_scotts_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.641, abs=TOL)
        _assert_brackets_point(result)


# ── BCa, end to end ──────────────────────────────────────────────────────── #
#
# Every test above uses percentile, for the reason in the module docstring.
# This is the one place the BCa path runs through the whole pipeline on real
# data, so that a break in it cannot hide behind the unit tests.


class TestBcaEndToEnd:
    """BCa at B=1000 on the Fleiss 1971 data (30 subjects, 6 raters)."""

    def test_bca_runs_without_falling_back(self, fleiss1971_nij, bca_fleiss_config):
        """The interval must be BCa's own work, not percentile wearing its name.

        ``method_used`` is the assertion that matters. BCa steps down to
        percentile when it cannot compute a correction, and the point
        estimate and bounds would look perfectly reasonable either way --
        only the interval's own account distinguishes them.
        """
        analyzer = FleissKappaAnalyzer.build_from_nij_matrix(
            fleiss1971_nij, bca_fleiss_config,
        )
        result = analyzer.estimate_from_nij_matrix()
        ci = result.confidence_interval

        assert result.point == pytest.approx(0.430, abs=TOL)
        assert ci.method_intended == "bca"
        assert ci.method_used == "bca"
        assert ci.degraded is False

    def test_bca_carries_its_diagnostics(self, fleiss1971_nij, bca_fleiss_config):
        """A real BCa interval reports the two constants it was built from."""
        analyzer = FleissKappaAnalyzer.build_from_nij_matrix(
            fleiss1971_nij, bca_fleiss_config,
        )
        ci = analyzer.estimate_from_nij_matrix().confidence_interval

        assert ci.acceleration is not None
        assert ci.bias_correction is not None
        assert abs(ci.acceleration) < 1.0

    def test_bca_is_given_enough_resamples(self, fleiss1971_nij, bca_fleiss_config):
        """B=1000 is Efron & Tibshirani's minimum for BCa; the interval agrees."""
        analyzer = FleissKappaAnalyzer.build_from_nij_matrix(
            fleiss1971_nij, bca_fleiss_config,
        )
        ci = analyzer.estimate_from_nij_matrix().confidence_interval

        assert ci.n_resamples == 1000
        assert ci.recommended_resamples == 1000
        assert ci.below_recommended_resamples is False

    # Containment is deliberately not asserted here -- see the module
    # docstring. It happens to hold on this dataset, and it is not a
    # property BCa guarantees.
