"""
Integration tests against datasets with externally published reference values.

Each test asserts the point estimate to 3 decimal places (abs=5e-4) and
checks CI structure (lower < point < upper) as a sanity guard with
n_resamples=50.

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
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

    def test_gwet_ac1(self, fleiss1971_nij, fast_gwet_config):
        """Gwet AC1 = 0.448 (V2 validated)."""
        analyzer = GwetAC1Analyzer.build_from_nij_matrix(
            fleiss1971_nij, fast_gwet_config,
        )
        result = analyzer.estimate_from_nij_matrix()

        assert result.point == pytest.approx(0.448, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

    def test_brennan_prediger(self, fleiss1971_nij, fast_bp_config):
        """Brennan-Prediger = 0.444 (V2 validated, Pe = 1/5)."""
        analyzer = BrennanPredigerAnalyzer.build_from_nij_matrix(
            fleiss1971_nij, fast_bp_config,
        )
        result = analyzer.estimate_from_nij_matrix()

        assert result.point == pytest.approx(0.444, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]


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
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

    def test_gwet_ac1_ratings(self, gwet_ratings, fast_gwet_config):
        """Gwet AC1 = 0.849 via ratings pathway (irrCAC: 0.84933)."""
        analyzer = GwetAC1Analyzer(gwet_ratings, fast_gwet_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.849, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

    def test_krippendorff_alpha_nij(self, gwet_nij, fast_krippendorff_config):
        """Krippendorff α = 0.797 via nij pathway (V2 validated)."""
        analyzer = KrippendorffAlphaAnalyzer.build_from_nij_matrix(
            gwet_nij, fast_krippendorff_config,
        )
        result = analyzer.estimate_from_nij_matrix()

        assert result.point == pytest.approx(0.797, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]


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
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

    def test_scotts_pi(self, synthetic_ratings, fast_scotts_config):
        """Scott's π = 0.641 (V2 validated, = Fleiss at n=2)."""
        analyzer = ScottsAnalyzer(synthetic_ratings, fast_scotts_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.641, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]
