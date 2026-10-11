"""
Integration tests against the original systematic review datasets.

Each SR gets one or two coefficients chosen to exercise a different
analyzer and avoid redundancy with test_reference_data.py.

Point estimates are asserted to 3 decimal places (abs=5e-4) and CI
structure (lower < point < upper) is checked as a sanity guard with
n_resamples=50.

Datasets:
    SR1 prelim  -- ratings, 4 raters, K=3, balanced, ~83% skew
    SR1 main    -- ratings, 4 raters, K=3, unbalanced, ~90% skew
    SR2         -- ratings, 4 raters, K=3, nearly balanced (80/81)
    SR3         -- ratings, 3 raters, K=2, balanced, ~90% skew
    SR4         -- ratings, 4 raters, K=2, balanced, ~82% skew
    SR5         -- ratings, 3 raters, K=2, balanced, ~88% skew
"""

import pytest

from interrater import (
    FleissKappaAnalyzer,
    GwetAC1Analyzer,
    BrennanPredigerAnalyzer,
    KrippendorffAlphaAnalyzer,
    CongerKappaAnalyzer,
    LightKappaAnalyzer,
)
from interrater.analyzers.configs import FleissKappaConfig

# Tolerance: 3 decimal places.
TOL = 5e-4


# ── SR1 prelim (ratings, 4 raters, K=3, balanced) ─────────────────────── #
#
# Source: Systematic Review 1, preliminary screening (199 subjects).
# Reference: example_data/README.md (V2 validated values).
# Also exercises Conger and Light (individual-marginal analyzers).


class TestSR1Prelim:
    """SR1 preliminary data (199 subjects, 4 raters, K=3, balanced)."""

    def test_fleiss_kappa(self, sr1_prelim_ratings, fast_fleiss_config):
        """Fleiss' kappa = 0.493 (V2 validated, group coefficient)."""
        analyzer = FleissKappaAnalyzer(sr1_prelim_ratings, fast_fleiss_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.493, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

    def test_conger_kappa(self, sr1_prelim_ratings, fast_conger_config):
        """Conger's kappa = 0.494 (group-level, individual marginals).

        Conger's group-level estimate (0.494) differs from the pairwise
        mean Cohen kappa (0.496) because it uses each rater's full
        marginal distribution, not pair-by-pair marginals.
        """
        analyzer = CongerKappaAnalyzer(sr1_prelim_ratings, fast_conger_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.494, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

    def test_light_kappa(self, sr1_prelim_ratings, fast_light_config):
        """Light's kappa = 0.496 (mean of C(4,2)=6 pairwise Cohens).

        Light's kappa is by definition the mean of all pairwise Cohen
        kappas, so 0.496 is exact (matches pairwise mean in README).
        """
        analyzer = LightKappaAnalyzer(sr1_prelim_ratings, fast_light_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.496, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]
        assert result.n_pairs == 6  # C(4, 2)


# ── SR1 main (ratings, 4 raters, K=3, unbalanced) ─────────────────────── #
#
# Source: Systematic Review 1, main screening (1617 subjects).
# Most subjects rated by exactly 2 of 4 raters (3 have 3 raters).
# Tests the unbalanced-data approximation (Fleiss) vs exact (Krippendorff).


class TestSR1Main:
    """SR1 main data (1617 subjects, 4 raters, K=3, unbalanced)."""

    def test_fleiss_kappa(self, sr1_main_ratings):
        """Fleiss' kappa = 0.619 (unbalanced approximation).

        FleissKappaConfig defaults to balanced_data_only=True, which
        rejects unbalanced data.  This test explicitly sets it to False
        to exercise the unbalanced-data approximation path.
        """
        config = FleissKappaConfig(
            n_resamples=50, random_state=42, ci_method="bca",
            balanced_data_only=False,
        )
        analyzer = FleissKappaAnalyzer(sr1_main_ratings, config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.619, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]

    def test_krippendorff_alpha(self, sr1_main_ratings, fast_krippendorff_config):
        """Krippendorff alpha = 0.618 (exact on unbalanced data)."""
        analyzer = KrippendorffAlphaAnalyzer(sr1_main_ratings, fast_krippendorff_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.618, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]


# ── SR2 (ratings, 4 raters, K=3, nearly balanced) ─────────────────────── #
#
# Source: Systematic Review 2 (81 subjects, 80 fully rated).


class TestSR2:
    """SR2 data (81 subjects, 4 raters, K=3, nearly balanced)."""

    def test_gwet_ac1(self, sr2_ratings, fast_gwet_config):
        """Gwet AC1 = 0.782 (V2 validated)."""
        analyzer = GwetAC1Analyzer(sr2_ratings, fast_gwet_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.782, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]


# ── SR3 FIS Psych (ratings, 3 raters, K=2, balanced) ──────────────────── #
#
# Source: SR1 full-text screening, FIS + Psych databases (770 subjects).
# Binary include/exclude only (no uncertain category).


class TestSR3:
    """SR3 FIS Psych data (770 subjects, 3 raters, K=2, balanced)."""

    def test_brennan_prediger(self, sr3_ratings, fast_bp_config):
        """Brennan-Prediger = 0.886 (V2 validated, binary 3 raters)."""
        analyzer = BrennanPredigerAnalyzer(sr3_ratings, fast_bp_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.886, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]


# ── SR4 FIS CVD (ratings, 4 raters, K=2, balanced) ────────────────────── #
#
# Source: SR1 full-text screening, FIS CVD database (787 subjects).
# Binary include/exclude only (no uncertain category).


class TestSR4:
    """SR4 FIS CVD data (787 subjects, 4 raters, K=2, balanced)."""

    def test_gwet_ac1(self, sr4_ratings, fast_gwet_config):
        """Gwet AC1 = 0.879 (V2 validated, binary 4 raters)."""
        analyzer = GwetAC1Analyzer(sr4_ratings, fast_gwet_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.879, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]


# ── SR5 Cocaine Cardio (ratings, 3 raters, K=2, balanced) ─────────────── #
#
# Source: SR2 full-text screening, Cocaine-Cardio database (881 subjects).
# Binary include/exclude only (no uncertain category).


class TestSR5:
    """SR5 Cocaine Cardio data (881 subjects, 3 raters, K=2, balanced)."""

    def test_brennan_prediger(self, sr5_ratings, fast_bp_config):
        """Brennan-Prediger = 0.788 (V2 validated, binary 3 raters)."""
        analyzer = BrennanPredigerAnalyzer(sr5_ratings, fast_bp_config)
        result = analyzer.estimate()

        assert result.point == pytest.approx(0.788, abs=TOL)
        assert result.confidence_interval[0] < result.point < result.confidence_interval[1]
