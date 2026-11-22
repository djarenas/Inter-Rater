"""
Integration tests for weighted analysis (linear + quadratic weights).

Tests that distance-based weighting changes the point estimate, for the group
analyzers with a weighted form (Fleiss, Gwet's AC2, B-P, Krippendorff) and for
Cohen.

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

Weights are declared in the semantics, so each weighted test runs the
ordinary config on ratings whose categories are declared ordered
(``ordered_ratings``), and compares against the same data declared unordered.
"""

from pathlib import Path

import numpy as np
import pytest
import yaml

from interrater import (
    FleissKappaAnalyzer,
    GwetAC1Analyzer,
    KrippendorffAlphaAnalyzer,
    CohenKappaAnalyzer,
    BrennanPredigerAnalyzer,
    OrderedCategories,
    RatingsSemantics,
)
from interrater.analyzers.permutation_analysis import (
    PairwisePermutationAnalysis,
    SubsetPermutationAnalysis,
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
# Two weight types on the same data, with two different group
# analyzers.  Categories are ordinal (-1, 0, 1), so
# -1/0 is closer than -1/1 under distance-based weights.
#
# Unweighted reference values (from test_research_data.py):
#   Fleiss kappa  = 0.493
#   Krippendorff  = 0.493  (approx; 0.493269 from ratings pathway)


class TestSR1PrelimWeighted:
    """Weighted analysis on SR1 prelim (199 subjects, 4 raters, K=3)."""

    def test_fleiss_kappa_quadratic(self, sr1_prelim_ratings, ordered_ratings, fast_fleiss_config):
        """Fleiss' kappa with quadratic weights = 0.593.

        Quadratic weights (w_ij = 1 - (i-j)^2 / (K-1)^2) upweight partial
        agreements, increasing the estimate from 0.493 to 0.593.
        """
        ratings = ordered_ratings("sr1_prelim", "quadratic")
        analyzer = FleissKappaAnalyzer(ratings, fast_fleiss_config)
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

    def test_krippendorff_alpha_linear(self, sr1_prelim_ratings, ordered_ratings, fast_krippendorff_config):
        """Krippendorff's alpha with linear weights = 0.555.

        Linear weights (w_ij = 1 - |i-j| / (K-1)) are less aggressive
        than quadratic, but still upweight partial agreements.
        Different weight type on same data as Fleiss quadratic above.
        """
        ratings = ordered_ratings("sr1_prelim", "linear")
        analyzer = KrippendorffAlphaAnalyzer(ratings, fast_krippendorff_config)
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


# ── SR1 prelim, weights declared in the semantics ──────────────────── #
#
# The shipped weighted example gives its matrix two ways: in a separate
# weights file (SR1_prelim_weighted_semantics.yaml) and inline
# (SR1_prelim_weighted_inline_semantics.yaml). Same matrix, so every
# coefficient and every pair or subset must come out identical.

SR1_DIR = Path(__file__).resolve().parents[2] / "example_data" / "SR1_OSF2018"
SR1_MATRIX = [[1.0, 0.8, 0.0], [0.8, 1.0, 0.2], [0.0, 0.2, 1.0]]


@pytest.fixture(scope="module")
def from_file(sr1_prelim_with_semantics):
    """SR1 prelim with its matrix in a separate weights file."""
    return sr1_prelim_with_semantics(SR1_DIR / "SR1_prelim_weighted_semantics.yaml")


@pytest.fixture(scope="module")
def inline(sr1_prelim_with_semantics):
    """SR1 prelim with the same matrix written inline."""
    return sr1_prelim_with_semantics(SR1_DIR / "SR1_prelim_weighted_inline_semantics.yaml")


class TestSR1PrelimWeightsFromSemantics:
    """The weights come from the semantics, however they are written down."""

    @pytest.mark.parametrize("analyzer_cls, config_name", [
        (FleissKappaAnalyzer, "fast_fleiss_config"),
        (GwetAC1Analyzer, "fast_gwet_config"),
        (BrennanPredigerAnalyzer, "fast_bp_config"),
        (KrippendorffAlphaAnalyzer, "fast_krippendorff_config"),
    ], ids=["fleiss", "gwet_ac2", "brennan_prediger", "krippendorff"])
    def test_file_and_inline_give_identical_results(
        self, from_file, inline, analyzer_cls, config_name, request,
    ):
        config = request.getfixturevalue(config_name)
        a = analyzer_cls(from_file, config).estimate()
        b = analyzer_cls(inline, config).estimate()
        assert a.point == b.point
        assert a.confidence_interval.as_tuple() == b.confidence_interval.as_tuple()
        assert a.weights_type == b.weights_type == "custom"
        assert a.weight_matrix == b.weight_matrix == SR1_MATRIX

    @pytest.mark.parametrize("analyzer_cls, config_name, expected", [
        (FleissKappaAnalyzer, "fast_fleiss_config", 0.5929),
        (GwetAC1Analyzer, "fast_gwet_config", 0.8675),
        (BrennanPredigerAnalyzer, "fast_bp_config", 0.7180),
        (KrippendorffAlphaAnalyzer, "fast_krippendorff_config", 0.5934),
    ], ids=["fleiss", "gwet_ac2", "brennan_prediger", "krippendorff"])
    def test_quadratic_matrix_in_a_file_reproduces_ordered_quadratic(
        self, tmp_path, sr1_prelim_with_semantics, ordered_ratings,
        analyzer_cls, config_name, expected, request,
    ):
        """A file holding the quadratic matrix is the quadratic scheme, number for number."""
        quadratic = OrderedCategories(("1", "0", "-1"), weighting="quadratic").weight_matrix()
        (tmp_path / "quadratic_weights.yaml").write_text(yaml.safe_dump(
            {"names": ["1", "0", "-1"], "matrix": quadratic.tolist()}
        ))
        semantics = tmp_path / "semantics.yaml"
        semantics.write_text(
            "raters: 4\n"
            "categories:\n"
            "    kind: weighted\n"
            "    names: ['1', '0', '-1']\n"
            "    weights_file: quadratic_weights.yaml\n"
            "n_subjects: 199\n"
        )
        config = request.getfixturevalue(config_name)
        declared = analyzer_cls(sr1_prelim_with_semantics(semantics), config).estimate()
        ordered = analyzer_cls(ordered_ratings("sr1_prelim", "quadratic"), config).estimate()
        assert declared.point == pytest.approx(ordered.point, abs=1e-12)
        assert declared.point == pytest.approx(expected, abs=5e-5)
        assert declared.weight_matrix == ordered.weight_matrix

    def test_loads_from_any_working_directory(self, tmp_path, monkeypatch):
        """weights_file is read next to the semantics file, not from where Python runs."""
        monkeypatch.chdir(tmp_path)
        semantics = RatingsSemantics.build_from_yaml(SR1_DIR / "SR1_prelim_weighted_semantics.yaml")
        np.testing.assert_array_equal(semantics.categories.weight_matrix(), np.array(SR1_MATRIX))

    def test_pairs_use_the_semantics_weights(self, from_file, pairwise_config):
        result = PairwisePermutationAnalysis(from_file, pairwise_config()).run()
        for pair, unit in result.pair_results.items():
            assert unit.weights_type == "custom", pair
            assert unit.weight_matrix == SR1_MATRIX, pair

    def test_subsets_use_the_semantics_weights(self, from_file, subset_config):
        result = SubsetPermutationAnalysis(from_file, subset_config(3)).run()
        for subset, unit in result.per_k_results[3].subset_results.items():
            assert unit.weights_type == "custom", subset
            assert unit.weight_matrix == SR1_MATRIX, subset


# ── Gwet cont3x3 (ratings, 2 raters, K=3, balanced) ─────────────────── #
#
# Cohen is the individual-marginal analyzer with a weighted form. irrCAC
# gives 0.89216 for the quadratic case (analyses/external_validation/).
#
# Unweighted reference (from test_reference_data.py):
#   Cohen kappa = 0.796


class TestGwetCont3x3Weighted:
    """Weighted analysis on Gwet cont3x3 (100 subjects, 2 raters, K=3)."""

    def test_cohen_kappa_quadratic(self, gwet_ratings, ordered_ratings, fast_cohen_config):
        """Cohen's kappa with quadratic weights = 0.892.

        Quadratic weights increase the estimate from 0.796 to 0.892 because
        adjacent-category disagreements get partial credit.
        """
        ratings = ordered_ratings("gwet", "quadratic")
        analyzer = CohenKappaAnalyzer(ratings, fast_cohen_config)
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


# ── Gwet's AC2 against irrCAC ───────────────────────────────────────── #
#
# AC2 is AC1 with weights (Gwet 2021, Handbook, 5th ed., Chapter 4).
# The references are the values of Gwet's R package irrCAC 1.4
# (gwet.ac1.raw), which rounds to 5 decimals
# (analyses/external_validation/r_tools.csv).


class TestGwetAC2:
    """AC2 through the full pathway, against irrCAC 1.4."""

    @pytest.mark.parametrize("dataset, weighting, irrcac", [
        ("gwet", "linear", 0.90002),
        ("gwet", "quadratic", 0.94024),
        ("sr1_prelim", "linear", 0.85247),
        ("sr1_prelim", "quadratic", 0.86749),
    ])
    def test_matches_irrcac(
        self, ordered_ratings, fast_gwet_config, dataset, weighting, irrcac, request,
    ):
        result = GwetAC1Analyzer(ordered_ratings(dataset, weighting), fast_gwet_config).estimate()
        assert result.point == pytest.approx(irrcac, abs=5e-6)
        _assert_brackets_point(result)
        unweighted = request.getfixturevalue(f"{dataset}_ratings")
        assert result.point != pytest.approx(
            GwetAC1Analyzer(unweighted, fast_gwet_config).estimate().point, abs=1e-6
        )
        assert result.weights_type == weighting
        assert f"{weighting} weights (AC2)" in str(result)

    def test_unweighted_is_ac1(self, gwet_ratings, fast_gwet_config):
        """Unordered categories: no weights, AC1, and the data line says unweighted."""
        result = GwetAC1Analyzer(gwet_ratings, fast_gwet_config).estimate()
        assert result.point == pytest.approx(0.84933, abs=5e-6)
        assert result.weights_type == "identity"
        assert "AC2" not in str(result)


# ── Synthetic 2-rater (ratings, 2 raters, K=3, balanced, ordinal) ────── #
#
# A second group analyzer with weights, different dataset.
#
# Unweighted reference (from test_reference_data.py):
#   Brennan-Prediger = 0.670  (synthetic_2rater B-P via Scott's equivalence)


class TestSynthetic2RaterWeighted:
    """Weighted analysis on Synthetic 2-rater (80 subjects, 2 raters, K=3)."""

    def test_brennan_prediger_linear(self, synthetic_ratings, ordered_ratings, fast_bp_config):
        """Brennan-Prediger with linear weights = 0.651.

        Linear weights lower the estimate from 0.670 to 0.651 here because
        the B-P chance correction (Pe = 1/K) interacts differently with
        partial credit than the observed-agreement correction in kappa.
        """
        ratings = ordered_ratings("synthetic", "linear")
        analyzer = BrennanPredigerAnalyzer(ratings, fast_bp_config)
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
