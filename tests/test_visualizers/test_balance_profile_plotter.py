"""
BalanceProfilePlotter: what the plot says, not how it looks.

Each test reads the figure back and checks that it states the facts of the
balance profile: one bar for every count of raters from the fewest to the
most, zero bars for counts no subject had, each bar's count written on it,
and whether the data is balanced.
"""

import matplotlib

matplotlib.use("Agg")  # draw without a screen, as on CI

import matplotlib.pyplot as plt
import numpy as np
import pytest

from interrater import BalanceProfilePlotter
from interrater.analyzers import KrippendorffAlphaAnalyzer
from interrater.analyzers.configs import KrippendorffAlphaConfig
from interrater.analyzers.configs.cohen_kappa_config import CohenKappaConfig
from interrater.analyzers.configs.pairwise_permutation_config import (
    PairwisePermutationConfig,
)
from interrater.analyzers.permutation_analysis import PairwisePermutationAnalysis
from interrater.contracts import EMPTY
from interrater.data_objects.ratings import Ratings
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def gappy():
    """6 subjects, 4 raters: five subjects have 2 raters, one has all 4, none has 3."""
    index = np.array([
        [0, 0, 1, 0],
        [1, 1, EMPTY, EMPTY],
        [0, EMPTY, 0, EMPTY],
        [EMPTY, 1, 1, EMPTY],
        [0, EMPTY, EMPTY, 0],
        [EMPTY, EMPTY, 1, 1],
    ])
    semantics = RatingsSemantics(raters=4, categories=UnorderedCategories(("x", "y")))
    return Ratings(index, semantics)


@pytest.fixture
def complete():
    """20 subjects, 4 raters, every rater rated every subject."""
    rng = np.random.default_rng(42)
    codes = rng.integers(0, 3, (20, 4)).astype(np.int8)
    semantics = RatingsSemantics(raters=4, categories=UnorderedCategories(("A", "B", "C")))
    return Ratings(codes, semantics)


def _config():
    return KrippendorffAlphaConfig(n_resamples=20, random_state=42)


# -- Construction ------------------------------------------------------------ #


class TestConstruction:

    def test_from_ratings_and_from_counts_agree(self, gappy):
        from_ratings = BalanceProfilePlotter.build_from_ratings(gappy)
        from_counts = BalanceProfilePlotter.build_from_nij_matrix(gappy.to_nij_matrix())
        assert from_ratings.balance_profile.subject_counts == {2: 5, 4: 1}
        assert from_counts.balance_profile.subject_counts == {2: 5, 4: 1}

    def test_from_point_estimate_result(self, gappy):
        result = KrippendorffAlphaAnalyzer(gappy, _config()).estimate()
        plotter = BalanceProfilePlotter.build_from_result(result)
        assert plotter.balance_profile is result.balance_profile

    def test_from_result_built_from_counts(self, gappy):
        """Unlike the rater profile, a counts matrix has a balance profile."""
        result = KrippendorffAlphaAnalyzer.build_from_nij_matrix(
            gappy.to_nij_matrix(), _config(),
        ).estimate_from_nij_matrix()
        plotter = BalanceProfilePlotter.build_from_result(result)
        assert plotter.balance_profile is result.balance_profile

    def test_from_permutation_result(self, complete):
        config = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(n_resamples=20, random_state=42),
        )
        result = PairwisePermutationAnalysis(complete, config).run()
        plotter = BalanceProfilePlotter.build_from_result(result)
        assert plotter.balance_profile is result.balance_profile

    def test_not_a_result_refused(self, gappy):
        with pytest.raises(TypeError, match="PointEstimateResult or a PermutationResult"):
            BalanceProfilePlotter.build_from_result(gappy)

    def test_not_a_balance_profile_refused(self, gappy):
        with pytest.raises(TypeError, match="must be a BalanceProfile"):
            BalanceProfilePlotter(gappy)


# -- draw -------------------------------------------------------------------- #


class TestDraw:

    def test_returns_figure_and_uses_given_axes(self, gappy):
        fig, ax = plt.subplots()
        assert BalanceProfilePlotter.build_from_ratings(gappy).draw(ax=ax) is fig

    def test_every_count_from_fewest_to_most_with_zero_bars(self, gappy):
        ax = BalanceProfilePlotter.build_from_ratings(gappy).draw().axes[0]
        bars = [(p.get_x() + p.get_width() / 2, p.get_height()) for p in ax.patches]
        assert bars == [(2, 5), (3, 0), (4, 1)]

    def test_count_written_on_every_bar(self, gappy):
        ax = BalanceProfilePlotter.build_from_ratings(gappy).draw().axes[0]
        assert [t.get_text() for t in ax.texts] == ["5", "0", "1"]

    def test_unbalanced_title(self, gappy):
        title = BalanceProfilePlotter.build_from_ratings(gappy).draw().axes[0].get_title()
        assert title == "Unbalanced: subjects have 2 to 4 raters\n6 subjects"

    def test_balanced_title(self, complete):
        ax = BalanceProfilePlotter.build_from_ratings(complete).draw().axes[0]
        assert ax.get_title() == "Balanced: every subject has 4 raters\n20 subjects"
        assert [p.get_height() for p in ax.patches] == [20]
