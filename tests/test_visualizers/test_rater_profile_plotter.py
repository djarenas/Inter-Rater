"""
RaterProfilePlotter: what the plot says, not how it looks.

Each test reads the figure back and checks that it states the facts of the
rater profile: the shared counts, a white diagonal outside the colour
scale, hatching on exactly the pairs the analysis gives no coefficient,
and each rater's total against all subjects.
"""

import matplotlib

matplotlib.use("Agg")  # draw without a screen, as on CI

import matplotlib.pyplot as plt
import numpy as np
import pytest

from interrater import RaterProfilePlotter
from interrater.analyzers import KrippendorffAlphaAnalyzer
from interrater.analyzers.configs import KrippendorffAlphaConfig
from interrater.analyzers.configs.cohen_kappa_config import CohenKappaConfig
from interrater.analyzers.configs.pairwise_permutation_config import (
    PairwisePermutationConfig,
)
from interrater.analyzers.permutation_analysis import PairwisePermutationAnalysis
from interrater.contracts import EMPTY
from interrater.data_objects.ratings import Ratings
from interrater.intermediates.data_profiles import RaterProfile
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.visualizers.rater_profile_plotter import THIN_PAIR_LABEL

NAMES = ("ann", "bob", "cat", "dan")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def gappy():
    """
    6 subjects, 4 raters. ann rated 4, bob 3, cat 4, dan 3. Every pair
    shares 2 subjects except bob and dan, who share 1: the one thin pair.
    """
    index = np.array([
        [0, 0, 1, 0],
        [1, 1, EMPTY, EMPTY],
        [0, EMPTY, 0, EMPTY],
        [EMPTY, 1, 1, EMPTY],
        [0, EMPTY, EMPTY, 0],
        [EMPTY, EMPTY, 1, 1],
    ])
    semantics = RatingsSemantics(raters=list(NAMES), categories=UnorderedCategories(("x", "y")))
    return Ratings(index, semantics)


@pytest.fixture
def complete():
    """20 subjects, 4 raters, every rater rated every subject."""
    rng = np.random.default_rng(42)
    codes = rng.integers(0, 3, (20, 4)).astype(np.int8)
    semantics = RatingsSemantics(raters=list(NAMES), categories=UnorderedCategories(("A", "B", "C")))
    return Ratings(codes, semantics)


def _image(fig):
    return fig.axes[0].images[0]


def _hatched_cells(fig):
    """(row, col) of every hatched cell."""
    return {
        (int(round(p.get_y() + 0.5)), int(round(p.get_x() + 0.5)))
        for p in fig.axes[0].patches if p.get_hatch() == "///"
    }


# -- Construction ------------------------------------------------------------ #


class TestConstruction:

    def test_from_ratings(self, gappy):
        plotter = RaterProfilePlotter.build_from_ratings(gappy)
        assert plotter.rater_profile.rater_names == NAMES
        assert plotter.annotate is True

    def test_from_point_estimate_result(self, gappy):
        result = KrippendorffAlphaAnalyzer(
            gappy, KrippendorffAlphaConfig(n_resamples=20, random_state=42),
        ).estimate()
        plotter = RaterProfilePlotter.build_from_result(result)
        assert plotter.rater_profile is result.rater_profile

    def test_from_permutation_result(self, complete):
        config = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(n_resamples=20, random_state=42),
        )
        result = PairwisePermutationAnalysis(complete, config).run()
        plotter = RaterProfilePlotter.build_from_result(result)
        assert plotter.rater_profile is result.rater_profile

    def test_result_from_counts_refused_when_built(self, gappy):
        result = KrippendorffAlphaAnalyzer.build_from_nij_matrix(
            gappy.to_nij_matrix(), KrippendorffAlphaConfig(n_resamples=20, random_state=42),
        ).estimate_from_nij_matrix()
        with pytest.raises(ValueError, match="built from counts: no rater identity"):
            RaterProfilePlotter.build_from_result(result)

    def test_not_a_result_refused(self, gappy):
        with pytest.raises(TypeError, match="PointEstimateResult or a PermutationResult"):
            RaterProfilePlotter.build_from_result(gappy)

    def test_not_a_rater_profile_refused(self, gappy):
        with pytest.raises(TypeError, match="must be a RaterProfile"):
            RaterProfilePlotter(gappy)

    def test_annotate_must_be_bool(self, gappy):
        profile = RaterProfile.build_from_ratings(gappy)
        with pytest.raises(TypeError, match="annotate must be bool"):
            RaterProfilePlotter(profile, annotate=1)


# -- draw_overlap ------------------------------------------------------------ #


class TestDrawOverlap:

    def test_returns_figure_and_uses_given_axes(self, gappy):
        fig, ax = plt.subplots()
        assert RaterProfilePlotter.build_from_ratings(gappy).draw_overlap(ax=ax) is fig

    def test_raters_named_on_both_axes(self, gappy):
        ax = RaterProfilePlotter.build_from_ratings(gappy).draw_overlap().axes[0]
        assert tuple(t.get_text() for t in ax.get_xticklabels()) == NAMES
        assert tuple(t.get_text() for t in ax.get_yticklabels()) == NAMES

    def test_full_matrix_with_diagonal_masked(self, gappy):
        """Every off-diagonal cell is coloured, both triangles; the diagonal is not."""
        fig = RaterProfilePlotter.build_from_ratings(gappy).draw_overlap()
        cells = _image(fig).get_array()
        assert cells.shape == (4, 4)
        assert np.array_equal(np.ma.getmaskarray(cells), np.eye(4, dtype=bool))
        shared = RaterProfile.build_from_ratings(gappy).rater_overlap.shared
        off = ~np.eye(4, dtype=bool)
        assert np.array_equal(cells.data[off], shared[off])

    def test_diagonal_white_and_outside_the_scale(self, gappy):
        """The scale runs to the largest pair (2), not the largest total (4)."""
        image = _image(RaterProfilePlotter.build_from_ratings(gappy).draw_overlap())
        assert image.cmap.get_bad() == pytest.approx((1.0, 1.0, 1.0, 1.0))
        assert image.norm.vmin == 0
        assert image.norm.vmax == 2

    def test_thin_pair_hatched_with_legend(self, gappy):
        fig = RaterProfilePlotter.build_from_ratings(gappy).draw_overlap()
        assert _hatched_cells(fig) == {(1, 3), (3, 1)}  # bob and dan, both triangles
        legend = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
        assert legend == [THIN_PAIR_LABEL]

    def test_no_thin_pair_no_hatching_no_legend(self, complete):
        fig = RaterProfilePlotter.build_from_ratings(complete).draw_overlap()
        assert _hatched_cells(fig) == set()
        assert fig.axes[0].get_legend() is None

    def test_annotate_writes_every_cell(self, gappy):
        fig = RaterProfilePlotter.build_from_ratings(gappy).draw_overlap()
        shared = RaterProfile.build_from_ratings(gappy).rater_overlap.shared
        written = {(int(t.get_position()[1]), int(t.get_position()[0])): t.get_text()
                   for t in fig.axes[0].texts}
        assert written == {(i, j): str(shared[i, j]) for i in range(4) for j in range(4)}

    def test_annotate_off_writes_nothing(self, gappy):
        plotter = RaterProfilePlotter.build_from_ratings(gappy, annotate=False)
        assert len(plotter.draw_overlap().axes[0].texts) == 0


# -- draw_counts ------------------------------------------------------------- #


class TestDrawCounts:

    def test_returns_figure_and_uses_given_axes(self, gappy):
        fig, ax = plt.subplots()
        assert RaterProfilePlotter.build_from_ratings(gappy).draw_counts(ax=ax) is fig

    def test_one_bar_per_rater_in_declared_order(self, gappy):
        ax = RaterProfilePlotter.build_from_ratings(gappy).draw_counts().axes[0]
        assert [p.get_height() for p in ax.patches] == [4, 3, 4, 3]
        assert tuple(t.get_text() for t in ax.get_xticklabels()) == NAMES

    def test_dashed_line_at_all_subjects(self, gappy):
        ax = RaterProfilePlotter.build_from_ratings(gappy).draw_counts().axes[0]
        (line,) = ax.get_lines()
        assert line.get_linestyle() == "--"
        assert set(line.get_ydata()) == {6}
        assert [t.get_text() for t in ax.get_legend().get_texts()] == ["all 6 subjects"]
