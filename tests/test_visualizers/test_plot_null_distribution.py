"""
plot_null_distribution: what the plot says, not how it looks.

Each test reads the figure back and checks that it states the facts of the
test: which coefficient, the p-value with the counts it came from, and that
the coloured shuffles are exactly the ones counted toward p.
"""

import matplotlib

matplotlib.use("Agg")  # draw without a screen, as on CI

import matplotlib.pyplot as plt
import numpy as np
import pytest

from interrater import plot_null_distribution
from interrater.analyzers import FleissKappaAnalyzer, KrippendorffAlphaAnalyzer
from interrater.analyzers.configs import (
    FleissKappaConfig,
    HypothesisTestConfig,
    KrippendorffAlphaConfig,
)
from interrater.data_objects.ratings import Ratings
from interrater.intermediates.results import HypothesisTestResult
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _ratings(agreement: float) -> Ratings:
    """30 subjects, 3 raters, 3 categories; each rater copies a true answer with this probability."""
    rng = np.random.default_rng(42)
    truth = rng.integers(0, 3, 30)
    noise = rng.integers(0, 3, (30, 3))
    codes = np.where(rng.random((30, 3)) < agreement, truth[:, None], noise)
    cats = UnorderedCategories(names=("A", "B", "C"), labels={"A": "A", "B": "B", "C": "C"})
    return Ratings(category_index=codes.astype(np.int8),
                   semantics=RatingsSemantics(raters=3, categories=cats))


def _test(ratings, alternative="two-sided", analyzer=FleissKappaAnalyzer,
          config=FleissKappaConfig):
    ht = HypothesisTestConfig(n_permutations=500, alternative=alternative)
    return analyzer(ratings, config(hypothesis_test_config=ht)).hypothesis_test()


def _legend(fig):
    return [t.get_text() for t in fig.axes[0].get_legend().get_texts()]


def _box(fig):
    return " ".join(t.get_text() for t in fig.axes[0].texts)


def _red_count(fig):
    """Shuffles in the coloured bars: the second stacked series."""
    bars = fig.axes[0].containers[1]
    return int(round(sum(p.get_height() for p in bars)))


# -- The tail and the arithmetic -------------------------------------------- #


@pytest.mark.parametrize("alternative", ["two-sided", "greater", "less"])
def test_coloured_shuffles_are_the_counted_ones(alternative):
    result = _test(_ratings(0.0), alternative)
    fig = plot_null_distribution(result)
    side = "above" if any("at or above" in s for s in _legend(fig)) else "below"
    count = result.n_at_or_above if side == "above" else result.n_at_or_below
    assert _red_count(fig) == count
    assert f"({count:,} + 1) / ({result.n_usable:,} + 1)" in _box(fig)
    assert f"{alternative} p = {result.p_value:.4f}" in _box(fig)


def test_two_sided_shows_the_doubling():
    fig = plot_null_distribution(_test(_ratings(0.0), "two-sided"))
    assert "= 2 x (" in _box(fig)


def test_one_sided_has_no_doubling():
    fig = plot_null_distribution(_test(_ratings(0.0), "greater"))
    assert "2 x" not in _box(fig)
    assert any("at or above" in s for s in _legend(fig))


def test_strong_agreement_nothing_as_extreme():
    """Observed beyond every shuffle: the smallest p 500 shuffles can give."""
    result = _test(_ratings(0.9))
    fig = plot_null_distribution(result)
    assert "2 x (0 + 1) / (500 + 1)" in _box(fig)
    assert _red_count(fig) == 0
    assert "reject at alpha" in _box(fig)


def test_cap_at_one_is_stated():
    """Both tails hold 7 of 9 shuffles: 2 x 8/10 is capped at 1."""
    null = np.array([-0.2, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.2])
    result = HypothesisTestResult(
        p_value=1.0, alpha=0.05, reject_null=False, method="permutation",
        point=0.0, coefficient="fleiss_kappa", null_distribution=null,
        n_permutations=9, alternative="two-sided",
        n_usable=9, n_at_or_above=7, n_at_or_below=7,
    )
    assert "capped at 1" in _box(plot_null_distribution(result))


# -- Labels ----------------------------------------------------------------- #


def test_title_names_the_coefficient():
    result = _test(_ratings(0.0), analyzer=KrippendorffAlphaAnalyzer,
                   config=KrippendorffAlphaConfig)
    assert plot_null_distribution(result).axes[0].get_title() == (
        "Krippendorff's alpha: permutation test"
    )


def test_caption_counts_usable_shuffles():
    """Undefined shuffles are left out and counted, as the p-value left them out."""
    null = np.array([0.1, np.nan, 0.2, 0.3, np.nan, 0.0, 0.1, 0.2, 0.1, 0.0])
    result = HypothesisTestResult(
        p_value=1.0, alpha=0.05, reject_null=False, method="permutation",
        point=0.1, coefficient="fleiss_kappa", null_distribution=null,
        n_permutations=10, alternative="two-sided",
        n_usable=8, n_at_or_above=6, n_at_or_below=5,
    )
    assert "8 of 10 shuffles usable" in plot_null_distribution(result).axes[0].get_xlabel()


def test_draws_into_given_axes():
    fig, (left, right) = plt.subplots(1, 2)
    assert plot_null_distribution(_test(_ratings(0.0)), ax=right) is fig
    assert right.get_title() and not left.get_title()


# -- Refusals --------------------------------------------------------------- #


@pytest.mark.filterwarnings("ignore::UserWarning")  # parametric_null warns; not the subject
def test_parametric_null_refused():
    ht = HypothesisTestConfig(method="parametric_null")
    result = FleissKappaAnalyzer(
        _ratings(0.5), FleissKappaConfig(hypothesis_test_config=ht),
    ).hypothesis_test()
    with pytest.raises(ValueError, match="only a permutation test"):
        plot_null_distribution(result)


def test_result_without_counts_refused():
    result = HypothesisTestResult(
        p_value=0.5, alpha=0.05, reject_null=False, method="permutation",
        point=0.1, coefficient="fleiss_kappa",
        null_distribution=np.array([0.0, 0.1, 0.2]), n_permutations=3,
    )
    with pytest.raises(ValueError, match="does not record the counts"):
        plot_null_distribution(result)


def test_not_a_result_refused():
    with pytest.raises(TypeError, match="HypothesisTestResult"):
        plot_null_distribution(np.zeros(10))
