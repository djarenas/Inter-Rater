"""
plot_bootstrap_distribution: what the plot says, not how it looks.

Each test reads the figure back -- legend labels, the BCa box, the caption
under the x-axis -- and checks that it states the facts of the result: the
method that actually produced the interval, BCa's two corrections when BCa
ran, how many draws were usable, and a refusal when there is nothing to plot.
"""

import dataclasses

import matplotlib

matplotlib.use("Agg")  # draw without a screen, as on CI

import matplotlib.pyplot as plt
import numpy as np
import pytest

from interrater import plot_bootstrap_distribution
from interrater.analyzers import BrennanPredigerAnalyzer, FleissKappaAnalyzer
from interrater.analyzers.configs import BrennanPredigerConfig, FleissKappaConfig
from interrater.data_objects.ratings import Ratings
from interrater.intermediates.results.ci_results import BootstrapCIResult
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def ratings():
    """40 subjects, 3 raters, 3 categories; each rater copies a true answer 70% of the time."""
    rng = np.random.default_rng(42)
    truth = rng.integers(0, 3, 40)
    noise = rng.integers(0, 3, (40, 3))
    copy = rng.random((40, 3)) < 0.7
    codes = np.where(copy, truth[:, None], noise).astype(np.int8)
    cats = UnorderedCategories(names=("A", "B", "C"), labels={"A": "A", "B": "B", "C": "C"})
    return Ratings(category_index=codes, semantics=RatingsSemantics(raters=3, categories=cats))


def _fleiss(ratings, **config):
    return FleissKappaAnalyzer(ratings, FleissKappaConfig(random_state=42, **config)).estimate()


def _legend(fig):
    return [t.get_text() for t in fig.axes[0].get_legend().get_texts()]


def _texts(fig):
    return " ".join(t.get_text() for t in fig.axes[0].texts)


def _caption(fig):
    return fig.axes[0].get_xlabel()


# -- What is drawn ---------------------------------------------------------- #


def test_bca_shows_both_intervals_and_corrections(ratings):
    result = _fleiss(ratings, n_resamples=1000, ci_method="bca")
    ci = result.confidence_interval
    fig = plot_bootstrap_distribution(result)

    legend = _legend(fig)
    assert any(f"BCa interval [{ci.lower:.3f}, {ci.upper:.3f}]" in s for s in legend)
    assert any("percentile interval, for comparison" in s for s in legend)
    assert any(f"point estimate {result.point:.3f}" in s for s in legend)
    box = _texts(fig)
    assert f"z0 = {ci.bias_correction:+.3f}" in box
    assert f"a  = {ci.acceleration:+.4f}" in box


def test_percentile_shows_no_bca_view(ratings):
    fig = plot_bootstrap_distribution(_fleiss(ratings, n_resamples=300, ci_method="percentile"))
    assert any("percentile interval [" in s for s in _legend(fig))
    assert not any("for comparison" in s for s in _legend(fig))
    assert _texts(fig) == ""


def test_bca_fallback_says_so(ratings):
    """Asked for BCa, got percentile: the box says the interval is percentile."""
    result = _fleiss(ratings, n_resamples=1000, ci_method="bca")
    ci = result.confidence_interval
    fell_back = BootstrapCIResult(
        lower=ci.lower, upper=ci.upper, alpha=ci.alpha,
        n_resamples=ci.n_resamples, n_failed=ci.n_failed,
        method_intended="bca", method_used="percentile",
    )
    fig = plot_bootstrap_distribution(dataclasses.replace(result, confidence_interval=fell_back))
    assert "BCa could not be applied" in _texts(fig)
    assert any("percentile interval [" in s for s in _legend(fig))


# -- The caption ------------------------------------------------------------ #


def test_caption_counts_usable_draws(ratings):
    result = _fleiss(ratings, n_resamples=300, ci_method="percentile")
    ci = result.confidence_interval
    assert f"{ci.n_usable:,} of {ci.n_resamples:,} draws usable" in _caption(
        plot_bootstrap_distribution(result)
    )


def test_caption_flags_too_few_draws(ratings):
    fig = plot_bootstrap_distribution(_fleiss(ratings, n_resamples=100, ci_method="percentile"))
    assert "below the 250 recommended for percentile" in _caption(fig)


def test_discrete_draws_get_one_bar_per_value(ratings):
    """Brennan-Prediger depends only on raw agreement: few distinct values, many ties."""
    result = BrennanPredigerAnalyzer(
        ratings, BrennanPredigerConfig(n_resamples=300, random_state=42),
    ).estimate()
    n_distinct = np.unique(np.round(result.bootstrap_distribution, 10)).size
    fig = plot_bootstrap_distribution(result)
    assert f"{n_distinct} distinct values" in _caption(fig)
    assert len(fig.axes[0].patches) >= n_distinct  # one bar per value, plus the band


def test_last_bit_differences_are_one_value(ratings):
    """The same value reached by different routes can differ in the last bit."""
    result = _fleiss(ratings, n_resamples=300, ci_method="percentile")
    draws = np.repeat([0.4, 0.5, 0.6], 100)
    draws[::2] += 1e-16  # last-bit noise on half the draws
    noisy = dataclasses.replace(result, bootstrap_distribution=draws)
    assert "3 distinct values" in _caption(plot_bootstrap_distribution(noisy))


# -- Drawing into a given axes ---------------------------------------------- #


def test_draws_into_given_axes(ratings):
    fig, (left, right) = plt.subplots(1, 2)
    returned = plot_bootstrap_distribution(
        _fleiss(ratings, n_resamples=300, ci_method="percentile"), ax=right,
    )
    assert returned is fig
    assert right.get_title() and not left.get_title()


# -- Refusals --------------------------------------------------------------- #


# parametric_null warns twice about its width; expected here, not the subject.
@pytest.mark.filterwarnings("ignore:.*parametric_null:UserWarning")
@pytest.mark.filterwarnings("ignore:.*scaled as though:UserWarning")
def test_parametric_null_refused(ratings):
    with pytest.raises(ValueError, match="no bootstrap distribution"):
        plot_bootstrap_distribution(_fleiss(ratings, ci_method="parametric_null"))


def test_all_draws_undefined_refused(ratings):
    result = _fleiss(ratings, n_resamples=300, ci_method="percentile")
    empty = dataclasses.replace(
        result, bootstrap_distribution=np.full(300, np.nan),
    )
    with pytest.raises(ValueError, match="every one of the 300"):
        plot_bootstrap_distribution(empty)


def test_not_a_result_refused():
    with pytest.raises(TypeError, match="result must be"):
        plot_bootstrap_distribution(np.zeros(10))
