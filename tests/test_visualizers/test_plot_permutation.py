"""
Pairwise and subset plots: what the legend says about the mean's interval.

An "independent" mean interval ignores the correlation between combinations
that share raters, so it tends to be too narrow. Wherever it is drawn, the
legend says so; a "correlated" one carries no such note.
"""

import matplotlib

matplotlib.use("Agg")  # draw without a screen, as on CI

import matplotlib.pyplot as plt
import numpy as np
import pytest

from interrater import (
    plot_pairwise_permutation,
    plot_subset_sorted,
    plot_subset_strip,
)
from interrater.analyzers.configs.cohen_kappa_config import CohenKappaConfig
from interrater.analyzers.configs.fleiss_kappa_config import FleissKappaConfig
from interrater.analyzers.configs.pairwise_permutation_config import (
    PairwisePermutationConfig,
)
from interrater.analyzers.configs.permutation_bootstrap_mean_config import (
    PermutationBootstrapMeanConfig,
)
from interrater.analyzers.configs.subset_permutation_config import (
    SubsetPermutationConfig,
)
from interrater.analyzers.permutation_analysis import (
    PairwisePermutationAnalysis,
    SubsetPermutationAnalysis,
)
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics

from fixtures.through_cleaner import clean_ratings

# A few words of the caveat, as they appear once the legend wraps it
CAVEAT_WORDS = ("independent bootstrap", "too narrow")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def ratings():
    """40 subjects, 4 raters, 3 categories; each rater copies a true answer 70% of the time."""
    rng = np.random.default_rng(42)
    truth = rng.integers(0, 3, 40)
    noise = rng.integers(0, 3, (40, 4))
    copy = rng.random((40, 4)) < 0.7
    codes = np.where(copy, truth[:, None], noise).astype(np.int8)
    cats = UnorderedCategories(names=("A", "B", "C"))
    return clean_ratings(category_index=codes, semantics=RatingsSemantics(raters=4, categories=cats))


def _mean_config(method):
    return PermutationBootstrapMeanConfig(method=method)


def _pairwise(ratings, method):
    config = PairwisePermutationConfig(
        coefficient_to_compute="cohen_kappa",
        config_of_coefficient=CohenKappaConfig(n_resamples=50, random_state=42),
        bootstrap_mean_config=_mean_config(method),
    )
    return PairwisePermutationAnalysis(ratings, config).run()


def _subsets(ratings, method):
    config = SubsetPermutationConfig(
        coefficient_to_compute="fleiss_kappa",
        config_of_coefficient=FleissKappaConfig(n_resamples=50, random_state=42),
        subset_k=[2, 3],
        bootstrap_mean_config=_mean_config(method),
    )
    return SubsetPermutationAnalysis(ratings, config).run()


def _legend(fig):
    text = " ".join(t.get_text() for t in fig.axes[0].get_legend().get_texts())
    return " ".join(text.split())  # undo the wrapping


def _figures(ratings, method):
    subsets = _subsets(ratings, method)
    return {
        "pairwise": plot_pairwise_permutation(_pairwise(ratings, method)),
        "sorted": plot_subset_sorted(subsets[3]),
        "strip": plot_subset_strip(subsets),
    }


@pytest.mark.parametrize("plot", ["pairwise", "sorted", "strip"])
def test_independent_mean_carries_the_caveat(ratings, plot):
    legend = _legend(_figures(ratings, "independent")[plot])
    for words in CAVEAT_WORDS:
        assert words in legend


@pytest.mark.parametrize("plot", ["pairwise", "sorted", "strip"])
def test_correlated_mean_has_no_caveat(ratings, plot):
    legend = _legend(_figures(ratings, "correlated")[plot])
    assert "too narrow" not in legend


@pytest.mark.parametrize("plot, title, ylabel", [
    ("pairwise", "Pairwise Cohen's κ (correlated bootstrap)", "Cohen's κ"),
    ("sorted", "Subset k=3 Fleiss' κ — sorted by value (correlated bootstrap)", "Fleiss' κ"),
    ("strip", "Fleiss' κ across subset sizes (correlated bootstrap)", "Fleiss' κ"),
])
def test_titles_name_the_coefficient(ratings, plot, title, ylabel):
    """Named as the coefficient's result class names it, not 'cohen_kappa'."""
    ax = _figures(ratings, "correlated")[plot].axes[0]
    assert ax.get_title() == title
    assert ax.get_ylabel() == ylabel


def test_pairwise_labels_name_the_raters(ratings):
    """Rater columns and partner labels use the result's names, never "Rater 0"."""
    ax = _figures(ratings, "correlated")["pairwise"].axes[0]
    ticks = [t.get_text() for t in ax.get_xticklabels()]
    assert ticks == ["rater0", "rater1", "rater2", "rater3"]
    assert {t.get_text() for t in ax.texts} <= set(ticks)
