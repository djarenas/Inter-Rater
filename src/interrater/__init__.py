"""
interrater
==========
Agreement between raters, for rating tasks where not every rater sees every
subject.
"""

from interrater.contracts import EMPTY
from interrater.semantics.categories import (
    Categories,
    DistanceCategories,
    OrderedCategories,
    UnorderedCategories,
    WeightedCategories,
    build_categories,
)
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.data_objects.ratings import Ratings
from interrater.cleaners.cleaner_config import CleanerConfig
from interrater.cleaners.frame_layout import FrameLayout
from interrater.cleaners.insufficient_data_policy import InsufficientDataPolicy
from interrater.cleaners.ratings_cleaner import RatingsCleaner

from interrater.intermediates.estimator_parts import (Weights, NijMatrix)

from interrater.analyzers import FleissKappaAnalyzer, ScottsAnalyzer, GwetAC1Analyzer, BrennanPredigerAnalyzer, KrippendorffAlphaAnalyzer, CohenKappaAnalyzer, CongerKappaAnalyzer, LightKappaAnalyzer
from interrater.analyzers.permutation_analysis import PairwisePermutationAnalysis, SubsetPermutationAnalysis
from interrater.analyzers.configs import FleissKappaConfig, ScottsConfig, CustomWeightsInput, GwetAC1Config, BrennanPredigerConfig, KrippendorffAlphaConfig, CohenKappaConfig, CongerKappaConfig, LightKappaConfig, PairwisePermutationConfig, SubsetPermutationConfig, PermutationBootstrapMeanConfig, HypothesisTestConfig

from interrater.injectors import RaterInjector, WithinSubjectInjector, RaterBootstrapInjector, PooledInjector, DistributionInjector, ExtremalInjector
from interrater.injectors.configs import InjectorConfig, WithinSubjectInjectorConfig, RaterBootstrapInjectorConfig, PooledInjectorConfig, DistributionInjectorConfig, ExtremalInjectorConfig

from interrater.reports.cleaning_report import CleaningReport
from interrater.reports.matching_report import MatchingReport
from interrater.reports.drop_report import DropReport
from interrater.reports.category_report import CategoryReport

from importlib.metadata import version as _version, PackageNotFoundError

try:
    __version__ = _version("interrater")
except PackageNotFoundError:
    __version__ = "2.0.0"  # fallback for editable/dev installs


# -- Lazy plotting exports -------------------------------------------------- #
#
# The plotting functions import matplotlib, which costs roughly 2.5 seconds.
# Nobody computing a coefficient should pay that, so these three names are
# resolved on first use rather than at import time (PEP 562). The effect is
# that `from interrater import plot_pairwise_permutation` works exactly as if
# it were imported above, while `import interrater` never touches matplotlib.
#
# This also keeps matplotlib separable: if it ever moves to an optional
# extra, `import interrater` keeps working and only these names raise.

_LAZY_EXPORTS = {
    "plot_pairwise_permutation": "interrater.visualizers",
    "plot_subset_sorted": "interrater.visualizers",
    "plot_subset_strip": "interrater.visualizers",
}


def __getattr__(name: str):
    """Resolve a lazily exported name on first access, then cache it."""
    module_path = _LAZY_EXPORTS.get(name)
    if module_path is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    from importlib import import_module

    value = getattr(import_module(module_path), name)
    globals()[name] = value  # subsequent lookups skip __getattr__ entirely
    return value


def __dir__() -> list[str]:
    return sorted(__all__)


__all__ = [
    "EMPTY",
    "Categories",
    "UnorderedCategories",
    "OrderedCategories",
    "WeightedCategories",
    "DistanceCategories",
    "build_categories",
    "RatingsSemantics",
    "Ratings",
    "RatingsCleaner",
    "CleanerConfig",
    "FrameLayout",
    "InsufficientDataPolicy",
    "CleaningReport",
    "MatchingReport",
    "DropReport",
    "CategoryReport",
    # Analyzer objects
    "FleissKappaAnalyzer",
    "ScottsAnalyzer",
    "GwetAC1Analyzer",
    "BrennanPredigerAnalyzer",
    "KrippendorffAlphaAnalyzer",
    "CohenKappaAnalyzer",
    "CongerKappaAnalyzer",
    "LightKappaAnalyzer",
    # Permutation analyses
    "PairwisePermutationAnalysis",
    "SubsetPermutationAnalysis",
    # Configs for Analyzers
    "FleissKappaConfig",
    "ScottsConfig",
    "GwetAC1Config",
    "BrennanPredigerConfig",
    "KrippendorffAlphaConfig",
    "CohenKappaConfig",
    "CongerKappaConfig",
    "LightKappaConfig",
    "PairwisePermutationConfig",
    "SubsetPermutationConfig",
    "CustomWeightsInput",
    "Weights",
    "NijMatrix",
    # Bootstrap of the mean over rater combinations: "correlated" (default,
    # resamples subjects once per draw) or "independent"
    "PermutationBootstrapMeanConfig",
    # Hypothesis testing against H0: coefficient = 0
    "HypothesisTestConfig",
    # Synthetic rater injectors, ordered by how much real structure they keep:
    # Echo, Clone, Chance, Chance, Extremal (conformist/dissenter/contrarian/median)
    "WithinSubjectInjector",
    "RaterBootstrapInjector",
    "PooledInjector",
    "DistributionInjector",
    "ExtremalInjector",
    "RaterInjector",
    # Configs for injectors
    "WithinSubjectInjectorConfig",
    "RaterBootstrapInjectorConfig",
    "PooledInjectorConfig",
    "DistributionInjectorConfig",
    "ExtremalInjectorConfig",
    "InjectorConfig",
    # Plotting (imported lazily; see above)
    "plot_pairwise_permutation",
    "plot_subset_sorted",
    "plot_subset_strip",
]
