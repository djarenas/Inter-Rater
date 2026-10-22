"""Configuration for group coefficients that do not support weighting."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .group_coefficient_config import GroupCoefficientConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[NoWeightGroupConfig] Error"


@dataclass(frozen=True)
class NoWeightGroupConfig(GroupCoefficientConfig):
    """
    I am the base configuration for group coefficients that do not support weighting.

    These coefficients (Gwet's AC1, Scott's pi) have no weighted forms
    because their chance terms are defined independently of how raters use
    categories. Weighting would create an inconsistency: weighted observed
    agreement but unweighted chance term.

    I add no fields. The analyzers for these coefficients refuse ratings whose
    semantics declares ordered or weighted categories, with a clear error.

    Weights are not set here. They are declared in the semantics of the
    ratings (categories of kind 'ordered' or 'weighted'), and the analyzer
    reads them from there. Before 2.1.0 this class carried weights_type and
    custom_weights_path; both were removed.

    Leaf classes (GwetAC1Config, ScottsConfig) inherit from me without
    adding new fields, existing only for type clarity.

    Attributes
    ----------
    (all fields inherited from GroupCoefficientConfig and AnalyzerConfig)
    """

    # Attribute declaration

    # Inherited from parent class AnalyzerConfig
    n_resamples: int
    alpha: float
    random_state: int | None
    ci_method: str
    convergence_test: ConvergenceTest | None
    hypothesis_test_config: HypothesisTestConfig | None
    injection: tuple[InjectorConfig] | None
    # Inherited from parent class GroupCoefficientConfig
    balanced_data_only: bool
    rater_jackknife: bool

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Run the parent's validation, naming me in any error."""
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc
