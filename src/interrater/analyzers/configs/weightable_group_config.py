"""Configuration for group coefficients that have a weighted form."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .group_coefficient_config import GroupCoefficientConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[WeightableGroupConfig] Error"


@dataclass(frozen=True)
class WeightableGroupConfig(GroupCoefficientConfig):
    """
    I am the configuration for group coefficients that have a weighted form.

    Fleiss' kappa, Brennan-Prediger and Krippendorff's alpha can score a
    disagreement by how far apart the two categories are, which matters for
    ordinal data. I add no fields. I exist so that the analyzers for those
    coefficients, and only those, accept ratings whose semantics declares
    weights.

    Weights are not set here. They are declared in the semantics of the
    ratings (categories of kind 'ordered' or 'weighted'), and the analyzer
    reads them from there. Before 2.1.0 this class carried weights_type and
    custom_weights_path; both were removed.

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
