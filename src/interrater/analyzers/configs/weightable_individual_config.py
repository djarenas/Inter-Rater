"""Configuration for individual marginal coefficients that have a weighted form."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .individual_marginal_config import IndividualMarginalConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[WeightableIndividualConfig] Error"


@dataclass(frozen=True)
class WeightableIndividualConfig(IndividualMarginalConfig):
    """
    I am the configuration for individual marginal coefficients that have a weighted form.

    Cohen's kappa has a weighted form for ordinal data (Cohen 1968). I add no
    fields. I exist so that its analyzer accepts ratings whose semantics
    declares weights.

    Weights are not set here. They are declared in the semantics of the
    ratings (categories of kind 'ordered' or 'weighted'), and the analyzer
    reads them from there. Before 2.1.0 this class carried weights_type and
    custom_weights_path; both were removed.

    Leaf classes (CohenKappaConfig) inherit from me without adding new fields,
    existing only for type clarity.

    Attributes
    ----------
    (all fields inherited from IndividualMarginalConfig and AnalyzerConfig)
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
    # Inherited from parent class IndividualMarginalConfig
    rater_jackknife: bool

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Run the parent's validation, naming me in any error."""
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc
