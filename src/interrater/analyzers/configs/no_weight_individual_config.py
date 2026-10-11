"""Configuration for individual marginal coefficients that do not support weighting."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING
from .individual_marginal_config import IndividualMarginalConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[NoWeightIndividualConfig] Error"


@dataclass(frozen=True)
class NoWeightIndividualConfig(IndividualMarginalConfig):
    """
    I am the configuration for individual marginal coefficients without weighting.

    These coefficients (Conger's kappa) have no weighted forms. I validate that
    weights are not requested, matching the pattern of NoWeightGroupConfig.

    Leaf classes (CongerKappaConfig) inherit from me without adding new fields,
    existing only for type clarity.

    Attributes
    ----------
    weights_type : str, default "identity"
        Must be "identity". Any other value raises ValueError in __post_init__.
    custom_weights_path : str | None, default None
        Must be None. Any other value raises ValueError in __post_init__.
    (all other fields inherited from IndividualMarginalConfig and AnalyzerConfig)
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
    # Populated at this level 
    # These must stay at defaults for non-weightable individual marginals
    weights_type: str = "identity"
    custom_weights_path: str | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate that no weights are requested."""

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        name = "weights_type"
        if not isinstance(getattr(self, name), str):
            raise TypeError(
                f"{_ERROR_PREFIX}: {name} must be a string, "
                f"got {type(getattr(self, name)).__name__}"
            )

        # This coefficient does not support weighting
        if self.weights_type != "identity":
            raise ValueError(
                f"{_ERROR_PREFIX}: This coefficient does not support weighting. "
                f"Set weights_type='identity', got '{self.weights_type}'"
            )

        if self.custom_weights_path is not None:
            raise ValueError(
                f"{_ERROR_PREFIX}: This coefficient does not support weighting. "
                f"custom_weights_path must be None, got '{self.custom_weights_path}'"
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration (version added by parent)."""
        return super().summary()
