"""Configuration for individual marginal coefficients that support weighting."""

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
    I am the configuration for individual marginal coefficients that support weighting.

    These coefficients (Cohen's kappa) have weighted forms for ordinal data.
    I validate weights configuration and also ensure rater_jackknife stays False
    (Cohen is pairwise; you cannot drop one of two raters).

    Leaf classes (CohenKappaConfig) inherit from me without adding new fields,
    existing only for type clarity.

    Attributes
    ----------
    weights_type : str, default "identity"
        Type of weights to use. "identity", "linear", "quadratic", or "custom".
    custom_weights_path : str | None, default None
        Path to custom weights definition. Used only if weights_type="custom".
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
    weights_type: str = "identity"
    custom_weights_path: str | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate weight-specific and rater-jackknife configuration."""

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

        # Validate weights_type
        if self.weights_type not in {"identity", "linear", "quadratic", "custom"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: weights_type must be one of "
                f"'identity', 'linear', 'quadratic', 'custom', "
                f"got '{self.weights_type}'"
            )

        # Validate custom_weights_path type
        if self.custom_weights_path is not None:
            if not isinstance(self.custom_weights_path, str):
                raise TypeError(
                    f"{_ERROR_PREFIX}: custom_weights_path must be str or None, "
                    f"got {type(self.custom_weights_path).__name__}"
                )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration including weight-specific fields."""
        result = super().summary()
        result.update({
            "weights_type": self.weights_type,
            "custom_weights_path": self.custom_weights_path,
        })
        return result
