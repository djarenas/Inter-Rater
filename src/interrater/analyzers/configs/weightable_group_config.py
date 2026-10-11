"""Configuration for group coefficients that support weighting."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .group_coefficient_config import GroupCoefficientConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from .hypothesis_test_config import HypothesisTestConfig
    from ..convergence_test import ConvergenceTest

_ERROR_PREFIX = "[WeightableGroupConfig] Error"


@dataclass(frozen=True)
class WeightableGroupConfig(GroupCoefficientConfig):
    """
    I am the configuration for group coefficients that support weighting.

    These coefficients (Fleiss' kappa, Scott's pi) have weighted forms where
    disagreement can be scored by distance or custom metric. Weighting is
    especially useful for ordinal data where some disagreements matter more
    than others.

    I add weight configuration: which type of weights to use (identity, linear,
    quadratic, or custom) and where to load custom weights if needed. The
    analyzer enforces that if weights_type="custom", either a custom_weights_path
    is set in config or a Weights object is passed to the constructor.

    Attributes
    ----------
    weights_type : str, default "identity"
        Type of weights to use. "identity" means no weights (unweighted
        coefficient). "linear" applies linear distance weights for ordinal
        data. "quadratic" applies squared distance weights. "custom" means
        load from custom_weights_path or pass Weights object to analyzer.
    custom_weights_path : str | None, default None
        Path to custom weights definition (YAML or other format). Only used
        if weights_type="custom". If None but weights_type="custom", the
        analyzer will require a Weights object in constructor.
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
    # Populated at this level
    weights_type: str = "identity"
    custom_weights_path: str | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate weight-specific configuration."""

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
