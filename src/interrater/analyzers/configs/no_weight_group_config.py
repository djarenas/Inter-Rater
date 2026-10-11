"""Configuration for group coefficients that do not support weighting."""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .group_coefficient_config import GroupCoefficientConfig

if TYPE_CHECKING:
    from ...injectors import InjectorConfig
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

    I include weights_type and custom_weights_path fields (inherited from
    subclasses loading from YAML) but validate in __post_init__ that they are
    not used. If someone tries to set weights in a YAML config for a
    non-weightable coefficient, they get a clear error message.

    Leaf classes (GwetAC1Config, ScottsConfig) inherit from me without
    adding new fields, existing only for type clarity.

    Attributes
    ----------
    weights_type : str, default "identity"
        Must be "identity" for non-weightable coefficients. Any other value
        raises ValueError in __post_init__.
    custom_weights_path : str | None, default None
        Must be None for non-weightable coefficients. Any other value raises
        ValueError in __post_init__.
    (all other fields inherited from GroupCoefficientConfig and AnalyzerConfig)
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
    # These must stay at defaults (identity + None) for non-weightable coefficients
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
