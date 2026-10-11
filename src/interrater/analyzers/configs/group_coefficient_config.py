"""Configuration for group coefficients (Fleiss, AC1, BP, Scott, Krippendorff)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .analyzer_config import AnalyzerConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[GroupCoefficientConfig] Error"


@dataclass(frozen=True)
class GroupCoefficientConfig(AnalyzerConfig):
    """
    I am the base configuration for group coefficients.

    Group coefficients (Fleiss' kappa, Gwet AC1, Brennan-Prediger, Scott's pi,
    Krippendorff's alpha) read a counts matrix (NijMatrix) that discards rater
    identity. 


    I add two decisions specific to group analysis:
    - Whether to enforce balance (constant raters per subject). Fleiss requires
      this mathematically; others are more flexible.
    - Whether to run leave-one-rater-out diagnostics to see how much each
      rater matters to the agreement.

    Attributes
    ----------
    balanced_data_only : bool, default False
        If True, require that every subject was rated the same number of times.
        If False, allow unbalanced designs (different subjects rated by
        different numbers of raters). Fleiss' kappa requires balance; other
        group coefficients are flexible.
    rater_jackknife : bool, default False
        If True, compute leave-one-rater-out diagnostics: drop each rater in
        turn and recompute the coefficient. Useful for detecting which rater
        contributes most to agreement or causes instability.
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
    # Populated at this level
    balanced_data_only: bool = False
    rater_jackknife: bool = False

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate group-specific configuration."""

        # Call parent validation first: n_resamples, alpha,...
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate balanced_data_only
        if not isinstance(self.balanced_data_only, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: balanced_data_only must be bool, "
                f"got {type(self.balanced_data_only).__name__}"
            )

        # Validate rater_jackknife
        if not isinstance(self.rater_jackknife, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: rater_jackknife must be bool, "
                f"got {type(self.rater_jackknife).__name__}"
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration including group-specific fields."""
        # Overrides the parent class method, but still utilizes it
        result = super().summary() # Returns dictionary

        result.update({
            "balanced_data_only": self.balanced_data_only,
            "rater_jackknife": self.rater_jackknife,
        })
        return result
