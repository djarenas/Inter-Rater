"""Configuration for Fleiss' kappa analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .weightable_group_config import WeightableGroupConfig

if TYPE_CHECKING:
    from ...injectors import InjectorConfig
    from .hypothesis_test_config import HypothesisTestConfig
    from ..convergence_test import ConvergenceTest

_ERROR_PREFIX = "[FleissKappaConfig] Error"


@dataclass(frozen=True)
class FleissKappaConfig(WeightableGroupConfig):
    """
    I am the configuration for Fleiss' kappa analysis.

    Fleiss' kappa works exactly on balanced data (constant raters per subject).
    On unbalanced data, an approximation is available using per-subject
    averaging of marginals.

    I default balanced_data_only=True (prefer exact formula), but the user can
    set it to False to allow approximations. The analyzer then reports what
    was used and warns if approximation occurred.

    Attributes
    ----------
    balanced_data_only : bool, default True
        If True, the analyzer will throw an error on unbalanced data (exact formula only).
        If False, the analyzer will allow unbalanced data (use approximation + warning).
    Fleiss' kappa supports weighting for ordinal data. Weights are declared
    in the semantics of the ratings (categories kind 'ordered' or 'weighted'),
    not here.
    (all other fields inherited from WeightableGroupConfig and parents)
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
    rater_jackknife: bool
    # Override parent default: Fleiss prefers balance by default
    balanced_data_only: bool = True

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Fleiss' kappa configuration."""

        # Call validation by parent classes
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Fleiss has a published closed-form null variance (Fleiss, Nee &
        # Landis 1979, Eq. 12), so parametric_null is permitted here -- and
        # warned about, since the variance it uses is the right one for a
        # test and the wrong one for an interval. Called last, and from here
        # rather than from the base class: the five configs that reject this
        # method do so after super(), so a base-class warning would fire for
        # them too. See AnalyzerConfig._warn_parametric_null.
        if self.ci_method == "parametric_null":
            self._warn_parametric_null()

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration (inherits parent's summary)."""
        return super().summary()
