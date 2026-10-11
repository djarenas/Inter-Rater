"""Result for Fleiss' kappa analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[FleissKappaResult] Error"


@dataclass(frozen=True)
class FleissKappaResult(PointEstimateResult):
    """
    I am Fleiss' kappa, exact or approximate, weighted or unweighted,
    with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers whether this is the exact coefficient from Fleiss 1971
    or a generalization applied to unbalanced data.

    When is_exact is True, I am the exact formula (balanced data, constant
    raters per subject). When False, I am a valid generalization using
    per-subject averaging of p_j, but not the original coefficient.

    weights_type is on me (not just on the config) because weighted and
    unweighted Fleiss' kappa are different statistics with different
    interpretations. A result must be self-describing.


    Attributes
    ----------
    point : float
        Fleiss' kappa (inherited from PointEstimateResult).
    confidence_interval : tuple[float, float] | None
        Bootstrap confidence interval (lower, upper) (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Its is_balanced property
        determines whether this result is exact or approximate.
    ci_method : str
        Method used: "percentile", "bca", or "parametric_null".
    weights_type : str
        Weighting scheme: "identity" (unweighted), "linear", "quadratic",
        or "custom". Determines what statistic this result represents.
    """

    # Attribute declaration

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: tuple[float, float] | None
    bootstrap_distribution: np.ndarray | None
    balance_profile: BalanceProfile | None

    # Populated at this level
    ci_method: str
    weights_type: str

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Fleiss' kappa result attributes."""

        # Fleiss always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for Fleiss' kappa"
            )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca", "parametric_null"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile', 'bca', "
                f"or 'parametric_null', got '{self.ci_method}'"
            )

        # Validate weights_type
        if not isinstance(self.weights_type, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: weights_type must be str, "
                f"got {type(self.weights_type).__name__}"
            )
        if self.weights_type not in {"identity", "linear", "quadratic", "custom"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: weights_type must be 'identity', 'linear', "
                f"'quadratic', or 'custom', got '{self.weights_type}'"
            )

        # Call parent validation (validates point, confidence_interval,
        # bootstrap_distribution, balance_profile type)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Return a compact string representation with coefficient name and context.

        Format:
            Fleiss' κ = 0.523 (95% CI: [0.481, 0.565], n_bootstrap=1000)
            n_subjects=100; n_raters=5; unweighted; ci_method=bca; exact formula
        """
        parent_str = super().__str__()
        weight_label = (
            "unweighted" if self.weights_type == "identity"
            else f"weighted ({self.weights_type})"
        )
        exact_str = "exact formula" if self.is_exact else "approximation (unbalanced data)"
        return (
            f"Fleiss' κ {parent_str}\n"
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"{self.balance_profile.n_raters_display}; {weight_label}; "
            f"ci_method={self.ci_method}; {exact_str}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and Fleiss-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'is_exact', 'ci_method',
            'weights_type'.
        """
        result = super().summary()
        result.update({
            "is_exact": self.is_exact,
            "ci_method": self.ci_method,
            "weights_type": self.weights_type,
        })
        return result

    @property
    def is_exact(self) -> bool:
        """
        Whether this is the exact Fleiss formula or an approximation.

        Returns True if data was balanced (constant raters per subject).
        Returns False if an approximation was used for unbalanced data.
        """
        return self.balance_profile.is_balanced
