"""Result for Gwet's AC1 analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[GwetsAC1Result] Error"


@dataclass(frozen=True)
class GwetsAC1Result(PointEstimateResult):
    """
    I am Gwet's AC1, exact or approximate, with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers whether this is the exact coefficient or a generalization
    applied to unbalanced data.

    When is_exact is True, I am the exact formula (balanced data, constant
    raters per subject). When False, I am a valid generalization using
    per-subject averaging of p_j, but not the original coefficient.

    I exist as a separate class (rather than just using PointEstimateResult)
    for type clarity. This is especially useful for error chaining.

    Attributes
    ----------
    point : float
        Gwet's AC1 (inherited from PointEstimateResult).
    confidence_interval : tuple[float, float] | None
        Bootstrap confidence interval (lower, upper) (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Its is_balanced property
        determines whether this result is exact or approximate.
    ci_method : str
        Method used: "percentile" or "bca".
    """

    # Attribute declaration

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: tuple[float, float] | None
    bootstrap_distribution: np.ndarray | None
    balance_profile: BalanceProfile | None

    # Populated at this level
    ci_method: str

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Gwet's AC1 result attributes."""

        # Gwet's AC1 always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for Gwet's AC1"
            )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method 'parametric_null' is not implemented "
                f"for Gwet's AC1"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or 'bca', "
                f"got '{self.ci_method}'"
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
            Gwet's AC1 = 0.523 (95% CI: [0.481, 0.565], n_bootstrap=1000)
            n_subjects=100; n_raters=5; ci_method=bca; exact formula
        """
        parent_str = super().__str__()
        exact_str = "exact formula" if self.is_exact else "approximation (unbalanced data)"
        return (
            f"Gwet's AC1 {parent_str}\n"
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"{self.balance_profile.n_raters_display}; "
            f"ci_method={self.ci_method}; {exact_str}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and AC1-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'is_exact', 'ci_method'.
        """
        result = super().summary()
        result.update({
            "is_exact": self.is_exact,
            "ci_method": self.ci_method,
        })
        return result

    @property
    def is_exact(self) -> bool:
        """
        Whether this is the exact AC1 formula or an approximation.

        Returns True if data was balanced (constant raters per subject).
        Returns False if an approximation was used for unbalanced data.
        """
        return self.balance_profile.is_balanced
