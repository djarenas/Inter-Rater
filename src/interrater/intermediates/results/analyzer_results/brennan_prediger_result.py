"""Result for Brennan-Prediger analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING
import numpy as np

from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile


_ERROR_PREFIX = "[BrennanPredigerResult] Error"


@dataclass(frozen=True)
class BrennanPredigerResult(PointEstimateResult):
    """
    I am Brennan-Prediger's coefficient, exact or approximate, weighted or
    unweighted, with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers whether this is the exact coefficient or a generalization
    applied to unbalanced data.

    When is_exact is True, the data was balanced and the formula is the
    standard Brennan-Prediger with Pe = 1/K. When False, the per-subject
    averaging generalization was used -- a reasonable number, but not the
    original coefficient.

    weights_type is on me (not just on the config) because weighted and
    unweighted Brennan-Prediger are different statistics with different
    interpretations. A result must be self-describing.

    Attributes
    ----------
    point : float
        Brennan-Prediger coefficient (inherited from PointEstimateResult).
    confidence_interval : tuple[float, float] | None
        Bootstrap confidence interval (lower, upper) (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Its is_balanced property
        determines whether this result is exact or approximate.
    ci_method : str
        Method used: "percentile" or "bca".
    weights_type : str
        Weighting scheme: "identity" (unweighted), "linear", "quadratic",
        or "custom". Determines what statistic this result represents.
    weight_matrix : list[list[float]] | None, default None
        The weight matrix that was used, as nested lists (row i, column j is
        the credit for categories i and j). None when unweighted. Recorded so
        a report shows the numbers, not only the name of the scheme.
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
    weight_matrix: list[list[float]] | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Brennan-Prediger result attributes."""

        # Brennan-Prediger always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for "
                f"Brennan-Prediger"
            )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or 'bca', "
                f"got '{self.ci_method}'"
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

        # Validate weight_matrix: none when unweighted, square when given
        if self.weight_matrix is not None:
            if self.weights_type == "identity":
                raise ValueError(
                    f"{_ERROR_PREFIX}: weight_matrix must be None when "
                    f"weights_type is 'identity' (unweighted)"
                )
            matrix = np.asarray(self.weight_matrix, dtype=float)
            if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
                raise ValueError(
                    f"{_ERROR_PREFIX}: weight_matrix must be square, "
                    f"got shape {matrix.shape}"
                )

        # Call parent validation (validates point, confidence_interval,
        # bootstrap_distribution, balance_profile type)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API -------------------------------------------------------------  #

    def __str__(self) -> str:
        """
        Return a compact string representation with coefficient name and context.

        Format:
            Brennan-Prediger = 0.835 (95% CI: [0.764, 0.919], n_bootstrap=4000)
            n_subjects=100; n_raters=5; unweighted; ci_method=bca; exact formula
        """
        parent_str = super().__str__()
        weight_label = (
            "unweighted" if self.weights_type == "identity"
            else f"weighted ({self.weights_type})"
        )
        exact_str = "exact formula" if self.is_exact else "approximation (unbalanced data)"
        return (
            f"Brennan-Prediger {parent_str}\n"
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"{self.balance_profile.n_raters_display}; {weight_label}; "
            f"ci_method={self.ci_method}; {exact_str}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and BP-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'is_exact', 'ci_method',
            'weights_type', 'weight_matrix'.
        """
        result = super().summary()
        result.update({
            "is_exact": self.is_exact,
            "ci_method": self.ci_method,
            "weights_type": self.weights_type,
            "weight_matrix": self.weight_matrix,
        })
        return result

    @property
    def is_exact(self) -> bool:
        """
        Whether this is the exact formula or an approximation.

        Returns True if data was balanced (constant raters per subject).
        Returns False if an approximation was used for unbalanced data.
        """
        return self.balance_profile.is_balanced
