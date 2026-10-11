"""Result for Krippendorff's alpha analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING
import numpy as np

from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile


_ERROR_PREFIX = "[KrippendorffAlphaResult] Error"


@dataclass(frozen=True)
class KrippendorffAlphaResult(PointEstimateResult):
    """
    I am Krippendorff's alpha, weighted or unweighted, with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers about the balance structure of the data.

    Unlike the kappa family, alpha handles unbalanced data natively via
    pair-pooling -- a subject rated by six raters contributes fifteen pairs
    and one rated by two contributes a single pair. There is no
    approximation involved: the unbalanced formula *is* the formula. So
    is_exact is always True, and the balance profile is informational
    context rather than a correctness flag.

    Alpha measures disagreement rather than agreement:

        alpha = 1 - Do / De

    This is not the same expression as (Po - Pe) / (1 - Pe) rearranged.
    The two use different units (pairs pooled vs. subjects averaged) and
    different chance models (finite-sample permutation vs. population
    rates).

    weights_type is on me (not just on the config) because weighted and
    unweighted alpha are different statistics with different
    interpretations. A result must be self-describing.

    Attributes
    ----------
    point : float
        Krippendorff's alpha (inherited from PointEstimateResult).
    confidence_interval : tuple[float, float] | None
        Bootstrap confidence interval (lower, upper) (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Informational only --
        alpha handles both balanced and unbalanced data exactly.
    ci_method : str
        Method used: "percentile" or "bca".
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
        """Validate Krippendorff's alpha result attributes."""

        # Krippendorff's alpha always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for "
                f"Krippendorff's alpha"
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
            Krippendorff's alpha = 0.625 (95% CI: [0.551, 0.681], n_bootstrap=4000)
            n_subjects=100; n_raters=5; unweighted; ci_method=bca; exact formula
        or (unbalanced):
            Krippendorff's alpha = 0.625 (95% CI: [0.551, 0.681], n_bootstrap=4000)
            n_subjects=100; n_raters=5; unweighted; ci_method=bca; exact formula (unbalanced data)
        """
        parent_str = super().__str__()
        weight_label = (
            "unweighted" if self.weights_type == "identity"
            else f"weighted ({self.weights_type})"
        )
        if self.balance_profile.is_balanced:
            balance_str = "exact formula"
        else:
            balance_str = "exact formula (unbalanced data)"
        return (
            f"Krippendorff's alpha {parent_str}\n"
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"{self.balance_profile.n_raters_display}; {weight_label}; "
            f"ci_method={self.ci_method}; {balance_str}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and alpha-specific metadata.

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
        Whether this is the exact formula or an approximation.

        Always returns True for Krippendorff's alpha. Alpha handles
        unbalanced data natively via pair-pooling -- the unbalanced
        formula is the formula, not an approximation. The balance
        profile is available for informational context but does not
        affect the correctness of the result.
        """
        return True
