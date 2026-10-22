"""Result for Cohen's kappa analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .individual_marginal_result import IndividualMarginalResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[CohenKappaResult] Error"


@dataclass(frozen=True)
class CohenKappaResult(IndividualMarginalResult):
    """
    I am Cohen's kappa, weighted or unweighted, for a single pair of raters.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. I add CI method and weights_type
    so anyone reading me knows what statistic I am.

    weights_type is on me (not just on the config) because weighted Cohen's
    kappa and unweighted Cohen's kappa are different statistics with different
    interpretations. A result must be self-describing.

    Cohen's kappa operates on paired data (two raters, fully crossed), so
    balance is guaranteed by construction. My balance_profile is always
    populated and always balanced. n_subjects is available via
    balance_profile.n_subjects.

    I carry no rater names. Who the raters are is input metadata, available
    on the Ratings object the caller passed in. I describe what the data said,
    not who said it.

    References
    ----------
    Cohen, J. (1960). A coefficient of agreement for nominal scales.
        Educational and Psychological Measurement, 20(1), 37--46.
    Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision
        for scaled disagreement or partial credit. Psychological Bulletin,
        70(4), 213--220.

    Attributes
    ----------
    point : float
        Cohen's kappa (inherited from PointEstimateResult).
    confidence_interval : tuple[float, float] | None
        Bootstrap or parametric confidence interval (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Always populated for
        Cohen's kappa. n_subjects available via balance_profile.n_subjects.
    ci_method : str
        Method used: "percentile", "bca", or "parametric_null".
    weights_type : str
        Weighting scheme: "identity" (unweighted), "linear", "quadratic",
        or "custom". Determines what statistic this result represents.
    weight_matrix : list[list[float]] | None, default None
        The weight matrix that was used, as nested lists (row i, column j is
        the credit for categories i and j). None when unweighted. Recorded so
        a report shows the numbers, not only the name of the scheme.
    """

    # Attribute declaration

    # Inherited from IndividualMarginalResult (which validates
    # balance_profile is not None and n_raters is not None)
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
        """Validate Cohen's kappa result attributes."""

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

        # Call parent validation (IndividualMarginalResult validates
        # balance_profile and n_raters; PointEstimateResult validates
        # point, confidence_interval, bootstrap_distribution)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Return a compact string representation with coefficient name.

        Format:
            Cohen's κ = 0.723 (95% CI: [0.681, 0.765], n_bootstrap=1000)
            n_subjects=50; unweighted; ci_method=percentile
        """
        parent_str = super().__str__()
        weight_label = (
            "unweighted" if self.weights_type == "identity"
            else f"weighted ({self.weights_type})"
        )
        return (
            f"Cohen's κ {parent_str}\n"
            f"n_subjects={self.balance_profile.n_subjects}; {weight_label}; "
            f"ci_method={self.ci_method}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and Cohen-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'ci_method', 'weights_type', 'weight_matrix'.
        """
        result = super().summary()
        result.update({
            "ci_method": self.ci_method,
            "weights_type": self.weights_type,
            "weight_matrix": self.weight_matrix,
        })
        return result
