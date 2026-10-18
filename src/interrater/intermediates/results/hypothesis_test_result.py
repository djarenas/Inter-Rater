"""Result of a hypothesis test for inter-rater agreement coefficients."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_ERROR_PREFIX = "[HypothesisTestResult] Error"


@dataclass(frozen=True)
class HypothesisTestResult:
    """
    I hold the result of a hypothesis test against H₀: coefficient = 0.

    I store the test outcome and enough metadata to interpret it.
    Different methods produce different statistics -- z_statistic
    is only meaningful for parametric_null, while null_distribution
    is only meaningful for permutation.

    I am not on the AnalysisResult tree. A hypothesis test is a
    fundamentally different kind of output from a point estimate
    with confidence interval.

    Attributes
    ----------
    p_value : float
        The p-value of the test (two-sided). Must be in [0, 1].
    alpha : float
        The significance level used. Must be in (0, 1).
    reject_null : bool
        Whether H₀ is rejected at the given alpha. Must be consistent
        with p_value and alpha.
    method : str
        The method used: "permutation" or "parametric_null".
    point : float
        The observed point estimate being tested.
    z_statistic : float | None
        Z-score (parametric_null only, None for permutation).
    null_variance : float | None
        The null variance used (parametric_null only). Must be >= 0.
    null_distribution : np.ndarray | None
        The permutation distribution (permutation only, None for parametric).
    n_permutations : int | None
        Number of permutations performed (permutation only). Must be >= 1.
    """

    # Required fields
    p_value: float
    alpha: float
    reject_null: bool
    method: str
    point: float

    # Parametric-null fields (None for permutation)
    z_statistic: float | None = None
    null_variance: float | None = None

    # Permutation fields (None for parametric)
    null_distribution: np.ndarray | None = None
    n_permutations: int | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate hypothesis test result attributes."""

        # Validate p_value
        if not isinstance(self.p_value, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: p_value must be float, "
                f"got {type(self.p_value).__name__}"
            )
        if not (0.0 <= self.p_value <= 1.0):
            raise ValueError(
                f"{_ERROR_PREFIX}: p_value must be in [0, 1], "
                f"got {self.p_value}"
            )

        # Validate alpha
        if not isinstance(self.alpha, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: alpha must be float, "
                f"got {type(self.alpha).__name__}"
            )
        if not (0.0 < self.alpha < 1.0):
            raise ValueError(
                f"{_ERROR_PREFIX}: alpha must be in (0, 1), "
                f"got {self.alpha}"
            )

        # Validate reject_null
        if not isinstance(self.reject_null, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: reject_null must be bool, "
                f"got {type(self.reject_null).__name__}"
            )

        # Validate reject_null consistency with p_value and alpha
        # (boundary case p_value == alpha accepts either decision)
        if self.p_value < self.alpha and not self.reject_null:
            raise ValueError(
                f"{_ERROR_PREFIX}: reject_null is False but "
                f"p_value ({self.p_value}) < alpha ({self.alpha})"
            )
        if self.p_value > self.alpha and self.reject_null:
            raise ValueError(
                f"{_ERROR_PREFIX}: reject_null is True but "
                f"p_value ({self.p_value}) > alpha ({self.alpha})"
            )

        # Validate method
        if not isinstance(self.method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: method must be str, "
                f"got {type(self.method).__name__}"
            )
        if self.method not in {"permutation", "parametric_null"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: method must be 'permutation' or "
                f"'parametric_null', got '{self.method}'"
            )

        # Validate point
        if not isinstance(self.point, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: point must be float, "
                f"got {type(self.point).__name__}"
            )

        # Cross-validate method-specific fields
        if self.method == "parametric_null":
            if self.z_statistic is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: z_statistic must not be None "
                    f"for method='parametric_null'"
                )
            if self.null_variance is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: null_variance must not be None "
                    f"for method='parametric_null'"
                )
        if self.method == "permutation":
            if self.n_permutations is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: n_permutations must not be None "
                    f"for method='permutation'"
                )

        # Validate optional field values when present
        if self.z_statistic is not None:
            if not isinstance(self.z_statistic, (int, float)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: z_statistic must be float, "
                    f"got {type(self.z_statistic).__name__}"
                )
        if self.null_variance is not None:
            if not isinstance(self.null_variance, (int, float)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: null_variance must be float, "
                    f"got {type(self.null_variance).__name__}"
                )
            if self.null_variance < 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: null_variance must be >= 0, "
                    f"got {self.null_variance}"
                )
        if self.null_distribution is not None:
            if not isinstance(self.null_distribution, np.ndarray):
                raise TypeError(
                    f"{_ERROR_PREFIX}: null_distribution must be np.ndarray, "
                    f"got {type(self.null_distribution).__name__}"
                )
        if self.n_permutations is not None:
            if not isinstance(self.n_permutations, int):
                raise TypeError(
                    f"{_ERROR_PREFIX}: n_permutations must be int, "
                    f"got {type(self.n_permutations).__name__}"
                )
            if self.n_permutations < 1:
                raise ValueError(
                    f"{_ERROR_PREFIX}: n_permutations must be >= 1, "
                    f"got {self.n_permutations}"
                )

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Human-readable summary of the hypothesis test result.

        Examples
        --------
        Parametric:
            H₀: point estimate = 0, z = 3.42, p = 0.0006, reject at α = 0.05
        Permutation:
            H₀: point estimate = 0, p = 0.0030 (5000 permutations),
            reject at α = 0.05
        """
        # Format p-value: use scientific notation for very small values
        if self.p_value < 0.000001:
            p_str = "p < 0.000001"
        else:
            p_str = f"p = {self.p_value:.6f}"

        # Decision
        decision = "reject" if self.reject_null else "fail to reject"

        if self.method == "parametric_null" and self.z_statistic is not None:
            return (
                f"H₀: point estimate = 0, z = {self.z_statistic:.2f}, "
                f"{p_str}, {decision} at α = {self.alpha}"
            )
        elif self.method == "permutation" and self.n_permutations is not None:
            return (
                f"H₀: point estimate = 0, {p_str} "
                f"({self.n_permutations} permutations), "
                f"{decision} at α = {self.alpha}"
            )
        else:
            return (
                f"H₀: point estimate = 0, {p_str}, method={self.method}, "
                f"{decision} at α = {self.alpha}"
            )

    def summary(self) -> dict:
        """
        Return a plain-type summary of the hypothesis test result.

        Returns
        -------
        dict
            Keys: p_value, alpha, reject_null, method, point,
            and method-specific fields.
        """
        result = {
            "p_value": self.p_value,
            "alpha": self.alpha,
            "reject_null": self.reject_null,
            "method": self.method,
            "point": self.point,
        }

        if self.z_statistic is not None:
            result["z_statistic"] = self.z_statistic
        if self.null_variance is not None:
            result["null_variance"] = self.null_variance
        if self.n_permutations is not None:
            result["n_permutations"] = self.n_permutations

        return result
