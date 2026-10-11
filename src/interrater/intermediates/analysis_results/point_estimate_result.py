"""Result for any coefficient that produces a single point estimate."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .analysis_result import AnalysisResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[PointEstimateResult] Error"

@dataclass(frozen=True)
class PointEstimateResult(AnalysisResult):
    """
    I hold a point estimate, optionally a confidence interval, and optionally
    a bootstrap distribution.

    I am the result for any coefficient that produces a single number: Fleiss'
    kappa, Cohen's kappa, Conger's kappa, Gwet's AC1, Brennan-Prediger,
    Krippendorff's alpha, etc...    
    
    I am immutable (frozen dataclass) because a result is a snapshot of what the data 
    said at one moment and should never change.

    I do not compute my own confidence interval. That is the analyzer's job.
    If confidence_interval is None, no interval was computed (perhaps the user
    asked only for a point estimate, or the data was too sparse to resample).

    If bootstrap_distribution is provided, it contains all coefficient values
    computed during resampling. This is useful for visualization, diagnostics,
    and convergence checking.

    I carry a balance_profile so every result is self-describing: the number
    of subjects, the number of raters, whether the data was balanced, and the
    quality of rater-level information are all available without needing to
    go back to the input data. When balance_profile is None, this structural
    metadata was not available at construction time.

    Attributes
    ----------
    point : float
        The coefficient value.
    confidence_interval : tuple[float, float] | None
        (lower, upper) bounds of the confidence interval at the chosen alpha
        (conventionally 95%). None if no interval was computed.
    bootstrap_distribution : np.ndarray | None
        Array of shape (n_resamples,) containing coefficient values from each
        bootstrap resample. None if bootstrap was not performed or not saved.
        Useful for plotting and diagnostics.
    balance_profile : BalanceProfile | None
        Structural metadata about the data: n_subjects, n_raters, balance
        structure. None if not available at construction time.
    """

    point: float
    confidence_interval: tuple[float, float] | None
    bootstrap_distribution: np.ndarray | None
    balance_profile: BalanceProfile | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction, before I freeze."""

        # point must be a float
        if not isinstance(self.point, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: point must be float, got {type(self.point).__name__}"
            )

        # If confidence_interval is provided, validate it fully
        if self.confidence_interval is not None:
            # Must be a tuple of exactly 2 elements
            if not isinstance(self.confidence_interval, tuple) or len(self.confidence_interval) != 2:
                raise ValueError(
                    f"{_ERROR_PREFIX}: confidence_interval must be (float, float) tuple, "
                    f"got {type(self.confidence_interval).__name__}"
                )

            lower, upper = self.confidence_interval

            # Both bounds must be floats
            if not isinstance(lower, (int, float)) or not isinstance(upper, (int, float)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: CI bounds must be float, "
                    f"got ({type(lower).__name__}, {type(upper).__name__})"
                )

            # Lower bound must not exceed upper bound
            if lower > upper:
                raise ValueError(
                    f"{_ERROR_PREFIX}: CI lower bound {lower} > upper bound {upper}"
                )

            # Point estimate must fall within the confidence interval
            if not (lower <= self.point <= upper):
                raise ValueError(
                    f"{_ERROR_PREFIX}: point estimate {self.point} outside CI [{lower}, {upper}]"
                )

        # If bootstrap_distribution is provided, validate it
        if self.bootstrap_distribution is not None:
            if not isinstance(self.bootstrap_distribution, np.ndarray):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution must be np.ndarray or None, "
                    f"got {type(self.bootstrap_distribution).__name__}"
                )

            if self.bootstrap_distribution.ndim != 1:
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution must be 1-dimensional, "
                    f"got shape {self.bootstrap_distribution.shape}"
                )

            if len(self.bootstrap_distribution) == 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution must have at least 1 value"
                )

        # If balance_profile is provided, validate its type
        if self.balance_profile is not None:
            from ..data_profiles.balance_profile import BalanceProfile

            if not isinstance(self.balance_profile, BalanceProfile):
                raise TypeError(
                    f"{_ERROR_PREFIX}: balance_profile must be BalanceProfile or None, "
                    f"got {type(self.balance_profile).__name__}"
                )

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Return a compact string representation of the result.

        Format: = 0.523 (95% CI: [0.481, 0.565], n_bootstrap=1000)
        """
        meta = []
        if self.confidence_interval is not None:
            lower, upper = self.confidence_interval
            meta.append(f"95% CI: [{lower:.3f}, {upper:.3f}]")
        if self.bootstrap_distribution is not None:
            meta.append(f"n_bootstrap={len(self.bootstrap_distribution)}")

        suffix = f" ({', '.join(meta)})" if meta else ""
        return f"= {self.point:.3f}{suffix}"

    def summary(self) -> dict:
        """
        Return my point estimate, interval, bootstrap diagnostics, and
        balance metadata in plain types.

        Returns
        -------
        dict
            Keys: 'point' (float), optionally 'ci_lower', 'ci_upper' if
            confidence_interval is not None, optionally 'bootstrap_n',
            'bootstrap_mean', 'bootstrap_std' if bootstrap_distribution is
            not None, and optionally 'balance' if balance_profile is not None.
        """
        result = {"point": self.point}
        if self.confidence_interval is not None:
            result["ci_lower"], result["ci_upper"] = self.confidence_interval
        if self.bootstrap_distribution is not None:
            result["bootstrap_n"] = len(self.bootstrap_distribution)
            result["bootstrap_mean"] = float(np.mean(self.bootstrap_distribution))
            result["bootstrap_std"] = float(np.std(self.bootstrap_distribution))
        if self.balance_profile is not None:
            result["balance"] = self.balance_profile.summary()
        return result
