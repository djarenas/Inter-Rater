"""Result for any coefficient that produces a single point estimate."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..analysis_result import AnalysisResult
from ..ci_results.ci_result import CIResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[PointEstimateResult] Error"
_WARNING_PREFIX = "[PointEstimateResult] Warning"

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
    When one was computed it arrives as a BootstrapCIResult, which carries not
    only the bounds but which method was asked for, which one ran, and how many
    resamples failed.

    I do not require my point estimate to fall inside my interval. A strongly
    bias-corrected BCa interval can sit entirely to one side of it -- that is
    correct behaviour, not a contradiction, and refusing to exist would be the
    wrong response to it. I warn instead, loudly and three ways, because it is
    also the signature of a degenerate bootstrap and a reader should not have
    to notice it for themselves.

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
    confidence_interval : BootstrapCIResult | None
        The interval and an account of how it was produced: bounds, alpha,
        the method asked for against the one that ran, and the number of
        resamples that failed. None if no interval was computed.
    bootstrap_distribution : np.ndarray | None
        Array of shape (n_resamples,) containing coefficient values from each
        bootstrap resample. None if bootstrap was not performed or not saved.
        Useful for plotting and diagnostics.
    balance_profile : BalanceProfile | None
        Structural metadata about the data: n_subjects, n_raters, balance
        structure. None if not available at construction time.
    """

    point: float
    confidence_interval: BootstrapCIResult | None
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

        # If confidence_interval is provided, it must be a BootstrapCIResult.
        # Its own bounds, counts and methods were validated when it was built,
        # so there is nothing here to re-check -- only the relationship between
        # it and my point estimate, which it cannot see.
        if self.confidence_interval is not None:
            if not isinstance(self.confidence_interval, CIResult):
                raise TypeError(
                    f"{_ERROR_PREFIX}: confidence_interval must be a CIResult "
                    f"or None, got {type(self.confidence_interval).__name__}"
                )

            # The point estimate falling outside the interval is permitted --
            # strong bias correction does that legitimately -- but it is also
            # what a degenerate bootstrap looks like, so it is never silent.
            # The interval describes itself here rather than being inspected,
            # so this message stays right whatever kind of interval it is.
            if self.point_outside_interval:
                warnings.warn(
                    f"{_WARNING_PREFIX}: the point estimate {self.point:.4f} "
                    f"falls outside its own interval "
                    f"{self.confidence_interval}. This is legitimate when "
                    f"bias correction is strong, and it is also what a "
                    f"degenerate bootstrap looks like. Check which before "
                    f"interpreting the interval.",
                    UserWarning,
                    stacklevel=2,
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
            from ...data_profiles.balance_profile import BalanceProfile

            if not isinstance(self.balance_profile, BalanceProfile):
                raise TypeError(
                    f"{_ERROR_PREFIX}: balance_profile must be BalanceProfile or None, "
                    f"got {type(self.balance_profile).__name__}"
                )

    # Public API ------------------------------------------------------------- #

    @property
    def point_outside_interval(self) -> bool:
        """
        Whether my point estimate falls outside my own interval.

        False when there is no interval, and when the interval exists but has
        no bounds because nothing ran. Legitimate under strong bias
        correction; also the signature of a degenerate bootstrap.
        """
        ci = self.confidence_interval
        if ci is None or not ci.is_defined:
            return False
        return not (ci.lower <= self.point <= ci.upper)

    def __str__(self) -> str:
        """
        Return a compact string representation of the result.

        Format: = 0.523 [0.481, 0.565] (95% CI, bca, n=1000)

        My interval describes itself -- I do not inspect it, so whatever
        qualifies it travels with it unchanged, whether that is a fallback,
        a count of failed resamples, or a caveat about the variance it was
        built from. I add only what it cannot know, which is where my point
        estimate sits relative to it:

            = 1.000 [0.812, 0.961] (95% CI, bca, n=1000, 61 failed) POINT OUTSIDE INTERVAL
            = 0.523 [0.401, 0.639] (95% CI, null variance -- width assumes kappa = 0)
            = 0.523 (n_bootstrap=1000)
        """
        parts = [f"= {self.point:.3f}"]

        if self.confidence_interval is not None:
            parts.append(str(self.confidence_interval))
            if self.point_outside_interval:
                parts.append("POINT OUTSIDE INTERVAL")
        elif self.bootstrap_distribution is not None:
            # Only worth saying when no interval came of it
            parts.append(f"(n_bootstrap={len(self.bootstrap_distribution)})")

        return " ".join(parts)

    def summary(self) -> dict:
        """
        Return my point estimate, interval, bootstrap diagnostics, and
        balance metadata in plain types.

        Returns
        -------
        dict
            Keys: 'point' (float); if confidence_interval is not None, 'ci'
            (the interval's own summary dict) plus 'ci_lower', 'ci_upper' and
            'point_outside_interval' hoisted to the top level, since those are
            what a reader scans for; optionally 'bootstrap_n',
            'bootstrap_mean', 'bootstrap_std' if bootstrap_distribution is
            not None, and optionally 'balance' if balance_profile is not None.
        """
        result = {"point": self.point}
        if self.confidence_interval is not None:
            result["ci"] = self.confidence_interval.summary()
            result["ci_lower"] = self.confidence_interval.lower
            result["ci_upper"] = self.confidence_interval.upper
            result["point_outside_interval"] = self.point_outside_interval
        if self.bootstrap_distribution is not None:
            result["bootstrap_n"] = len(self.bootstrap_distribution)
            result["bootstrap_mean"] = float(np.mean(self.bootstrap_distribution))
            result["bootstrap_std"] = float(np.std(self.bootstrap_distribution))
        if self.balance_profile is not None:
            result["balance"] = self.balance_profile.summary()
        return result
