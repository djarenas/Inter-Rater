"""Result for a percentile bootstrap confidence interval."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ...not_computed import NotComputed
from .bootstrap_ci_result import BootstrapCIResult

_ERROR_PREFIX = "[PercentileBootstrapCIResult] Error"

# Bounds recomputed from the draws must match mine to this
_TOL = 1e-12


@dataclass(frozen=True)
class PercentileBootstrapCIResult(BootstrapCIResult):
    """
    I am a percentile bootstrap interval: my bounds are the alpha/2 and
    1 - alpha/2 quantiles of my usable draws.

    I am also what BCa returns when it cannot apply its corrections: then
    method_intended is 'bca', method_used is 'percentile', and my bounds are
    still plain quantiles, with the bias correction dropped too.

    When my draws are kept, I recompute my bounds from them and refuse to be
    built if they differ. When they were not kept, I was checked when I was
    first built, before they were dropped.

    Attributes
    ----------
    (lower, upper, alpha, method_intended and method_used are inherited from
    CIResult; n_resamples, n_failed and draws from BootstrapCIResult.)
    """

    # Inherited from and populated at CIResult parent class
    lower: float
    upper: float
    alpha: float
    method_intended: str
    method_used: str | None

    # Inherited from and populated at BootstrapCIResult parent class
    n_resamples: int
    n_failed: int
    draws: np.ndarray | NotComputed = field(compare=False, repr=False)

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction, before I freeze."""
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    def _check_method(self) -> None:
        """Percentile ran, on some usable draw, and my bounds are its quantiles."""
        if self.method_used != "percentile":
            raise ValueError(
                f"method_used must be 'percentile', got {self.method_used!r}"
            )
        if self.n_usable == 0:
            raise ValueError(
                "a percentile interval needs some usable draw; with none, "
                "the interval is an EmptyBootstrapCIResult"
            )
        if isinstance(self.draws, NotComputed):
            return
        lower, upper = self.compute_bounds(self.draws, self.alpha)
        if not (
            np.isclose(self.lower, lower, rtol=0, atol=_TOL)
            and np.isclose(self.upper, upper, rtol=0, atol=_TOL)
        ):
            raise ValueError(
                f"bounds ({self.lower}, {self.upper}) are not the percentile "
                f"bounds of my draws, ({lower}, {upper})"
            )

    # Public API ------------------------------------------------------------- #

    @staticmethod
    def compute_bounds(draws: np.ndarray, alpha: float) -> tuple[float, float]:
        """
        Compute the alpha/2 and 1 - alpha/2 quantiles of the usable draws.

        Parameters
        ----------
        draws : np.ndarray
            The draws; non-finite ones are left out.
        alpha : float
            Significance level, e.g. 0.05 for a 95% interval.

        Returns
        -------
        tuple[float, float]
            (lower, upper).
        """
        usable = BootstrapCIResult.usable_draws(draws)
        if usable.size == 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: no usable draw to take quantiles of"
            )
        return (
            float(np.percentile(usable, (alpha / 2) * 100)),
            float(np.percentile(usable, (1 - alpha / 2) * 100)),
        )
