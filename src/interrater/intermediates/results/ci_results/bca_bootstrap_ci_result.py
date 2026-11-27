"""Result for a bias-corrected and accelerated (BCa) bootstrap interval."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ...not_computed import NotComputed
from .bootstrap_ci_result import BootstrapCIResult

_ERROR_PREFIX = "[BCABootstrapCIResult] Error"

# Values recomputed from my inputs must match mine to this
_TOL = 1e-12


@dataclass(frozen=True)
class BCABootstrapCIResult(BootstrapCIResult):
    """
    I am a BCa bootstrap interval (Efron & Tibshirani 1986, eq. 7.15): the
    percentile interval with two corrections, the bias correction z0 and the
    acceleration a.

    I hold what those corrections are computed from: the point estimate z0
    is measured against, and the leave-one-subject-out estimates a is read
    from. With those and my draws I recompute a, z0 and my bounds, and refuse
    to be built if any of them differs. When my draws or jackknife estimates
    were not kept, I was checked when I was first built, before they were
    dropped.

    I always ran as BCa: when BCa cannot apply its corrections the interval
    is a PercentileBootstrapCIResult, and when no draw was usable it is an
    EmptyBootstrapCIResult.

    Attributes
    ----------
    point : float
        The coefficient on the full data, which z0 is measured against.
    jackknife_estimates : np.ndarray | NotComputed
        One leave-one-subject-out value of the coefficient per subject.
        NotComputed when they were not kept, saying why.
    acceleration : float
        BCa's acceleration constant a.
    bias_correction : float
        BCa's bias correction z0.
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

    # Populated at this level
    point: float
    # Left out of equality and repr, as the draws are
    jackknife_estimates: np.ndarray | NotComputed = field(compare=False, repr=False)
    acceleration: float
    bias_correction: float

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction, before I freeze."""
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    def _check_method(self) -> None:
        """BCa ran, and a, z0 and my bounds follow from my inputs."""
        if self.method_intended != "bca" or self.method_used != "bca":
            raise ValueError(
                f"method_intended and method_used must both be 'bca', got "
                f"{self.method_intended!r} and {self.method_used!r}"
            )
        if self.n_usable == 0:
            raise ValueError(
                "a BCa interval needs some usable draw; with none, the "
                "interval is an EmptyBootstrapCIResult"
            )

        # point, acceleration and bias_correction must be numbers; the two
        # corrections must also be finite, or BCa could not have run
        for name, value in (
            ("point", self.point),
            ("acceleration", self.acceleration),
            ("bias_correction", self.bias_correction),
        ):
            if isinstance(value, bool) or not isinstance(
                value, (int, float, np.floating)
            ):
                raise TypeError(
                    f"{name} must be float, got {type(value).__name__}"
                )
        for name, value in (
            ("acceleration", self.acceleration),
            ("bias_correction", self.bias_correction),
        ):
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite, got {value}")
        object.__setattr__(self, "point", float(self.point))
        object.__setattr__(self, "acceleration", float(self.acceleration))
        object.__setattr__(self, "bias_correction", float(self.bias_correction))

        # The acceleration is the one my jackknife estimates give
        if not isinstance(self.jackknife_estimates, NotComputed):
            if not isinstance(self.jackknife_estimates, np.ndarray):
                raise TypeError(
                    f"jackknife_estimates must be np.ndarray or NotComputed, "
                    f"got {type(self.jackknife_estimates).__name__}"
                )
            acceleration = self.compute_acceleration(self.jackknife_estimates)
            if acceleration is None or not np.isclose(
                self.acceleration, acceleration, rtol=0, atol=_TOL
            ):
                raise ValueError(
                    f"acceleration {self.acceleration} is not the one my "
                    f"jackknife estimates give, {acceleration}"
                )

        # z0 and the bounds are the ones my draws give
        if isinstance(self.draws, NotComputed):
            return
        bias_correction = self.compute_bias_correction(self.draws, self.point)
        if not np.isclose(self.bias_correction, bias_correction, rtol=0, atol=_TOL):
            raise ValueError(
                f"bias_correction {self.bias_correction} is not the one my "
                f"draws and point give, {bias_correction}"
            )
        bounds = self.compute_bounds(
            self.draws, self.alpha, self.bias_correction, self.acceleration
        )
        if bounds is None or not (
            np.isclose(self.lower, bounds[0], rtol=0, atol=_TOL)
            and np.isclose(self.upper, bounds[1], rtol=0, atol=_TOL)
        ):
            raise ValueError(
                f"bounds ({self.lower}, {self.upper}) are not the BCa bounds "
                f"of my draws, {bounds}"
            )

    # Public API ------------------------------------------------------------- #

    @staticmethod
    def compute_acceleration(jackknife_estimates: np.ndarray) -> float | None:
        """
        Compute the acceleration: Efron & Tibshirani (1986) eq. (8.5), in its
        jackknife form.

            a = sum(d_i^3) / (6 * [sum(d_i^2)]^(3/2)),  d_i = mean(theta_(.)) - theta_(i)

        where ``theta_(i)`` is the coefficient recomputed with subject i left
        out. The published formula is written in empirical influence values
        ``U_i^0``, which differ from these ``d_i`` by a factor of (n-1); the
        factor cancels between numerator and denominator, so the two forms
        agree exactly.

        Parameters
        ----------
        jackknife_estimates : np.ndarray
            Leave-one-subject-out values of the coefficient, one per subject.
            Non-finite entries are dropped -- a leave-one-out dataset can be
            degenerate for the same reason a resample can.

        Returns
        -------
        float | None
            The acceleration, or None when it cannot be computed: fewer than
            two usable values, or every usable value identical, which leaves
            no skewness to correct.

        Notes
        -----
        Identical values are caught by their range, not by their deviations
        from the mean. The mean of identical floats can be one ulp off
        (twenty copies of 0.6 average to 0.5999999999999999), leaving
        deviations near 1e-16. ``a`` divides one such tiny number by another,
        so it would come out ordinary-looking and BCa would correct for skew
        that is not there. The range of identical values is exactly zero, and
        the cut-off is a few ulps of the values' size.
        """
        values = BootstrapCIResult.usable_draws(jackknife_estimates)
        if values.size < 2:
            return None

        # See Notes: the gate is on the spread of the values themselves.
        scale = float(np.max(np.abs(values)))
        noise_floor = float(np.finfo(np.float64).eps) * max(scale, 1.0) * 8.0
        if float(np.ptp(values)) <= noise_floor:
            return None

        deviations = values.mean() - values
        denominator = 6.0 * float(np.sum(deviations ** 2) ** 1.5)
        if denominator == 0.0 or not np.isfinite(denominator):
            return None

        acceleration = float(np.sum(deviations ** 3) / denominator)
        return acceleration if np.isfinite(acceleration) else None

    @staticmethod
    def compute_bias_correction(draws: np.ndarray, point: float) -> float:
        """
        Compute the bias correction z0 over the usable draws.

        Taking it over every draw would count each failed draw as lying above
        the point estimate, since ``np.nan <= point`` is False -- a bias
        rather than an error, and so one nothing downstream would catch.

        Draws equal to the point estimate count as half below and half above.
        Counting them all as below (``<=``) shifts the interval up, counting
        them all as above (Efron & Tibshirani's strict ``<``) shifts it down,
        and on small discrete data ties are common, e.g. Brennan-Prediger,
        which depends only on raw agreement. This is the 'mean' rule
        scipy.stats.bootstrap uses.

        Parameters
        ----------
        draws : np.ndarray
            The draws; at least one must be usable.
        point : float
            The coefficient on the full data.

        Returns
        -------
        float
            z0, from a share clipped to [1e-10, 1 - 1e-10] so it is finite.
        """
        from scipy.stats import norm

        usable = BootstrapCIResult.usable_draws(draws)
        if usable.size == 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: no usable draw to measure the bias against"
            )
        n_below = np.sum(usable < point)
        n_at_or_below = np.sum(usable <= point)
        share_below = float((n_below + n_at_or_below) / (2 * usable.size))
        return float(norm.ppf(np.clip(share_below, 1e-10, 1 - 1e-10)))

    @staticmethod
    def compute_bounds(
        draws: np.ndarray,
        alpha: float,
        bias_correction: float,
        acceleration: float,
    ) -> tuple[float, float] | None:
        """
        Compute the BCa bounds: quantiles of the usable draws at levels the
        two corrections adjust (Efron & Tibshirani 1986, eq. 7.15).

        Parameters
        ----------
        draws : np.ndarray
            The draws; at least one must be usable.
        alpha : float
            Significance level, e.g. 0.05 for a 95% interval.
        bias_correction : float
            z0.
        acceleration : float
            a.

        Returns
        -------
        tuple[float, float] | None
            (lower, upper), or None when 1 - a(z0 + z) is not positive at
            either end. At zero the adjusted level runs off to infinity;
            below zero the bounds come out on the wrong side. The data are
            fine, the correction is not, and BCa does not apply.
        """
        from scipy.stats import norm

        usable = BootstrapCIResult.usable_draws(draws)
        if usable.size == 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: no usable draw to take quantiles of"
            )
        z0, a = bias_correction, acceleration
        z_lower = float(norm.ppf(alpha / 2))
        z_upper = float(norm.ppf(1 - alpha / 2))
        denom_lower = 1 - a * (z0 + z_lower)
        denom_upper = 1 - a * (z0 + z_upper)
        if not (denom_lower > 0 and denom_upper > 0):
            return None
        p_lower = float(np.clip(norm.cdf(z0 + (z0 + z_lower) / denom_lower), 0.0, 1.0))
        p_upper = float(np.clip(norm.cdf(z0 + (z0 + z_upper) / denom_upper), 0.0, 1.0))
        return (
            float(np.quantile(usable, p_lower)),
            float(np.quantile(usable, p_upper)),
        )

    def summary(self) -> dict:
        """
        Return my parent's summary plus my two corrections, in plain types.

        Returns
        -------
        dict
            My parent's keys, plus 'acceleration' and 'bias_correction'
            (float).
        """
        result = super().summary()
        result["acceleration"] = self.acceleration
        result["bias_correction"] = self.bias_correction
        return result
