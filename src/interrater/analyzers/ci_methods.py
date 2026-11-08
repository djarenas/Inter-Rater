"""Confidence interval computation methods for bootstrap distributions.

Both methods return a :class:`BootstrapCIResult`, which records the method
asked for, the method that ran, and how many draws failed.

**Failed draws.** A coefficient of the form ``(p_o - p_e) / (1 - p_e)`` is
undefined when a resample has a single category (p_e = 1). The agreement
functions return NaN there. Both methods take their bounds from the usable
draws and count the rest. Neither thresholds on that count: whether failures
matter depends on whether failing correlates with the coefficient's value.

**Acceleration needs the data.** BCa's ``a`` measures how the coefficient's
standard error changes with its value (agreement near 1 varies less than
agreement near 0.5). A bootstrap distribution describes the spread at one
value only, so ``a`` comes from the leave-one-subject-out estimates instead
(Efron & Tibshirani 1986, eq. 8.5). These functions never see the ratings, so
:func:`bca_ci` takes those estimates as a required argument and raises
without them. (Versions before 2.0.1 took the skewness of the draws, which is
not the acceleration.)

**When BCa steps down:**

- ``a`` undefined (fewer than two jackknife values, or all identical):
  percentile.
- ``1 - a(z0 + z)`` not positive at either end, which would send a bound
  to infinity or to the wrong side: percentile.
- no draw usable: nothing runs, and the method is None.

A fallback to percentile drops ``z0`` as well as ``a``, so ``method_used``
'percentile' means plain quantiles. Each method in Efron & Tibshirani's
Table 6 is correct under weaker assumptions than the one above it, so the
fallback costs accuracy, not validity. ``method_used`` of None is an absence:
there are no bounds.

References
----------
Efron, B., & Tibshirani, R. (1986). "Bootstrap Methods for Standard Errors,
Confidence Intervals, and Other Measures of Statistical Accuracy."
*Statistical Science*, 1(1), 54-77. Table 6 (p. 68) for the nested sequence
of intervals; eq. (8.5) (p. 71) for the acceleration; section 9 (p. 72) for
the resample counts each method needs.

Efron, B. (1987). "Better Bootstrap Confidence Intervals." *Journal of the
American Statistical Association*, 82(397), 171-185.
"""

from __future__ import annotations

import numpy as np

from ..intermediates.results.ci_results.bootstrap_ci_result import (
    BootstrapCIResult,
    MIN_RESAMPLES_BCA,
    MIN_RESAMPLES_PERCENTILE,
)

_ERROR_PREFIX = "[ci_methods] Error"

# The resample minima live on BootstrapCIResult, because the interval is what
# reports against them. Re-exported here so that `from ci_methods import
# MIN_RESAMPLES_BCA` keeps working and so that a reader of this module finds
# them next to the functions they describe.
__all__ = [
    "percentile_ci",
    "bca_ci",
    "MIN_RESAMPLES_BCA",
    "MIN_RESAMPLES_PERCENTILE",
]


# -- Shared machinery ------------------------------------------------------- #


def _checked(bootstrap_distribution: np.ndarray, alpha: float) -> np.ndarray:
    """Validate the two arguments every method here takes, and return floats."""
    if len(bootstrap_distribution) == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: bootstrap_distribution cannot be empty"
        )
    if not (0 < alpha < 1):
        raise ValueError(
            f"{_ERROR_PREFIX}: alpha must be in (0, 1), got {alpha}"
        )
    return np.asarray(bootstrap_distribution, dtype=np.float64)


def _usable(array: np.ndarray) -> tuple[np.ndarray, int]:
    """
    Split an array into its usable values and a count of the rest.

    Usable means finite. The agreement functions signal an undefined
    coefficient with NaN, which is what this is here for, but an infinite
    draw is equally unusable for a quantile and is counted the same way.
    """
    keep = np.isfinite(array)
    return array[keep], int(array.size - int(keep.sum()))


def _acceleration(jackknife_estimates: np.ndarray) -> float | None:
    """
    Efron & Tibshirani (1986) eq. (8.5), in its jackknife form.

        a = sum(d_i^3) / (6 * [sum(d_i^2)]^(3/2)),  d_i = mean(theta_(.)) - theta_(i)

    where ``theta_(i)`` is the coefficient recomputed with subject i left out.
    The published formula is written in empirical influence values ``U_i^0``,
    which differ from these ``d_i`` by a factor of (n-1); the factor cancels
    between numerator and denominator, so the two forms agree exactly.

    Parameters
    ----------
    jackknife_estimates : np.ndarray
        Leave-one-subject-out values of the coefficient, one per subject.
        Non-finite entries are dropped -- a leave-one-out dataset can be
        degenerate for the same reason a resample can.

    Returns
    -------
    float | None
        The acceleration, or None when it cannot be computed: fewer than two
        usable values, or every usable value identical, which leaves no
        skewness to correct.

    Notes
    -----
    Identical values are detected by their range, not by their deviations
    from the mean. The mean of identical floats is not always one of them:
    twenty copies of 0.6 average to a value one ulp away, leaving deviations
    near 1e-17. Because ``a`` is a standardised skewness, that noise divided
    by noise of the same order gives an ordinary-looking acceleration, and
    BCa would apply a correction to a coefficient no subject moves (0.5 and
    0.75 are exact in binary and did not show it; 0.6 did).

    The range of identical values is exactly zero. The tolerance is a few
    ulps of the values' magnitude; real leave-one-out spread is around 1e-3
    (the coefficient over the number of subjects), far above it.
    """
    values, _ = _usable(np.asarray(jackknife_estimates, dtype=np.float64))
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


def _nothing_ran(
    draws: np.ndarray, n_failed: int, alpha: float, intended: str
) -> BootstrapCIResult:
    """The result when no resample was usable, so no method could be applied."""
    return BootstrapCIResult(
        lower=float("nan"),
        upper=float("nan"),
        alpha=alpha,
        n_resamples=draws.size,
        n_failed=n_failed,
        method_intended=intended,
        method_used=None,
    )


# -- The methods ------------------------------------------------------------ #


def percentile_ci(
    bootstrap_distribution: np.ndarray, alpha: float
) -> BootstrapCIResult:
    """
    The alpha/2 and 1 - alpha/2 quantiles of the bootstrap distribution.

    The simpler of the two methods in Efron & Tibshirani's Table 6, and the
    floor of this module: BCa can decline and hand its work here, while this
    has nothing simpler to hand off to in turn. So this is where failed
    resamples are dealt with rather than passed along.

    Parameters
    ----------
    bootstrap_distribution : np.ndarray
        Shape (n_resamples,). Coefficient values from each resample. NaN
        entries mark resamples whose coefficient was undefined.
    alpha : float
        Significance level, e.g. 0.05 for a 95% interval.

    Returns
    -------
    BootstrapCIResult
        With ``method_intended`` of 'percentile', and ``method_used`` of
        'percentile', or None when no draw was usable.

    Raises
    ------
    ValueError
        If bootstrap_distribution is empty or alpha is outside (0, 1). An
        all-NaN distribution is not an error: it is a result.
    """
    draws = _checked(bootstrap_distribution, alpha)
    usable, n_failed = _usable(draws)

    if usable.size == 0:
        return _nothing_ran(draws, n_failed, alpha, "percentile")

    return BootstrapCIResult(
        lower=float(np.percentile(usable, (alpha / 2) * 100)),
        upper=float(np.percentile(usable, (1 - alpha / 2) * 100)),
        alpha=alpha,
        n_resamples=draws.size,
        n_failed=n_failed,
        method_intended="percentile",
        method_used="percentile",
    )


def bca_ci(
    bootstrap_distribution: np.ndarray,
    point: float,
    alpha: float,
    jackknife_estimates: np.ndarray,
) -> BootstrapCIResult:
    """
    Bias-corrected and accelerated interval (Efron & Tibshirani 1986, eq. 7.15).

    Two corrections sit on top of the percentile method. The bias correction
    ``z0`` asks whether the bootstrap distribution sits centred on the point
    estimate, and shifts the interval when it does not. The acceleration ``a``
    asks whether the coefficient's standard error changes as the coefficient
    itself changes, and lets the interval lean when it does; it is read from
    how much each individual subject moves the estimate.

    When either correction cannot be computed I fall back to the percentile
    interval, dropping ``z0`` as well as ``a``, so ``method_used`` of
    'percentile' means plain quantiles.

    Parameters
    ----------
    bootstrap_distribution : np.ndarray
        Shape (n_resamples,). NaN marks a resample whose coefficient was
        undefined; those draws are excluded and counted.
    point : float
        The coefficient on the full data, used for the bias correction.
    alpha : float
        Significance level, e.g. 0.05 for a 95% interval.
    jackknife_estimates : np.ndarray
        One leave-one-subject-out value of the coefficient per subject.
        ``resampling_functions.subject_jackknife`` yields the datasets these
        are computed from. Required: asking for BCa without the means to
        compute it is a mistake in the calling analyzer, not a property of
        the data, so it raises rather than quietly producing something else.

    Returns
    -------
    BootstrapCIResult
        With ``method_intended`` of 'bca', and ``method_used`` of 'bca',
        'percentile' when a correction was unavailable, or None when no draw
        was usable.

    Raises
    ------
    ValueError
        If bootstrap_distribution is empty, alpha is outside (0, 1), or
        jackknife_estimates is None or empty.

    Notes
    -----
    Efron & Tibshirani (1986) section 9 puts the rough minimum for BCa at
    B = 1000 resamples, against 250 for a useful percentile interval. This
    function does not enforce that; ``n_usable`` is reported so the caller can.
    """
    from scipy.stats import norm

    draws = _checked(bootstrap_distribution, alpha)

    if jackknife_estimates is None or len(jackknife_estimates) == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: bca_ci requires jackknife_estimates -- "
            f"leave-one-subject-out values of the coefficient. Acceleration "
            f"is a property of the data and cannot be recovered from the "
            f"bootstrap distribution."
        )

    usable, n_failed = _usable(draws)

    # No usable draw means no z0, so BCa cannot run whatever the jackknife
    # says, and neither can the percentile method it would fall back to.
    if usable.size == 0:
        return _nothing_ran(draws, n_failed, alpha, "bca")

    def _fall_back_to_percentile() -> BootstrapCIResult:
        """A true percentile interval, bias correction discarded with the rest."""
        return BootstrapCIResult(
            lower=float(np.percentile(usable, (alpha / 2) * 100)),
            upper=float(np.percentile(usable, (1 - alpha / 2) * 100)),
            alpha=alpha,
            n_resamples=draws.size,
            n_failed=n_failed,
            method_intended="bca",
            method_used="percentile",
        )

    acceleration = _acceleration(jackknife_estimates)
    if acceleration is None:
        return _fall_back_to_percentile()

    # Bias correction, over the usable draws only. Taking it over the whole
    # array would count every failed draw as lying above the point estimate,
    # since `np.nan <= point` is False -- a bias rather than an error, and so
    # one that nothing downstream would catch.
    #
    # Draws equal to the point estimate count as half below and half above.
    # Counting them all as below (`<=`) shifts the interval up, counting them
    # all as above (Efron & Tibshirani's strict `<`) shifts it down, and on
    # small discrete data ties are common: Brennan-Prediger, which depends
    # only on raw agreement, ties the point in up to 12% of draws in the
    # example data. This is the 'mean' rule scipy.stats.bootstrap uses.
    n_below = np.sum(usable < point)
    n_at_or_below = np.sum(usable <= point)
    ecdf_point = float((n_below + n_at_or_below) / (2 * usable.size))
    z0 = float(norm.ppf(np.clip(ecdf_point, 1e-10, 1 - 1e-10)))

    z_lower = float(norm.ppf(alpha / 2))
    z_upper = float(norm.ppf(1 - alpha / 2))

    # The correction needs 1 - a*(z0 + z) > 0 at both ends. At zero the
    # adjusted quantile runs off to infinity; below zero it lands on the wrong
    # side and the bounds come out inverted (possible when the point lies
    # outside nearly all the draws, so z0 sits at its clip, and |a| is large).
    # The data are fine, the correction is not. Table 6 makes stepping down
    # the right answer: each method there is correct under weaker assumptions
    # than the one below it, so the fallback costs accuracy, not validity.
    denom_lower = 1 - acceleration * (z0 + z_lower)
    denom_upper = 1 - acceleration * (z0 + z_upper)
    if not (denom_lower > 0 and denom_upper > 0):
        return _fall_back_to_percentile()

    z_lower_adj = z0 + (z0 + z_lower) / denom_lower
    z_upper_adj = z0 + (z0 + z_upper) / denom_upper

    p_lower = float(np.clip(norm.cdf(z_lower_adj), 0.0, 1.0))
    p_upper = float(np.clip(norm.cdf(z_upper_adj), 0.0, 1.0))

    return BootstrapCIResult(
        lower=float(np.quantile(usable, p_lower)),
        upper=float(np.quantile(usable, p_upper)),
        alpha=alpha,
        n_resamples=draws.size,
        n_failed=n_failed,
        method_intended="bca",
        method_used="bca",
        acceleration=acceleration,
        bias_correction=z0,
    )
