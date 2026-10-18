"""Confidence interval computation methods for bootstrap distributions.

Three facts shape everything here.

**A resample can fail.** An agreement coefficient of the form
``(p_o - p_e) / (1 - p_e)`` is undefined when chance agreement reaches 1,
which happens whenever a resample contains only one category. The agreement
functions return NaN there, so a bootstrap distribution may arrive with gaps
in it. Both methods below compute their bounds from the usable draws and
report how many failed. Neither thresholds on that count: how much the
failures matter depends on whether failing correlates with the coefficient's
value, which no count can tell you. Interpreting it belongs upstairs.

**Acceleration needs the data, not just the resamples.** The percentile
method assumes a coefficient's standard error is the same whatever the
coefficient's value happens to be. Often it is not: agreement near 1.0 varies
less between samples than agreement near 0.5, because there is less room to
move. BCa's ``a`` corrects for that, and it measures how fast the standard
error changes as the coefficient changes -- which cannot be read off a
bootstrap distribution, since that describes the spread at the single value
the data actually produced. Measuring it means disturbing the data and
watching the coefficient move: drop one subject, recompute, repeat.

Nothing here can do that. These functions receive an array of numbers and
never see the ratings or the estimator. So :func:`bca_ci` takes the
leave-one-out estimates as a **required** argument: asking for BCa without
the means to compute it is a mistake in the calling analyzer rather than a
property of the data, and it raises.

**BCa can decline.** When either correction cannot be computed, it falls all
the way back to the percentile interval -- ``z0`` discarded along with ``a``,
so a reported method of 'percentile' is literally what was computed and not a
relabelled half-corrected interval. When no draw is usable at all, nothing
runs and the method is ``None``.

Every one of those outcomes is recorded on the :class:`BootstrapCIResult`:
which method was asked for, which one ran, and how many draws failed. A
caller never has to infer from the numbers what happened to produce them.

References
----------
Efron, B., & Tibshirani, R. (1986). "Bootstrap Methods for Standard Errors,
Confidence Intervals, and Other Measures of Statistical Accuracy."
*Statistical Science*, 1(1), 54-77. Table 6 (p. 68) defines the percentile
and BCa intervals as a nested sequence, each correct under weaker assumptions
than the one before it; eq. (8.5) (p. 71) gives the acceleration; section 9
(p. 72) gives the resample counts each method needs.

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
    The identical-values case is gated on the spread of the values, not on
    the deviations from their mean, because the mean of identical floats is
    not reliably one of them: summing twenty copies of 0.6 and dividing by
    twenty lands an ulp away, leaving deviations around 1e-17 and a sum of
    squares around 1e-31 rather than the exact zero one would expect. That
    is not a harmless rounding error here. ``a`` is a *standardised*
    skewness, so dividing noise by noise of the same order returns a value
    in the ordinary range for a real acceleration -- a coefficient that no
    single subject moves would get a fabricated correction, and the interval
    would report ``method_used='bca'`` as though it meant something. Whether
    it happened at all depended on the bit pattern of the value: 0.5 and
    0.75 are exact in binary and came out clean, 0.6 did not.

    Identical values have a range of exactly zero whatever the value is,
    since the noise enters through the mean rather than the values, so the
    range is what this tests. The tolerance is a few ulps of the values'
    own magnitude. Real leave-one-out spread is on the order of the
    coefficient divided by the number of subjects -- around 1e-3 against a
    floor near 1e-15 -- so nothing legitimate comes close to the gate.
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

    When either correction cannot be computed I fall all the way back to the
    percentile interval -- ``z0`` discarded, not merely ``a`` -- so that
    ``method_used`` of 'percentile' is literally true rather than a
    convenient relabelling of a half-corrected interval.

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
    ecdf_point = float(np.mean(usable <= point))
    z0 = float(norm.ppf(np.clip(ecdf_point, 1e-10, 1 - 1e-10)))

    z_lower = float(norm.ppf(alpha / 2))
    z_upper = float(norm.ppf(1 - alpha / 2))

    z_lower_adj = z0 + (z0 + z_lower) / (1 - acceleration * (z0 + z_lower))
    z_upper_adj = z0 + (z0 + z_upper) / (1 - acceleration * (z0 + z_upper))

    # When 1 - a*(z0 + z) crosses zero the adjusted quantile runs off to
    # infinity: the data are fine, the correction is not. Table 6 makes
    # stepping down the right answer -- each method there is correct under
    # weaker assumptions than the one below it, so the fallback costs
    # accuracy, not validity.
    if not (np.isfinite(z_lower_adj) and np.isfinite(z_upper_adj)):
        return _fall_back_to_percentile()

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
