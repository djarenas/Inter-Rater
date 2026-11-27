"""Confidence interval computation methods for bootstrap distributions.

Two ci methods: percentile_ci and bca_ci

Both of the ci methods return a :class:`BootstrapCIResult`. This result object records the
method asked for, the method that ran, and how many draws failed.

**Failed draws.** Coefficients of the form ``(p_o - p_e) / (1 - p_e)`` are
undefined when a resample has a single category (p_e = 1). The agreement
functions return NaN there. Both methods take their bounds from the usable
draws and count the rest. Neither thresholds on that count: whether failures
matter depends on whether failing correlates with the coefficient's value.

**Acceleration needs the data.** BCa's ``a`` measures how the coefficient's
standard error changes with its value (for example, agreement near 1 varies less than
agreement near 0.5). ``a`` comes from the leave-one-subject-out estimates
(Efron & Tibshirani 1986, eq. 8.5). These functions never see the ratings, so
:func:`bca_ci` takes those estimates as a required argument and raises
without them.

**When BCa steps down:**
If the method requested is bca_ci, this package handles three scenarios at
which the method could fail:

- ``a`` undefined (fewer than two jackknife values, or all identical) ->
  percentile method is used instead.
- ``1 - a(z0 + z)`` not positive at either end, which would send a bound
  to infinity or to the wrong side -> percentile method is used instead.
- no draw usable: nothing runs, and the method is None.

A fallback to percentile drops ``z0`` as well as ``a``, so ``method_used``
'percentile' means plain quantiles. With z0 = a = 0, BCa reduces to the
percentile interval, so falling back drops the two corrections and keeps the
method they correct.

``method_used`` is None only when no draw was usable: there are no bounds.

References
----------
Efron, B., & Tibshirani, R. (1986). "Bootstrap Methods for Standard Errors,
Confidence Intervals, and Other Measures of Statistical Accuracy."
*Statistical Science*, 1(1), 54-77. Table 6 (p. 68) for the four methods,
each correct under more general assumptions than the one before; eq. (8.5)
(p. 71) for the acceleration; section 9 (p. 72) for the resample counts each
method needs.

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
from ..intermediates.results.ci_results.bca_bootstrap_ci_result import (
    BCABootstrapCIResult,
)
from ..intermediates.results.ci_results.empty_bootstrap_ci_result import (
    EmptyBootstrapCIResult,
)
from ..intermediates.results.ci_results.percentile_bootstrap_ci_result import (
    PercentileBootstrapCIResult,
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


def _nothing_ran(draws: np.ndarray, alpha: float, intended: str) -> EmptyBootstrapCIResult:
    """The result when no resample was usable, so no method could be applied."""
    return EmptyBootstrapCIResult(
        lower=float("nan"),
        upper=float("nan"),
        alpha=alpha,
        method_intended=intended,
        method_used=None,
        n_resamples=draws.size,
        n_failed=draws.size,
        draws=draws,
    )


# -- The methods ------------------------------------------------------------ #


def percentile_ci(
    bootstrap_distribution: np.ndarray, alpha: float
) -> BootstrapCIResult:
    """
    The alpha/2 and 1 - alpha/2 quantiles of the bootstrap distribution.

    The percentile method, line 2 of Efron & Tibshirani's Table 6 (1986,
    p. 68), and the floor of this module: BCa can decline and hand its work
    here, while this has nothing simpler to hand off to in turn.

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
    usable = BootstrapCIResult.usable_draws(draws)

    if usable.size == 0:
        return _nothing_ran(draws, alpha, "percentile")

    lower, upper = PercentileBootstrapCIResult.compute_bounds(draws, alpha)
    return PercentileBootstrapCIResult(
        lower=lower,
        upper=upper,
        alpha=alpha,
        method_intended="percentile",
        method_used="percentile",
        n_resamples=draws.size,
        n_failed=draws.size - usable.size,
        draws=draws,
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

    When either correction cannot be computed, falls back to the percentile
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
    draws = _checked(bootstrap_distribution, alpha)

    if jackknife_estimates is None or len(jackknife_estimates) == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: bca_ci requires jackknife_estimates -- "
            f"leave-one-subject-out values of the coefficient. Acceleration "
            f"is a property of the data and cannot be recovered from the "
            f"bootstrap distribution."
        )

    usable = BootstrapCIResult.usable_draws(draws)
    n_failed = draws.size - usable.size

    # No usable draw means no z0, so BCa cannot run whatever the jackknife
    # says, and neither can the percentile method it would fall back to.
    if usable.size == 0:
        return _nothing_ran(draws, alpha, "bca")

    def _fall_back_to_percentile() -> PercentileBootstrapCIResult:
        """A true percentile interval, bias correction discarded with the rest."""
        lower, upper = PercentileBootstrapCIResult.compute_bounds(draws, alpha)
        return PercentileBootstrapCIResult(
            lower=lower,
            upper=upper,
            alpha=alpha,
            method_intended="bca",
            method_used="percentile",
            n_resamples=draws.size,
            n_failed=n_failed,
            draws=draws,
        )

    jackknife = np.asarray(jackknife_estimates, dtype=np.float64)
    acceleration = BCABootstrapCIResult.compute_acceleration(jackknife)
    if acceleration is None:
        return _fall_back_to_percentile()

    z0 = BCABootstrapCIResult.compute_bias_correction(draws, point)

    # The correction needs 1 - a*(z0 + z) > 0 at both ends; when it does not
    # hold, compute_bounds returns None. With z0 = a = 0, BCa reduces to the
    # percentile interval, so falling back drops the two corrections and
    # keeps the method they correct.
    bounds = BCABootstrapCIResult.compute_bounds(draws, alpha, z0, acceleration)
    if bounds is None:
        return _fall_back_to_percentile()

    return BCABootstrapCIResult(
        lower=bounds[0],
        upper=bounds[1],
        alpha=alpha,
        method_intended="bca",
        method_used="bca",
        n_resamples=draws.size,
        n_failed=n_failed,
        draws=draws,
        point=point,
        jackknife_estimates=jackknife,
        acceleration=acceleration,
        bias_correction=z0,
    )
