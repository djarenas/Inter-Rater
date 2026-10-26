"""
Bootstrap draws and the interval built from them.

A result that holds raw draws also holds the interval made from those draws.
They are two views of one bootstrap, so they must agree, and anything that
summarises the draws must count them the way the interval did: a draw is
usable when it is finite, as in ``analyzers.ci_methods``.

Callers wrap the errors raised here with their own prefix.
"""

from __future__ import annotations

import numpy as np

from .ci_results.bootstrap_ci_result import BootstrapCIResult


def check_draws_match_interval(
    draws: np.ndarray | None,
    interval: object,
    *,
    draws_name: str,
    interval_name: str,
) -> None:
    """
    Refuse draws and an interval that do not describe the same bootstrap.

    Draws are present exactly when the interval is a ``BootstrapCIResult``;
    there are as many draws as it resampled, and as many unusable draws as
    it counted failed.

    Raises
    ------
    ValueError
        If any of the three does not hold.
    """
    has_draws = draws is not None
    is_bootstrap = isinstance(interval, BootstrapCIResult)
    if has_draws and not is_bootstrap:
        raise ValueError(
            f"{draws_name} is present but {interval_name} is "
            f"{type(interval).__name__}, not a BootstrapCIResult; draws "
            f"come only with the bootstrap interval made from them"
        )
    if is_bootstrap and not has_draws:
        raise ValueError(
            f"{interval_name} is a BootstrapCIResult but {draws_name} is "
            f"None; a bootstrap interval comes with the draws it was made from"
        )
    if not has_draws:
        return
    if len(draws) != interval.n_resamples:
        raise ValueError(
            f"{draws_name} has {len(draws)} draws but {interval_name} "
            f"resampled {interval.n_resamples}"
        )
    n_unusable = int(np.count_nonzero(~np.isfinite(draws)))
    if n_unusable != interval.n_failed:
        raise ValueError(
            f"{draws_name} has {n_unusable} unusable (non-finite) draws but "
            f"{interval_name} counted {interval.n_failed} failed"
        )


def draw_moments(draws: np.ndarray) -> tuple[float | None, float | None]:
    """
    Mean and bootstrap standard error of the usable draws.

    The standard error divides by n - 1 (Efron & Tibshirani 1993, eq. 6.6).
    The mean is None with no usable draw; the standard error is None with
    fewer than two.

    Returns
    -------
    tuple[float | None, float | None]
        (mean, standard error)
    """
    usable = draws[np.isfinite(draws)]
    mean = float(np.mean(usable)) if usable.size >= 1 else None
    se = float(np.std(usable, ddof=1)) if usable.size >= 2 else None
    return mean, se
