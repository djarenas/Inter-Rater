"""Confidence interval computation methods for bootstrap distributions."""

from __future__ import annotations

import numpy as np

_ERROR_PREFIX = "[ci_methods] Error"


def percentile_ci(bootstrap_distribution: np.ndarray, alpha: float) -> tuple[float, float]:
    """
    Compute percentile confidence interval from bootstrap distribution.

    The percentile method is the simplest bootstrap CI: take the alpha/2 and
    1-alpha/2 quantiles of the bootstrap distribution directly.

    Parameters
    ----------
    bootstrap_distribution : np.ndarray
        Array of shape (n_resamples,) with coefficient values from each resample.
    alpha : float
        Significance level (e.g., 0.05 for 95% CI).

    Returns
    -------
    tuple[float, float]
        (lower, upper) bounds of the confidence interval.

    Raises
    ------
    ValueError
        If bootstrap_distribution is empty or alpha out of range.
    """
    if len(bootstrap_distribution) == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: bootstrap_distribution cannot be empty"
        )

    if not (0 < alpha < 1):
        raise ValueError(
            f"{_ERROR_PREFIX}: alpha must be in (0, 1), got {alpha}"
        )

    # Quantiles for CI bounds
    lower_percentile = (alpha / 2) * 100
    upper_percentile = (1 - alpha / 2) * 100

    lower = float(np.percentile(bootstrap_distribution, lower_percentile))
    upper = float(np.percentile(bootstrap_distribution, upper_percentile))

    return (lower, upper)


def bca_ci(
    bootstrap_distribution: np.ndarray, point: float, alpha: float
) -> tuple[float, float]:
    """
    Compute bias-corrected and accelerated (BCa) confidence interval.

    BCa adjusts the percentile method for bias and skewness in the bootstrap
    distribution. It is more accurate than the percentile method, especially
    for skewed distributions.

    Algorithm (Efron 1987):
    1. Compute bias correction z0 = Phi^-1(F(point))
       where F is the ECDF of bootstrap_distribution and Phi^-1 is inverse normal CDF
    2. Compute acceleration a using jackknife influence values
    3. Adjust lower and upper percentiles using z0 and a
    4. Return quantiles at adjusted percentiles

    Parameters
    ----------
    bootstrap_distribution : np.ndarray
        Array of shape (n_resamples,) with coefficient values from each resample.
    point : float
        The point estimate (e.g., point kappa from the original data).
    alpha : float
        Significance level (e.g., 0.05 for 95% CI).

    Returns
    -------
    tuple[float, float]
        (lower, upper) bounds of the BCa confidence interval.

    Raises
    ------
    ValueError
        If bootstrap_distribution is empty, alpha out of range, or acceleration
        calculation fails.

    References
    ----------
    Efron, B. (1987). Better bootstrap confidence intervals.
    Journal of the American Statistical Association, 82(397), 171-185.
    """
    from scipy.stats import norm

    if len(bootstrap_distribution) == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: bootstrap_distribution cannot be empty"
        )

    if not (0 < alpha < 1):
        raise ValueError(
            f"{_ERROR_PREFIX}: alpha must be in (0, 1), got {alpha}"
        )

    n = len(bootstrap_distribution)

    # Step 1: Bias correction z0
    # z0 = Phi^-1(F_hat(theta_hat))
    # F_hat is ECDF: proportion of bootstrap values <= point
    ecdf_point = np.mean(bootstrap_distribution <= point)
    # Clamp to avoid infinities at exactly 0 or 1
    ecdf_point = np.clip(ecdf_point, 1e-10, 1 - 1e-10)
    z0 = norm.ppf(ecdf_point)

    # Step 2: Acceleration a
    # Computed from jackknife influence values
    # a = sum(L_i^3) / (6 * (sum(L_i^2))^1.5)
    # where L_i is the ith jackknife influence value
    
    # Jackknife: leave-one-out replicates
    jackknife_values = []
    for i in range(n):
        # Leave out the ith value
        mask = np.arange(n) != i
        jackknife_values.append(np.mean(bootstrap_distribution[mask]))

    jackknife_values = np.array(jackknife_values)
    
    # Compute influence values L_i
    # L_i = (mean(boot) - boot_i) for each leave-one-out
    jackknife_mean = np.mean(jackknife_values)
    influence_values = jackknife_mean - jackknife_values

    # Acceleration
    numerator = np.sum(influence_values**3)
    denominator = 6 * (np.sum(influence_values**2) ** 1.5)

    if denominator == 0:
        # No acceleration (symmetric distribution)
        acceleration = 0.0
    else:
        acceleration = numerator / denominator

    # Guard: NaN in bootstrap_distribution propagates through the
    # jackknife into acceleration. If acceleration is not finite,
    # the BCa correction is unreliable -- fall back to percentile.
    if not np.isfinite(acceleration):
        return percentile_ci(bootstrap_distribution, alpha)

    # Step 3: Adjusted percentiles
    # p_lower = Phi(z0 + (z0 + z_alpha/2) / (1 - a*(z0 + z_alpha/2)))
    # p_upper = Phi(z0 + (z0 + z_{1-alpha/2}) / (1 - a*(z0 + z_{1-alpha/2})))
    
    z_lower = norm.ppf(alpha / 2)
    z_upper = norm.ppf(1 - alpha / 2)

    # Adjusted z values
    z_lower_adj = z0 + (z0 + z_lower) / (1 - acceleration * (z0 + z_lower))
    z_upper_adj = z0 + (z0 + z_upper) / (1 - acceleration * (z0 + z_upper))

    # Guard: when the denominator 1 - a*(z0 + z) crosses zero, the
    # adjusted z-value becomes inf or NaN. norm.cdf(inf) = 0 or 1
    # (degenerate CI), norm.cdf(NaN) = NaN (crashes np.quantile).
    # Either way the BCa correction is unstable -- fall back to
    # the simpler percentile method.
    if not (np.isfinite(z_lower_adj) and np.isfinite(z_upper_adj)):
        return percentile_ci(bootstrap_distribution, alpha)

    # Convert back to percentiles
    p_lower = norm.cdf(z_lower_adj)
    p_upper = norm.cdf(z_upper_adj)

    # Clamp to [0, 1]
    p_lower = np.clip(p_lower, 0, 1)
    p_upper = np.clip(p_upper, 0, 1)

    # Belt-and-suspenders: np.clip does not catch NaN
    # (clip(NaN, 0, 1) returns NaN). If NaN survived to here
    # despite the checks above, fall back to percentile.
    if np.isnan(p_lower) or np.isnan(p_upper):
        return percentile_ci(bootstrap_distribution, alpha)

    # Step 4: Quantiles at adjusted percentiles
    lower = float(np.quantile(bootstrap_distribution, p_lower))
    upper = float(np.quantile(bootstrap_distribution, p_upper))

    return (lower, upper)


def scipy_bca_ci(
    bootstrap_distribution: np.ndarray, point: float, alpha: float
) -> tuple[float, float]:
    """
    Compute BCa CI using scipy.stats (for validation/testing only).

    This is a reference implementation for testing and comparing against the
    custom bca_ci() implementation. It uses scipy.stats.bootstrap for BCa
    computation.

    Parameters
    ----------
    bootstrap_distribution : np.ndarray
        Array of shape (n_resamples,) with coefficient values from each resample.
    point : float
        The point estimate.
    alpha : float
        Significance level.

    Returns
    -------
    tuple[float, float]
        (lower, upper) bounds of the BCa CI.

    Raises
    ------
    ValueError
        If scipy.stats does not have the required bootstrap function.
    """
    from scipy.stats import bootstrap

    if len(bootstrap_distribution) == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: bootstrap_distribution cannot be empty"
        )

    if not (0 < alpha < 1):
        raise ValueError(
            f"{_ERROR_PREFIX}: alpha must be in (0, 1), got {alpha}"
        )

    # scipy.stats.bootstrap computes CI from bootstrap replicates
    # We pass our bootstrap distribution directly and use bca method
    try:
        # Create a callable that returns the bootstrap distribution
        def statistic(x):
            return x.mean()

        # Use scipy's BCa
        result = bootstrap(
            (bootstrap_distribution,),
            statistic,
            method="bca",
            n_resamples=len(bootstrap_distribution),
            random_state=None,
            confidence_level=1 - alpha,
        )

        lower, upper = result.confidence_interval

        return (float(lower), float(upper))

    except Exception as e:
        raise ValueError(
            f"{_ERROR_PREFIX}: scipy.stats.bootstrap failed: {e}"
        ) from e
