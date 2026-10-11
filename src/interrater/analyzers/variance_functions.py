"""Analytical variance formulas for inter-rater agreement coefficients.

This module holds closed-form variance expressions used for parametric
confidence intervals and hypothesis tests. Each function computes the
variance of a specific coefficient under a specific assumption (null or
nonnull). The formulas come from the original statistical literature and
are cited in each docstring.

Bootstrap and permutation approaches do not need these formulas -- they
estimate variance empirically. These exist for:
- Parametric null CIs (ci_method="parametric_null")
- Parametric hypothesis tests (method="parametric_null")
- Reproducing published results that used these formulas

Note on bias correction:
    Fleiss (1979) Eq. 17 includes a bias correction of 1/(N*(n-1)) added
    to kappa before computing the z-statistic. This corrects for the slight
    negative bias in kappa under the null. The correction is O(1/N) and
    negligible for large N, but is included in our hypothesis_test() for
    exactness. The variance formula itself (Eq. 12) does not include the
    bias correction -- it is applied at the z-test level.
"""

from __future__ import annotations

import numpy as np

_ERROR_PREFIX = "[variance_functions] Error"


def fleiss_null_variance(
    pj: np.ndarray,
    n_subjects: int,
    n_raters: int,
) -> float:
    """
    Null variance of Fleiss' kappa under H₀: κ = 0.

    This is Equation 12 from Fleiss, Nee & Landis (1979), which corrected
    the incorrect variance formulas published in Fleiss (1971) and confirmed
    by Cicchetti (1978) via Monte Carlo simulation.

    This variance is valid for:
    - Testing H₀: κ = 0 (hypothesis testing)
    - Large N (asymptotic approximation)
    - Constant n raters per subject (balanced data)
    - Null hypothesis is true

    This variance is NOT valid for:
    - Confidence interval construction (use bootstrap instead)
    - Small N (asymptotic approximation breaks down)
    - Nonnull κ (underestimates variance when κ > 0)

    Parameters
    ----------
    pj : np.ndarray
        Category proportions, shape (n_categories,). Must sum to 1.
    n_subjects : int
        Number of subjects (N in Fleiss' notation).
    n_raters : int
        Number of raters per subject (n in Fleiss' notation).

    Returns
    -------
    float
        Var(κ) under the null hypothesis.

    Raises
    ------
    ValueError
        If pj does not sum to 1 (within tolerance).
        If n_subjects or n_raters are not positive.

    References
    ----------
    Fleiss, J. L., Nee, J. C. M., & Landis, J. R. (1979). "Large sample
        variance of kappa in the case of different sets of raters."
        *Psychological Bulletin*, 86(5), 974-977. Equation 12.

    Fleiss, J. L. (1971). "Measuring nominal scale agreement among many
        raters." *Psychological Bulletin*, 76(5), 378-382.

    Cicchetti, D. V. (1978). Confirmed the corrected formulas via
        Monte Carlo simulation.
    """
    # Validate inputs
    if n_subjects <= 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_subjects must be positive, got {n_subjects}."
        )
    if n_raters <= 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_raters must be positive, got {n_raters}."
        )
    if not np.isclose(pj.sum(), 1.0, atol=1e-6):
        raise ValueError(
            f"{_ERROR_PREFIX}: pj must sum to 1, got {pj.sum():.6f}."
        )

    # qj = 1 - pj
    qj = 1.0 - pj

    # Fleiss (1979) Eq. 12:
    # Var(κ) = [2 / (N * n * (n-1) * (Σ pj*qj)²)]
    #        × [(Σ pj*qj)² - Σ pj*qj*(qj - pj)]

    pq = pj * qj                    # pj * qj for each category
    sum_pq = np.sum(pq)             # Σ pj * qj
    sum_pq_sq = sum_pq ** 2         # (Σ pj * qj)²
    sum_pq_diff = np.sum(pq * (qj - pj))  # Σ pj * qj * (qj - pj)

    N = n_subjects
    n = n_raters

    numerator = 2.0 * (sum_pq_sq - sum_pq_diff)
    denominator = N * n * (n - 1) * sum_pq_sq

    if denominator == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: denominator is zero. This happens when all "
            f"ratings fall into a single category (perfect marginal agreement)."
        )

    return numerator / denominator


def cohen_null_variance(
    p_a: np.ndarray,
    p_b: np.ndarray,
    n_subjects: int,
    weight_matrix: np.ndarray | None = None,
) -> float:
    """
    Null variance of Cohen's (weighted) kappa under H₀: κ_w = 0.

    I am the large-sample null variance for two raters with individual
    marginals, supporting both unweighted and weighted kappa. Under the null,
    the two raters are independent, so the joint probability factorises as
    p_ij = p_a_i * p_b_j.

    The formula (weighted case):

        σ₀²(κ_w) = [1 / (N * (1 - Pe_w)²)]
                  × [Σ_ij p_a_i * p_b_j * (w_ij - w̄_i· - w̄·_j)² - Pe_w²]

    where:
        w̄_i·  = Σ_j  w_ij * p_b_j     (row-weighted mean of weights)
        w̄·_j  = Σ_i  w_ij * p_a_i     (column-weighted mean of weights)
        Pe_w   = Σ_ij w_ij * p_a_i * p_b_j   (chance agreement)

    For unweighted kappa (w_ij = δ_ij), this reduces to Equation 14:

        σ₀²(κ) = [Pe + Pe² - Σ_i p_a_i * p_b_i * (p_a_i + p_b_i)]
               / [N * (1 - Pe)²]

    Provenance chain -- the formula has been independently verified:

        1. Cohen (1960, 1968) published standard error formulas for kappa
           and weighted kappa. These were shown to be INCORRECT.
        2. Fleiss, Cohen & Everitt (1969) derived CORRECTED formulas
           (Equations 9 and 14). This is the source of the formula above.
        3. Landis & Koch (1977a) CONFIRMED the 1969 corrections analytically.
        4. Cicchetti & Fleiss (1977) and Fleiss & Cicchetti (1978) CONFIRMED
           them by Monte Carlo simulation.
        5. Fleiss, Nee & Landis (1979) explicitly states that the 1969 two-
           rater formulas "have been confirmed analytically by Landis and
           Koch (1977a) and have been confirmed by means of Monte Carlo
           simulation" (p. 974). The 1979 paper corrects a DIFFERENT formula
           (Fleiss' 1971 multi-rater kappa variance), not the 1969 two-rater
           formula implemented here.

    This variance is valid for:
    - Testing H₀: κ_w = 0 (hypothesis testing)
    - Parametric null confidence intervals
    - Large N (asymptotic approximation)
    - Exactly two raters with individual marginals
    - Both weighted and unweighted kappa

    This variance is NOT valid for:
    - Nonnull κ (underestimates variance when κ > 0)
    - Small N (asymptotic approximation breaks down)
    - More than two raters (use Conger/Fleiss instead)

    Parameters
    ----------
    p_a : np.ndarray
        Rater A's marginal category proportions, shape (K,). Must sum to 1.
    p_b : np.ndarray
        Rater B's marginal category proportions, shape (K,). Must sum to 1.
    n_subjects : int
        Number of subjects both raters rated (N in the paper's notation).
    weight_matrix : np.ndarray or None
        K×K weight matrix. None means identity (unweighted kappa).
        Must satisfy w_ij = w_ji (symmetric) and w_ii = 1 (diagonal).

    Returns
    -------
    float
        Var(κ_w) under the null hypothesis.

    Raises
    ------
    ValueError
        If p_a or p_b does not sum to 1 (within tolerance).
        If p_a and p_b have different shapes.
        If n_subjects is not positive.
        If weight_matrix has wrong shape.
        If denominator is zero (perfect marginal concentration).

    References
    ----------
    Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). "Large sample
        standard errors of kappa and weighted kappa." *Psychological
        Bulletin*, 72(5), 323-327. Equations 9 (weighted) and 14
        (unweighted).

    Cohen, J. (1960). "A coefficient of agreement for nominal scales."
        *Educational and Psychological Measurement*, 20(1), 37-46.

    Cohen, J. (1968). "Weighted kappa: nominal scale agreement with
        provision for scaled disagreement or partial credit."
        *Psychological Bulletin*, 70(4), 213-220.

    Fleiss, J. L., Nee, J. C. M., & Landis, J. R. (1979). "Large sample
        variance of kappa in the case of different sets of raters."
        *Psychological Bulletin*, 86(5), 974-977. Confirms the 1969
        two-rater formulas on p. 974.
    """
    # Validate inputs
    if n_subjects <= 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_subjects must be positive, got {n_subjects}."
        )
    if not np.isclose(p_a.sum(), 1.0, atol=1e-6):
        raise ValueError(
            f"{_ERROR_PREFIX}: p_a must sum to 1, got {p_a.sum():.6f}."
        )
    if not np.isclose(p_b.sum(), 1.0, atol=1e-6):
        raise ValueError(
            f"{_ERROR_PREFIX}: p_b must sum to 1, got {p_b.sum():.6f}."
        )
    if p_a.shape != p_b.shape:
        raise ValueError(
            f"{_ERROR_PREFIX}: p_a and p_b must have the same shape, "
            f"got {p_a.shape} and {p_b.shape}."
        )

    K = p_a.shape[0]

    # Build or validate weight matrix
    if weight_matrix is None:
        W = np.eye(K)
    else:
        W = weight_matrix
        if W.shape != (K, K):
            raise ValueError(
                f"{_ERROR_PREFIX}: weight_matrix must be ({K}, {K}), "
                f"got {W.shape}."
            )

    # Fleiss, Cohen & Everitt (1969):
    #
    # w̄_i·  = Σ_j  w_ij * p_b_j       row-weighted average of weights
    # w̄·_j  = Σ_i  w_ij * p_a_i       column-weighted average of weights
    # Pe_w   = Σ_ij w_ij * p_a_i * p_b_j   chance agreement under weighting

    w_bar_row = W @ p_b            # shape (K,) -- w̄_i· for each i
    w_bar_col = W.T @ p_a          # shape (K,) -- w̄·_j for each j
    Pe_w = float(p_a @ W @ p_b)    # scalar

    # diff_ij = w_ij - w̄_i· - w̄·_j
    diff = W - w_bar_row[:, np.newaxis] - w_bar_col[np.newaxis, :]

    # Σ_ij p_a_i * p_b_j * diff_ij²
    joint_null = p_a[:, np.newaxis] * p_b[np.newaxis, :]  # (K, K)
    sum_term = float(np.sum(joint_null * diff ** 2))

    numerator = sum_term - Pe_w ** 2
    denominator = n_subjects * (1.0 - Pe_w) ** 2

    if denominator == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: denominator is zero. This happens when "
            f"Pe_w = 1, i.e. chance agreement is perfect because all "
            f"probability mass concentrates on one category."
        )

    return numerator / denominator
