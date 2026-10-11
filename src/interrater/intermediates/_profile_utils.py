"""
Arithmetic over what raters actually chose.

Every function here reads category values, which is what separates them from
``interrater.array_utils`` at the root: those read only whether a cell holds a
rating, and are needed by data objects and extractors as well as by me. These
are reporting arithmetic and belong with the profile that reports them.

Free functions: numpy in, numpy out. Nothing here imports a data object, so a
bootstrapper or an estimator can call these with a bare array. Nothing here
validates either.
"""

from __future__ import annotations

import numpy as np

from ..array_utils import category_counts, n_per_rater, present_mask


def prevalence(category_index: np.ndarray, k: int) -> np.ndarray:
    """
    Share of present ratings in each category.

    Empties are excluded from numerator and denominator alike, so this says how
    raters used the categories, not how much of the data was covered.

    Skew here is the thing to watch. When one category dominates, expected
    agreement is already high and kappa is left dividing by a small remainder,
    so it falls even when raters agreed on almost everything.
    """
    present = category_index[present_mask(category_index)]
    counts  = np.bincount(present.astype(np.int64), minlength=k)
    total   = counts.sum()
    if total == 0:
        return np.zeros(k, dtype=np.float64)
    return counts / total


def prevalence_per_rater(category_index: np.ndarray, k: int) -> np.ndarray:
    """
    Per-rater category shares, (n_raters, k).

    A row unlike the others means that rater used the categories differently:
    systematic bias rather than noise, which kappa cannot see but a confusion
    matrix can.
    """
    n_raters = category_index.shape[1]
    out      = np.zeros((n_raters, k), dtype=np.float64)
    for r in range(n_raters):
        out[r] = prevalence(category_index[:, r : r + 1], k)
    return out


def coverage_per_rater(category_index: np.ndarray) -> np.ndarray:
    """Fraction of subject rows each rater rated."""
    return n_per_rater(category_index) / category_index.shape[0]
