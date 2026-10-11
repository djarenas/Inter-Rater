"""Rater-level resampling: jackknife and subset."""

from __future__ import annotations

from itertools import combinations
from collections.abc import Iterator

import numpy as np

_ERROR_PREFIX = "[rater_resampling] Error"


def rater_jackknife(category_index: np.ndarray) -> Iterator[np.ndarray]:
    """
    Leave-one-rater-out: yield data with each rater removed in turn.

    I iterate through all raters, yielding the data with one rater removed
    each time. This gives n_raters datasets, each missing one rater's column.

    Useful for diagnostics: how much does agreement depend on each rater?
    Which rater is most critical? Does removing one rater reveal a problem?

    Parameters
    ----------
    category_index : np.ndarray
        Shape (n_subjects, n_raters), dtype int8. The ratings array.

    Yields
    ------
    np.ndarray
        Each iteration yields data with one rater's column removed. Shape is
        (n_subjects, n_raters - 1), dtype int8.

    Examples
    --------
    >>> import numpy as np
    >>> ratings = np.array([[0, 1, 0], [1, 0, 1]], dtype=np.int8)
    >>> for i, jackknife_draw in enumerate(rater_jackknife(ratings)):
    ...     print(f"Omit rater {i}: shape {jackknife_draw.shape}")
    Raises
    ------
    TypeError
        If category_index is not np.ndarray.
    ValueError
        If category_index is not 2D or has fewer than 2 raters.

    Omit rater 0: shape (2, 2)
    Omit rater 1: shape (2, 2)
    Omit rater 2: shape (2, 2)
    """
    # Type checks
    if not isinstance(category_index, np.ndarray):
        raise TypeError(
            f"{_ERROR_PREFIX}: category_index must be np.ndarray, "
            f"got {type(category_index).__name__}"
        )

    # Value checks
    if category_index.ndim != 2:
        raise ValueError(
            f"{_ERROR_PREFIX}: category_index must be 2D, "
            f"got {category_index.ndim}D"
        )
    if category_index.shape[1] < 2:
        raise ValueError(
            f"{_ERROR_PREFIX}: rater_jackknife requires >= 2 "
            f"raters, got {category_index.shape[1]}"
        )

    n_raters = category_index.shape[1]

    for omit_idx in range(n_raters):
        # Keep all raters except omit_idx
        keep_idx = np.concatenate([
            np.arange(0, omit_idx),
            np.arange(omit_idx + 1, n_raters)
        ])
        yield category_index[:, keep_idx]


def rater_subset(category_index: np.ndarray, k: int) -> Iterator[np.ndarray]:
    """
    All combinations of k raters from n: yield each subset in turn.

    I generate all C(n_raters, k) combinations of k raters and yield the data
    restricted to each combination. This is useful for understanding how
    agreement scales with panel size or for analyzing subgroups.

    Parameters
    ----------
    category_index : np.ndarray
        Shape (n_subjects, n_raters), dtype int8. The ratings array.
    k : int
        Size of each rater subset. Must be >= 1 and <= n_raters.

    Yields
    ------
    np.ndarray
        Each iteration yields data with only k raters' columns. Shape is
        (n_subjects, k), dtype int8.

    Raises
    ------
    TypeError
        If category_index is not np.ndarray or k is not int.
    ValueError
        If category_index is not 2D, or k < 1 or k > n_raters.

    Examples
    --------
    >>> import numpy as np
    >>> ratings = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=np.int8)
    >>> for i, subset in enumerate(rater_subset(ratings, k=2)):
    ...     print(f"Subset {i}: shape {subset.shape}")
    Subset 0: shape (3, 2)
    Subset 1: shape (3, 2)
    Subset 2: shape (3, 2)
    """
    # Type checks
    if not isinstance(category_index, np.ndarray):
        raise TypeError(
            f"{_ERROR_PREFIX}: category_index must be np.ndarray, "
            f"got {type(category_index).__name__}"
        )
    if not isinstance(k, int):
        raise TypeError(
            f"{_ERROR_PREFIX}: k must be int, "
            f"got {type(k).__name__}"
        )

    # Value checks
    if category_index.ndim != 2:
        raise ValueError(
            f"{_ERROR_PREFIX}: category_index must be 2D, "
            f"got {category_index.ndim}D"
        )

    n_raters = category_index.shape[1]

    if k < 1 or k > n_raters:
        raise ValueError(
            f"{_ERROR_PREFIX}: k must be between 1 and "
            f"n_raters ({n_raters}), got {k}"
        )

    # Generate all combinations of k rater indices
    for rater_indices in combinations(range(n_raters), k):
        rater_indices = np.array(rater_indices)
        yield category_index[:, rater_indices]
