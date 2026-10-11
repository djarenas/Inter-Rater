"""Subject-level resampling: bootstrap and jackknife."""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np

_ERROR_PREFIX = "[subject_resampling] Error"


def subject_bootstrap(
    category_index: np.ndarray,
    n_resamples: int = 1000,
    random_state: int | None = 42,
) -> Iterator[np.ndarray]:
    """
    Resample subjects (rows) with replacement, yielding one at a time.

    I draw subject rows (with replacement) from the category_index array,
    yielding n_resamples new datasets one at a time. This is memory-efficient:
    only one resample exists in memory at any moment, even for 10,000 draws.

    The resampling is subject-level, never rater-level. A rater drawn twice is
    one rater agreeing perfectly with themselves, which inflates every
    coefficient. Subject resampling is the one that works.

    Parameters
    ----------
    category_index : np.ndarray
        Shape (n_subjects, n_raters), dtype int8. The ratings array where each
        entry is a category index or EMPTY. This is the raw data, validated
        already by Ratings.
    n_resamples : int, default 1000
        Number of bootstrap draws to yield.
    random_state : int | None, default 42
        Random seed for reproducibility. Default of 42 ensures reproducible
        results unless explicitly overridden. Pass None for unseeded randomness.

    Yields
    ------
    np.ndarray
        Each iteration yields one resampled array of shape (n_subjects, n_raters)
        and dtype int8.

    Raises
    ------
    TypeError
        If category_index is not np.ndarray, n_resamples is not int, or
        random_state is not int or None.
    ValueError
        If category_index is not 2D or n_resamples < 1.

    Examples
    --------
    >>> import numpy as np
    >>> ratings = np.array([[0, 1], [1, 0], [1, 1]], dtype=np.int8)
    >>> for i, resample in enumerate(subject_bootstrap(ratings, n_resamples=5)):
    ...     print(f"Draw {i}: shape {resample.shape}")
    Draw 0: shape (3, 2)
    Draw 1: shape (3, 2)
    ...
    """
    # Type checks
    if not isinstance(category_index, np.ndarray):
        raise TypeError(
            f"{_ERROR_PREFIX}: category_index must be np.ndarray, "
            f"got {type(category_index).__name__}"
        )
    if not isinstance(n_resamples, int):
        raise TypeError(
            f"{_ERROR_PREFIX}: n_resamples must be int, "
            f"got {type(n_resamples).__name__}"
        )
    if random_state is not None and not isinstance(random_state, int):
        raise TypeError(
            f"{_ERROR_PREFIX}: random_state must be int or None, "
            f"got {type(random_state).__name__}"
        )

    # Value checks
    if category_index.ndim != 2:
        raise ValueError(
            f"{_ERROR_PREFIX}: category_index must be 2D, "
            f"got {category_index.ndim}D"
        )
    if n_resamples < 1:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_resamples must be >= 1, "
            f"got {n_resamples}"
        )

    rng = np.random.default_rng(random_state)
    n_subjects = category_index.shape[0]

    for _ in range(n_resamples):
        # Draw subject indices with replacement
        indices = rng.integers(0, n_subjects, size=n_subjects)
        # Yield this resample (numpy copies, never views)
        yield category_index[indices]


def subject_jackknife(category_index: np.ndarray) -> Iterator[np.ndarray]:
    """
    Leave-one-subject-out: yield data with each subject removed in turn.

    I iterate through all subjects, yielding the data with one subject removed
    each time. This gives n_subjects datasets, each missing one observation.

    Useful for bias diagnostics and stability analysis: how much does the
    estimate change when we drop each subject?

    Parameters
    ----------
    category_index : np.ndarray
        Shape (n_subjects, n_raters), dtype int8. The ratings array.

    Yields
    ------
    np.ndarray
        Each iteration yields data with one subject removed. Shape is
        (n_subjects - 1, n_raters), dtype int8.

    Examples
    --------
    >>> import numpy as np
    >>> ratings = np.array([[0, 1], [1, 0], [1, 1]], dtype=np.int8)
    Raises
    ------
    TypeError
        If category_index is not np.ndarray.
    ValueError
        If category_index is not 2D or has fewer than 2 subjects.

    >>> for i, jackknife_draw in enumerate(subject_jackknife(ratings)):
    ...     print(f"Omit subject {i}: shape {jackknife_draw.shape}")
    Omit subject 0: shape (2, 2)
    Omit subject 1: shape (2, 2)
    Omit subject 2: shape (2, 2)
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
    if category_index.shape[0] < 2:
        raise ValueError(
            f"{_ERROR_PREFIX}: subject_jackknife requires >= 2 "
            f"subjects, got {category_index.shape[0]}"
        )

    n_subjects = category_index.shape[0]

    for omit_idx in range(n_subjects):
        # Keep all subjects except omit_idx
        keep_idx = np.concatenate([
            np.arange(0, omit_idx),
            np.arange(omit_idx + 1, n_subjects)
        ])
        yield category_index[keep_idx]


def subject_bootstrap_from_nij_matrix(
    nij_matrix: np.ndarray,
    n_resamples: int = 1000,
    random_state: int | None = 42,
) -> Iterator[np.ndarray]:
    """
    Resample subjects (rows) from nij_matrix with replacement, yielding one at a time.

    I draw subject rows (with replacement) from the nij_matrix, yielding n_resamples
    new matrices one at a time. This is memory-efficient: only one resample exists
    in memory at any moment, even for 10,000 draws.

    The resampling is subject-level, never category-level. A subject drawn twice
    is one subject appearing twice, which is correct for bootstrap inference on
    group coefficients.

    Parameters
    ----------
    nij_matrix : np.ndarray
        Shape (n_subjects, n_categories), dtype int. The counts matrix where each
        row is a subject and each column is a category count. This is the raw data,
        validated already by NijMatrix.
    n_resamples : int, default 1000
        Number of bootstrap draws to yield.
    random_state : int | None, default 42
        Random seed for reproducibility. Default of 42 ensures reproducible
        results unless explicitly overridden. Pass None for unseeded randomness.

    Yields
    ------
    np.ndarray
        Shape (n_subjects, n_categories). A bootstrap resample of the nij_matrix
        with subjects drawn with replacement.

    Raises
    ------
    TypeError
        If nij_matrix is not np.ndarray, n_resamples is not int, or
        random_state is not int or None.
    ValueError
        If nij_matrix is not 2D or n_resamples < 1.

    References
    ----------
    Efron, Bradley. 1979. "Bootstrap Methods: Another Look at the Jackknife."
    *The Annals of Statistics*, 7(1), 1-26.
    """
    # Type checks
    if not isinstance(nij_matrix, np.ndarray):
        raise TypeError(
            f"{_ERROR_PREFIX}: nij_matrix must be np.ndarray, "
            f"got {type(nij_matrix).__name__}"
        )
    if not isinstance(n_resamples, int):
        raise TypeError(
            f"{_ERROR_PREFIX}: n_resamples must be int, "
            f"got {type(n_resamples).__name__}"
        )
    if random_state is not None and not isinstance(random_state, int):
        raise TypeError(
            f"{_ERROR_PREFIX}: random_state must be int or None, "
            f"got {type(random_state).__name__}"
        )

    # Value checks
    if nij_matrix.ndim != 2:
        raise ValueError(
            f"{_ERROR_PREFIX}: nij_matrix must be 2D, "
            f"got {nij_matrix.ndim}D"
        )
    if n_resamples < 1:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_resamples must be >= 1, "
            f"got {n_resamples}"
        )

    rng = np.random.default_rng(random_state)
    n_subjects = nij_matrix.shape[0]

    for _ in range(n_resamples):
        # Sample subject indices with replacement
        resampled_indices = rng.choice(n_subjects, size=n_subjects, replace=True)
        # Yield the resampled nij_matrix (rows in new order)
        yield nij_matrix[resampled_indices]
