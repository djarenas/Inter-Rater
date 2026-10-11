"""Permutation resampling for null hypothesis testing.

This module holds generators that produce permuted versions of the ratings
data under the null hypothesis of no agreement. Each scheme shuffles the
data differently and therefore tests a subtly different null hypothesis.
See the docstring of each function for what it tests.

These are conceptually different from bootstrap resampling:
- Bootstrap: resample WITH replacement from observed data →
  sampling distribution under the ALTERNATIVE (for CIs)
- Permutation: shuffle WITHOUT replacement to destroy structure →
  distribution under the NULL (for hypothesis testing)

All three schemes operate on category_index (the raw ratings array,
shape n_subjects × n_raters, integer-coded). They yield one permuted
array at a time for memory efficiency.

References
----------
Efron, B. & Tibshirani, R. J. (1993). *An Introduction to the Bootstrap*.
    Chapman & Hall. Chapter 15: Permutation tests.
"""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np

_ERROR_PREFIX = "[permutation_resampling] Error"


def column_shuffle(
    category_index: np.ndarray,
    n_permutations: int = 5000,
    random_state: int | None = 42,
) -> Iterator[np.ndarray]:
    """
    Shuffle each rater's column independently.

    For each rater (column), I randomly reassign their ratings to different
    subjects. This preserves each rater's marginal distribution -- how often
    they use each category -- but breaks the subject-rater pairing.

    **Null hypothesis tested:** "Raters agree no more than expected given
    their individual category-usage biases." This is the strongest test
    that respects rater-level marginals.

    Parameters
    ----------
    category_index : np.ndarray
        Shape (n_subjects, n_raters), dtype int. The raw ratings array
        where -1 indicates missing. Validated already by Ratings.
    n_permutations : int, default 5000
        Number of permuted arrays to yield.
    random_state : int | None, default 42
        Random seed for reproducibility.

    Yields
    ------
    np.ndarray
        Shape (n_subjects, n_raters). A permuted ratings array with each
        column shuffled independently.

    Raises
    ------
    TypeError
        If category_index is not np.ndarray, n_permutations is not int,
        or random_state is not int or None.
    ValueError
        If category_index is not 2D or n_permutations < 1.
    """
    # Type checks
    if not isinstance(category_index, np.ndarray):
        raise TypeError(
            f"{_ERROR_PREFIX}: category_index must be np.ndarray, "
            f"got {type(category_index).__name__}"
        )
    if not isinstance(n_permutations, int):
        raise TypeError(
            f"{_ERROR_PREFIX}: n_permutations must be int, "
            f"got {type(n_permutations).__name__}"
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
    if n_permutations < 1:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_permutations must be >= 1, "
            f"got {n_permutations}"
        )

    rng = np.random.default_rng(random_state)
    n_subjects, n_raters = category_index.shape

    for _ in range(n_permutations):
        permuted = category_index.copy()
        for col in range(n_raters):
            # Find non-missing entries for this rater
            valid_mask = permuted[:, col] >= 0
            valid_values = permuted[valid_mask, col]
            # Shuffle only the valid entries
            rng.shuffle(valid_values)
            permuted[valid_mask, col] = valid_values
        yield permuted


def within_subject_shuffle(
    category_index: np.ndarray,
    n_permutations: int = 5000,
    random_state: int | None = 42,
) -> Iterator[np.ndarray]:
    """
    Shuffle category assignments within each subject's row.

    For each subject (row), I randomly permute the ratings across rater
    positions. This preserves the number of raters per subject and the
    overall count of each category within that subject, but breaks the
    association between rater identity and category choice.

    **Null hypothesis tested:** "The pattern of which rater chose which
    category within a subject is random." This tests whether rater
    identity matters for the category chosen.

    Note: for group coefficients that discard rater identity (Fleiss,
    Scott, AC1), this shuffle does not change the nij matrix and therefore
    does not change the coefficient. It is more meaningful for pairwise
    coefficients (Cohen's kappa) where rater identity matters.

    Parameters
    ----------
    category_index : np.ndarray
        Shape (n_subjects, n_raters), dtype int. The raw ratings array
        where -1 indicates missing. Validated already by Ratings.
    n_permutations : int, default 5000
        Number of permuted arrays to yield.
    random_state : int | None, default 42
        Random seed for reproducibility.

    Yields
    ------
    np.ndarray
        Shape (n_subjects, n_raters). A permuted ratings array with each
        row shuffled independently.

    Raises
    ------
    TypeError
        If category_index is not np.ndarray, n_permutations is not int,
        or random_state is not int or None.
    ValueError
        If category_index is not 2D or n_permutations < 1.
    """
    # Type checks
    if not isinstance(category_index, np.ndarray):
        raise TypeError(
            f"{_ERROR_PREFIX}: category_index must be np.ndarray, "
            f"got {type(category_index).__name__}"
        )
    if not isinstance(n_permutations, int):
        raise TypeError(
            f"{_ERROR_PREFIX}: n_permutations must be int, "
            f"got {type(n_permutations).__name__}"
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
    if n_permutations < 1:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_permutations must be >= 1, "
            f"got {n_permutations}"
        )

    rng = np.random.default_rng(random_state)
    n_subjects, n_raters = category_index.shape

    for _ in range(n_permutations):
        permuted = category_index.copy()
        for row in range(n_subjects):
            # Find non-missing entries for this subject
            valid_mask = permuted[row, :] >= 0
            valid_values = permuted[row, valid_mask]
            # Shuffle only the valid entries
            rng.shuffle(valid_values)
            permuted[row, valid_mask] = valid_values
        yield permuted


def pool_and_redistribute(
    category_index: np.ndarray,
    n_permutations: int = 5000,
    random_state: int | None = 42,
) -> Iterator[np.ndarray]:
    """
    Pool all ratings and redistribute randomly across the matrix.

    I collect all non-missing ratings into a single pool, shuffle them,
    and redistribute them back into the original matrix positions (respecting
    which cells are missing). This preserves only the global category
    proportions -- all other structure is destroyed.

    **Null hypothesis tested:** "There is no structure in the data beyond
    the overall category proportions." This is the strongest null -- it
    breaks both subject-rater pairing AND rater-level marginals.

    Parameters
    ----------
    category_index : np.ndarray
        Shape (n_subjects, n_raters), dtype int. The raw ratings array
        where -1 indicates missing. Validated already by Ratings.
    n_permutations : int, default 5000
        Number of permuted arrays to yield.
    random_state : int | None, default 42
        Random seed for reproducibility.

    Yields
    ------
    np.ndarray
        Shape (n_subjects, n_raters). A permuted ratings array with all
        non-missing values pooled and redistributed randomly.

    Raises
    ------
    TypeError
        If category_index is not np.ndarray, n_permutations is not int,
        or random_state is not int or None.
    ValueError
        If category_index is not 2D or n_permutations < 1.
    """
    # Type checks
    if not isinstance(category_index, np.ndarray):
        raise TypeError(
            f"{_ERROR_PREFIX}: category_index must be np.ndarray, "
            f"got {type(category_index).__name__}"
        )
    if not isinstance(n_permutations, int):
        raise TypeError(
            f"{_ERROR_PREFIX}: n_permutations must be int, "
            f"got {type(n_permutations).__name__}"
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
    if n_permutations < 1:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_permutations must be >= 1, "
            f"got {n_permutations}"
        )

    rng = np.random.default_rng(random_state)

    # Find all non-missing positions and their values
    valid_mask = category_index >= 0
    valid_values = category_index[valid_mask].copy()

    for _ in range(n_permutations):
        # Shuffle all valid values
        rng.shuffle(valid_values)
        # Place back into original positions
        permuted = category_index.copy()
        permuted[valid_mask] = valid_values
        yield permuted
