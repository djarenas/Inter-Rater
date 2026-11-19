"""Subject-level resampling: bootstrap and jackknife."""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np

from ...intermediates.estimator_parts.estimator_input import EstimatorInput

_ERROR_PREFIX = "[subject_resampling] Error"


def subject_bootstrap(
    data: EstimatorInput,
    n_resamples: int = 1000,
    random_state: int | None = 42,
) -> Iterator[EstimatorInput]:
    """
    Resample subjects (rows) with replacement, yielding one at a time.

    I draw subject rows (with replacement) from data,
    yielding n_resamples new datasets one at a time. Each keeps data's K and
    weights, so a draw that misses a category is still computed against the
    rating task's K. This is memory-efficient:
    only one resample exists in memory at any moment, even for 10,000 draws.

    The resampling is subject-level, never rater-level. A rater drawn twice is
    one rater agreeing perfectly with themselves, which inflates every
    coefficient. Subject resampling is the one that works.

    Parameters
    ----------
    data : EstimatorInput
        The ratings array, K and weights, validated already by Ratings.
    n_resamples : int, default 1000
        Number of bootstrap draws to yield.
    random_state : int | None, default 42
        Random seed for reproducibility. Default of 42 ensures reproducible
        results unless explicitly overridden. Pass None for unseeded randomness.

    Yields
    ------
    EstimatorInput
        Each iteration yields one resample with n_subjects rows.

    Raises
    ------
    TypeError
        If data is not an EstimatorInput, n_resamples is not int, or
        random_state is not int or None.
    ValueError
        If n_resamples < 1.

    Examples
    --------
    >>> import numpy as np
    >>> data = EstimatorInput(np.array([[0, 1], [1, 0], [1, 1]], dtype=np.int8), 2)
    >>> for i, resample in enumerate(subject_bootstrap(data, n_resamples=5)):
    ...     print(f"Draw {i}: shape {resample.category_index.shape}")
    Draw 0: shape (3, 2)
    Draw 1: shape (3, 2)
    Draw 2: shape (3, 2)
    Draw 3: shape (3, 2)
    Draw 4: shape (3, 2)
    """
    # Type checks
    _check_data(data)
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
    if n_resamples < 1:
        raise ValueError(
            f"{_ERROR_PREFIX}: n_resamples must be >= 1, "
            f"got {n_resamples}"
        )

    rng = np.random.default_rng(random_state)
    n_subjects = data.n_subjects

    for _ in range(n_resamples):
        # Draw subject indices with replacement
        indices = rng.integers(0, n_subjects, size=n_subjects)
        # Yield this resample (numpy copies, never views)
        yield data.take_subjects(indices)


def subject_jackknife(data: EstimatorInput) -> Iterator[EstimatorInput]:
    """
    Leave-one-subject-out: yield data with each subject removed in turn.

    I iterate through all subjects, yielding the data with one subject removed
    each time. This gives n_subjects datasets, each missing one observation.

    Useful for bias diagnostics and stability analysis: how much does the
    estimate change when we drop each subject?

    Parameters
    ----------
    data : EstimatorInput
        The ratings array, K and weights.

    Yields
    ------
    EstimatorInput
        Each iteration yields data with one subject removed: n_subjects - 1
        rows, same K and weights.

    Raises
    ------
    TypeError
        If data is not an EstimatorInput.
    ValueError
        If data has fewer than 2 subjects.

    Examples
    --------
    >>> import numpy as np
    >>> data = EstimatorInput(np.array([[0, 1], [1, 0], [1, 1]], dtype=np.int8), 2)
    >>> for i, jackknife_draw in enumerate(subject_jackknife(data)):
    ...     print(f"Omit subject {i}: shape {jackknife_draw.category_index.shape}")
    Omit subject 0: shape (2, 2)
    Omit subject 1: shape (2, 2)
    Omit subject 2: shape (2, 2)
    """
    # Type checks
    _check_data(data)

    # Value checks
    if data.n_subjects < 2:
        raise ValueError(
            f"{_ERROR_PREFIX}: subject_jackknife requires >= 2 "
            f"subjects, got {data.n_subjects}"
        )

    n_subjects = data.n_subjects

    for omit_idx in range(n_subjects):
        # Keep all subjects except omit_idx
        keep_idx = np.concatenate([
            np.arange(0, omit_idx),
            np.arange(omit_idx + 1, n_subjects)
        ])
        yield data.take_subjects(keep_idx)


def _check_data(data: object) -> None:
    if not isinstance(data, EstimatorInput):
        raise TypeError(
            f"{_ERROR_PREFIX}: data must be EstimatorInput, "
            f"got {type(data).__name__}"
        )
