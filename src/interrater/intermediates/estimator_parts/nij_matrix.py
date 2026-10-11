"""The counts matrix: how many raters chose each category per subject."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_NIJ_ERROR = "[NijMatrix] Error"


@dataclass(frozen=True, eq=False)
class NijMatrix:
    """
    I am the counts matrix: category tallies per subject.

    One row per subject, one column per category. Entry (i, j) is how many
    raters chose category j for subject i. Rater identity is discarded; only
    counts remain.

    I am the sole input to group coefficients. Pairwise coefficients read
    EstimatorInput directly and compute me on demand, which is why I exist
    separately: group coefficients are written around me, and they need me
    to be the thing they read.

    Attributes
    ----------
    matrix : np.ndarray
        A numpy array of shape (n_subjects, n_categories), dtype int, non-negative.
    """

    matrix: np.ndarray

    def __post_init__(self) -> None:
        """Validate the counts matrix."""
        m = self.matrix
        if m.ndim != 2:
            raise ValueError(
                f"{_NIJ_ERROR}: matrix must be 2-dimensional "
                f"(subjects × categories), got shape {m.shape}."
            )

        if m.dtype not in (np.int8, np.int16, np.int32, np.int64, int):
            raise ValueError(
                f"{_NIJ_ERROR}: matrix must hold integers (counts), got dtype {m.dtype}."
            )

        if (m < 0).any():
            raise ValueError(f"{_NIJ_ERROR}: counts cannot be negative.")

    def to_dataframe(self):
        """
        Return the counts matrix as a bare pandas DataFrame. No column names, no index names -- just the raw integer counts
        so the caller gets a clean numeric frame without metadata.
        Returns
        -------
        pd.DataFrame
            Shape (n_subjects, n_categories), values are counts.
        """
        import pandas as pd

        return pd.DataFrame(self.matrix)

    @property
    def n_subjects(self) -> int:
        """Number of subjects."""
        return int(self.matrix.shape[0])

    @property
    def n_categories(self) -> int:
        """Number of categories."""
        return int(self.matrix.shape[1])

    @property
    def n_per_subject(self) -> np.ndarray:
        """
        Number of raters per subject.

        Returns array of shape (n_subjects,) where element i is the sum of
        counts across categories for subject i (total raters for that subject).
        """
        return np.sum(self.matrix, axis=1)

    @classmethod
    def build_from_dataframe(cls, df) -> NijMatrix:
        """
        Build NijMatrix from a pandas DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Shape (n_subjects, n_categories), values are counts.
            Row index: subject identifiers (for reference only).
            Column names: category labels (for reference only).
            All values must be non-negative integers.

        Returns
        -------
        NijMatrix
            Instance created from df.values (the underlying numpy array).

        Raises
        ------
        TypeError
            If df is not a pandas DataFrame.
        ValueError
            If df contains non-integer or negative values.
        """
        import pandas as pd

        if not isinstance(df, pd.DataFrame):
            raise TypeError(
                f"{_NIJ_ERROR}: df must be a pandas DataFrame, "
                f"got {type(df).__name__}."
            )

        # Extract numpy array from DataFrame
        matrix = df.values.astype(np.int64)

        # Create and return NijMatrix (validation happens in __post_init__)
        return cls(matrix=matrix)
