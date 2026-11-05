"""Profile describing the raters: who they are and how much each one rated."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ...array_utils import n_per_rater as get_n_per_rater
from .rater_shared_matrix import RaterSharedMatrix

if TYPE_CHECKING:
    from ...data_objects.ratings import Ratings

_ERROR_PREFIX = "[RaterProfile] Error"


@dataclass(frozen=True, eq=False)
class RaterProfile:
    """
    I describe the raters: who they are, how much each one rated, and how
    much each pair of them overlapped.

    I need rater identity, so I am built only from Ratings. A counts matrix
    keeps how many raters chose each category, not who they were, so a
    result built from one holds NotComputed in my place. What a counts
    matrix does know -- how many ratings each subject got -- is
    BalanceProfile.

    I am the place a result's rater names come from: the printed tables,
    to_frame, summaries and plots of pairs and subsets name raters as I do.

    Attributes
    ----------
    rater_names : tuple[str, ...]
        The raters, in position order, as the semantics name them (injected
        raters included). Position i is column i of the ratings.
    n_subjects : int
        How many subjects there were to rate.
    n_per_rater : ndarray
        (n_raters,) -- how many subjects each rater rated.
    rater_overlap : RaterSharedMatrix
        How many subjects each pair of raters both rated.
    """

    rater_names: tuple[str, ...]
    n_subjects: int
    n_per_rater: np.ndarray
    rater_overlap: RaterSharedMatrix

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes, then keep my own copy of the counts."""
        names = self.rater_names
        if not isinstance(names, tuple) or not all(isinstance(n, str) for n in names):
            raise TypeError(
                f"{_ERROR_PREFIX}: rater_names must be a tuple of str, got {names!r}"
            )
        if len(set(names)) != len(names):
            raise ValueError(f"{_ERROR_PREFIX}: rater_names must not repeat, got {names!r}")

        if not isinstance(self.n_subjects, int) or isinstance(self.n_subjects, bool) or self.n_subjects < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_subjects must be a positive int, got {self.n_subjects!r}"
            )

        if not isinstance(self.n_per_rater, np.ndarray):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_per_rater must be numpy array, "
                f"got {type(self.n_per_rater).__name__}"
            )
        if self.n_per_rater.shape != (len(names),):
            raise ValueError(
                f"{_ERROR_PREFIX}: n_per_rater has shape {self.n_per_rater.shape} "
                f"but there are {len(names)} raters"
            )
        if np.any(self.n_per_rater < 0) or np.any(self.n_per_rater > self.n_subjects):
            raise ValueError(
                f"{_ERROR_PREFIX}: each rater rated between 0 and the "
                f"{self.n_subjects} subjects, but n_per_rater is "
                f"{self.n_per_rater.tolist()}"
            )

        if not isinstance(self.rater_overlap, RaterSharedMatrix):
            raise TypeError(
                f"{_ERROR_PREFIX}: rater_overlap must be RaterSharedMatrix, "
                f"got {type(self.rater_overlap).__name__}"
            )
        # The overlap's diagonal is each rater's own total
        if not np.array_equal(self.rater_overlap.per_rater, self.n_per_rater):
            raise ValueError(
                f"{_ERROR_PREFIX}: rater_overlap says the raters rated "
                f"{self.rater_overlap.per_rater.tolist()} subjects, but "
                f"n_per_rater says {self.n_per_rater.tolist()}"
            )

        object.__setattr__(self, "n_per_rater", np.array(self.n_per_rater, copy=True))

    # Properties ------------------------------------------------------------- #

    @property
    def n_raters(self) -> int:
        """How many raters there are, exactly."""
        return len(self.rater_names)

    @property
    def fill_rate(self) -> float:
        """
        Fraction of cells filled in the (n_subjects, n_raters) grid; 1.0
        means every rater rated every subject. Derived from the counts, so
        it cannot disagree with them.
        """
        return float(self.n_per_rater.sum()) / (self.n_subjects * self.n_raters)

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_ratings(cls, ratings: Ratings) -> RaterProfile:
        """
        Construct a RaterProfile from a Ratings object.

        Parameters
        ----------
        ratings : Ratings
            The ratings to describe.

        Returns
        -------
        RaterProfile
            Every field exact.
        """
        return cls(
            rater_names=tuple(str(name) for name in ratings.semantics.raters),
            n_subjects=ratings.n_subjects,
            n_per_rater=get_n_per_rater(ratings.category_index),
            rater_overlap=RaterSharedMatrix.build_from_ratings(ratings),
        )

    def summary(self) -> dict:
        """
        Return a plain-type summary of the raters.

        Returns
        -------
        dict
            rater_names, n_raters, n_subjects, n_per_rater, fill_rate, and
            the overlap matrix as nested lists.
        """
        return {
            "rater_names": list(self.rater_names),
            "n_raters": self.n_raters,
            "n_subjects": self.n_subjects,
            "n_per_rater": self.n_per_rater.tolist(),
            "fill_rate": round(self.fill_rate, 4),
            "rater_overlap": self.rater_overlap.shared.tolist(),
        }
