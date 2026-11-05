"""Profile describing the balance structure of ratings."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ...array_utils import (
    is_balanced as check_balanced,
    n_per_subject as get_n_per_subject,
    subject_count_distribution,
)

if TYPE_CHECKING:
    from ...data_objects.ratings import Ratings
    from ...intermediates.estimator_parts import NijMatrix

_ERROR_PREFIX = "[BalanceProfile] Error"


@dataclass(frozen=True)
class BalanceProfile:
    """
    I describe whether ratings are balanced and how they are distributed.

    I answer structural questions about data: Are all subjects rated the same
    number of times? If not, what's the distribution? How many subjects did
    each rater see? I exist because visualization and interpretation depend on
    understanding the data shape, not just the coefficient value.

    If is_balanced is True, the data is uniform and I report it simply. If
    False, I hold the granular distribution for diagnostics and visualization.

    I describe subjects only: how many ratings each one got. That much is in
    a counts matrix as well as in Ratings, so I am always built, from
    either. Who the raters are, how much each one rated and who overlapped
    with whom need rater identity, which a counts matrix does not have;
    that is RaterProfile, built only from Ratings.

    Attributes
    ----------
    is_balanced : bool
        Every subject was rated the same number of times.
    n_per_subject : ndarray
        (n_subjects,) -- how many raters rated each subject.
    subject_counts : dict[int, int]
        Distribution: {n_raters: count}. Example: {2: 1614, 3: 3} means 1614
        subjects had 2 raters, 3 subjects had 3 raters.
    raters_per_subject : int | None
        If is_balanced, the constant number of raters per subject. None if
        unbalanced.
    """

    is_balanced: bool
    n_per_subject: np.ndarray
    subject_counts: dict[int, int]
    raters_per_subject: int | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction."""

        # n_per_subject must be a numpy array of positive integers
        if not isinstance(self.n_per_subject, np.ndarray):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_per_subject must be numpy array, "
                f"got {type(self.n_per_subject).__name__}"
            )
        if self.n_per_subject.size == 0:
            raise ValueError(f"{_ERROR_PREFIX}: n_per_subject is empty")
        if np.any(self.n_per_subject <= 0):
            raise ValueError(
                f"{_ERROR_PREFIX}: all values in n_per_subject must be > 0, "
                f"got min {self.n_per_subject.min()}"
            )

        # subject_counts must be a dict with positive int keys and values
        if not isinstance(self.subject_counts, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: subject_counts must be dict, "
                f"got {type(self.subject_counts).__name__}"
            )
        for k, v in self.subject_counts.items():
            if not isinstance(k, int) or k <= 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subject_counts keys must be positive ints, got {k}"
                )
            if not isinstance(v, int) or v <= 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subject_counts values must be positive ints, got {v}"
                )

        # is_balanced consistency: if balanced, raters_per_subject should match
        if self.is_balanced:
            if self.raters_per_subject is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: if is_balanced=True, raters_per_subject must not be None"
                )
            if not isinstance(self.raters_per_subject, int) or self.raters_per_subject <= 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: raters_per_subject must be positive int, "
                    f"got {self.raters_per_subject}"
                )
            # All n_per_subject should equal raters_per_subject
            if not np.all(self.n_per_subject == self.raters_per_subject):
                raise ValueError(
                    f"{_ERROR_PREFIX}: is_balanced=True but n_per_subject is not constant "
                    f"(got values {np.unique(self.n_per_subject)})"
                )
        else:
            # If unbalanced, raters_per_subject should be None
            if self.raters_per_subject is not None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: if is_balanced=False, raters_per_subject must be None, "
                    f"got {self.raters_per_subject}"
                )
            # subject_counts should reflect the distribution of n_per_subject
            unique_counts, count_frequencies = np.unique(self.n_per_subject, return_counts=True)
            observed_counts = {int(k): int(v) for k, v in zip(unique_counts, count_frequencies)}
            if observed_counts != self.subject_counts:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subject_counts does not match n_per_subject distribution "
                    f"(expected {observed_counts}, got {self.subject_counts})"
                )

    # Properties ------------------------------------------------------------- #

    @property
    def n_subjects(self) -> int:
        """Return the number of subjects."""
        return len(self.n_per_subject)

    @property
    def n_raters_lower_bound(self) -> int:
        """
        At least this many raters -- max of n_per_subject.

        All a counts matrix can say about the raters: a rater who rated zero
        subjects is invisible in the counts, so the true number is at least
        this but may be higher. The exact number is RaterProfile.n_raters,
        when there is rater identity.
        """
        return int(self.n_per_subject.max())

    # Public API ------------------------------------------------------------- #

    @staticmethod
    def from_ratings(ratings: Ratings) -> BalanceProfile:
        """
        Construct a BalanceProfile from a Ratings object.

        Parameters
        ----------
        ratings : Ratings
            The ratings to analyze.

        Returns
        -------
        BalanceProfile
            A profile describing the balance structure.
        """
        index = ratings.category_index

        is_balanced_flag = check_balanced(index)
        n_per_subject_arr = get_n_per_subject(index)
        subject_counts_dict = subject_count_distribution(index)

        raters_per_subject = int(n_per_subject_arr[0]) if is_balanced_flag else None

        return BalanceProfile(
            is_balanced=is_balanced_flag,
            n_per_subject=n_per_subject_arr,
            subject_counts=subject_counts_dict,
            raters_per_subject=raters_per_subject,
        )

    @staticmethod
    def from_nij_matrix(nij: NijMatrix) -> BalanceProfile:
        """
        Construct a BalanceProfile from an NijMatrix (counts matrix).

        Everything I hold is subject-level, so it is all in the counts.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).

        Returns
        -------
        BalanceProfile
            A profile describing the balance structure.
        """
        n_per_subject_arr = nij.n_per_subject

        if int(np.max(n_per_subject_arr)) == 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: all subjects have zero raters "
                f"(every row of the nij matrix sums to 0)"
            )

        is_balanced_flag = len(np.unique(n_per_subject_arr)) == 1

        unique_counts, count_frequencies = np.unique(n_per_subject_arr, return_counts=True)
        subject_counts_dict = {int(k): int(v) for k, v in zip(unique_counts, count_frequencies)}

        raters_per_subject = int(n_per_subject_arr[0]) if is_balanced_flag else None

        return BalanceProfile(
            is_balanced=is_balanced_flag,
            n_per_subject=n_per_subject_arr,
            subject_counts=subject_counts_dict,
            raters_per_subject=raters_per_subject,
        )

    def summary(self) -> dict:
        """
        Return a plain-type summary of the balance structure.

        Returns
        -------
        dict
            Keys always include 'is_balanced', 'n_subjects',
            'n_raters_lower_bound'. If balanced, adds 'raters_per_subject'.
            If unbalanced, adds 'subject_counts'.
        """
        result: dict = {
            "is_balanced": self.is_balanced,
            "n_subjects": self.n_subjects,
            "n_raters_lower_bound": self.n_raters_lower_bound,
        }

        if self.is_balanced:
            result["raters_per_subject"] = self.raters_per_subject
        else:
            result["subject_counts"] = dict(self.subject_counts)

        return result
