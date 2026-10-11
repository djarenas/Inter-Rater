"""Profile describing the balance structure of ratings."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ...array_utils import (
    is_balanced as check_balanced,
    n_per_subject as get_n_per_subject,
    n_per_rater as get_n_per_rater,
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

    Some of my fields can be None. When constructed from a Ratings object,
    everything is known exactly -- Ratings knows every rater by name. When
    constructed from an nij matrix, rater identity is gone: I cannot know
    how many raters there were (only a lower bound), which rater gave which
    rating, or what fraction of the grid was filled. Those fields are None
    rather than estimates, because honest absence is better than fabricated
    data.

    Attributes
    ----------
    is_balanced : bool
        Every subject was rated the same number of times.
    n_raters : int | None
        Number of raters. Exact when built from Ratings, None when built
        from an nij matrix (use n_raters_lower_bound for what is known).
    n_raters_exact : bool
        True if n_raters is known exactly (from Ratings). False if n_raters
        is None or a lower bound.
    n_per_subject : ndarray
        (n_subjects,) -- how many raters rated each subject.
    n_per_rater : ndarray | None
        (n_raters,) -- how many subjects each rater rated. None when built
        from an nij matrix, because rater identity is not in the counts.
    subject_counts : dict[int, int]
        Distribution: {n_raters: count}. Example: {2: 1614, 3: 3} means 1614
        subjects had 2 raters, 3 subjects had 3 raters.
    fill_rate : float | None
        Fraction of cells filled in the (n_subjects, n_raters) grid. None
        when n_raters is unknown.
    raters_per_subject : int | None
        If is_balanced, the constant number of raters per subject. None if
        unbalanced.
    """

    is_balanced: bool
    n_raters: int | None
    n_raters_exact: bool
    n_per_subject: np.ndarray
    n_per_rater: np.ndarray | None
    subject_counts: dict[int, int]
    fill_rate: float | None
    raters_per_subject: int | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction."""

        # n_raters must be a positive integer or None
        if self.n_raters is not None:
            if not isinstance(self.n_raters, int) or self.n_raters <= 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: n_raters must be positive int or None, "
                    f"got {self.n_raters}"
                )

        # n_raters_exact must be bool
        if not isinstance(self.n_raters_exact, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_raters_exact must be bool, "
                f"got {type(self.n_raters_exact).__name__}"
            )

        # n_raters_exact=True requires n_raters to be set
        if self.n_raters_exact and self.n_raters is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_raters_exact=True but n_raters is None. "
                f"Cannot claim exact knowledge without a value."
            )

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

        # n_per_rater: numpy array of non-negative integers, or None
        if self.n_per_rater is not None:
            if not isinstance(self.n_per_rater, np.ndarray):
                raise TypeError(
                    f"{_ERROR_PREFIX}: n_per_rater must be numpy array or None, "
                    f"got {type(self.n_per_rater).__name__}"
                )
            if self.n_per_rater.size == 0:
                raise ValueError(f"{_ERROR_PREFIX}: n_per_rater is empty")
            if np.any(self.n_per_rater < 0):
                raise ValueError(
                    f"{_ERROR_PREFIX}: all values in n_per_rater must be >= 0, "
                    f"got min {self.n_per_rater.min()}"
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

        # fill_rate: float in [0, 1], or None
        if self.fill_rate is not None:
            if not isinstance(self.fill_rate, (int, float)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: fill_rate must be float or None, "
                    f"got {type(self.fill_rate).__name__}"
                )
            if not (0.0 <= self.fill_rate <= 1.0):
                raise ValueError(
                    f"{_ERROR_PREFIX}: fill_rate must be in [0, 1], got {self.fill_rate}"
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

        Always available, even when n_raters is None. From an nij matrix this
        is the best we can do: a rater who rated zero subjects is invisible in
        the counts, so the true number is at least this but may be higher.
        When n_raters is known exactly, this equals n_raters or is less (if
        some raters rated nothing, but Ratings forbids that so in practice
        they are equal).
        """
        return int(self.n_per_subject.max())

    @property
    def n_raters_display(self) -> str:
        """Return 'n_raters=X' if exact, 'n_raters>=X' if lower bound."""
        if self.n_raters is not None and self.n_raters_exact:
            return f"n_raters={self.n_raters}"
        return f"n_raters≥{self.n_raters_lower_bound}"

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
            A profile describing the balance structure. All fields are exact.
        """
        index = ratings.category_index
        n_subjects = ratings.n_subjects
        n_raters = ratings.n_raters

        is_balanced_flag = check_balanced(index)
        n_per_subject_arr = get_n_per_subject(index)
        n_per_rater_arr = get_n_per_rater(index)
        subject_counts_dict = subject_count_distribution(index)

        total_ratings = n_per_rater_arr.sum()
        fill_rate = float(total_ratings) / (n_subjects * n_raters)

        raters_per_subject = int(n_per_subject_arr[0]) if is_balanced_flag else None

        return BalanceProfile(
            is_balanced=is_balanced_flag,
            n_raters=n_raters,
            n_raters_exact=True,
            n_per_subject=n_per_subject_arr,
            n_per_rater=n_per_rater_arr,
            subject_counts=subject_counts_dict,
            fill_rate=fill_rate,
            raters_per_subject=raters_per_subject,
        )

    @staticmethod
    def from_nij_matrix(nij: NijMatrix) -> BalanceProfile:
        """
        Construct a BalanceProfile from an NijMatrix (counts matrix).

        Rater identity is gone in an nij matrix, so three fields are None:
        n_raters (use n_raters_lower_bound for what is known), n_per_rater
        (unknowable), and fill_rate (meaningless without n_raters). Everything
        subject-level is exact.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).

        Returns
        -------
        BalanceProfile
            A profile with subject-level fields exact and rater-level fields
            set to None.
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
            n_raters=None,
            n_raters_exact=False,
            n_per_subject=n_per_subject_arr,
            n_per_rater=None,
            subject_counts=subject_counts_dict,
            fill_rate=None,
            raters_per_subject=raters_per_subject,
        )

    def summary(self) -> dict:
        """
        Return a plain-type summary of the balance structure.

        Returns
        -------
        dict
            Keys always include 'is_balanced', 'n_subjects', 'n_raters_exact',
            'n_raters_lower_bound'. Adds 'n_raters' and 'fill_rate' when known.
            If balanced, adds 'raters_per_subject'. If unbalanced, adds
            'subject_counts'.
        """
        result: dict = {
            "is_balanced": self.is_balanced,
            "n_subjects": self.n_subjects,
            "n_raters_exact": self.n_raters_exact,
            "n_raters_lower_bound": self.n_raters_lower_bound,
        }

        if self.n_raters is not None:
            result["n_raters"] = self.n_raters
        if self.fill_rate is not None:
            result["fill_rate"] = round(self.fill_rate, 4)

        if self.is_balanced:
            result["raters_per_subject"] = self.raters_per_subject
        else:
            result["subject_counts"] = dict(self.subject_counts)

        return result
