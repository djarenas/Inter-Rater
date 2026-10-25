"""Injector that generates deterministic raters based on prevalence rank."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from ..data_objects.ratings import Ratings
from ..contracts import EMPTY, INDEX_DTYPE
from .rater_injector import RaterInjector
from .configs import ExtremalInjectorConfig

_ERROR_PREFIX = "[ExtremalInjector] Error"


class ExtremalInjector(RaterInjector):
    """
    I inject raters who pick a category based on prevalence rank.

    At each subject, I look at how the real raters voted and pick a
    category deterministically -- with ties broken randomly so that
    multiple injected raters can disagree at tied subjects.

    Four strategies:

    - ``"majority"`` (Conformist): picks the most common category.
      Always goes with the crowd.
    - ``"minority"`` (Dissenter): picks the least common category
      that was actually used. The unpopular option someone did choose.
    - ``"median"``: picks the category with middle-ranked prevalence
      among those that appeared. Ranks by count, not category index,
      so it works for nominal (unordered) data.
    - ``"contrarian"``: picks a category that no real rater used at
      this subject. The true contrarian -- goes where nobody went.
      Falls back to minority when all categories were used.

    I validate that:
    - Config is ExtremalInjectorConfig

    Attributes
    ----------
    ratings : Ratings
        The original ratings data (from RaterInjector).
    config : ExtremalInjectorConfig
        Injection configuration including strategy.
    """

    # Attribute declarations
    ratings: Ratings
    config: ExtremalInjectorConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: ExtremalInjectorConfig, *, step_position: int = 0
    ) -> None:
        """
        Initialize with ratings and extremal config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to inject into.
        config : ExtremalInjectorConfig
            Must include strategy and n_new_raters.

        Raises
        ------
        TypeError
            If config is not ExtremalInjectorConfig.
        """
        try:
            super().__init__(ratings, config, step_position=step_position)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._validate_extremal_config()

    # Public API ------------------------------------------------------------- #

    def generate(self) -> np.ndarray:
        """
        Generate rater columns by picking categories at prevalence extremes.

        At each subject, I count how many raters chose each category,
        find the candidates for the configured strategy, and let each
        injected rater independently pick from the candidate set
        (relevant only when there are ties or multiple valid choices).

        Returns
        -------
        np.ndarray
            Shape ``(n_subjects, n_new_raters)``, dtype INDEX_DTYPE.
        """
        category_index = self.ratings.category_index
        n_subjects = category_index.shape[0]
        n_injected = self.config.n_new_raters
        n_categories = self.ratings.semantics.n_categories
        strategy = self.config.strategy

        rng = self._rng()
        result = np.empty((n_subjects, n_injected), dtype=INDEX_DTYPE)

        # Dispatch to the right candidate finder
        candidate_fn = self._get_candidate_fn(strategy)

        for i in range(n_subjects):
            row = category_index[i, :]
            valid = row[row != EMPTY]

            # Ratings guarantees at least one non-empty per subject,
            # but guard defensively.
            if len(valid) == 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subject {i} has no valid "
                    f"ratings. Cannot determine extremal category."
                )

            counts = np.bincount(valid, minlength=n_categories)
            candidates = candidate_fn(counts, n_categories)
            result[i, :] = rng.choice(candidates, size=n_injected,
                                      replace=True)

        return result

    # Helpers ---------------------------------------------------------------- #

    def _get_candidate_fn(
        self, strategy: str
    ) -> Callable[[np.ndarray, int], np.ndarray]:
        """
        Return the candidate-finding function for a strategy.

        Parameters
        ----------
        strategy : str
            One of "majority", "minority", "median", "contrarian".

        Returns
        -------
        callable
            A function ``(counts, n_categories) -> np.ndarray`` that
            returns the candidate category indices for one subject.
        """
        dispatch = {
            "majority": self._candidates_majority,
            "minority": self._candidates_minority,
            "median": self._candidates_median,
            "contrarian": self._candidates_contrarian,
        }
        return dispatch[strategy]

    @staticmethod
    def _candidates_majority(
        counts: np.ndarray, n_categories: int
    ) -> np.ndarray:
        """
        Find categories tied for the highest count.

        Parameters
        ----------
        counts : np.ndarray
            Length n_categories, count per category at this subject.
        n_categories : int
            Total number of categories (unused here, kept for
            uniform signature).

        Returns
        -------
        np.ndarray
            Category indices tied for the maximum count.
        """
        max_count = counts[counts > 0].max()
        return np.where(counts == max_count)[0].astype(INDEX_DTYPE)

    @staticmethod
    def _candidates_minority(
        counts: np.ndarray, n_categories: int
    ) -> np.ndarray:
        """
        Find categories tied for the lowest nonzero count.

        Parameters
        ----------
        counts : np.ndarray
            Length n_categories, count per category at this subject.
        n_categories : int
            Total number of categories (unused here, kept for
            uniform signature).

        Returns
        -------
        np.ndarray
            Category indices tied for the minimum nonzero count.
        """
        used_mask = counts > 0
        min_count = counts[used_mask].min()
        return np.where(used_mask & (counts == min_count))[0].astype(
            INDEX_DTYPE
        )

    @staticmethod
    def _candidates_median(
        counts: np.ndarray, n_categories: int
    ) -> np.ndarray:
        """
        Find categories at the median prevalence rank.

        I rank used categories by count, then pick those at the
        middle position. With an even number of used categories, I
        take the lower-middle (floor). Multiple categories can share
        the median count, so all tied at that rank are candidates.

        Parameters
        ----------
        counts : np.ndarray
            Length n_categories, count per category at this subject.
        n_categories : int
            Total number of categories (unused here, kept for
            uniform signature).

        Returns
        -------
        np.ndarray
            Category indices whose count equals the median count.
        """
        used_mask = counts > 0
        used_counts = counts[used_mask]

        # Sort the distinct used counts to find the median count value
        sorted_counts = np.sort(used_counts)
        median_idx = (len(sorted_counts) - 1) // 2
        median_count = sorted_counts[median_idx]

        return np.where(
            used_mask & (counts == median_count)
        )[0].astype(INDEX_DTYPE)

    @staticmethod
    def _candidates_contrarian(
        counts: np.ndarray, n_categories: int
    ) -> np.ndarray:
        """
        Find categories that no rater used at this subject.

        If all categories were used (no unused available), I fall
        back to minority behavior -- the least common used category.

        Parameters
        ----------
        counts : np.ndarray
            Length n_categories, count per category at this subject.
        n_categories : int
            Total number of categories.

        Returns
        -------
        np.ndarray
            Unused category indices, or minority candidates if none.
        """
        unused = np.where(counts == 0)[0].astype(INDEX_DTYPE)

        if len(unused) > 0:
            return unused

        # Fallback: minority behavior
        used_mask = counts > 0
        min_count = counts[used_mask].min()
        return np.where(
            used_mask & (counts == min_count)
        )[0].astype(INDEX_DTYPE)

    def _validate_extremal_config(self) -> None:
        """
        Validate that config is ExtremalInjectorConfig.

        Raises
        ------
        TypeError
            If config is not ExtremalInjectorConfig.
        """
        if not isinstance(self.config, ExtremalInjectorConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be "
                f"ExtremalInjectorConfig, "
                f"got {type(self.config).__name__}"
            )
