"""Injector that generates raters by copying real answers within each subject."""

from __future__ import annotations

import numpy as np

from ..data_objects.ratings import Ratings
from ..contracts import EMPTY, INDEX_DTYPE
from .rater_injector import RaterInjector
from .configs import WithinSubjectInjectorConfig

_ERROR_PREFIX = "[WithinSubjectInjector] Error"


class WithinSubjectInjector(RaterInjector):
    """
    I inject raters who copy real answers at each subject.

    At every subject, each injected rater randomly selects one of the
    real raters who actually rated that subject (i.e. not EMPTY) and
    copies their answer. This is the most structure-preserving strategy:
    where the panel agreed, I agree; where they disagreed, I pick a side.

    With ``noise_rate > 0``, at each subject there is a probability of
    instead picking a category that no real rater used for that subject.
    If all categories were already used, I fall back to copying normally.

    I simulate "a rater who tried" -- someone reading from the same
    signal as the panel, just noisily. Statistically, this is a
    within-subject bootstrap of rater identity, with optional injected
    noise.

    I validate that:
    - Config is WithinSubjectInjectorConfig

    Attributes
    ----------
    ratings : Ratings
        The original ratings data (from RaterInjector).
    config : WithinSubjectInjectorConfig
        Injection configuration including noise_rate.
    """

    # Attribute declarations
    
    # Inherited from parent RaterInjector
    ratings: Ratings
    config: WithinSubjectInjectorConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: WithinSubjectInjectorConfig
    ) -> None:
        """
        Initialize with ratings and within-subject config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to inject into.
        config : WithinSubjectInjectorConfig
            Must include noise_rate and n_new_raters.

        Raises
        ------
        TypeError
            If config is not WithinSubjectInjectorConfig.
        """
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._validate_within_subject_config()

    # Public API ------------------------------------------------------------- #

    def generate(self) -> np.ndarray:
        """
        Generate rater columns by copying real ratings within each subject.

        For each subject, I find the raters who actually rated it
        (non-EMPTY) and randomly select one of their answers for each
        injected rater. With noise_rate > 0, some selections are
        replaced by an unused category.

        Returns
        -------
        np.ndarray
            Shape ``(n_subjects, n_new_raters)``, dtype INDEX_DTYPE.
        """
        category_index = self.ratings.category_index
        n_subjects = category_index.shape[0]
        n_injected = self.config.n_new_raters
        n_categories = self.ratings.semantics.n_categories
        noise_rate = self.config.noise_rate

        rng = np.random.default_rng(self.config.random_state)

        result = np.empty((n_subjects, n_injected), dtype=INDEX_DTYPE)

        for i in range(n_subjects):
            row = category_index[i, :]
            valid_ratings = row[row != EMPTY]

            # Ratings guarantees at least one non-empty per subject
            # (validate_no_unrated_subjects), but guard defensively.
            if len(valid_ratings) == 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subject {i} has no valid "
                    f"ratings. Cannot sample from an empty set."
                )

            if noise_rate == 0.0:
                # Fast path: pure copy, no noise check needed
                result[i, :] = rng.choice(
                    valid_ratings, size=n_injected, replace=True
                )
            else:
                self._generate_subject_with_noise(
                    result, i, valid_ratings, n_injected,
                    n_categories, noise_rate, rng,
                )

        return result

    # Helpers ---------------------------------------------------------------- #

    def _generate_subject_with_noise(
        self,
        result: np.ndarray,
        subject_idx: int,
        valid_ratings: np.ndarray,
        n_injected: int,
        n_categories: int,
        noise_rate: float,
        rng: np.random.Generator,
    ) -> None:
        """
        Generate ratings for one subject with noise injection.

        I decide which injected raters get noise vs copy, then fill
        the result array in place.

        Parameters
        ----------
        result : np.ndarray
            The output array to fill at row ``subject_idx``.
        subject_idx : int
            Which subject row to fill.
        valid_ratings : np.ndarray
            The non-EMPTY ratings at this subject.
        n_injected : int
            Number of injected raters.
        n_categories : int
            Total number of categories in the rating scheme.
        noise_rate : float
            Probability of noise path per rater.
        rng : np.random.Generator
            The seeded random number generator.
        """
        noise_flags = rng.random(n_injected) < noise_rate
        n_noise = int(noise_flags.sum())
        n_copy = n_injected - n_noise

        # Copy path
        if n_copy > 0:
            result[subject_idx, ~noise_flags] = rng.choice(
                valid_ratings, size=n_copy, replace=True
            )

        # Noise path
        if n_noise > 0:
            used_set = set(valid_ratings)
            unused = np.array(
                [c for c in range(n_categories) if c not in used_set],
                dtype=INDEX_DTYPE,
            )

            if len(unused) > 0:
                result[subject_idx, noise_flags] = rng.choice(
                    unused, size=n_noise, replace=True
                )
            else:
                # All categories used -- fall back to copy
                result[subject_idx, noise_flags] = rng.choice(
                    valid_ratings, size=n_noise, replace=True
                )

    def _validate_within_subject_config(self) -> None:
        """
        Validate that config is WithinSubjectInjectorConfig.

        Raises
        ------
        TypeError
            If config is not WithinSubjectInjectorConfig.
        """
        if not isinstance(self.config, WithinSubjectInjectorConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be "
                f"WithinSubjectInjectorConfig, "
                f"got {type(self.config).__name__}"
            )
