"""Injector that generates raters from a categorical distribution."""

from __future__ import annotations

import numpy as np

from ..data_objects.ratings import Ratings
from ..contracts import INDEX_DTYPE
from .rater_injector import RaterInjector
from .configs import DistributionInjectorConfig

_ERROR_PREFIX = "[DistributionInjector] Error"


class DistributionInjector(RaterInjector):
    """
    I inject raters drawn from a categorical distribution.

    Each injected rater's ratings are independent draws from the
    prevalence vector. The rater has no relationship to the real
    data -- they rate based on base rates alone.

    This is the null model injector. "What does agreement look like
    if we add someone who is just guessing at these base rates?"

    When ``prevalences`` is None, I draw uniformly across all
    categories. When provided, each draw is categorical (multinomial
    with n=1). For binary categories with prevalences (p, 1-p), each
    draw is Bernoulli.

    I validate that:
    - Config is DistributionInjectorConfig
    - If prevalences is provided, its length matches n_categories

    Attributes
    ----------
    ratings : Ratings
        The original ratings data (from RaterInjector).
    config : DistributionInjectorConfig
        Injection configuration including prevalences.
    """

    # Attribute declarations
    ratings: Ratings
    config: DistributionInjectorConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: DistributionInjectorConfig
    ) -> None:
        """
        Initialize with ratings and distribution config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to inject into.
        config : DistributionInjectorConfig
            Must include prevalences (or None for uniform) and
            n_new_raters.

        Raises
        ------
        TypeError
            If config is not DistributionInjectorConfig.
        ValueError
            If prevalences length does not match n_categories.
        """
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._validate_distribution_config()
        self._validate_prevalences_length()

    # Public API ------------------------------------------------------------- #

    def generate(self) -> np.ndarray:
        """
        Generate rater columns from the prevalence distribution.

        When prevalences is None, I draw uniformly. Otherwise I draw
        from the specified categorical distribution. All draws are
        independent -- no subject-level structure.

        Returns
        -------
        np.ndarray
            Shape ``(n_subjects, n_new_raters)``, dtype INDEX_DTYPE.
        """
        n_subjects = self.ratings.category_index.shape[0]
        n_injected = self.config.n_new_raters
        n_categories = self.ratings.semantics.n_categories

        probs = self._build_prevalences(n_categories)
        rng = np.random.default_rng(self.config.random_state)

        result = rng.choice(
            n_categories,
            size=(n_subjects, n_injected),
            p=probs,
        )

        return result.astype(INDEX_DTYPE, copy=False)

    # Helpers ---------------------------------------------------------------- #

    def _build_prevalences(self, n_categories: int) -> np.ndarray:
        """
        Build the probability vector for sampling.

        Parameters
        ----------
        n_categories : int
            Number of categories in the rating scheme.

        Returns
        -------
        np.ndarray
            Probability vector of length n_categories, sums to 1.
        """
        if self.config.prevalences is None:
            return np.full(n_categories, 1.0 / n_categories)
        return np.array(self.config.prevalences, dtype=float)

    def _validate_distribution_config(self) -> None:
        """
        Validate that config is DistributionInjectorConfig.

        Raises
        ------
        TypeError
            If config is not DistributionInjectorConfig.
        """
        if not isinstance(self.config, DistributionInjectorConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be "
                f"DistributionInjectorConfig, "
                f"got {type(self.config).__name__}"
            )

    def _validate_prevalences_length(self) -> None:
        """
        Validate that prevalences length matches n_categories.

        Only checked when prevalences is not None. None means uniform,
        which is built at generate time from n_categories.

        Raises
        ------
        ValueError
            If len(prevalences) != n_categories.
        """
        if self.config.prevalences is None:
            return

        n_categories = self.ratings.semantics.n_categories
        n_prevalences = len(self.config.prevalences)
        if n_prevalences != n_categories:
            raise ValueError(
                f"{_ERROR_PREFIX}: prevalences has {n_prevalences} "
                f"entries, but ratings has {n_categories} categories. "
                f"Provide one probability per category, or None "
                f"for uniform."
            )
