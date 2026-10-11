"""Injector that generates raters by bootstrapping existing raters' columns."""

from __future__ import annotations

import numpy as np

from ..data_objects.ratings import Ratings
from ..contracts import EMPTY, INDEX_DTYPE
from .rater_injector import RaterInjector
from .configs import RaterBootstrapInjectorConfig

_ERROR_PREFIX = "[RaterBootstrapInjector] Error"


class RaterBootstrapInjector(RaterInjector):
    """
    I inject raters cloned from one or more existing raters' columns.

    I pool the non-EMPTY ratings from the specified source raters and
    resample with replacement across all subjects. The injected rater
    has the same marginal distribution as the source subset -- their
    habits are cloned -- but no subject-specific information. Their
    judgments are scrambled.

    With a single source rater, this clones one person's biases. With
    all raters as sources, this collapses into a PooledInjector
    (without smoothing or exclusion knobs).

    "What if we added someone with rater k's biases but who can't
    actually distinguish between subjects?"

    I validate that:
    - Config is RaterBootstrapInjectorConfig
    - source_raters indices are in bounds
    - At least one non-EMPTY rating exists in the source pool

    Attributes
    ----------
    ratings : Ratings
        The original ratings data (from RaterInjector).
    config : RaterBootstrapInjectorConfig
        Injection configuration including source_raters.
    """

    # Attribute declarations
    ratings: Ratings
    config: RaterBootstrapInjectorConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: RaterBootstrapInjectorConfig
    ) -> None:
        """
        Initialize with ratings and rater bootstrap config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to inject into.
        config : RaterBootstrapInjectorConfig
            Must include source_raters and n_new_raters.

        Raises
        ------
        TypeError
            If config is not RaterBootstrapInjectorConfig.
        ValueError
            If source_raters indices are out of bounds or the source
            pool contains no valid ratings.
        """
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._validate_rater_bootstrap_config()
        self._validate_source_raters_in_bounds()
        self._validate_source_pool_not_empty()

    # Public API ------------------------------------------------------------- #

    def generate(self) -> np.ndarray:
        """
        Generate rater columns by bootstrapping the source raters.

        I pool the non-EMPTY ratings from all source raters and
        resample with replacement to fill ``(n_subjects, n_new_raters)``
        positions.

        Returns
        -------
        np.ndarray
            Shape ``(n_subjects, n_new_raters)``, dtype INDEX_DTYPE.
        """
        n_subjects = self.ratings.category_index.shape[0]
        n_injected = self.config.n_new_raters

        pool = self._build_source_pool()
        rng = np.random.default_rng(self.config.random_state)

        result = rng.choice(
            pool, size=(n_subjects, n_injected), replace=True
        )

        return result.astype(INDEX_DTYPE, copy=False)

    # Helpers ---------------------------------------------------------------- #

    def _build_source_pool(self) -> np.ndarray:
        """
        Build the flat pool of valid ratings from the source raters.

        Returns
        -------
        np.ndarray
            1-D array of non-EMPTY category indices from the source
            raters' columns.
        """
        category_index = self.ratings.category_index
        source_cols = list(self.config.source_raters)
        flat = category_index[:, source_cols].ravel()
        return flat[flat != EMPTY]

    def _validate_rater_bootstrap_config(self) -> None:
        """
        Validate that config is RaterBootstrapInjectorConfig.

        Raises
        ------
        TypeError
            If config is not RaterBootstrapInjectorConfig.
        """
        if not isinstance(self.config, RaterBootstrapInjectorConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be "
                f"RaterBootstrapInjectorConfig, "
                f"got {type(self.config).__name__}"
            )

    def _validate_source_raters_in_bounds(self) -> None:
        """
        Validate that all source_raters indices are valid rater positions.

        Raises
        ------
        ValueError
            If any index in source_raters >= n_raters.
        """
        n_raters = self.ratings.n_raters
        for idx in self.config.source_raters:
            if idx >= n_raters:
                raise ValueError(
                    f"{_ERROR_PREFIX}: source_raters contains index "
                    f"{idx}, but ratings has only {n_raters} raters "
                    f"(valid indices: 0..{n_raters - 1})"
                )

    def _validate_source_pool_not_empty(self) -> None:
        """
        Validate that the source raters have at least one non-EMPTY rating.

        Raises
        ------
        ValueError
            If the source pool contains no valid ratings.
        """
        pool = self._build_source_pool()
        if len(pool) == 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: source raters "
                f"{self.config.source_raters} have no valid ratings. "
                f"Cannot bootstrap from an empty pool."
            )
