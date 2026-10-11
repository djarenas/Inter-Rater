"""Injector that generates raters by sampling from the global rating pool."""

from __future__ import annotations

import numpy as np

from ..data_objects.ratings import Ratings
from ..contracts import EMPTY, INDEX_DTYPE
from .rater_injector import RaterInjector
from .configs import PooledInjectorConfig

_ERROR_PREFIX = "[PooledInjector] Error"


class PooledInjector(RaterInjector):
    """
    I inject raters who sample from the global pool of all ratings.

    At every subject, each injected rater draws from a single global
    distribution derived from the entire dataset. This preserves the
    empirical category prevalences but destroys all subject-level
    structure. The injected rater "knows" how often each category
    appears overall, but not which subjects got which answers.

    With ``smoothing > 0``, I blend the empirical distribution with a
    uniform distribution over categories. At smoothing=0.0 you get the
    pure empirical marginals; at smoothing=1.0 the rater draws
    uniformly regardless of the data.

    With ``exclude_raters``, I leave specific raters out of the pool
    before computing the distribution. Useful for leave-one-out
    sensitivity analysis: "what if rater 2 hadn't been there, and we
    replaced them with a chance-level rater?"

    I validate that:
    - Config is PooledInjectorConfig
    - exclude_raters indices are in bounds
    - At least one rater remains after exclusion
    - At least one non-EMPTY rating exists after exclusion

    Attributes
    ----------
    ratings : Ratings
        The original ratings data (from RaterInjector).
    config : PooledInjectorConfig
        Injection configuration including smoothing and exclude_raters.
    """

    # Attribute declarations
    ratings: Ratings
    config: PooledInjectorConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: PooledInjectorConfig
    ) -> None:
        """
        Initialize with ratings and pooled config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to inject into.
        config : PooledInjectorConfig
            Must include smoothing, exclude_raters, and n_new_raters.

        Raises
        ------
        TypeError
            If config is not PooledInjectorConfig.
        ValueError
            If exclude_raters indices are out of bounds, all raters
            are excluded, or no valid ratings remain after exclusion.
        """
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._validate_pooled_config()
        self._validate_exclude_raters_in_bounds()
        self._validate_raters_remain_after_exclusion()
        self._validate_pool_not_empty()

    # Public API ------------------------------------------------------------- #

    def generate(self) -> np.ndarray:
        """
        Generate rater columns by sampling from the global pool.

        I collect all non-EMPTY ratings (excluding any raters in
        ``exclude_raters``), then draw from that pool. With
        ``smoothing=0.0`` I resample directly from the pool values.
        With ``smoothing > 0`` I compute the empirical category
        distribution, blend it toward uniform, and draw from the
        blended categorical distribution.

        Returns
        -------
        np.ndarray
            Shape ``(n_subjects, n_new_raters)``, dtype INDEX_DTYPE.
        """
        n_subjects = self.ratings.category_index.shape[0]
        n_injected = self.config.n_new_raters
        smoothing = self.config.smoothing

        pool = self._build_pool()
        rng = np.random.default_rng(self.config.random_state)

        if smoothing == 0.0:
            # Fast path: resample directly from pool values
            result = rng.choice(
                pool, size=(n_subjects, n_injected), replace=True
            )
        else:
            # Probability path: blend empirical with uniform
            probs = self._build_blended_distribution(pool, smoothing)
            categories = np.arange(
                self.ratings.semantics.n_categories, dtype=INDEX_DTYPE
            )
            result = rng.choice(
                categories, size=(n_subjects, n_injected),
                replace=True, p=probs,
            )

        return result.astype(INDEX_DTYPE, copy=False)

    # Helpers ---------------------------------------------------------------- #

    def _build_pool(self) -> np.ndarray:
        """
        Build the flat pool of valid ratings, excluding specified raters.

        Returns
        -------
        np.ndarray
            1-D array of non-EMPTY category indices from the included
            raters.
        """
        category_index = self.ratings.category_index
        exclude = set(self.config.exclude_raters)

        if len(exclude) == 0:
            flat = category_index.ravel()
        else:
            n_raters = category_index.shape[1]
            included_cols = [
                c for c in range(n_raters) if c not in exclude
            ]
            flat = category_index[:, included_cols].ravel()

        return flat[flat != EMPTY]

    def _build_blended_distribution(
        self, pool: np.ndarray, smoothing: float
    ) -> np.ndarray:
        """
        Compute the blended empirical-uniform distribution.

        Parameters
        ----------
        pool : np.ndarray
            1-D array of valid category indices.
        smoothing : float
            Blend weight toward uniform. 0.0 = pure empirical,
            1.0 = pure uniform.

        Returns
        -------
        np.ndarray
            Probability vector of length n_categories, sums to 1.
        """
        n_categories = self.ratings.semantics.n_categories
        counts = np.bincount(pool, minlength=n_categories).astype(float)
        empirical = counts / counts.sum()
        uniform = np.full(n_categories, 1.0 / n_categories)
        return (1.0 - smoothing) * empirical + smoothing * uniform

    def _validate_pooled_config(self) -> None:
        """
        Validate that config is PooledInjectorConfig.

        Raises
        ------
        TypeError
            If config is not PooledInjectorConfig.
        """
        if not isinstance(self.config, PooledInjectorConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be PooledInjectorConfig, "
                f"got {type(self.config).__name__}"
            )

    def _validate_exclude_raters_in_bounds(self) -> None:
        """
        Validate that all exclude_raters indices are valid rater positions.

        Raises
        ------
        ValueError
            If any index in exclude_raters >= n_raters.
        """
        n_raters = self.ratings.n_raters
        for idx in self.config.exclude_raters:
            if idx >= n_raters:
                raise ValueError(
                    f"{_ERROR_PREFIX}: exclude_raters contains index "
                    f"{idx}, but ratings has only {n_raters} raters "
                    f"(valid indices: 0..{n_raters - 1})"
                )

    def _validate_raters_remain_after_exclusion(self) -> None:
        """
        Validate that at least one rater remains after exclusion.

        Raises
        ------
        ValueError
            If all raters are excluded.
        """
        n_remaining = (
            self.ratings.n_raters - len(self.config.exclude_raters)
        )
        if n_remaining <= 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: exclude_raters removes all "
                f"{self.ratings.n_raters} raters. At least one "
                f"rater must remain to build the pool."
            )

    def _validate_pool_not_empty(self) -> None:
        """
        Validate that after exclusion, at least one non-EMPTY rating exists.

        Raises
        ------
        ValueError
            If the pool contains no valid ratings.
        """
        pool = self._build_pool()
        if len(pool) == 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: no valid ratings remain after "
                f"excluding raters {self.config.exclude_raters}. "
                f"Cannot sample from an empty pool."
            )
