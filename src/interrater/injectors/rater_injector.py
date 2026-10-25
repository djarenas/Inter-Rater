"""Base class for rater injection strategies."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from ..data_objects.ratings import Ratings
from ..semantics.ratings_semantics import RatingsSemantics
from .configs import InjectorConfig

_ERROR_PREFIX = "[RaterInjector] Error"


class RaterInjector(ABC):
    """
    I am the base class for rater injection strategies.

    I take a Ratings object and an InjectorConfig, generate one or more
    synthetic rater columns, and return a new Ratings with those columns
    appended. The original Ratings is never modified -- I build a new one.

    My subclasses implement ``generate()``, which produces the raw columns.
    I handle the mechanical work of appending them to the Ratings and
    building a valid new Ratings object with updated semantics.

    I validate that:
    - ratings is a Ratings instance (not None)
    - config is an InjectorConfig instance

    Attributes
    ----------
    ratings : Ratings
        The original ratings data. Never modified.
    config : InjectorConfig
        The injection configuration.
    """

    # Attribute declarations
    ratings: Ratings
    config: InjectorConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: InjectorConfig, *, step_position: int = 0
    ) -> None:
        """
        Initialize the injector with ratings and config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to inject into. Must not be None.
        config : InjectorConfig
            The injection configuration.
        step_position : int, default 0
            Where this step stands in an injection pipeline, counting from
            0. Mixed into the seed (see ``_rng``) so that chained steps
            with the same ``random_state`` do not produce the same raters.

        Raises
        ------
        TypeError
            If ratings is not Ratings or config is not InjectorConfig, or
            step_position is not an int.
        ValueError
            If ratings is None or step_position is negative.
        """
        self._validate_base_inputs(ratings, config)
        if isinstance(step_position, bool) or not isinstance(step_position, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: step_position must be int, "
                f"got {type(step_position).__name__}"
            )
        if step_position < 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: step_position must be >= 0, got {step_position}"
            )
        self.ratings = ratings
        self.config = config
        self.step_position = step_position

    def _rng(self) -> np.random.Generator:
        """
        The random generator for this step, seeded from its seed and position.

        Seeding with the pair ``[random_state, step_position]`` lets numpy mix
        both into one seed, so two steps given the same ``random_state`` get
        unrelated streams, and different pairs cannot collide the way
        ``seed + position`` can. One rule for every step, the first included.
        ``random_state=None`` means fresh randomness, as before.
        """
        if self.config.random_state is None:
            return np.random.default_rng()
        return np.random.default_rng([self.config.random_state, self.step_position])

    # Public API ------------------------------------------------------------- #

    @abstractmethod
    def generate(self) -> np.ndarray:
        """
        Generate synthetic rater columns.

        Returns
        -------
        np.ndarray
            Shape ``(n_subjects,)`` for a single rater, or
            ``(n_subjects, k)`` for k raters. Values are integer
            category indices matching the Ratings' category scheme.

        Notes
        -----
        Subclasses implement the generation strategy. The base class
        calls this from ``inject()`` and handles appending and Ratings
        construction.
        """
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: {type(self).__name__} must implement "
            f"generate(). This abstract method should not be called "
            f"directly."
        )

    def inject(self) -> Ratings:
        """
        Generate synthetic raters and return a new Ratings with them appended.

        Calls ``generate()`` to produce columns, validates them, and
        delegates to ``_build_new_ratings`` for construction.

        Returns
        -------
        Ratings
            New Ratings with the original raters plus the injected ones.
            The original Ratings is not modified.
        """
        new_columns = self.generate()

        # Normalise to 2D
        if new_columns.ndim == 1:
            new_columns = new_columns.reshape(-1, 1)

        self._validate_generated_columns(new_columns)

        return self._build_new_ratings(new_columns)

    # Helpers ---------------------------------------------------------------- #

    def _validate_base_inputs(
        self, ratings: Ratings, config: InjectorConfig
    ) -> None:
        """
        Validate constructor inputs.

        Raises
        ------
        TypeError
            If ratings is not Ratings or config is not InjectorConfig.
        ValueError
            If ratings or config is None.
        """
        if ratings is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: ratings cannot be None. "
                f"Injection requires data to inject into."
            )

        if not isinstance(ratings, Ratings):
            raise TypeError(
                f"{_ERROR_PREFIX}: ratings must be Ratings, "
                f"got {type(ratings).__name__}"
            )

        if config is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: config cannot be None."
            )

        if not isinstance(config, InjectorConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be InjectorConfig "
                f"(or a subclass), got {type(config).__name__}"
            )

    def _validate_generated_columns(self, columns: np.ndarray) -> None:
        """
        Validate that generated columns are compatible with the Ratings.

        Parameters
        ----------
        columns : np.ndarray
            Shape ``(n_subjects, k)``. Already normalised to 2D.

        Raises
        ------
        ValueError
            If shape or values are incompatible.
        """
        n_subjects = self.ratings.category_index.shape[0]
        n_categories = self.ratings.semantics.n_categories

        if columns.shape[0] != n_subjects:
            raise ValueError(
                f"{_ERROR_PREFIX}: generated columns have "
                f"{columns.shape[0]} rows, expected {n_subjects} "
                f"(one per subject)."
            )

        # Check values are valid category indices
        min_val = columns.min()
        max_val = columns.max()

        if min_val < 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: generated columns contain negative "
                f"values (min={min_val}). Values must be category "
                f"indices in [0, {n_categories - 1}]."
            )

        if max_val >= n_categories:
            raise ValueError(
                f"{_ERROR_PREFIX}: generated columns contain values "
                f">= n_categories={n_categories} (max={max_val}). "
                f"Values must be category indices in "
                f"[0, {n_categories - 1}]."
            )

    def _build_rater_names(self, n_injected: int) -> tuple[str, ...]:
        """
        Build the full rater name tuple: original names + injected names.

        Injected raters are named ``injected_rater0``, ``injected_rater1``,
        etc. When chaining injections, I scan the existing rater names
        for the highest ``injected_raterN`` index and start numbering
        from N+1, so names never collide.

        Parameters
        ----------
        n_injected : int
            Number of injected raters.

        Returns
        -------
        tuple of str
            Original rater names followed by injected rater names.
        """
        original_names = self.ratings.semantics.raters

        # Find the highest existing injected_raterN index to avoid
        # collisions when chaining injections
        start_index = 0
        prefix = "injected_rater"
        for name in original_names:
            if name.startswith(prefix):
                suffix = name[len(prefix):]
                if suffix.isdigit():
                    start_index = max(start_index, int(suffix) + 1)

        injected_names = tuple(
            f"injected_rater{start_index + i}"
            for i in range(n_injected)
        )
        return original_names + injected_names

    def _build_new_semantics(self, n_injected: int) -> RatingsSemantics:
        """
        Build a new RatingsSemantics with injected rater names appended.

        Parameters
        ----------
        n_injected : int
            Number of injected raters.

        Returns
        -------
        RatingsSemantics
            New semantics with original + injected rater names.
        """
        new_rater_names = self._build_rater_names(n_injected)
        return RatingsSemantics(
            raters=list(new_rater_names),
            categories=self.ratings.semantics.categories,
            n_subjects=None,  # drop expectation -- we know the count
            unique_subjects=self.ratings.semantics.unique_subjects,
        )

    def _build_new_ratings(self, new_columns: np.ndarray) -> Ratings:
        """
        Build a new Ratings with the injected columns appended.

        Parameters
        ----------
        new_columns : np.ndarray
            Shape ``(n_subjects, k)``. Already validated by
            ``_validate_generated_columns``.

        Returns
        -------
        Ratings
            New Ratings with original + injected rater columns.
        """
        n_injected = new_columns.shape[1]

        new_category_index = np.column_stack([
            self.ratings.category_index, new_columns
        ])

        new_semantics = self._build_new_semantics(n_injected)

        return Ratings(
            category_index=new_category_index,
            semantics=new_semantics,
            subject_ids=self.ratings.subject_ids,
        )
