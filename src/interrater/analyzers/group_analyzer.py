"""Base analyzer for group coefficients."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

from .analyzer import Analyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
    from .configs import GroupCoefficientConfig

_ERROR_PREFIX = "[GroupAnalyzer] Error"


class GroupAnalyzer(Analyzer, ABC):
    """
    I am the abstract base for group coefficients (Fleiss, AC1, BP, Scott, Krippendorff).

    Group coefficients read a counts matrix (NijMatrix) that discards rater
    identity. I add validation specific to group analysis:
    - Config must be GroupCoefficientConfig
    - If config.balanced_data_only=True, ratings must have constant raters per subject
    - If config.balanced_data_only=False, allow unbalanced data

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    config : GroupCoefficientConfig
        Must be an instance of GroupCoefficientConfig or a subclass
        (NoWeightGroupConfig, WeightableGroupConfig, FleissKappaConfig, etc.)
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    # Checked at this level
    config: GroupCoefficientConfig


    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: GroupCoefficientConfig) -> None:
        """
        Initialize group analyzer: validate inputs and config type.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : GroupCoefficientConfig
            The analysis configuration for a group coefficient.

        Raises
        ------
        TypeError
            If config is not GroupCoefficientConfig.
        ValueError
            If config.balanced_data_only=True but data is unbalanced.
        """
        # Call parent __init__ (validates ratings and config, computes BalanceProfile)
        # Wrap to propagate errors with GroupAnalyzer context
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type and balance constraint
        self._validate_group_config()
        self._validate_balance_constraint()

        # Weights arrive with the semantics; refuse them if unusable here
        self._refuse_weights_if_unsupported()

    # Public API ------------------------------------------------------------- #

    # Resampling ------------------------------------------------------------- #

    @abstractmethod
    def _estimate_from_nij(self, nij_array: np.ndarray) -> float:
        """
        Compute my coefficient from a counts matrix, and nothing else.

        This is the one step that differs between group coefficients. My
        subclasses implement this and nothing else about resampling: the
        bootstrap and the jackknife are written once in terms of it.

        Parameters
        ----------
        nij_array : np.ndarray
            Shape (n_subjects, n_categories). Entry (i, j) is how many raters
            assigned category j to subject i.

        Returns
        -------
        float
            The coefficient, or NaN when this matrix cannot define it.
        """
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: _estimate_from_nij() must be implemented "
            f"by subclass"
        )

    def _estimate_from_array(self, category_index: np.ndarray) -> float:
        """
        Satisfy Analyzer's contract by counting, then deferring to nij.

        A group coefficient reads a counts matrix, so counting is the one step
        in front of it. Writing that step here means the five group analyzers
        implement one method rather than two.
        """
        from ..array_utils import category_counts

        return self._estimate_from_nij(
            category_counts(category_index, self._n_categories())
        )

    # Private helpers -------------------------------------------------------- #

    def _nij_array(self) -> np.ndarray:
        """
        The counts matrix of my ratings, for ``_result`` and
        ``_parametric_null_interval``.
        """
        from ..array_utils import category_counts

        return category_counts(self.ratings.category_index, self._n_categories())

    def _agreement_profile(self):
        """
        My agreement profile, from my ratings.
        """
        from ..intermediates.data_profiles import AgreementProfile

        return AgreementProfile.build_from_ratings(self.ratings)

    def _n_categories(self) -> int:
        """
        K, from the semantics.

        Never inferred from the values: a matrix whose largest observed
        category is 2 may still have been declared over five.
        """
        return self.ratings.semantics.n_categories

    def _validate_group_config(self) -> None:
        """
        Validate that self.config is a GroupCoefficientConfig or subclass.

        Raises
        ------
        TypeError
            If config is not GroupCoefficientConfig.
        """
        from .configs import GroupCoefficientConfig

        if not isinstance(self.config, GroupCoefficientConfig):
            expected_configs = GroupCoefficientConfig.__subclasses__()
            expected_names = ", ".join([c.__name__ for c in expected_configs])
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be GroupCoefficientConfig or subclass. "
                f"Expected one of: {expected_names}. "
                f"Got {type(self.config).__name__}"
            )

    def _validate_balance_constraint(self) -> None:
        """
        If config.balanced_data_only=True, validate data is balanced.

        Raises
        ------
        ValueError
            If balanced_data_only=True but data is unbalanced.
        """
        if self.config.balanced_data_only and not self.balance_profile.is_balanced:
            raise ValueError(
                f"{_ERROR_PREFIX}: config.balanced_data_only=True but data is unbalanced. "
                f"Subjects have varying numbers of raters: {self.balance_profile.subject_counts}. "
                f"Either set balanced_data_only=False to accept unbalanced data, "
                f"or provide balanced data (constant raters per subject)."
            )
