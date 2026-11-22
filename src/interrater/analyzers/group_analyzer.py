"""Base analyzer for group coefficients."""

from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

from .analyzer import Analyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
    from .configs import GroupCoefficientConfig

_ERROR_PREFIX = "[GroupAnalyzer] Error"


class GroupAnalyzer(Analyzer, ABC):
    """
    I am the abstract base for analyzers of group coefficients (Fleiss, AC1,
    BP, Scott, Krippendorff). These coefficients use how many raters chose
    each category for each subject, not which rater chose it.

    My inheritors' coefficient functions take an EstimatorInput and count it
    into an NijMatrix (``data.nij``), which discards rater identity.

    At this inheritance level, I check:
    - Config must be GroupCoefficientConfig
    - If config.balanced_data_only=True, ratings must have constant raters per subject
    - If config.balanced_data_only=False, allow unbalanced data
    - Weights in the semantics are refused if my coefficient has no weighted form

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    weights : Weights | None
        From the semantics (from Analyzer).
    estimator_input : EstimatorInput
        What my coefficient is computed from (from Analyzer).
    config : GroupCoefficientConfig
        Must be an instance of GroupCoefficientConfig or a subclass
        (FleissKappaConfig, GwetAC1Config, etc.)
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
        Run the parent's checks, then check the config type, the balance
        constraint and the weights.

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
            If config.balanced_data_only=True but data is unbalanced, or the
            semantics declares weights and my coefficient has no weighted form.
        """
        # Call parent __init__ (validates ratings and config, applies the
        # injection, builds the profiles, the weights and the EstimatorInput)
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

    # Private helpers -------------------------------------------------------- #

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
