"""Base analyzer for group coefficients."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from .analyzer import Analyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates import AnalysisResult
    from ..intermediates.data_profiles import BalanceProfile
    from ..intermediates.estimator_parts import NijMatrix
    from .configs import GroupCoefficientConfig

_ERROR_PREFIX = "[GroupAnalyzer] Error"


class GroupAnalyzer(Analyzer, ABC):
    """
    I am the abstract base for group coefficients (Fleiss, AC1, BP, Scott, Krippendorff).

    Group coefficients read a counts matrix (NijMatrix) that discards rater
    identity. I add validation specific to group analysis and support two
    initialization paths:
    - Ratings-based: __init__(ratings, config) computes nij internally
    - NijMatrix-based: build_from_nij_matrix(nij, config) uses nij directly

    I add validation specific to group analysis:
    - Config must be GroupCoefficientConfig
    - If config.balanced_data_only=True, ratings must have constant raters per subject
    - If config.balanced_data_only=False, allow unbalanced data (use approximation)

    Attributes
    ----------
    ratings : Ratings | None
        The ratings data (from Analyzer, can be None if initialized from nij_matrix).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    config : GroupCoefficientConfig
        Must be an instance of GroupCoefficientConfig or a subclass
        (NoWeightGroupConfig, WeightableGroupConfig, FleissKappaConfig, etc.)
    nij_matrix : NijMatrix | None
        The counts matrix (populated at this level if built from nij_matrix).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings | None
    balance_profile: BalanceProfile
    # Checked at this level
    config: GroupCoefficientConfig
    # Populated at this level
    nij_matrix: NijMatrix | None


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

        # Set nij_matrix to None (Ratings-based path)
        self.nij_matrix = None

        # Validate config type and balance constraint
        self._validate_group_config()
        self._validate_balance_constraint()

    @classmethod
    def build_from_nij_matrix(
        cls, nij: NijMatrix, config: GroupCoefficientConfig
    ) -> GroupAnalyzer:
        """
        Build group analyzer from NijMatrix (counts matrix) directly.

        This alternate constructor is useful for analyzing published data that only
        provides the counts matrix (nij) without the original ratings.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).
        config : GroupCoefficientConfig
            The analysis configuration for a group coefficient.

        Returns
        -------
        GroupAnalyzer
            An initialized analyzer ready to call estimate_from_nij_matrix().

        Raises
        ------
        TypeError
            If config is not GroupCoefficientConfig.
        ValueError
            If config.balanced_data_only=True but nij indicates unbalanced data.
        """
        from ..intermediates.data_profiles import BalanceProfile

        # Create instance via __new__ (bypass normal __init__)
        instance = cls.__new__(cls)

        # Call parent __init__ with None for ratings
        # Wrap to propagate errors with GroupAnalyzer context
        try:
            Analyzer.__init__(instance, None, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Store nij_matrix
        instance.nij_matrix = nij

        # Compute BalanceProfile from nij_matrix
        instance.balance_profile = BalanceProfile.from_nij_matrix(nij)

        # Validate config type and balance constraint
        instance._validate_group_config()
        instance._validate_balance_constraint()

        return instance

    # Public API ------------------------------------------------------------- #

    @abstractmethod
    def estimate(self) -> AnalysisResult:
        """
        Compute the group coefficient from ratings (Ratings-based path).

        Must be implemented by subclass. Called when analyzer was initialized
        via __init__(ratings, config).

        Returns
        -------
        AnalysisResult
            The coefficient with point estimate, CI, and metadata.

        Raises
        ------
        ValueError
            If self.ratings is None (should use estimate_from_nij_matrix instead).
        NotImplementedError
            If not implemented by subclass.
        """
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: estimate() must be implemented by subclass, "
            f"got {type(self).__name__}"
        )

    @abstractmethod
    def estimate_from_nij_matrix(self) -> AnalysisResult:
        """
        Compute the group coefficient from nij_matrix (NijMatrix-based path).

        Must be implemented by subclass. Called when analyzer was initialized
        via build_from_nij_matrix(nij, config).

        Returns
        -------
        AnalysisResult
            The coefficient with point estimate, CI, and metadata.

        Raises
        ------
        ValueError
            If self.nij_matrix is None (should use estimate instead).
        NotImplementedError
            If not implemented by subclass.
        """
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: estimate_from_nij_matrix() must be implemented by subclass, "
            f"got {type(self).__name__}"
        )

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
                f"Either set balanced_data_only=False to use the approximation, "
                f"or provide balanced data (constant raters per subject)."
            )
