"""Base analyzer for weightable group coefficients."""

from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

from .group_analyzer import GroupAnalyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.estimator_parts import Weights, NijMatrix
    from ..intermediates.data_profiles import BalanceProfile
    from .configs import WeightableGroupConfig

_ERROR_PREFIX = "[WeightableGroupAnalyzer] Error"


class WeightableGroupAnalyzer(GroupAnalyzer, ABC):
    """
    I am the abstract base for group coefficients that support weighting.

    Coefficients like Fleiss' kappa have weighted forms for ordinal data.
    I accept weights, and take them from the semantics of the ratings:
    - Config must be WeightableGroupConfig
    - Weights come from ``ratings.semantics.categories`` (unordered: none;
      ordered: linear or quadratic; weighted: the declared matrix)
    - Built from a NijMatrix there are no semantics, so no weights

    Attributes
    ----------
    ratings : Ratings | None
        The ratings data (from Analyzer, can be None if initialized from nij_matrix).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    config : WeightableGroupConfig
        Must be an instance of WeightableGroupConfig or a subclass
        (FleissKappaConfig, etc.)
    nij_matrix : NijMatrix | None
        The counts matrix (from GroupAnalyzer, populated if built from nij_matrix).
    weights : Weights | None
        The weights declared by the semantics. None for unordered categories
        and for the NijMatrix path (unweighted).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings | None
    balance_profile: BalanceProfile
    # Populated at parent class GroupAnalyzer
    nij_matrix: NijMatrix | None
    # Checked at this level
    config: WeightableGroupConfig
    # Populated at this level
    weights: Weights | None

    _ACCEPTS_WEIGHTS = True

    # Constructors ----------------------------------------------------------- #

    def __init__(
        self, ratings: Ratings, config: WeightableGroupConfig
    ) -> None:
        """
        Initialize weightable group analyzer: validate config and build Weights.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : WeightableGroupConfig
            The analysis configuration for a weightable group coefficient.

        Raises
        ------
        TypeError
            If config is not WeightableGroupConfig.
        ValueError
        """
        # Call parent __init__ (validates group config, balance constraint)
        # Wrap to propagate errors with WeightableGroupAnalyzer context
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type; weights come from the semantics
        self._validate_weightable_config()
        from ..intermediates.estimator_parts import Weights
        self.weights = Weights.from_categories(self.ratings.semantics.categories)

    @classmethod
    def build_from_nij_matrix(
        cls, nij: NijMatrix, config: WeightableGroupConfig
    ) -> WeightableGroupAnalyzer:
        """
        Build weightable group analyzer from NijMatrix directly.

        Calls parent build_from_nij_matrix() to set ratings, config,
        balance_profile, and nij_matrix, then resolves Weights from config.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).
        config : WeightableGroupConfig
            The analysis configuration for a weightable group coefficient.

        Returns
        -------
        WeightableGroupAnalyzer
            An initialized analyzer with weights resolved.

        Raises
        ------
        TypeError
            If config is not WeightableGroupConfig.
        ValueError
        """
        # Call parent build_from_nij_matrix (sets ratings, config, balance_profile, nij_matrix)
        # Wrap to propagate errors with WeightableGroupAnalyzer context
        try:
            instance = super().build_from_nij_matrix(nij, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type. Counts carry no semantics, so no weights.
        instance._validate_weightable_config()
        instance.weights = None

        return instance

    # Public API ------------------------------------------------------------- #

    # estimate() and estimate_from_nij_matrix() inherited from GroupAnalyzer
    # (abstract, implemented by leaf classes)

    # Private helpers -------------------------------------------------------- #

    def _validate_weightable_config(self) -> None:
        """
        Validate that self.config is a WeightableGroupConfig or subclass.

        Raises
        ------
        TypeError
            If config is not WeightableGroupConfig.
        """
        from .configs import WeightableGroupConfig

        if not isinstance(self.config, WeightableGroupConfig):
            expected_configs = WeightableGroupConfig.__subclasses__()
            expected_names = ", ".join([c.__name__ for c in expected_configs])
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be WeightableGroupConfig or subclass. "
                f"Expected one of: {expected_names}. "
                f"Got {type(self.config).__name__}"
            )
