"""Base analyzer for weightable group coefficients."""

from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

from .group_analyzer import GroupAnalyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.estimator_parts import Weights
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
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

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    config : WeightableGroupConfig
        Must be an instance of WeightableGroupConfig or a subclass
        (FleissKappaConfig, etc.)
    weights : Weights | None
        The weights declared by the semantics. None for unordered categories.
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
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
        self.weights = Weights.build_from_categories(self.ratings.semantics.categories)

    # Public API ------------------------------------------------------------- #

    # estimate() inherited from GroupAnalyzer
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
