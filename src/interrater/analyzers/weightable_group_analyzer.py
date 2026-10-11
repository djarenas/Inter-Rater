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
    I add validation and construction of the Weights object from the config:
    - Config must be WeightableGroupConfig
    - Build Weights from config.weights_type and custom_weights_input
    - Validate weights_type and custom_weights_input are consistent

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
        The weights object built from config.weights_type.
        None if weights_type="identity" (unweighted).
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
            If weights_type or custom_weights_input are invalid.
        """
        # Call parent __init__ (validates group config, balance constraint)
        # Wrap to propagate errors with WeightableGroupAnalyzer context
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type and resolve weights
        self._validate_weightable_config()
        self.weights = self._resolve_weights(self.ratings.semantics.n_categories)

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
            If weights_type or custom_weights_input are invalid.
        """
        # Call parent build_from_nij_matrix (sets ratings, config, balance_profile, nij_matrix)
        # Wrap to propagate errors with WeightableGroupAnalyzer context
        try:
            instance = super().build_from_nij_matrix(nij, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type and resolve weights
        instance._validate_weightable_config()
        instance.weights = instance._resolve_weights(nij.n_categories)

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

    def _resolve_weights(self, n_categories: int) -> Weights | None:
        """
        Resolve Weights object from config.weights_type.

        Parameters
        ----------
        n_categories : int
            Number of categories in the data. Needed to build weight matrices.

        Returns
        -------
        Weights | None
            The resolved Weights object, or None if weights_type="identity".

        Raises
        ------
        ValueError
            If weights_type is "custom" but custom_weights_input is None.
            If weights_type is unknown.
        """
        if self.config.weights_type == "identity":
            return None
        elif self.config.weights_type in {"linear", "quadratic"}:
            from ..intermediates.estimator_parts import Weights

            return Weights.from_type(self.config.weights_type, n_categories)
        elif self.config.weights_type == "custom":
            if self.config.custom_weights_input is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: weights_type='custom' but custom_weights_input is None. "
                    f"This should not happen after config validation."
                )
            return self.config.custom_weights_input.build_weights()
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: Unknown weights_type '{self.config.weights_type}'"
            )
