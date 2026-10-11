"""Base analyzer for weightable individual marginal coefficients."""

from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

from .individual_marginal_analyzer import IndividualMarginalAnalyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.data_profiles import BalanceProfile
    from ..intermediates.estimator_parts import Weights
    from ..intermediates.estimator_parts.estimator_input import EstimatorInput
    from .configs import WeightableIndividualConfig

_ERROR_PREFIX = "[WeightableIndividualAnalyzer] Error"


class WeightableIndividualAnalyzer(IndividualMarginalAnalyzer, ABC):
    """
    I am the abstract base for individual marginal coefficients that support weighting.

    Coefficients like Cohen's kappa have weighted forms for ordinal data
    (Cohen 1968). I add validation and construction of the Weights object from
    the config, mirroring what WeightableGroupAnalyzer does on the group side.

    I resolve the Weights object from config.weights_type and rebuild the
    EstimatorInput that my parent built with weights=None. After my __init__,
    self.estimator_input carries the resolved weights and is ready for the
    agreement functions.

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer, always required).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    config : WeightableIndividualConfig
        Must be an instance of WeightableIndividualConfig or a subclass
        (CohenKappaConfig).
    estimator_input : EstimatorInput
        Rebuilt from parent's version with resolved weights.
    weights : Weights | None
        The weights object built from config.weights_type.
        None if weights_type="identity" (unweighted).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    # Populated at parent class IndividualMarginalAnalyzer
    estimator_input: EstimatorInput
    # Checked at this level
    config: WeightableIndividualConfig
    # Populated at this level
    weights: Weights | None

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: WeightableIndividualConfig
    ) -> None:
        """
        Initialize weightable individual analyzer: validate config and build
        Weights.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : WeightableIndividualConfig
            The analysis configuration for a weightable individual marginal
            coefficient.

        Raises
        ------
        TypeError
            If config is not WeightableIndividualConfig.
        ValueError
            If weights_type or custom_weights_path are invalid.
        """
        # Call parent __init__ (validates ratings, config, computes BalanceProfile,
        # builds EstimatorInput with weights=None)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate config type
        self._validate_weightable_config()

        # Resolve weights and rebuild EstimatorInput with them
        self.weights = self._resolve_weights(
            self.ratings.semantics.n_categories
        )
        if self.weights is not None:
            self.estimator_input = self._rebuild_estimator_input_with_weights()

    # Public API ------------------------------------------------------------- #

    # estimate() and hypothesis_test() inherited from IndividualMarginalAnalyzer
    # (abstract, implemented by leaf classes)

    # Private helpers -------------------------------------------------------- #

    def _validate_weightable_config(self) -> None:
        """
        Validate that self.config is a WeightableIndividualConfig or subclass.

        Raises
        ------
        TypeError
            If config is not WeightableIndividualConfig.
        """
        from .configs import WeightableIndividualConfig

        if not isinstance(self.config, WeightableIndividualConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be "
                f"WeightableIndividualConfig or subclass "
                f"(CohenKappaConfig). "
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
            If weights_type is "custom" but custom_weights_path is None.
            If weights_type is unknown.
        """
        if self.config.weights_type == "identity":
            return None
        elif self.config.weights_type in {"linear", "quadratic"}:
            from ..intermediates.estimator_parts import Weights

            return Weights.from_type(self.config.weights_type, n_categories)
        elif self.config.weights_type == "custom":
            if self.config.custom_weights_path is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: weights_type='custom' but "
                    f"custom_weights_path is None. This should not happen "
                    f"after config validation."
                )
            from ..intermediates.estimator_parts import Weights

            return Weights.from_file(self.config.custom_weights_path)
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: Unknown weights_type "
                f"'{self.config.weights_type}'"
            )

    def _rebuild_estimator_input_with_weights(self) -> EstimatorInput:
        """
        Rebuild EstimatorInput with resolved weights.

        Called after _resolve_weights when weights are not None. Replaces the
        EstimatorInput that IndividualMarginalAnalyzer built with weights=None.

        Returns
        -------
        EstimatorInput
            Carrying category_index, n_categories, and resolved weights.
        """
        from ..intermediates.estimator_parts.estimator_input import EstimatorInput

        return EstimatorInput(
            category_index=self.ratings.category_index,
            n_categories=self.ratings.semantics.n_categories,
            weights=self.weights,
        )
