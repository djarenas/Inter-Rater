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
    (Cohen 1968). I accept weights and take them from the semantics of the
    ratings, mirroring what WeightableGroupAnalyzer does on the group side.

    I read the weights from ``ratings.semantics.categories`` and rebuild the
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
        The weights declared by the semantics. None for unordered categories
        (unweighted).
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

    _ACCEPTS_WEIGHTS = True

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
        """
        # Call parent __init__ (validates ratings, config, computes BalanceProfile,
        # builds EstimatorInput with weights=None)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate config type
        self._validate_weightable_config()

        # Weights come from the semantics; rebuild EstimatorInput with them
        from ..intermediates.estimator_parts import Weights
        self.weights = Weights.from_categories(self.ratings.semantics.categories)
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

    def _rebuild_estimator_input_with_weights(self) -> EstimatorInput:
        """
        Rebuild EstimatorInput with resolved weights.

        Called when the semantics declares weights. Replaces the
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
