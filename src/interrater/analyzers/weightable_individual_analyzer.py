"""Base analyzer for weightable individual marginal coefficients."""

from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

from .individual_marginal_analyzer import IndividualMarginalAnalyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
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

    Analyzer reads the weights from ``ratings.semantics.categories``
    because I say I accept them, so self.estimator_input carries them.

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer, always required).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped. Always present here:
        this coefficient needs rater identity.
    config : WeightableIndividualConfig
        Must be an instance of WeightableIndividualConfig or a subclass
        (CohenKappaConfig).
    estimator_input : EstimatorInput
        Carrying the weights (from Analyzer).
    weights : Weights | None
        The weights declared by the semantics (from Analyzer). None for
        unordered categories (unweighted).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    estimator_input: EstimatorInput
    weights: Weights | None
    # Checked at this level
    config: WeightableIndividualConfig

    _ACCEPTS_WEIGHTS = True

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: WeightableIndividualConfig
    ) -> None:
        """
        Initialize weightable individual analyzer: validate config.

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
        # builds EstimatorInput with the weights)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate config type
        self._validate_weightable_config()

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
