"""Analyzer for Brennan-Prediger's coefficient."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .agreement_functions.group_agreement import (
    chance_brennan_prediger,
    observed_agreement,
)
from .group_analyzer import GroupAnalyzer
from ..intermediates.estimator_parts import EstimatorInput, NijMatrix
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.brennan_prediger_result import (
    KAPPA_MAX_NOT_COMPUTED,
)

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.results import HypothesisTestResult, BrennanPredigerResult
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
    from ..intermediates.estimator_parts import Weights
    from .configs import BrennanPredigerConfig

_ERROR_PREFIX = "[BrennanPredigerAnalyzer] Error"


class BrennanPredigerAnalyzer(GroupAnalyzer):
    """
    I am the analyzer for Brennan-Prediger's coefficient.

    Brennan and Prediger (1981) defined chance agreement as a uniform guess
    among K categories: Pe = 1/K, a chance term independent of
    how the raters actually used the categories; unlike Fleiss' kappa, which
    derives chance from the pooled category proportions, and unlike Gwet's
    AC1, which derives it from category heterogeneity.

    The weighted form (as in Gwet's irrCAC 1.4, ``bp.coeff.raw``) replaces
    1/K with sum(W) / K^2, where W is the weight matrix. For identity
    weights sum(I) = K, so Pe = K/K^2 = 1/K -- the unweighted coefficient is
    the special case.

    The package implements no parametric null for Brennan-Prediger, so
    confidence intervals use bootstrap (percentile or BCa) and hypothesis
    tests use permutation only.

    I validate that:
    - Ratings and config pass my parent's checks
    - Config is BrennanPredigerConfig

    I compute:
    - Point estimate of Brennan-Prediger's coefficient
    - Bootstrap confidence interval (percentile or BCa)
    - Return a BrennanPredigerResult

    References
    ----------
    Brennan, R. L., & Prediger, D. J. (1981). Coefficient kappa: Some
        uses, misuses, and alternatives. Educational and Psychological
        Measurement, 41(3), 687-699.

    Attributes
    ----------
    ratings : Ratings
        The validated ratings data (from Analyzer).
    balance_profile : BalanceProfile
        Describes the balance structure of the data (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    weights : Weights | None
        The weights object (from Analyzer).
    estimator_input : EstimatorInput
        What my coefficient is computed from (from Analyzer).
    config : BrennanPredigerConfig
        The analysis configuration for Brennan-Prediger.
    """

    # Has a weighted form: Analyzer reads the weights from the semantics
    _ACCEPTS_WEIGHTS = True

    # Label in hypothesis test results
    _COEFFICIENT_KEY = "brennan_prediger"

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    weights: Weights | None
    estimator_input: EstimatorInput
    # Checked at this level
    config: BrennanPredigerConfig

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: BrennanPredigerConfig) -> None:
        """
        Run the parent's checks, then check the config type.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : BrennanPredigerConfig
            The analysis configuration for Brennan-Prediger.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object, or config is not
            BrennanPredigerConfig.
        ValueError
            If config.balanced_data_only=True but data is unbalanced.
        """
        # Call parent __init__ (validates ratings and config, computes
        # BalanceProfile, validates group config and balance constraint)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config is BrennanPredigerConfig
        self._validate_bp_config()

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: Brennan-Prediger = 0 by permutation (``_permutation_test``).

        The package implements no parametric null for Brennan-Prediger.
        """
        return self._permutation_test()

    # Private helpers -------------------------------------------------------- #

    def _validate_bp_config(self) -> None:
        """
        Validate that self.config is a BrennanPredigerConfig.

        Raises
        ------
        TypeError
            If config is not BrennanPredigerConfig.
        """
        from .configs import BrennanPredigerConfig

        if not isinstance(self.config, BrennanPredigerConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be BrennanPredigerConfig, "
                f"got {type(self.config).__name__}"
            )

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Brennan-Prediger's coefficient from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples, bootstrap or jackknife, is written once on my base
        classes in terms of this.
        """
        from .agreement_functions.group_agreement import brennan_prediger

        return brennan_prediger(data)

    def _agreement_parts(self, nij: NijMatrix) -> AgreementParts:
        """
        Po and Pe on the counts the point estimate was computed from.

        The same two functions brennan_prediger combines, so the result can
        check that they rebuild the point estimate.
        """
        return AgreementParts(
            observed=observed_agreement(nij, self.weights),
            chance=chance_brennan_prediger(nij, self.weights),
            kind="agreement",
            weighted=self.weights is not None and not self.weights.is_identity,
        )

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> BrennanPredigerResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.
        """
        from ..intermediates.results import BrennanPredigerResult

        return BrennanPredigerResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=self._agreement_parts(self.estimator_input.nij),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
            weights_type=self._weights_kind(),
            weight_matrix=self._weight_matrix(),
        )
