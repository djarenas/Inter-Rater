"""Analyzer for Gwet's AC1 coefficient."""

from __future__ import annotations

from typing import TYPE_CHECKING

from .agreement_functions.group_agreement import chance_gwet_ac1, observed_agreement
from .group_analyzer import GroupAnalyzer
from ..intermediates.estimator_parts import EstimatorInput, NijMatrix, Weights
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.gwets_ac1_result import (
    KAPPA_MAX_NOT_COMPUTED,
)

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.results import HypothesisTestResult, GwetsAC1Result
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
    from .configs import GwetAC1Config

_ERROR_PREFIX = "[GwetAC1Analyzer] Error"


class GwetAC1Analyzer(GroupAnalyzer):
    """
    I am the analyzer for Gwet's AC1, the prevalence-robust agreement
    coefficient, and for AC2, its weighted form.

    AC1 measures chance-corrected agreement among two or more raters
    on categorical data. It addresses the low values of Fleiss' and Cohen's
    kappa under skewed prevalence
    by replacing kappa's marginal-based chance term with one derived from
    category heterogeneity (Gwet 2008, Equation 38):

        Pe = sum_j p_j (1 - p_j) / (K - 1)

    When prevalence is skewed, Pe is small, so AC1 stays close to Po.

    With ordered or weighted categories in the semantics I compute AC2: Po
    credits near misses through the weights, and Pe is multiplied by
    T_w / K, where T_w is the sum of the weights (Gwet 2021, Chapter 4);
    see ``chance_gwet_ac1``.

    I validate that:
    - Ratings and config pass my parent's checks
    - Config is GwetAC1Config

    I compute:
    - Point estimate of AC1
    - Bootstrap confidence interval (percentile or BCa)
    - Return a GwetsAC1Result

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
    config : GwetAC1Config
        The analysis configuration for Gwet's AC1.
    """

    # Has a weighted form (AC2): Analyzer reads the weights from the semantics
    _ACCEPTS_WEIGHTS = True

    # Label in hypothesis test results
    _COEFFICIENT_KEY = "gwet_ac1"

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    weights: Weights | None
    estimator_input: EstimatorInput
    # Checked at this level
    config: GwetAC1Config

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: GwetAC1Config) -> None:
        """
        Run the parent's checks, then check the config type.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : GwetAC1Config
            The analysis configuration for Gwet's AC1.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object, or config is not GwetAC1Config.
        ValueError
            If config.balanced_data_only=True but data is unbalanced.
        """
        # Call parent __init__ (validates ratings and config, computes
        # BalanceProfile, validates group config and balance constraint)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config is GwetAC1Config
        self._validate_gwet_config()

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: AC1 = 0 by permutation (``_permutation_test``).

        The package does not implement a parametric null for AC1.
        """
        return self._permutation_test()

    # Private helpers -------------------------------------------------------- #

    def _validate_gwet_config(self) -> None:
        """
        Validate that self.config is a GwetAC1Config.

        Raises
        ------
        TypeError
            If config is not GwetAC1Config.
        """
        from .configs import GwetAC1Config

        if not isinstance(self.config, GwetAC1Config):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be GwetAC1Config, "
                f"got {type(self.config).__name__}"
            )

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Gwet's AC1 (AC2 when weighted) from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples, bootstrap or jackknife, is written once on my base
        classes in terms of this.
        """
        from .agreement_functions.group_agreement import gwet_ac1

        return gwet_ac1(data)

    def _agreement_parts(self, nij: NijMatrix) -> AgreementParts:
        """
        Po and Pe on the counts the point estimate was computed from.

        The same two functions gwet_ac1 combines, so the result can check
        that they rebuild the point estimate.
        """
        return AgreementParts(
            observed=observed_agreement(nij, self.weights),
            chance=chance_gwet_ac1(nij, self.weights),
            kind="agreement",
            weighted=self.weights is not None and not self.weights.is_identity,
        )

    def _result(self, point: float, ci) -> GwetsAC1Result:
        """
        My result: the point and interval ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.
        """
        from ..intermediates.results import GwetsAC1Result

        return GwetsAC1Result(
            point=point,
            confidence_interval=ci,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=self._agreement_parts(self.estimator_input.nij),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
            weights_type=self._weights_kind(),
            weight_matrix=self._weight_matrix(),
        )
