"""Analyzer for Krippendorff's alpha."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .agreement_functions.krippendorff import (
    expected_disagreement,
    observed_disagreement,
)
from .group_analyzer import GroupAnalyzer
from ..intermediates.estimator_parts import EstimatorInput, NijMatrix
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.krippendorff_alpha_result import (
    KAPPA_MAX_NOT_COMPUTED,
)

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.results import HypothesisTestResult, KrippendorffAlphaResult
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
    from ..intermediates.estimator_parts import Weights
    from .configs import KrippendorffAlphaConfig

_ERROR_PREFIX = "[KrippendorffAlphaAnalyzer] Error"


class KrippendorffAlphaAnalyzer(GroupAnalyzer):
    """
    I am the analyzer for Krippendorff's alpha.

    Krippendorff's alpha measures *disagreement* rather than agreement:

        alpha = 1 - Do / De

    where Do is the observed disagreement pooled across all pairable rating
    pairs, and De is the expected disagreement under random permutation of
    the same ratings. Written with Po = 1 - Do and Pe = 1 - De, this is the
    kappa family's ``(Po - Pe) / (1 - Pe)``, but the two compute Po and Pe
    differently: different units (pairs pooled vs. subjects averaged) and
    different chance models (drawing without vs. with replacement).

    Alpha handles missing data natively. There is no assumption that every
    subject was rated the same number of times or by the same raters. A
    subject rated once contributes no pairs and therefore contributes
    nothing.

    The weighted form uses a cost function d_ab = 1 - w_ab entering both
    the observed and expected disagreement terms identically. Identity
    weights give the nominal (unweighted) case.

    The package does not implement a parametric null for Krippendorff's alpha,
    so confidence intervals use bootstrap (percentile or BCa) and hypothesis
    tests use permutation only.

    I validate that:
    - Ratings and config pass my parent's checks
    - Config is KrippendorffAlphaConfig

    I compute:
    - Point estimate of Krippendorff's alpha
    - Bootstrap confidence interval (percentile or BCa)
    - Return a KrippendorffAlphaResult

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
    config : KrippendorffAlphaConfig
        The analysis configuration for Krippendorff's alpha.
    """

    # Has a weighted form: Analyzer reads the weights from the semantics
    _ACCEPTS_WEIGHTS = True

    # Label in hypothesis test results
    _COEFFICIENT_KEY = "krippendorff_alpha"

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    weights: Weights | None
    estimator_input: EstimatorInput
    # Checked at this level
    config: KrippendorffAlphaConfig

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: KrippendorffAlphaConfig) -> None:
        """
        Run the parent's checks, then check the config type.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : KrippendorffAlphaConfig
            The analysis configuration for Krippendorff's alpha.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object, or config is not
            KrippendorffAlphaConfig.
        ValueError
            If config.balanced_data_only=True but data is unbalanced.
        """
        # Call parent __init__ (validates ratings and config, computes
        # BalanceProfile, validates group config and balance constraint)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config is KrippendorffAlphaConfig
        self._validate_ka_config()

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: alpha = 0 by permutation (``_permutation_test``).

        The package does not implement a parametric null for Krippendorff's
        alpha.
        """
        return self._permutation_test()

    # Private helpers -------------------------------------------------------- #

    def _validate_ka_config(self) -> None:
        """
        Validate that self.config is a KrippendorffAlphaConfig.

        Raises
        ------
        TypeError
            If config is not KrippendorffAlphaConfig.
        """
        from .configs import KrippendorffAlphaConfig

        if not isinstance(self.config, KrippendorffAlphaConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be KrippendorffAlphaConfig, "
                f"got {type(self.config).__name__}"
            )

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Krippendorff's alpha from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples, bootstrap or jackknife, is written once on my base
        classes in terms of this.
        """
        from .agreement_functions.krippendorff import krippendorff_alpha

        return krippendorff_alpha(data)

    def _disagreement_parts(self, nij: NijMatrix) -> AgreementParts:
        """
        Do and De on the counts the point estimate was computed from.

        The same two functions krippendorff_alpha combines, so the result can
        check that they rebuild the point estimate.
        """
        return AgreementParts(
            observed=observed_disagreement(nij, self.weights),
            chance=expected_disagreement(nij, self.weights),
            kind="disagreement",
            weighted=self.weights is not None and not self.weights.is_identity,
        )

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> KrippendorffAlphaResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.
        """
        from ..intermediates.results import KrippendorffAlphaResult

        return KrippendorffAlphaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=self._disagreement_parts(self.estimator_input.nij),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
            weights_type=self._weights_kind(),
            weight_matrix=self._weight_matrix(),
        )
