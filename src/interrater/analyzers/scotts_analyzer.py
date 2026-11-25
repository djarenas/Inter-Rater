"""Analyzer for Scott's pi coefficient."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .group_analyzer import GroupAnalyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.estimator_parts import Weights
    from ..intermediates.results import ScottsPiResult, HypothesisTestResult
    from .configs import ScottsConfig
    from ..intermediates.data_profiles.balance_profile import BalanceProfile
    from ..intermediates.data_profiles.rater_profile import RaterProfile

from ..intermediates.estimator_parts import EstimatorInput, NijMatrix
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.scotts_pi_result import (
    KAPPA_MAX_NOT_COMPUTED,
)
from .agreement_functions.group_agreement import (
    category_proportions,
    chance_fleiss,
    observed_agreement,
    scotts_pi,
)
from ..intermediates.results.ci_results import ParametricCIResult

_ERROR_PREFIX = "[ScottsAnalyzer] Error"


class ScottsAnalyzer(GroupAnalyzer):
    """
    I am the analyzer for Scott's pi coefficient.

    Scott's pi (Scott 1955) measures agreement between exactly two raters on
    nominal (unweighted) categorical data, corrected for chance using pooled
    marginal proportions.

    At exactly two raters, Fleiss' kappa equals Scott's pi. This class exists
    to honor Scott's contribution, provide a clear two-rater-only entry
    point, and refuse weighted semantics: Scott's pi has no weighted form.

    For more than two raters, use FleissKappaAnalyzer.

    My constructor checks that the config is a ScottsConfig, that the ratings
    have exactly two raters after any injection, and that parametric_null is
    asked for only on balanced data. I supply the coefficient
    (``_estimate_from_input``), my result (``_result``), and the parametric
    null interval and test, with the null variance of Fleiss' kappa (Fleiss,
    Nee & Landis 1979). Analyzer does the bootstrap and the permutation test.

    References
    ----------
    Fleiss, J. L., Nee, J. C. M., & Landis, J. R. (1979). Large sample
        variance of kappa in the case of different sets of raters.
        Psychological Bulletin, 86(5), 974-977.
    Scott, W. A. (1955). Reliability of content analysis: the case of
        nominal scale coding. Public Opinion Quarterly, 19(3), 321-325.

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
        Always None: Scott's pi has no weighted form (from Analyzer).
    estimator_input : EstimatorInput
        What my coefficient is computed from (from Analyzer).
    config : ScottsConfig
        The analysis configuration for Scott's pi.
    """

    # Label in hypothesis test results
    _COEFFICIENT_KEY = "scotts_pi"

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    weights: Weights | None
    estimator_input: EstimatorInput
    # Checked at this level
    config: ScottsConfig

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: ScottsConfig) -> None:
        """
        Run the parent's checks, then check the config type, the two raters
        and parametric_null.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze. Must have exactly 2 raters.
        config : ScottsConfig
            The analysis configuration for Scott's pi.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object, or config is not ScottsConfig.
        ValueError
            If the ratings, after any injection, do not have exactly 2 raters,
            or parametric_null is asked for on unbalanced data.
        """
        # Call parent __init__ (validates ratings and config, applies the
        # injection, checks the group config and balance constraint)
        # Wrap to propagate errors with ScottsAnalyzer context
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type
        self._validate_scotts_config()

        # Exactly 2 raters, checked after the injection, which can add raters
        self._validate_two_raters_from_ratings(self.ratings)
        self._refuse_parametric_null_if_invalid(allow_weights=False)

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H₀: π = 0 by the method in config.hypothesis_test_config:
        "parametric_null" (``_parametric_null_test``) or "permutation"
        (``_permutation_test``).
        """
        if self._hypothesis_test_config().method == "parametric_null":
            return self._parametric_null_test()
        return self._permutation_test()

    def _parametric_null_test(self) -> HypothesisTestResult:
        """
        Test H₀: π = 0 with Fleiss' null variance.

        Scott's pi equals Fleiss' kappa at two raters, so it uses the same
        null variance: Fleiss, Nee & Landis (1979) Eq. 12, with the bias
        correction of Eq. 17. The result names the source.
        """
        from ..intermediates.results import HypothesisTestResult
        from .variance_functions import (
            FLEISS_NULL_VARIANCE_SOURCE,
            fleiss_null_variance,
        )

        ht_config = self._hypothesis_test_config()
        nij = self.estimator_input.nij
        point = scotts_pi(self.estimator_input)
        pj = category_proportions(nij)

        n_subjects = nij.n_subjects
        # Balanced: _refuse_parametric_null_if_invalid refused anything else
        n_raters = int(nij.n_per_subject[0])

        # Fleiss, Nee & Landis (1979), Eq. 12
        var_null = fleiss_null_variance(pj, n_subjects, n_raters)
        se_null = np.sqrt(var_null)

        # Bias correction from Fleiss, Nee & Landis (1979), Eq. 17
        bias_correction = 1.0 / (n_subjects * (n_raters - 1))
        z = (point + bias_correction) / se_null
        p_value = self._parametric_p_value(z, ht_config.alternative)

        return HypothesisTestResult(
            p_value=p_value,
            alpha=ht_config.alpha,
            reject_null=self._reject_null(p_value, ht_config.alpha),
            method="parametric_null",
            coefficient=self._COEFFICIENT_KEY,
            point=point,
            z_statistic=z,
            null_variance=var_null,
            null_variance_source=FLEISS_NULL_VARIANCE_SOURCE,
            alternative=ht_config.alternative,
        )

    # Private helpers -------------------------------------------------------- #

    def _validate_scotts_config(self) -> None:
        """Validate that self.config is a ScottsConfig."""
        from .configs import ScottsConfig

        if not isinstance(self.config, ScottsConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be ScottsConfig, "
                f"got {type(self.config).__name__}"
            )

    def _validate_two_raters_from_ratings(self, ratings) -> None:
        """Validate exactly 2 raters."""
        if ratings.n_raters != 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: Scott's pi is defined for exactly 2 raters, "
                f"got {ratings.n_raters}. For {ratings.n_raters} raters, use "
                f"FleissKappaAnalyzer instead (Fleiss reduces to Scott's pi at n=2)."
            )

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Scott's pi from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples, bootstrap or jackknife, is written once on my base
        classes in terms of this.
        """
        return scotts_pi(data)

    def _agreement_parts(self, nij: NijMatrix) -> AgreementParts:
        """
        Po and Pe on the counts the point estimate was computed from.

        The same two functions scotts_pi combines, so the result can check
        that they rebuild the point estimate.
        """
        return AgreementParts(
            observed=observed_agreement(nij, None),
            chance=chance_fleiss(nij, None),
            kind="agreement",
            weighted=False,
        )

    def _parametric_null_interval(self, point: float) -> ParametricCIResult:
        """
        The interval from the Fleiss, Nee & Landis (1979) null variance
        (Eq. 12), for ``ci_method="parametric_null"``: no bootstrap, so no
        draws.
        """
        from scipy import stats
        from .variance_functions import (
            FLEISS_NULL_VARIANCE_SOURCE,
            fleiss_null_variance,
        )

        nij = self.estimator_input.nij
        n_subjects = nij.n_subjects
        # Balanced: _refuse_parametric_null_if_invalid refused anything else
        n_raters = int(nij.n_per_subject[0])

        var_null = fleiss_null_variance(
            category_proportions(nij), n_subjects, n_raters
        )
        se_null = np.sqrt(var_null)
        z_crit = stats.norm.ppf(1.0 - self.config.alpha / 2.0)

        return ParametricCIResult(
            lower=point - z_crit * se_null,
            upper=point + z_crit * se_null,
            alpha=self.config.alpha,
            variance_source="null",
            variance_reference=FLEISS_NULL_VARIANCE_SOURCE,
        )

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> ScottsPiResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.
        """
        from ..intermediates.results import ScottsPiResult

        return ScottsPiResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=self._agreement_parts(self.estimator_input.nij),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
            # No weighted form, so always unweighted
            weights_type="identity",
            weight_matrix=None,
        )
