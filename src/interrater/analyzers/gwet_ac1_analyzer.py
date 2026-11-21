"""Analyzer for Gwet's AC1 coefficient."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .agreement_functions.group_agreement import chance_gwet_ac1, observed_agreement
from .group_analyzer import GroupAnalyzer
from ..intermediates.estimator_parts import EstimatorInput, NijMatrix
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
    I am the analyzer for Gwet's AC1, the prevalence-robust agreement coefficient.

    AC1 measures chance-corrected agreement among any number of raters (n >= 2)
    on nominal categorical data. It addresses the "kappa paradox" -- the
    collapse of Fleiss/Cohen kappa under skewed prevalence --
    by replacing kappa's marginal-based chance term with one derived from
    category heterogeneity (also known as Simpson's Index of Diversity):

        Pe(gamma) = Tw / (K - 1),  where Tw = 1 - sum(pj^2)

    When prevalence is skewed, Tw is small, Pe is small, and AC1 stays high --
    matching the intuition that raters who mostly agree really do agree.

    AC1 is unweighted by definition. Gwet's AC2, the weighted form, is not
    implemented. I refuse ratings whose semantics declares ordered or
    weighted categories, and gwet_ac1() raises ValueError on non-identity
    weights.

    I validate that:
    - Config is GwetAC1Config
    - Data structure is appropriate for AC1

    I compute:
    - Point estimate of AC1
    - Bootstrap confidence interval (percentile or BCa)
    - Balance profile
    - Return GwetsAC1Result with all metadata

    Attributes
    ----------
    ratings : Ratings
        The validated ratings data (inherited from Analyzer).
    config : GwetAC1Config
        The analysis configuration for Gwet's AC1.
    balance_profile : BalanceProfile
        Describes the balance structure of the data (inherited from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    # Checked at this level
    config: GwetAC1Config

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: GwetAC1Config) -> None:
        """
        Initialize Gwet's AC1 analyzer: validate config and ratings.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : GwetAC1Config
            The analysis configuration for Gwet's AC1.

        Raises
        ------
        TypeError
            If config is not GwetAC1Config.
            If ratings is not a Ratings object.
        ValueError
            If config.balanced_data_only=True but data is unbalanced.
        """
        from ..data_objects.ratings import Ratings as RatingsBase

        # Validate ratings is at least a Ratings object
        if not isinstance(ratings, RatingsBase):
            raise TypeError(
                f"{_ERROR_PREFIX}: ratings must be a Ratings object, "
                f"got {type(ratings).__name__}"
            )

        # Call parent __init__ (validates config, computes BalanceProfile,
        # validates group config and balance constraint)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config is GwetAC1Config
        self._validate_gwet_config()

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test the null hypothesis that AC1 = 0 (agreement is due to chance).

        Reads method and parameters from self.config.hypothesis_test_config.
        AC1 supports permutation testing only -- there is no parametric null
        variance formula analogous to Fleiss (1979). Gwet's large-sample
        variance is for the standard error of the point estimate, not for a
        null-distribution variance.

        Returns
        -------
        HypothesisTestResult
            p-value, test decision, null distribution, and metadata.

        Raises
        ------
        ValueError
            If hypothesis_test_config is None.
            If method is not "permutation".
        """
        from ..intermediates.results import HypothesisTestResult
        from .agreement_functions.group_agreement import gwet_ac1
        from .resampling_functions.permutation_resampling import (
            within_rater_shuffle,
            pool_and_redistribute,
        )

        # Guard: hypothesis_test_config must be set
        if self.config.hypothesis_test_config is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis_test() requires "
                f"config.hypothesis_test_config to be set."
            )

        ht_config = self.config.hypothesis_test_config

        # AC1 supports permutation only
        if ht_config.method != "permutation":
            raise ValueError(
                f"{_ERROR_PREFIX}: AC1 supports only method='permutation' "
                f"for hypothesis testing. Got '{ht_config.method}'. "
                f"There is no parametric null variance formula for AC1."
            )

        # Compute observed AC1
        observed_ac1 = gwet_ac1(self.estimator_input)

        # Select permutation scheme
        scheme_map = {
            "within_rater_shuffle": within_rater_shuffle,
            "pool_and_redistribute": pool_and_redistribute,
        }
        permute_fn = scheme_map.get(ht_config.scheme)
        if permute_fn is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown permutation scheme "
                f"'{ht_config.scheme}'. Expected one of: "
                f"{set(scheme_map.keys())}."
            )

        # Permutation loop
        null_values = []
        for permuted_array in permute_fn(
            self.estimator_input.category_index,
            n_permutations=ht_config.n_permutations,
            random_state=ht_config.random_state,
        ):
            null_values.append(self._estimate_from_input(
                self._with_category_index(permuted_array)
            ))

        null_distribution = np.array(null_values)

        counts = self._permutation_p_value(
            null_distribution, observed_ac1, ht_config.alternative
        )

        return HypothesisTestResult(
            alpha=ht_config.alpha,
            reject_null=self._reject_null(counts["p_value"], ht_config.alpha),
            method="permutation",
            coefficient="gwet_ac1",
            point=observed_ac1,
            z_statistic=None,
            null_variance=None,
            null_distribution=null_distribution,
            n_permutations=ht_config.n_permutations,
            **counts,
        )

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
        Gwet's AC1 from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples, bootstrap or jackknife, is written once on my base
        classes in terms of this.

        AC1 takes no weights: Gwet (2008) defines it as unweighted.
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
            observed=observed_agreement(nij),
            chance=chance_gwet_ac1(nij),
            kind="agreement",
            weighted=False,
        )

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> GwetsAC1Result:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.
        """
        from ..intermediates.results import GwetsAC1Result

        return GwetsAC1Result(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=self._agreement_parts(self.estimator_input.nij),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
        )
