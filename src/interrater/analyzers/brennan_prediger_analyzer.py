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
    among K categories: Pe = 1/K -- a chance term independent of
    how the raters actually used the categories; unlike Fleiss' kappa, which
    derives chance from the observed marginals, and unlike Gwet's AC1, which
    derives it from category heterogeneity.

    The weighted form (as in Gwet's irrCAC 1.4, ``bp.coeff.raw``) replaces
    1/K with sum(W) / K^2, where W is the weight matrix. For identity
    weights sum(I) = K, so Pe = K/K^2 = 1/K -- the unweighted coefficient is
    the special case.

    There is no parametric null variance formula for Brennan-Prediger, so
    confidence intervals use bootstrap (percentile or BCa) and hypothesis
    tests use permutation only.

    I validate that:
    - Config is BrennanPredigerConfig
    - Data structure is appropriate for Brennan-Prediger

    I compute:
    - Point estimate of Brennan-Prediger's coefficient
    - Bootstrap confidence interval (percentile or BCa)
    - Balance profile
    - Return BrennanPredigerResult with all metadata

    Attributes
    ----------
    ratings : Ratings
        The validated ratings data (inherited from Analyzer).
    config : BrennanPredigerConfig
        The analysis configuration for Brennan-Prediger.
    balance_profile : BalanceProfile
        Describes the balance structure of the data (inherited from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    weights : Weights | None
        The weights object (from Analyzer).
    """

    # Has a weighted form: Analyzer reads the weights from the semantics
    _ACCEPTS_WEIGHTS = True

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    weights: Weights | None
    # Checked at this level
    config: BrennanPredigerConfig

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: BrennanPredigerConfig) -> None:
        """
        Initialize Brennan-Prediger analyzer: validate config and ratings.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : BrennanPredigerConfig
            The analysis configuration for Brennan-Prediger.

        Raises
        ------
        TypeError
            If config is not BrennanPredigerConfig.
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

        # Call parent __init__ (validates group config; Analyzer builds the Weights)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config is BrennanPredigerConfig
        self._validate_bp_config()

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test the null hypothesis that BP = 0 (agreement is due to chance).

        Reads method and parameters from self.config.hypothesis_test_config.
        Brennan-Prediger supports permutation testing only -- there is no
        parametric null variance formula for this coefficient.

        Uses a one-sided p-value: proportion of null values >= observed,
        testing whether agreement is significantly better than uniform
        random guessing.

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
        from .agreement_functions.group_agreement import brennan_prediger
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

        # Brennan-Prediger supports permutation only
        if ht_config.method != "permutation":
            raise ValueError(
                f"{_ERROR_PREFIX}: Brennan-Prediger supports only "
                f"method='permutation' for hypothesis testing. "
                f"Got '{ht_config.method}'. There is no parametric null "
                f"variance formula for this coefficient."
            )

        # Compute observed Brennan-Prediger
        observed_bp = brennan_prediger(self.estimator_input)

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
            null_distribution, observed_bp, ht_config.alternative
        )

        return HypothesisTestResult(
            alpha=ht_config.alpha,
            reject_null=self._reject_null(counts["p_value"], ht_config.alpha),
            method="permutation",
            coefficient="brennan_prediger",
            point=observed_bp,
            z_statistic=None,
            null_variance=None,
            null_distribution=null_distribution,
            n_permutations=ht_config.n_permutations,
            **counts,
        )

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
