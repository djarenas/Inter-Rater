"""Analyzer for Brennan-Prediger's coefficient."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .agreement_functions.group_agreement import (
    chance_brennan_prediger,
    observed_agreement,
)
from .weightable_group_analyzer import WeightableGroupAnalyzer
from ..intermediates.estimator_parts import NijMatrix
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.brennan_prediger_result import (
    KAPPA_MAX_NOT_COMPUTED,
)

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.results import HypothesisTestResult, BrennanPredigerResult
    from ..intermediates.data_profiles import BalanceProfile
    from ..intermediates.estimator_parts import Weights
    from .configs import BrennanPredigerConfig

_ERROR_PREFIX = "[BrennanPredigerAnalyzer] Error"


class BrennanPredigerAnalyzer(WeightableGroupAnalyzer):
    """
    I am the analyzer for Brennan-Prediger's coefficient.

    Brennan and Prediger (1981) defined chance agreement as a uniform guess
    among K categories: Pe = 1/K -- a chance term independent of
    how the raters actually used the categories; unlike Fleiss' kappa, which
    derives chance from the observed marginals, and unlike Gwet's AC1, which
    derives it from category heterogeneity.

    The weighted form (Gwet 2014, implemented in irrCAC as ``bp.coeff.raw``)
    replaces 1/K with sum(W) / K^2, where W is the weight matrix. For
    identity weights sum(I) = K, so Pe = K/K^2 = 1/K -- the unweighted
    coefficient is the special case.

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
    ratings : Ratings | None
        The validated ratings data (inherited from Analyzer).
        None if initialized from NijMatrix via build_from_nij_matrix().
    config : BrennanPredigerConfig
        The analysis configuration for Brennan-Prediger.
    balance_profile : BalanceProfile
        Describes the balance structure of the data (inherited from Analyzer).
    nij_matrix : NijMatrix | None
        The counts matrix (from GroupAnalyzer, populated if built from nij_matrix).
    weights : Weights | None
        The weights object (from WeightableGroupAnalyzer).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings | None
    balance_profile: BalanceProfile
    # Populated at parent class GroupAnalyzer
    nij_matrix: NijMatrix | None
    # Populated at parent class WeightableGroupAnalyzer
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

        # Call parent __init__ (validates weightable config, builds Weights)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config is BrennanPredigerConfig
        self._validate_bp_config()

    @classmethod
    def build_from_nij_matrix(
        cls, nij: NijMatrix, config: BrennanPredigerConfig
    ) -> BrennanPredigerAnalyzer:
        """
        Build Brennan-Prediger analyzer from NijMatrix (counts matrix) directly.

        This alternate constructor is useful for analyzing published data that
        only provides the counts matrix without the original ratings.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).
        config : BrennanPredigerConfig
            The analysis configuration for Brennan-Prediger.

        Returns
        -------
        BrennanPredigerAnalyzer
            An initialized analyzer ready to call estimate_from_nij_matrix().

        Raises
        ------
        TypeError
            If config is not BrennanPredigerConfig.
        ValueError
            If config.balanced_data_only=True but nij indicates unbalanced data.
        """
        # Delegate to parent's build_from_nij_matrix (handles __new__,
        # Analyzer.__init__ with None, stores nij, computes BalanceProfile,
        # validates group config, balance constraint, and resolves weights)
        try:
            instance = super().build_from_nij_matrix(nij, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config is BrennanPredigerConfig
        instance._validate_bp_config()

        return instance

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
            If ratings is None (permutation requires raw ratings).
        """
        from ..intermediates.results import HypothesisTestResult
        from ..intermediates.estimator_parts import NijMatrix
        from .agreement_functions.group_agreement import brennan_prediger
        from .resampling_functions.permutation_resampling import (
            within_rater_shuffle,
            pool_and_redistribute,
        )
        from ..array_utils import category_counts

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

        # Permutation requires raw ratings
        self._check_permutation_available()

        n_categories = self.ratings.semantics.n_categories

        # Compute observed Brennan-Prediger
        nij_array = category_counts(self.ratings.category_index, n_categories)
        nij = NijMatrix(matrix=nij_array)
        observed_bp = brennan_prediger(nij, weights=self.weights)

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
            self.ratings.category_index,
            n_permutations=ht_config.n_permutations,
            random_state=ht_config.random_state,
        ):
            nij_perm = category_counts(permuted_array, n_categories)
            nij_perm_obj = NijMatrix(matrix=nij_perm)
            bp_perm = brennan_prediger(nij_perm_obj, weights=self.weights)
            null_values.append(bp_perm)

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

    def _estimate_from_nij(self, nij_array: np.ndarray) -> float:
        """
        Brennan-Prediger's coefficient from a counts matrix.

        The only coefficient-specific step in my resampling. Everything that
        draws samples -- bootstrap or jackknife, from ratings or from counts
        -- is written once on my base classes in terms of this.
        """
        from .agreement_functions.group_agreement import brennan_prediger
        from ..intermediates.estimator_parts import NijMatrix

        return brennan_prediger(
            NijMatrix(matrix=nij_array), weights=self.weights
        )

    def _agreement_parts(self, nij_array: np.ndarray) -> AgreementParts:
        """
        Po and Pe on the counts the point estimate was computed from.

        The same two functions brennan_prediger combines, so the result can
        check that they rebuild the point estimate.
        """
        nij = NijMatrix(matrix=nij_array)
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

        nij_array = self._nij_array()
        return BrennanPredigerResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            parts=self._agreement_parts(nij_array),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
            weights_type=self._weights_kind(),
            weight_matrix=self._weight_matrix(),
        )
