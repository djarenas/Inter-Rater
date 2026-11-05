"""Analyzer for Krippendorff's alpha."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .agreement_functions.krippendorff import (
    expected_disagreement,
    observed_disagreement,
)
from .weightable_group_analyzer import WeightableGroupAnalyzer
from ..intermediates.estimator_parts import NijMatrix
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


class KrippendorffAlphaAnalyzer(WeightableGroupAnalyzer):
    """
    I am the analyzer for Krippendorff's alpha.

    Krippendorff's alpha measures *disagreement* rather than agreement:

        alpha = 1 - Do / De

    where Do is the observed disagreement pooled across all pairable rating
    pairs, and De is the expected disagreement under random permutation of
    the same ratings. This is not the same expression as the kappa family's
    ``(Po - Pe) / (1 - Pe)`` rearranged -- the two use different units
    (pairs pooled vs. subjects averaged) and different chance models
    (finite-sample permutation vs. population rates).

    Alpha handles missing data natively. There is no assumption that every
    subject was rated the same number of times or by the same raters. A
    subject rated once contributes no pairs and therefore contributes
    nothing. Therefore, the alpha coefficient is frequently utilized for
    unbalanced data.

    The weighted form uses a cost function d_ab = 1 - w_ab entering both
    the observed and expected disagreement terms identically. Identity
    weights give the nominal (unweighted) case.

    There is no parametric null variance formula for Krippendorff's alpha,
    so confidence intervals use bootstrap (percentile or BCa) and hypothesis
    tests use permutation only (Hayes & Krippendorff 2007).

    I validate that:
    - Config is KrippendorffAlphaConfig
    - Data structure is appropriate for Krippendorff's alpha

    I compute:
    - Point estimate of Krippendorff's alpha
    - Bootstrap confidence interval (percentile or BCa)
    - Balance profile
    - Return KrippendorffAlphaResult with all metadata

    Attributes
    ----------
    ratings : Ratings | None
        The validated ratings data (inherited from Analyzer).
        None if initialized from NijMatrix via build_from_nij_matrix().
    config : KrippendorffAlphaConfig
        The analysis configuration for Krippendorff's alpha.
    balance_profile : BalanceProfile
        Describes the balance structure of the data (inherited from Analyzer).
    rater_profile : RaterProfile | NotComputed
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped. NotComputed when built
        from an nij_matrix, which carries no rater identity.
    nij_matrix : NijMatrix | None
        The counts matrix (from GroupAnalyzer, populated if built from nij_matrix).
    weights : Weights | None
        The weights object (from WeightableGroupAnalyzer).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings | None
    balance_profile: BalanceProfile
    rater_profile: RaterProfile | NotComputed
    # Populated at parent class GroupAnalyzer
    nij_matrix: NijMatrix | None
    # Populated at parent class WeightableGroupAnalyzer
    weights: Weights | None
    # Checked at this level
    config: KrippendorffAlphaConfig

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: KrippendorffAlphaConfig) -> None:
        """
        Initialize Krippendorff's alpha analyzer: validate config and ratings.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : KrippendorffAlphaConfig
            The analysis configuration for Krippendorff's alpha.

        Raises
        ------
        TypeError
            If config is not KrippendorffAlphaConfig.
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

        # Validate config is KrippendorffAlphaConfig
        self._validate_ka_config()

    @classmethod
    def build_from_nij_matrix(
        cls, nij: NijMatrix, config: KrippendorffAlphaConfig
    ) -> KrippendorffAlphaAnalyzer:
        """
        Build Krippendorff's alpha analyzer from NijMatrix (counts matrix) directly.

        This alternate constructor is useful for analyzing published data that
        only provides the counts matrix without the original ratings.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).
        config : KrippendorffAlphaConfig
            The analysis configuration for Krippendorff's alpha.

        Returns
        -------
        KrippendorffAlphaAnalyzer
            An initialized analyzer ready to call estimate_from_nij_matrix().

        Raises
        ------
        TypeError
            If config is not KrippendorffAlphaConfig.
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

        # Validate config is KrippendorffAlphaConfig
        instance._validate_ka_config()

        return instance

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test the null hypothesis that alpha = 0 (agreement is due to chance).

        Reads method and parameters from self.config.hypothesis_test_config.
        Krippendorff's alpha supports permutation testing only -- there is no
        parametric null variance formula for this coefficient (Hayes &
        Krippendorff 2007).

        Uses a one-sided p-value: proportion of null values >= observed,
        testing whether agreement is significantly better than what shuffling
        the same ratings would produce.

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
        from .agreement_functions.krippendorff import krippendorff_alpha
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

        # Krippendorff's alpha supports permutation only
        if ht_config.method != "permutation":
            raise ValueError(
                f"{_ERROR_PREFIX}: Krippendorff's alpha supports only "
                f"method='permutation' for hypothesis testing. "
                f"Got '{ht_config.method}'. There is no parametric null "
                f"variance formula for this coefficient "
                f"(Hayes & Krippendorff 2007)."
            )

        # Permutation requires raw ratings
        self._check_permutation_available()

        n_categories = self.ratings.semantics.n_categories

        # Compute observed Krippendorff's alpha
        nij_array = category_counts(self.ratings.category_index, n_categories)
        nij = NijMatrix(matrix=nij_array)
        observed_ka = krippendorff_alpha(nij, weights=self.weights)

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
            ka_perm = krippendorff_alpha(nij_perm_obj, weights=self.weights)
            null_values.append(ka_perm)

        null_distribution = np.array(null_values)

        counts = self._permutation_p_value(
            null_distribution, observed_ka, ht_config.alternative
        )

        return HypothesisTestResult(
            alpha=ht_config.alpha,
            reject_null=self._reject_null(counts["p_value"], ht_config.alpha),
            method="permutation",
            coefficient="krippendorff_alpha",
            point=observed_ka,
            z_statistic=None,
            null_variance=None,
            null_distribution=null_distribution,
            n_permutations=ht_config.n_permutations,
            **counts,
        )

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

    def _estimate_from_nij(self, nij_array: np.ndarray) -> float:
        """
        Krippendorff's alpha from a counts matrix.

        The only coefficient-specific step in my resampling. Everything that
        draws samples -- bootstrap or jackknife, from ratings or from counts
        -- is written once on my base classes in terms of this.
        """
        from .agreement_functions.krippendorff import krippendorff_alpha
        from ..intermediates.estimator_parts import NijMatrix

        return krippendorff_alpha(
            NijMatrix(matrix=nij_array), weights=self.weights
        )

    def _disagreement_parts(self, nij_array: np.ndarray) -> AgreementParts:
        """
        Do and De on the counts the point estimate was computed from.

        The same two functions krippendorff_alpha combines, so the result can
        check that they rebuild the point estimate.
        """
        nij = NijMatrix(matrix=nij_array)
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

        nij_array = self._nij_array()
        return KrippendorffAlphaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=self._disagreement_parts(nij_array),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
            weights_type=self._weights_kind(),
            weight_matrix=self._weight_matrix(),
        )
