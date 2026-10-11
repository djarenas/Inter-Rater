"""Analyzer for Gwet's AC1 coefficient."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .group_analyzer import GroupAnalyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.analysis_results import HypothesisTestResult, GwetsAC1Result
    from ..intermediates.data_profiles import BalanceProfile
    from ..intermediates.estimator_parts import NijMatrix
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

    AC1 is unweighted by definition. For ordinal (weighted) analysis, use
    Gwet's AC2 instead, which is a separate coefficient with its own published
    definition. The agreement function gwet_ac1() takes no weights argument
    -- passing weights is a TypeError at the call site, not a silent
    ignore.

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
    ratings : Ratings | None
        The validated ratings data (inherited from Analyzer).
        None if initialized from NijMatrix via build_from_nij_matrix().
    config : GwetAC1Config
        The analysis configuration for Gwet's AC1.
    balance_profile : BalanceProfile
        Describes the balance structure of the data (inherited from Analyzer).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings | None
    balance_profile: BalanceProfile
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

    @classmethod
    def build_from_nij_matrix(
        cls, nij: NijMatrix, config: GwetAC1Config
    ) -> GwetAC1Analyzer:
        """
        Build AC1 analyzer from NijMatrix (counts matrix) directly.

        This alternate constructor is useful for analyzing published data that
        only provides the counts matrix without the original ratings.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).
        config : GwetAC1Config
            The analysis configuration for Gwet's AC1.

        Returns
        -------
        GwetAC1Analyzer
            An initialized analyzer ready to call estimate_from_nij_matrix().

        Raises
        ------
        TypeError
            If config is not GwetAC1Config.
        ValueError
            If config.balanced_data_only=True but nij indicates unbalanced data.
        """
        # Delegate to parent's build_from_nij_matrix (handles __new__,
        # Analyzer.__init__ with None, stores nij, computes BalanceProfile,
        # validates group config and balance constraint)
        try:
            instance = super().build_from_nij_matrix(nij, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config is GwetAC1Config
        instance._validate_gwet_config()

        return instance

    # Public API ------------------------------------------------------------- #

    def estimate(self, use_scipy_bca: bool = False) -> GwetsAC1Result:
        """
        Compute AC1 from ratings and return the result.

        Performs point estimation and bootstrap confidence interval computation.
        The bootstrap loop resamples subjects with replacement, recomputes AC1
        on each resample, and builds a distribution for CI bounds.

        Parameters
        ----------
        use_scipy_bca : bool, default False
            If True and config.ci_method == "bca", use scipy.stats for BCa
            computation (for testing/validation only).

        Returns
        -------
        GwetsAC1Result
            Point estimate, confidence interval, and bootstrap distribution.

        Raises
        ------
        ValueError
            If self.ratings is None (should use estimate_from_nij_matrix instead).
        NotImplementedError
            If config.convergence_test is set (not yet implemented).
        """
        from ..intermediates.estimator_parts import NijMatrix
        from .agreement_functions.group_agreement import gwet_ac1
        from .resampling_functions.subject_resampling import subject_bootstrap
        from ..array_utils import category_counts
        from ..intermediates.analysis_results import GwetsAC1Result

        # Guard: ratings must be available for this path
        if self.ratings is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: estimate() requires ratings data. "
                f"This analyzer was built from NijMatrix -- use "
                f"estimate_from_nij_matrix() instead."
            )

        # Guard convergence_test (not yet implemented)
        if self.config.convergence_test is not None:
            raise NotImplementedError(
                f"{_ERROR_PREFIX}: convergence_test not yet implemented"
            )

        # Get number of categories
        n_categories = self.ratings.semantics.n_categories

        # Compute point estimate from original data
        nij_array = category_counts(self.ratings.category_index, n_categories)
        nij = NijMatrix(matrix=nij_array)
        point = gwet_ac1(nij)  # No weights -- AC1 is unweighted

        # Bootstrap loop
        bootstrap_values = []

        for resampled_array in subject_bootstrap(
            self.ratings.category_index,
            n_resamples=self.config.n_resamples,
            random_state=self.config.random_state,
        ):
            nij_array_resamp = category_counts(resampled_array, n_categories)
            nij_resampled = NijMatrix(matrix=nij_array_resamp)
            ac1_resample = gwet_ac1(nij_resampled)
            bootstrap_values.append(ac1_resample)

        bootstrap_distribution = np.array(bootstrap_values)

        # Compute confidence interval
        ci = self._compute_ci(bootstrap_distribution, point, use_scipy_bca)

        # Build and return result
        return GwetsAC1Result(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            ci_method=self.config.ci_method,
        )

    def estimate_from_nij_matrix(
        self, use_scipy_bca: bool = False
    ) -> GwetsAC1Result:
        """
        Compute AC1 from NijMatrix (counts matrix) and return the result.

        Uses subject bootstrap resampling on the counts matrix directly.

        Parameters
        ----------
        use_scipy_bca : bool, default False
            If True and config.ci_method == "bca", use scipy.stats for BCa.

        Returns
        -------
        GwetsAC1Result
            Point estimate, confidence interval, and bootstrap distribution.

        Raises
        ------
        ValueError
            If self.nij_matrix is None (should use estimate instead).
        NotImplementedError
            If config.convergence_test is set (not yet implemented).
        """
        from .agreement_functions.group_agreement import gwet_ac1
        from .resampling_functions.subject_resampling import subject_bootstrap_from_nij_matrix
        from ..intermediates.estimator_parts import NijMatrix
        from ..intermediates.analysis_results import GwetsAC1Result

        # Guard: nij_matrix must be available for this path
        if self.nij_matrix is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: estimate_from_nij_matrix() requires nij_matrix. "
                f"This analyzer was built from Ratings -- use estimate() instead."
            )

        # Guard convergence_test (not yet implemented)
        if self.config.convergence_test is not None:
            raise NotImplementedError(
                f"{_ERROR_PREFIX}: convergence_test not yet implemented"
            )

        # Compute point estimate from nij_matrix
        point = gwet_ac1(self.nij_matrix)

        # Bootstrap loop on nij_matrix
        bootstrap_values = []
        n_categories = self.nij_matrix.n_categories

        for resampled_nij_array in subject_bootstrap_from_nij_matrix(
            self.nij_matrix.matrix,
            n_resamples=self.config.n_resamples,
            random_state=self.config.random_state,
        ):
            nij_resampled = NijMatrix(matrix=resampled_nij_array)
            ac1_resample = gwet_ac1(nij_resampled)
            bootstrap_values.append(ac1_resample)

        bootstrap_distribution = np.array(bootstrap_values)

        # Compute confidence interval
        ci = self._compute_ci(bootstrap_distribution, point, use_scipy_bca)

        # Build and return result
        return GwetsAC1Result(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            ci_method=self.config.ci_method,
        )
    
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
            If ratings is None (permutation requires raw ratings).
        """
        from ..intermediates.analysis_results import HypothesisTestResult
        from ..intermediates.estimator_parts import NijMatrix
        from .agreement_functions.group_agreement import gwet_ac1
        from .resampling_functions.permutation_resampling import (
            column_shuffle,
            within_subject_shuffle,
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

        # AC1 supports permutation only
        if ht_config.method != "permutation":
            raise ValueError(
                f"{_ERROR_PREFIX}: AC1 supports only method='permutation' "
                f"for hypothesis testing. Got '{ht_config.method}'. "
                f"There is no parametric null variance formula for AC1."
            )

        # Permutation requires raw ratings
        self._check_permutation_available()

        n_categories = self.ratings.semantics.n_categories

        # Compute observed AC1
        nij_array = category_counts(self.ratings.category_index, n_categories)
        nij = NijMatrix(matrix=nij_array)
        observed_ac1 = gwet_ac1(nij)

        # Select permutation scheme
        scheme_map = {
            "column_shuffle": column_shuffle,
            "within_subject": within_subject_shuffle,
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
            ac1_perm = gwet_ac1(nij_perm_obj)
            null_values.append(ac1_perm)

        null_distribution = np.array(null_values)

        # p-value: proportion of null values >= observed
        p_value = float(np.mean(null_distribution >= observed_ac1))

        return HypothesisTestResult(
            p_value=p_value,
            alpha=ht_config.alpha,
            reject_null=(p_value < ht_config.alpha),
            method="permutation",
            point=observed_ac1,
            z_statistic=None,
            null_variance=None,
            null_distribution=null_distribution,
            n_permutations=ht_config.n_permutations,
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

    def _compute_ci(
        self,
        bootstrap_distribution: np.ndarray,
        point: float,
        use_scipy_bca: bool,
    ) -> tuple[float, float]:
        """
        Compute confidence interval from bootstrap distribution.

        Parameters
        ----------
        bootstrap_distribution : np.ndarray
            Array of bootstrap coefficient values.
        point : float
            Point estimate of AC1.
        use_scipy_bca : bool
            If True and ci_method == "bca", use scipy BCa.

        Returns
        -------
        tuple[float, float]
            Lower and upper bounds of the confidence interval.
        """
        from .ci_methods import percentile_ci, bca_ci, scipy_bca_ci

        if self.config.ci_method == "percentile":
            return percentile_ci(bootstrap_distribution, self.config.alpha)
        elif self.config.ci_method == "bca":
            if use_scipy_bca:
                return scipy_bca_ci(bootstrap_distribution, point, self.config.alpha)
            else:
                return bca_ci(bootstrap_distribution, point, self.config.alpha)
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown ci_method '{self.config.ci_method}'. "
                f"Expected 'percentile' or 'bca'."
            )
