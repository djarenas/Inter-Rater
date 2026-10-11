"""Analyzer for Fleiss' kappa coefficient."""

from __future__ import annotations

import numpy as np

from .weightable_group_analyzer import WeightableGroupAnalyzer

from ..data_objects.ratings import Ratings
from ..intermediates.analysis_results import FleissKappaResult, HypothesisTestResult
from .configs import FleissKappaConfig

from .agreement_functions.group_agreement import fleiss_kappa
from .resampling_functions.subject_resampling import subject_bootstrap
from .ci_methods import percentile_ci, bca_ci, scipy_bca_ci
from ..array_utils import category_counts

from ..intermediates.estimator_parts import NijMatrix, Weights
from ..intermediates.data_profiles import BalanceProfile

_ERROR_PREFIX = "[FleissKappaAnalyzer] Error"


class FleissKappaAnalyzer(WeightableGroupAnalyzer):
    """
    I am the analyzer for Fleiss' kappa, the multi-rater agreement coefficient.

    Fleiss' kappa measures agreement among many raters (n >= 2) on categorical
    data. It discards rater identity and reads a counts matrix. It supports
    weighting for ordinal data.

    I validate that:
    - Config is FleissKappaConfig
    - Ratings are balanced (or auto-convert if needed)
    - Data structure is appropriate for Fleiss' kappa

    I compute:
    - Point estimate of Fleiss' kappa
    - Bootstrap confidence interval (percentile or BCa)
    - Balance profile and exact/approximate flag
    - Return FleissKappaResult with all metadata

    Attributes
    ----------
    ratings : Ratings | None
        The ratings data (from Analyzer, can be None if initialized from nij_matrix).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    config : FleissKappaConfig
        The analysis configuration for Fleiss' kappa.
    nij_matrix : NijMatrix | None
        The counts matrix (from GroupAnalyzer, populated if built from nij_matrix).
    weights : Weights | None
        The weights object (from WeightableGroupAnalyzer).
    is_exact : bool
        Whether the exact formula applies (data is perfectly balanced).
        Computed during initialization based on balance_profile.
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
    config: FleissKappaConfig
    # Populated at this level
    is_exact: bool

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: FleissKappaConfig) -> None:
        """
        Initialize Fleiss' kappa analyzer: validate config and convert ratings if needed.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : FleissKappaConfig
            The analysis configuration for Fleiss' kappa.

        Raises
        ------
        TypeError
            If config is not FleissKappaConfig.
            If ratings is not a Ratings object.
        ValueError
            If config.balanced_data_only=True but data is not balanced.
        """
        from ..data_objects.ratings import Ratings as RatingsBase
        from ..data_objects.balanced_ratings import BalancedRatings

        # Validate ratings is at least a Ratings object
        if not isinstance(ratings, RatingsBase):
            raise TypeError(
                f"{_ERROR_PREFIX}: ratings must be a Ratings object, "
                f"got {type(ratings).__name__}"
            )

        # Auto-convert Ratings to BalancedRatings if config requires it
        if config.balanced_data_only and not isinstance(ratings, BalancedRatings):
            try:
                ratings = BalancedRatings.build_from_ratings(ratings)
            except ValueError as e:
                raise ValueError(
                    f"{_ERROR_PREFIX}: config.balanced_data_only=True but data is not balanced. "
                    f"Details: {str(e)}"
                ) from e

        # Call parent __init__ (validates weightable config, builds Weights)
        # Wrap to propagate errors with FleissKappaAnalyzer context
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type and compute is_exact flag
        self._validate_fleiss_config()
        self.is_exact = self.balance_profile.is_balanced

    @classmethod
    def build_from_nij_matrix(
        cls, nij, config: FleissKappaConfig
    ) -> FleissKappaAnalyzer:
        """
        Build Fleiss' kappa analyzer from NijMatrix directly.

        Calls parent build_from_nij_matrix() to set ratings, config,
        balance_profile, nij_matrix, and weights, then sets is_exact flag.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).
        config : FleissKappaConfig
            The analysis configuration for Fleiss' kappa.

        Returns
        -------
        FleissKappaAnalyzer
            An initialized analyzer ready to call estimate_from_nij_matrix().

        Raises
        ------
        TypeError
            If config is not FleissKappaConfig.
        ValueError
            If config is invalid.
        """
        # Call parent build_from_nij_matrix (sets ratings, config, balance_profile, nij_matrix, weights)
        # Wrap to propagate errors with FleissKappaAnalyzer context
        try:
            instance = super().build_from_nij_matrix(nij, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type and compute is_exact flag
        instance._validate_fleiss_config()
        instance.is_exact = instance.balance_profile.is_balanced

        return instance

    # Public API ------------------------------------------------------------- #

    def estimate(self, use_scipy_bca: bool = False) -> FleissKappaResult:
        """
        Compute Fleiss' kappa and return the result.
        Performs point estimation and confidence interval computation.

        When ci_method is "percentile" or "bca", runs the bootstrap loop.
        When ci_method is "parametric_null", uses Fleiss (1979) null variance
        -- no bootstrap is performed and bootstrap_distribution is None.

        Parameters
        ----------
        use_scipy_bca : bool, default False
            If True and config.ci_method == "bca", use scipy.stats for BCa
            computation (for testing/validation only).
            If False, use custom BCa implementation.

        Returns
        -------
        FleissKappaResult
        """

        # Guard: ratings must be available for this path
        if self.ratings is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: ratings is None. "
                f"Use estimate_from_nij_matrix() when initialized from NijMatrix."
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
        point = fleiss_kappa(nij, weights=self.weights)

        # CI path depends on ci_method
        if self.config.ci_method == "parametric_null":
            # Parametric null CI: no bootstrap needed
            pj = self._compute_pj(nij)
            ci = self._parametric_null_ci(point, pj)
            bootstrap_distribution = None
        else:
            # Bootstrap CI: run the loop
            bootstrap_values = []

            for resampled_array in subject_bootstrap(
                self.ratings.category_index,
                n_resamples=self.config.n_resamples,
                random_state=self.config.random_state,
            ):
                nij_array_resamp = category_counts(resampled_array, n_categories)
                nij_resampled = NijMatrix(matrix=nij_array_resamp)
                kappa_resample = fleiss_kappa(nij_resampled, weights=self.weights)
                bootstrap_values.append(kappa_resample)

            bootstrap_distribution = np.array(bootstrap_values)
            ci = self._compute_ci(bootstrap_distribution, point, use_scipy_bca)

        # Build and return result
        return FleissKappaResult(
            point=point,
            confidence_interval=ci,
            balance_profile=self.balance_profile,
            bootstrap_distribution=bootstrap_distribution,
            ci_method=self.config.ci_method,
            weights_type=self.config.weights_type,
        )

    def estimate_from_nij_matrix(self, use_scipy_bca: bool = False) -> FleissKappaResult:
        """
        Compute Fleiss' kappa from nij_matrix (NijMatrix-based path).

        Called when analyzer was initialized via build_from_nij_matrix(nij, config).
        Performs point estimation and confidence interval computation
        directly from the counts matrix.

        When ci_method is "percentile" or "bca", runs the bootstrap loop.
        When ci_method is "parametric_null", uses Fleiss (1979) null variance
        -- no bootstrap is performed and bootstrap_distribution is None.

        Parameters
        ----------
        use_scipy_bca : bool, default False
            If True and config.ci_method == "bca", use scipy.stats for BCa
            computation (for testing/validation only).
            If False, use custom BCa implementation.

        Returns
        -------
        FleissKappaResult
            Point estimate, CI, bootstrap distribution, balance profile, is_exact flag.

        Raises
        ------
        ValueError
            If self.nij_matrix is None.
        NotImplementedError
            If config.convergence_test is set.
        """
        from .resampling_functions import subject_bootstrap_from_nij_matrix

        # Guard: nij_matrix must be available for this path
        if self.nij_matrix is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: nij_matrix is None. "
                f"Use estimate() when initialized from Ratings."
            )

        # Guard convergence_test (not yet implemented)
        if self.config.convergence_test is not None:
            raise NotImplementedError(
                f"{_ERROR_PREFIX}: convergence_test not yet implemented"
            )

        # Compute point estimate from original nij_matrix
        point = fleiss_kappa(self.nij_matrix, weights=self.weights)

        # CI path depends on ci_method
        if self.config.ci_method == "parametric_null":
            # Parametric null CI: no bootstrap needed
            pj = self._compute_pj(self.nij_matrix)
            ci = self._parametric_null_ci(point, pj)
            bootstrap_distribution = None
        else:
            # Bootstrap CI: run the loop
            bootstrap_values = []

            for resampled_nij_array in subject_bootstrap_from_nij_matrix(
                self.nij_matrix.matrix,
                n_resamples=self.config.n_resamples,
                random_state=self.config.random_state,
            ):
                nij_resampled = NijMatrix(matrix=resampled_nij_array)
                kappa_resample = fleiss_kappa(nij_resampled, weights=self.weights)
                bootstrap_values.append(kappa_resample)

            bootstrap_distribution = np.array(bootstrap_values)
            ci = self._compute_ci(bootstrap_distribution, point, use_scipy_bca)

        # Build and return result
        return FleissKappaResult(
            point=point,
            confidence_interval=ci,
            balance_profile=self.balance_profile,
            bootstrap_distribution=bootstrap_distribution,
            ci_method=self.config.ci_method,
            weights_type=self.config.weights_type,
        )

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H₀: κ = 0 using the method specified in config.hypothesis_test_config.

        Two methods are supported:
        - "parametric_null": uses Fleiss (1979) null variance formula.
          Works with both Ratings and NijMatrix.
        - "permutation": builds null distribution by shuffling ratings.
          Requires Ratings data (not NijMatrix).

        Returns
        -------
        HypothesisTestResult
            Test statistic, p-value, and metadata.

        Raises
        ------
        ValueError
            If config.hypothesis_test_config is None.
            If method="permutation" but ratings is None.
        """
        from scipy import stats
        from ..intermediates.analysis_results import HypothesisTestResult

        # Guard: hypothesis_test_config must be set
        if self.config.hypothesis_test_config is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis_test_config is None. "
                f"Set hypothesis_test_config in your config to use hypothesis testing."
            )

        htc = self.config.hypothesis_test_config

        # Compute point estimate from whichever data source is available
        nij = self._get_nij()
        point = fleiss_kappa(nij, weights=self.weights)
        pj = self._compute_pj(nij)

        if htc.method == "parametric_null":
            # Fleiss (1979) Eq. 12: null variance, Eq. 17: z-test with bias correction
            from .variance_functions import fleiss_null_variance

            n_subjects = nij.n_subjects
            n_raters = int(nij.n_per_subject[0])  # Assumes balanced

            var_null = fleiss_null_variance(pj, n_subjects, n_raters)
            se_null = np.sqrt(var_null)

            # Bias correction from Fleiss (1979) Eq. 17:
            # Kappa is slightly negatively biased by -1/(N(n-1)).
            # Adding 1/(N(n-1)) to the observed kappa before computing z
            # corrects for this bias.
            bias_correction = 1.0 / (n_subjects * (n_raters - 1))
            z = (point + bias_correction) / se_null

            p_value = 2.0 * (1.0 - stats.norm.cdf(abs(z)))

            return HypothesisTestResult(
                p_value=p_value,
                alpha=htc.alpha,
                reject_null=p_value < htc.alpha,
                method="parametric_null",
                point=point,
                z_statistic=z,
                null_variance=var_null,
            )

        elif htc.method == "permutation":
            # Guard: permutation requires Ratings
            self._check_permutation_available()

            # Select permutation scheme
            from .resampling_functions.permutation_resampling import (
                column_shuffle,
                within_subject_shuffle,
                pool_and_redistribute,
            )

            scheme_map = {
                "column_shuffle": column_shuffle,
                "within_subject": within_subject_shuffle,
                "pool_and_redistribute": pool_and_redistribute,
            }

            scheme_fn = scheme_map.get(htc.scheme)
            if scheme_fn is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: unknown permutation scheme '{htc.scheme}'. "
                    f"Expected one of: {set(scheme_map.keys())}."
                )

            # Build null distribution
            n_categories = nij.n_categories
            null_values = []

            for permuted_array in scheme_fn(
                self.ratings.category_index,
                n_permutations=htc.n_permutations,
                random_state=htc.random_state,
            ):
                nij_array_perm = category_counts(permuted_array, n_categories)
                nij_permuted = NijMatrix(matrix=nij_array_perm)
                kappa_perm = fleiss_kappa(nij_permuted, weights=self.weights)
                null_values.append(kappa_perm)

            null_distribution = np.array(null_values)

            # p-value: proportion of null values >= observed (two-sided)
            p_value = float(np.mean(np.abs(null_distribution) >= abs(point)))

            return HypothesisTestResult(
                p_value=p_value,
                alpha=htc.alpha,
                reject_null=p_value < htc.alpha,
                method="permutation",
                point=point,
                null_distribution=null_distribution,
                n_permutations=htc.n_permutations,
            )

        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown hypothesis test method '{htc.method}'. "
                f"Expected 'parametric_null' or 'permutation'."
            )

    # Private helpers -------------------------------------------------------- #

    def _validate_fleiss_config(self) -> None:
        """
        Validate that self.config is a FleissKappaConfig.

        Raises
        ------
        TypeError
            If config is not FleissKappaConfig.
        """
        if not isinstance(self.config, FleissKappaConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be FleissKappaConfig, "
                f"got {type(self.config).__name__}"
            )

    def _get_nij(self) -> NijMatrix:
        """
        Get the NijMatrix from whichever data source is available.

        Returns NijMatrix from self.nij_matrix if available,
        otherwise computes it from self.ratings.

        Returns
        -------
        NijMatrix

        Raises
        ------
        ValueError
            If neither ratings nor nij_matrix is available.
        """
        if self.nij_matrix is not None:
            return self.nij_matrix
        elif self.ratings is not None:
            n_categories = self.ratings.semantics.n_categories
            nij_array = category_counts(self.ratings.category_index, n_categories)
            return NijMatrix(matrix=nij_array)
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: neither ratings nor nij_matrix is available."
            )

    def _compute_pj(self, nij: NijMatrix) -> np.ndarray:
        """
        Compute category proportions (pj) from NijMatrix.

        Each subject contributes equally: pj is the mean of per-subject
        proportions, not the global pooled proportion.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix.

        Returns
        -------
        np.ndarray
            Category proportions, shape (n_categories,), summing to 1.
        """
        n_per_subject = nij.n_per_subject
        proportions_per_subject = nij.matrix / n_per_subject[:, np.newaxis]
        return np.mean(proportions_per_subject, axis=0)

    def _parametric_null_ci(
        self, point: float, pj: np.ndarray
    ) -> tuple[float, float]:
        """
        Compute confidence interval using Fleiss (1979) null variance.

        Parameters
        ----------
        point : float
            The point estimate of kappa.
        pj : np.ndarray
            Category proportions.

        Returns
        -------
        tuple[float, float]
            (lower, upper) confidence interval bounds.
        """
        from scipy import stats
        from .variance_functions import fleiss_null_variance

        nij = self._get_nij()
        n_subjects = nij.n_subjects
        n_raters = int(nij.n_per_subject[0])  # Assumes balanced

        var_null = fleiss_null_variance(pj, n_subjects, n_raters)
        se_null = np.sqrt(var_null)
        z_crit = stats.norm.ppf(1.0 - self.config.alpha / 2.0)

        return (point - z_crit * se_null, point + z_crit * se_null)

    def _compute_ci(
        self,
        bootstrap_distribution: np.ndarray,
        point: float,
        use_scipy_bca: bool = False,
    ) -> tuple[float, float]:
        """
        Compute confidence interval from bootstrap distribution.

        Parameters
        ----------
        bootstrap_distribution : np.ndarray
            Array of bootstrap coefficient values.
        point : float
            The point estimate of the coefficient.
        use_scipy_bca : bool, default False
            If True and ci_method == "bca", use scipy for BCa.

        Returns
        -------
        tuple[float, float]
            (lower, upper) confidence interval bounds.

        Raises
        ------
        ValueError
            If ci_method is unknown.
        """
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
                f"Expected 'percentile', 'bca', or 'parametric_null'."
            )
