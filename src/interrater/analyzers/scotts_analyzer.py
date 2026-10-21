"""Analyzer for Scott's pi coefficient."""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from .group_analyzer import GroupAnalyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.results import ScottsPiResult, HypothesisTestResult
    from .configs import ScottsConfig
    from ..intermediates.data_profiles.balance_profile import BalanceProfile

from ..intermediates.estimator_parts import NijMatrix
from .agreement_functions.group_agreement import scotts_pi
from .ci_methods import percentile_ci, bca_ci
from ..intermediates.results.ci_results import (
    BootstrapCIResult,
    ParametricCIResult,
)
from ..array_utils import category_counts

_ERROR_PREFIX = "[ScottsAnalyzer] Error"


class ScottsAnalyzer(GroupAnalyzer):
    """
    I am the analyzer for Scott's pi coefficient.

    Scott's pi (Scott 1955) measures agreement between exactly two raters on
    nominal (unweighted) categorical data, corrected for chance using pooled
    marginal proportions.

    At exactly two raters, Fleiss' kappa equals Scott's pi -- they are
    algebraically identical. This class exists to honor Scott's priority,
    provide a clear two-rater-only entry point, and enforce the nominal-data
    constraint. For more than two raters, use FleissKappaAnalyzer.

    I validate that:
    - Config is ScottsConfig
    - Ratings have exactly two raters (or nij_matrix has row sums of 2)
    - Data structure is appropriate for Scott's pi

    I compute:
    - Point estimate of Scott's pi
    - Bootstrap or parametric null confidence interval
    - Hypothesis test (permutation or parametric null)
    - Return ScottsPiResult with all metadata

    Attributes
    ----------
    ratings : Ratings | None
        The validated ratings data (from Analyzer, can be None if from nij_matrix).
    balance_profile : BalanceProfile
        Describes the balance structure of the data (from Analyzer).
    config : ScottsConfig
        The analysis configuration for Scott's pi.
    nij_matrix : NijMatrix | None
        The counts matrix (from GroupAnalyzer, populated if built from nij_matrix).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings | None
    balance_profile: BalanceProfile
    # Populated at parent class GroupAnalyzer
    nij_matrix: NijMatrix | None
    # Checked at this level
    config: ScottsConfig

    # Constructors ----------------------------------------------------------- #

    def __init__(self, ratings: Ratings, config: ScottsConfig) -> None:
        """
        Initialize Scott's pi analyzer: validate config and ratings type.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze. Must have exactly 2 raters.
        config : ScottsConfig
            The analysis configuration for Scott's pi.

        Raises
        ------
        TypeError
            If config is not ScottsConfig.
            If ratings is not a Ratings object.
        ValueError
            If ratings do not have exactly 2 raters.
        """
        from ..data_objects.ratings import Ratings as RatingsBase

        # Validate ratings is at least a Ratings object
        if not isinstance(ratings, RatingsBase):
            raise TypeError(
                f"{_ERROR_PREFIX}: ratings must be a Ratings object, "
                f"got {type(ratings).__name__}"
            )

        # Validate exactly 2 raters (Scott's pi is two-rater only)
        self._validate_two_raters_from_ratings(ratings)

        # Call parent __init__ (validates group config, balance constraint)
        # Wrap to propagate errors with ScottsAnalyzer context
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type
        self._validate_scotts_config()

    @classmethod
    def build_from_nij_matrix(
        cls, nij: NijMatrix, config: ScottsConfig
    ) -> ScottsAnalyzer:
        """
        Build Scott's pi analyzer from NijMatrix directly.

        Calls parent build_from_nij_matrix() to set ratings, config,
        balance_profile, and nij_matrix, then validates two-rater constraint.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories). Row sums must be 2.
        config : ScottsConfig
            The analysis configuration for Scott's pi.

        Returns
        -------
        ScottsAnalyzer
            An initialized analyzer ready to call estimate_from_nij_matrix().

        Raises
        ------
        TypeError
            If config is not ScottsConfig.
        ValueError
            If nij row sums are not all 2.
        """
        # Call parent build_from_nij_matrix (sets ratings, config, balance_profile, nij_matrix)
        # Wrap to propagate errors with ScottsAnalyzer context
        try:
            instance = super().build_from_nij_matrix(nij, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate two raters from nij row sums
        instance._validate_two_raters_from_nij(nij)

        # Validate config type
        instance._validate_scotts_config()

        return instance

    # Public API ------------------------------------------------------------- #

    def estimate(self) -> ScottsPiResult:
        """
        Compute Scott's pi and return the result.

        When ci_method is "percentile" or "bca", runs the bootstrap loop.
        When ci_method is "parametric_null", uses Fleiss (1979) null variance
        -- no bootstrap is performed and bootstrap_distribution is None.

        Returns
        -------
        ScottsPiResult
        """
        from ..intermediates.results import ScottsPiResult

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

        # Compute point estimate from original data
        nij_array = category_counts(
            self.ratings.category_index, self._n_categories()
        )
        point = self._estimate_from_nij(nij_array)

        # CI path depends on ci_method
        if self.config.ci_method == "parametric_null":
            pj = self._compute_pj(NijMatrix(matrix=nij_array))
            ci = self._parametric_null_ci(point, pj)
            bootstrap_distribution = None
        else:
            bootstrap_distribution = self._bootstrap_distribution()
            ci = self._compute_ci(
                bootstrap_distribution, point, self._jackknife_estimates
            )

        return ScottsPiResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            ci_method=self.config.ci_method
        )

    def estimate_from_nij_matrix(self) -> ScottsPiResult:
        """
        Compute Scott's pi from nij_matrix (NijMatrix-based path).

        Returns
        -------
        ScottsPiResult
        """
        from .resampling_functions import subject_bootstrap_from_nij_matrix
        from ..intermediates.results import ScottsPiResult

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
        point = self._estimate_from_nij(self.nij_matrix.matrix)

        # CI path depends on ci_method
        if self.config.ci_method == "parametric_null":
            pj = self._compute_pj(self.nij_matrix)
            ci = self._parametric_null_ci(point, pj)
            bootstrap_distribution = None
        else:
            nij_array = self.nij_matrix.matrix
            bootstrap_distribution = self._bootstrap_distribution_from_nij(
                nij_array
            )
            ci = self._compute_ci(
                bootstrap_distribution,
                point,
                lambda: self._jackknife_estimates_from_nij(nij_array),
            )

        return ScottsPiResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            ci_method=self.config.ci_method
        )

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H₀: π = 0 using the method specified in config.hypothesis_test_config.

        Scott's pi uses the same null variance as Fleiss' kappa (they are
        algebraically identical at two raters). The parametric_null method
        uses Fleiss (1979) Eq. 12 with bias correction from Eq. 17.

        Returns
        -------
        HypothesisTestResult
        """
        from scipy import stats
        from ..intermediates.results import HypothesisTestResult

        # Guard: hypothesis_test_config must be set
        if self.config.hypothesis_test_config is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis_test_config is None. "
                f"Set hypothesis_test_config in your config to use hypothesis testing."
            )

        ht_config = self.config.hypothesis_test_config

        # Compute point estimate from whichever data source is available
        nij = self._get_nij()
        point = scotts_pi(nij, weights=None)
        pj = self._compute_pj(nij)

        if ht_config.method == "parametric_null":
            # Scott's pi uses Fleiss' null variance (identical at n=2)
            from .variance_functions import fleiss_null_variance

            n_subjects = nij.n_subjects
            n_raters = int(nij.n_per_subject[0])

            var_null = fleiss_null_variance(pj, n_subjects, n_raters)
            se_null = np.sqrt(var_null)

            # Bias correction from Fleiss (1979) Eq. 17
            bias_correction = 1.0 / (n_subjects * (n_raters - 1))
            z = (point + bias_correction) / se_null
            p_value = 2.0 * (1.0 - stats.norm.cdf(abs(z)))

            return HypothesisTestResult(
                p_value=p_value,
                alpha=ht_config.alpha,
                reject_null=self._reject_null(p_value, ht_config.alpha),
                method="parametric_null",
                point=point,
                z_statistic=z,
                null_variance=var_null,
            )

        elif ht_config.method == "permutation":
            # Guard: permutation requires Ratings
            self._check_permutation_available()

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

            scheme_fn = scheme_map.get(ht_config.scheme)
            if scheme_fn is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: unknown permutation scheme '{ht_config.scheme}'. "
                    f"Expected one of: {set(scheme_map.keys())}."
                )

            # Build null distribution
            n_categories = nij.n_categories
            null_values = []

            for permuted_array in scheme_fn(
                self.ratings.category_index,
                n_permutations=ht_config.n_permutations,
                random_state=ht_config.random_state,
            ):
                nij_array_perm = category_counts(permuted_array, n_categories)
                nij_permuted = NijMatrix(matrix=nij_array_perm)
                pi_perm = scotts_pi(nij_permuted, weights=None)
                null_values.append(pi_perm)

            null_distribution = np.array(null_values)
            p_value = float(np.mean(np.abs(null_distribution) >= abs(point)))

            return HypothesisTestResult(
                p_value=p_value,
                alpha=ht_config.alpha,
                reject_null=self._reject_null(p_value, ht_config.alpha),
                method="permutation",
                point=point,
                null_distribution=null_distribution,
                n_permutations=ht_config.n_permutations,
            )

        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown hypothesis test method '{ht_config.method}'. "
                f"Expected 'parametric_null' or 'permutation'."
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
        """Validate exactly 2 raters from Ratings object."""
        if ratings.n_raters != 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: Scott's pi is defined for exactly 2 raters, "
                f"got {ratings.n_raters}. For {ratings.n_raters} raters, use "
                f"FleissKappaAnalyzer instead (Fleiss reduces to Scott's pi at n=2)."
            )

    def _validate_two_raters_from_nij(self, nij) -> None:
        """Validate exactly 2 raters from NijMatrix row sums."""
        n_per_subject = nij.n_per_subject
        if not np.all(n_per_subject == 2):
            raise ValueError(
                f"{_ERROR_PREFIX}: Scott's pi is defined for exactly 2 raters, "
                f"but nij row sums are not all 2 (got unique values: "
                f"{np.unique(n_per_subject).tolist()}). For more raters, use "
                f"FleissKappaAnalyzer instead."
            )

    def _get_nij(self) -> NijMatrix:
        """Get NijMatrix from whichever data source is available."""
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
        """Compute category proportions (pj) from NijMatrix."""
        n_per_subject = nij.n_per_subject
        proportions_per_subject = nij.matrix / n_per_subject[:, np.newaxis]
        return np.mean(proportions_per_subject, axis=0)

    def _estimate_from_nij(self, nij_array: np.ndarray) -> float:
        """
        Scott's pi from a counts matrix.

        The only coefficient-specific step in my resampling. Everything that
        draws samples -- bootstrap or jackknife, from ratings or from counts
        -- is written once on my base classes in terms of this.
        """
        return scotts_pi(NijMatrix(matrix=nij_array), weights=None)

    def _parametric_null_ci(
        self, point: float, pj: np.ndarray
    ) -> tuple[float, float]:
        """Compute CI using Fleiss (1979) null variance."""
        from scipy import stats
        from .variance_functions import fleiss_null_variance

        nij = self._get_nij()
        n_subjects = nij.n_subjects
        n_raters = int(nij.n_per_subject[0])

        var_null = fleiss_null_variance(pj, n_subjects, n_raters)
        se_null = np.sqrt(var_null)
        z_crit = stats.norm.ppf(1.0 - self.config.alpha / 2.0)

        return ParametricCIResult(
            lower=point - z_crit * se_null,
            upper=point + z_crit * se_null,
            alpha=self.config.alpha,
            variance_source="null",
        )

    def _compute_ci(
        self,
        bootstrap_distribution: np.ndarray,
        point: float,
        jackknife_fn: Callable[[], np.ndarray],
    ) -> BootstrapCIResult:
        """
        Compute CI from bootstrap distribution.

        ``jackknife_fn`` is a callable rather than an array because only BCa
        needs leave-one-out estimates, and producing them costs one
        coefficient evaluation per subject. I am the only place that knows
        which methods need them, so I am the only place that pays for them.
        """
        if self.config.ci_method == "percentile":
            return percentile_ci(bootstrap_distribution, self.config.alpha)
        elif self.config.ci_method == "bca":
            return bca_ci(
                bootstrap_distribution,
                point,
                self.config.alpha,
                jackknife_estimates=jackknife_fn(),
            )
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown ci_method '{self.config.ci_method}'. "
                f"Expected 'percentile', 'bca', or 'parametric_null'."
            )
