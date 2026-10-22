"""Analyzer for Cohen's kappa coefficient."""

from __future__ import annotations

from typing import Callable

import numpy as np

from .weightable_individual_analyzer import WeightableIndividualAnalyzer

from ..data_objects.ratings import Ratings
from ..intermediates.results import CohenKappaResult, HypothesisTestResult
from ..intermediates.data_profiles.balance_profile import BalanceProfile
from .configs import CohenKappaConfig

from ..intermediates.estimator_parts import Weights
from ..intermediates.estimator_parts.estimator_input import EstimatorInput
from .agreement_functions.individual_marginals import (
    cohen_kappa,
    rater_marginals,
    shared_rows,
)
from .resampling_functions.subject_resampling import subject_bootstrap
from .ci_methods import percentile_ci, bca_ci
from ..intermediates.results.ci_results import (
    BootstrapCIResult,
    ParametricCIResult,
)

_ERROR_PREFIX = "[CohenKappaAnalyzer] Error"


class CohenKappaAnalyzer(WeightableIndividualAnalyzer):
    """
    I am the analyzer for Cohen's kappa, the pairwise agreement coefficient
    with individual marginals.

    Cohen's kappa (Cohen 1960) measures agreement between exactly two raters,
    correcting for chance using each rater's own category distribution. The
    weighted form (Cohen 1968) gives partial credit for near misses on ordinal
    data.

    I validate that:
    - Config is CohenKappaConfig
    - Data contains exactly two raters (n_raters == 2)

    I compute:
    - Point estimate of Cohen's (weighted) kappa
    - Confidence interval (bootstrap or parametric null)
    - Return CohenKappaResult with all metadata

    There is no NijMatrix path. Cohen's kappa reads rater identity -- who said
    what -- and counts discard that. If you have only counts, you want Fleiss'
    kappa or Scott's pi, not Cohen's.

    There is no balance_profile or is_exact distinction. Cohen operates on the
    subjects both raters rated (shared_rows), so the formula is always exact
    on the pair it was given.

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer, always required).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    config : CohenKappaConfig
        The analysis configuration for Cohen's kappa.
    estimator_input : EstimatorInput
        Built from ratings, carrying category_index, n_categories, and weights.
    weights : Weights | None
        The weights declared by the semantics (from WeightableIndividualAnalyzer).
        None for unordered categories (unweighted).

    References
    ----------
    Cohen, J. (1960). A coefficient of agreement for nominal scales.
        Educational and Psychological Measurement, 20(1), 37--46.
    Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision
        for scaled disagreement or partial credit. Psychological Bulletin,
        70(4), 213--220.
    Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). Large sample
        standard errors of kappa and weighted kappa. Psychological
        Bulletin, 72(5), 323--327.
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    # Populated at parent class IndividualMarginalAnalyzer
    estimator_input: EstimatorInput
    # Populated at parent class WeightableIndividualAnalyzer
    weights: Weights | None
    # Checked at this level
    config: CohenKappaConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(self, ratings: Ratings, config: CohenKappaConfig) -> None:
        """
        Initialize Cohen's kappa analyzer: validate config and rater count.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze. Must contain exactly two raters.
        config : CohenKappaConfig
            The analysis configuration for Cohen's kappa.

        Raises
        ------
        TypeError
            If config is not CohenKappaConfig.
        ValueError
            If data does not contain exactly two raters.
        """
        # Call parent __init__ (validates ratings, config, computes BalanceProfile,
        # builds EstimatorInput, resolves weights)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate config type and rater count
        self._validate_cohen_config()
        self._validate_n_raters()

    # Public API ------------------------------------------------------------- #

    def estimate(self) -> CohenKappaResult:
        """
        Compute Cohen's kappa and return the result.

        Performs point estimation and confidence interval computation.

        When ci_method is "percentile" or "bca", runs the bootstrap loop.
        When ci_method is "parametric_null", uses Fleiss, Cohen & Everitt (1969)
        null variance -- no bootstrap is performed and bootstrap_distribution
        is None.

        Returns
        -------
        CohenKappaResult
        """
        # Guard convergence_test (not yet implemented)
        if self.config.convergence_test is not None:
            raise NotImplementedError(
                f"{_ERROR_PREFIX}: convergence_test not yet implemented"
            )

        n_categories = self.ratings.semantics.n_categories
        n_subjects = self._get_n_subjects()

        # Compute point estimate from original data
        point = cohen_kappa(self.estimator_input, self.weights)

        # CI path depends on ci_method
        if self.config.ci_method == "parametric_null":
            # Parametric null CI: no bootstrap needed
            p_a, p_b = self._compute_marginals()
            ci = self._parametric_null_ci(point, p_a, p_b, n_subjects)
            bootstrap_distribution = None
        else:
            bootstrap_distribution = self._bootstrap_distribution()
            ci = self._compute_ci(
                bootstrap_distribution, point, self._jackknife_estimates
            )

        # Build and return result
        return CohenKappaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            ci_method=self.config.ci_method,
            weights_type=self._weights_kind(),
            weight_matrix=self._weight_matrix(),
        )

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: kappa_w = 0 using the method specified in config.

        Two methods are supported:
        - "parametric_null": uses Fleiss, Cohen & Everitt (1969) null variance.
          No bias correction -- the Fleiss (1979) bias correction +1/(N(n-1))
          was derived for the multi-rater case and does not apply here.
        - "permutation": builds null distribution by shuffling ratings.

        Returns
        -------
        HypothesisTestResult
            Test statistic, p-value, and metadata.

        Raises
        ------
        ValueError
            If config.hypothesis_test_config is None.
        """
        from scipy import stats

        # Guard: hypothesis_test_config must be set
        if self.config.hypothesis_test_config is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis_test_config is None. "
                f"Set hypothesis_test_config in your config to use "
                f"hypothesis testing."
            )

        ht_config = self.config.hypothesis_test_config
        n_categories = self.ratings.semantics.n_categories

        # Compute point estimate
        point = cohen_kappa(self.estimator_input, self.weights)

        if ht_config.method == "parametric_null":
            # Fleiss, Cohen & Everitt (1969): null variance for two raters
            from .variance_functions import cohen_null_variance

            p_a, p_b = self._compute_marginals()
            n_subjects = self._get_n_subjects()
            weight_matrix = (
                self.weights.matrix if self.weights is not None else None
            )

            var_null = cohen_null_variance(
                p_a, p_b, n_subjects, weight_matrix
            )
            se_null = np.sqrt(var_null)

            # No bias correction for two-rater case.
            # The +1/(N(n-1)) correction in Fleiss (1979) Eq. 17 was derived
            # for multi-rater kappa and does not apply to Cohen's kappa.
            z = point / se_null

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
            # Permutation test: shuffle ratings to build null distribution
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
                    f"{_ERROR_PREFIX}: unknown permutation scheme "
                    f"'{ht_config.scheme}'. Expected one of: "
                    f"{set(scheme_map.keys())}."
                )

            # Build null distribution
            null_values = []

            for permuted_array in scheme_fn(
                self.ratings.category_index,
                n_permutations=ht_config.n_permutations,
                random_state=ht_config.random_state,
            ):
                ei = EstimatorInput(
                    category_index=permuted_array,
                    n_categories=n_categories,
                    weights=self.weights,
                )
                kappa_perm = cohen_kappa(ei, self.weights)
                null_values.append(kappa_perm)

            null_distribution = np.array(null_values)

            # p-value: proportion of null values >= observed (two-sided)
            p_value = float(
                np.mean(np.abs(null_distribution) >= abs(point))
            )

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
                f"{_ERROR_PREFIX}: unknown hypothesis test method "
                f"'{ht_config.method}'. Expected 'parametric_null' or "
                f"'permutation'."
            )

    # Private helpers -------------------------------------------------------- #

    def _validate_cohen_config(self) -> None:
        """
        Validate that self.config is a CohenKappaConfig.

        Raises
        ------
        TypeError
            If config is not CohenKappaConfig.
        """
        if not isinstance(self.config, CohenKappaConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be CohenKappaConfig, "
                f"got {type(self.config).__name__}"
            )

    def _validate_n_raters(self) -> None:
        """
        Validate that data contains exactly two raters.

        Cohen's kappa is defined for a single pair. The agreement function
        in individual_marginals.py assumes two columns and does not check;
        validation is my responsibility.

        Raises
        ------
        ValueError
            If the number of rater columns is not 2.
        """
        n_raters = self.ratings.category_index.shape[1]
        if n_raters != 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: Cohen's kappa requires exactly 2 raters, "
                f"got {n_raters}. For more than 2 raters, use "
                f"CongerKappaAnalyzer (individual marginals) or "
                f"FleissKappaAnalyzer (pooled marginals)."
            )

    def _get_n_subjects(self) -> int:
        """
        Count subjects both raters rated (the effective sample size).

        This is the N in the variance formula and in the result. Subjects
        only one rater saw contribute to neither agreement nor disagreement,
        so they are not counted.

        Returns
        -------
        int
            Number of shared subjects.
        """
        return int(shared_rows(self.estimator_input, 0, 1).size)

    def _compute_marginals(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute each rater's marginal category proportions on shared subjects.

        Marginals are computed on the subjects both raters rated, not on all
        subjects. A rater who was strict on the shared subjects and lenient
        elsewhere would otherwise bring the wrong distribution.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (p_a, p_b) -- each shape (K,), each summing to 1.
        """
        rows = shared_rows(self.estimator_input, 0, 1)
        marginals = rater_marginals(self.estimator_input, rows)
        return marginals[0], marginals[1]

    def _estimate_from_array(self, category_index: np.ndarray) -> float:
        """
        Cohen's kappa from a ratings array.

        The only coefficient-specific step in my resampling. Everything that
        draws samples is written once on Analyzer in terms of this.

        I build my own EstimatorInput because the array I am handed is a
        resample, not the original data: rater identity matters to me, and it
        is the column positions that carry it.
        """
        ei = EstimatorInput(
            category_index=category_index,
            n_categories=self.ratings.semantics.n_categories,
            weights=self.weights,
        )
        return cohen_kappa(ei, self.weights)

    def _parametric_null_ci(
        self,
        point: float,
        p_a: np.ndarray,
        p_b: np.ndarray,
        n_subjects: int,
    ) -> ParametricCIResult:
        """
        Compute confidence interval using Fleiss-Cohen-Everitt (1969) null
        variance.

        Parameters
        ----------
        point : float
            The point estimate of kappa.
        p_a : np.ndarray
            Rater A's marginal proportions.
        p_b : np.ndarray
            Rater B's marginal proportions.
        n_subjects : int
            Number of subjects both raters rated.

        Returns
        -------
        ParametricCIResult
            Bounds, the level, and the variance they were built from.
        """
        from scipy import stats
        from .variance_functions import cohen_null_variance

        weight_matrix = (
            self.weights.matrix if self.weights is not None else None
        )

        var_null = cohen_null_variance(p_a, p_b, n_subjects, weight_matrix)
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
        Compute confidence interval from bootstrap distribution.

        Parameters
        ----------
        bootstrap_distribution : np.ndarray
            Array of bootstrap kappa values.
        point : float
            The point estimate of kappa.
        jackknife_fn : Callable[[], np.ndarray]
            Produces leave-one-subject-out estimates when called. A callable
            rather than an array because only BCa needs them, and producing
            them costs one coefficient evaluation per subject.

        Returns
        -------
        BootstrapCIResult
            Bounds, the method asked for against the one that ran, and how
            many resamples failed.

        Raises
        ------
        ValueError
            If ci_method is unknown.
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
                f"{_ERROR_PREFIX}: unknown ci_method "
                f"'{self.config.ci_method}'. Expected 'percentile', "
                f"'bca', or 'parametric_null'."
            )
