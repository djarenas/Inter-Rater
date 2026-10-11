"""Analyzer for Light's kappa coefficient."""

from __future__ import annotations

from math import comb

import numpy as np

from .individual_marginal_analyzer import IndividualMarginalAnalyzer

from ..data_objects.ratings import Ratings
from ..intermediates.analysis_results import LightKappaResult, HypothesisTestResult
from ..intermediates.data_profiles.balance_profile import BalanceProfile
from .configs import LightKappaConfig

from ..intermediates.estimator_parts.estimator_input import EstimatorInput
from .agreement_functions.individual_marginals import light_kappa
from .resampling_functions.subject_resampling import subject_bootstrap
from .ci_methods import percentile_ci, bca_ci, scipy_bca_ci

_ERROR_PREFIX = "[LightKappaAnalyzer] Error"


class LightKappaAnalyzer(IndividualMarginalAnalyzer):
    """
    I am the analyzer for Light's kappa, the mean of all C(n, 2) pairwise
    Cohen's kappas.

    Light's kappa (Light 1971) keeps rater identity through Cohen's kappa:
    each pair uses its own marginals. The mean then compresses the C(n, 2)
    values into a single summary. For per-pair detail, spread, and per-rater
    diagnostics, use the pairwise permutation analyzer -- Light's kappa is
    the headline number, not the diagnostic.

    At exactly two raters (n_pairs=1), Light's kappa reduces to Cohen's
    unweighted kappa (there is only one pair).

    I validate that:
    - Config is LightKappaConfig

    I compute:
    - Point estimate of Light's kappa (always unweighted)
    - Bootstrap confidence interval (percentile or BCa -- no parametric null)
    - Return LightKappaResult with all metadata including n_pairs

    There is no NijMatrix path. Light's kappa reads rater identity -- who said
    what -- and counts discard that.

    There is no weighted form in this package. Light (1971) defined the coefficient as the
    mean of unweighted pairwise kappas.

    There is no parametric null variance. The C(n, 2) pairwise kappas are
    correlated (pairs sharing a rater are not independent), so confidence
    intervals and hypothesis tests use bootstrap and permutation methods
    respectively.

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer, always required).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    config : LightKappaConfig
        The analysis configuration for Light's kappa.
    estimator_input : EstimatorInput
        Built from ratings, carrying category_index, n_categories, and
        weights=None (from IndividualMarginalAnalyzer).

    References
    ----------
    Light, R. J. (1971). Measures of response agreement for qualitative
        data: Some generalizations and alternatives. Psychological
        Bulletin, 76(5), 365--377.
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    # Populated at parent class IndividualMarginalAnalyzer
    estimator_input: EstimatorInput
    # Checked at this level
    config: LightKappaConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(self, ratings: Ratings, config: LightKappaConfig) -> None:
        """
        Initialize Light's kappa analyzer: validate config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze. Must contain at least two raters.
        config : LightKappaConfig
            The analysis configuration for Light's kappa.

        Raises
        ------
        TypeError
            If config is not LightKappaConfig.
        ValueError
            If ratings is None.
        """
        # Call parent __init__ (validates ratings not None, validates config is
        # IndividualMarginalConfig, computes BalanceProfile, builds
        # EstimatorInput with weights=None)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate config type
        self._validate_light_config()

    # Public API ------------------------------------------------------------- #

    def estimate(self, use_scipy_bca: bool = False) -> LightKappaResult:
        """
        Compute Light's kappa and return the result.

        Performs point estimation and bootstrap confidence interval computation.
        Only "percentile" and "bca" CI methods are supported -- there is no
        closed-form null variance for Light's kappa.

        Parameters
        ----------
        use_scipy_bca : bool, default False
            If True and config.ci_method == "bca", use scipy.stats for BCa
            computation (for testing/validation only).

        Returns
        -------
        LightKappaResult
        """
        # Guard convergence_test (not yet implemented)
        if self.config.convergence_test is not None:
            raise NotImplementedError(
                f"{_ERROR_PREFIX}: convergence_test not yet implemented"
            )

        n_categories = self.ratings.semantics.n_categories

        # Compute point estimate from original data
        point = light_kappa(self.estimator_input)

        # Bootstrap CI: resample subjects, recompute kappa each time
        bootstrap_values = []

        for resampled_array in subject_bootstrap(
            self.ratings.category_index,
            n_resamples=self.config.n_resamples,
            random_state=self.config.random_state,
        ):
            ei = EstimatorInput(
                category_index=resampled_array,
                n_categories=n_categories,
            )
            kappa_resample = light_kappa(ei)
            bootstrap_values.append(kappa_resample)

        bootstrap_distribution = np.array(bootstrap_values)
        ci = self._compute_ci(bootstrap_distribution, point, use_scipy_bca)

        # Build and return result
        n_pairs = comb(self.balance_profile.n_raters, 2)

        return LightKappaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            n_pairs=n_pairs,
            ci_method=self.config.ci_method,
        )

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: kappa = 0 using permutation.

        Only "permutation" is supported -- there is no closed-form null
        variance for Light's kappa.

        Returns
        -------
        HypothesisTestResult
            Test statistic, p-value, and metadata.

        Raises
        ------
        ValueError
            If config.hypothesis_test_config is None.
            If method is not "permutation".
        """
        # Guard: hypothesis_test_config must be set
        if self.config.hypothesis_test_config is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis_test_config is None. "
                f"Set hypothesis_test_config in your config to use "
                f"hypothesis testing."
            )

        htc = self.config.hypothesis_test_config
        n_categories = self.ratings.semantics.n_categories

        # Compute point estimate
        point = light_kappa(self.estimator_input)

        if htc.method == "permutation":
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

            scheme_fn = scheme_map.get(htc.scheme)
            if scheme_fn is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: unknown permutation scheme "
                    f"'{htc.scheme}'. Expected one of: "
                    f"{set(scheme_map.keys())}."
                )

            # Build null distribution
            null_values = []

            for permuted_array in scheme_fn(
                self.ratings.category_index,
                n_permutations=htc.n_permutations,
                random_state=htc.random_state,
            ):
                ei = EstimatorInput(
                    category_index=permuted_array,
                    n_categories=n_categories,
                )
                kappa_perm = light_kappa(ei)
                null_values.append(kappa_perm)

            null_distribution = np.array(null_values)

            # p-value: proportion of null values >= observed (two-sided)
            p_value = float(
                np.mean(np.abs(null_distribution) >= abs(point))
            )

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
                f"{_ERROR_PREFIX}: unknown hypothesis test method "
                f"'{htc.method}'. Light's kappa has no closed-form null "
                f"variance, so only 'permutation' is available. "
                f"Got '{htc.method}'."
            )

    # Private helpers -------------------------------------------------------- #

    def _validate_light_config(self) -> None:
        """
        Validate that self.config is a LightKappaConfig.

        Raises
        ------
        TypeError
            If config is not LightKappaConfig.
        """
        if not isinstance(self.config, LightKappaConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be LightKappaConfig, "
                f"got {type(self.config).__name__}"
            )

    def _compute_ci(
        self,
        bootstrap_distribution: np.ndarray,
        point: float,
        use_scipy_bca: bool = False,
    ) -> tuple[float, float]:
        """
        Compute confidence interval from bootstrap distribution.

        Only "percentile" and "bca" are supported. There is no closed-form
        null variance for Light's kappa, so "parametric_null" is not
        available.

        Parameters
        ----------
        bootstrap_distribution : np.ndarray
            Array of bootstrap kappa values.
        point : float
            The point estimate of kappa.
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
                return scipy_bca_ci(
                    bootstrap_distribution, point, self.config.alpha
                )
            else:
                return bca_ci(
                    bootstrap_distribution, point, self.config.alpha
                )
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown ci_method "
                f"'{self.config.ci_method}'. Light's kappa has no "
                f"closed-form null variance, so only 'percentile' and "
                f"'bca' are available."
            )
