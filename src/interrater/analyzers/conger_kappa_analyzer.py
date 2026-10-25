"""Analyzer for Conger's kappa coefficient."""

from __future__ import annotations

from typing import Callable

import numpy as np

from .individual_marginal_analyzer import IndividualMarginalAnalyzer

from ..data_objects.ratings import Ratings
from ..intermediates.results import CongerKappaResult, HypothesisTestResult
from ..intermediates.data_profiles.balance_profile import BalanceProfile
from .configs import CongerKappaConfig

from ..intermediates.estimator_parts.estimator_input import EstimatorInput
from .agreement_functions.individual_marginals import conger_kappa
from .resampling_functions.subject_resampling import subject_bootstrap
from .ci_methods import percentile_ci, bca_ci
from ..intermediates.results.ci_results import BootstrapCIResult

_ERROR_PREFIX = "[CongerKappaAnalyzer] Error"


class CongerKappaAnalyzer(IndividualMarginalAnalyzer):
    """
    I am the analyzer for Conger's kappa, the multi-rater agreement coefficient
    with individual marginals.

    Conger's kappa (Conger 1980) generalises Cohen's kappa to any number of
    raters (n >= 2) while preserving each rater's own category distribution in
    the chance term. Where Fleiss' kappa pools marginal distributions across
    raters, Conger's kappa subtracts the variance of those distributions from
    the pooled chance term -- the same principle as Cohen's, extended beyond a
    single pair.

    At exactly two raters, Conger's kappa reduces to Cohen's unweighted kappa.

    I validate that:
    - Config is CongerKappaConfig

    I compute:
    - Point estimate of Conger's kappa (always unweighted)
    - Bootstrap confidence interval (percentile or BCa -- no parametric null)
    - Return CongerKappaResult with all metadata

    There is no NijMatrix path. Conger's kappa reads rater identity -- who said
    what -- and counts discard that. If you have only counts, you want Fleiss'
    kappa, not Conger's.

    There is no weighted form. Conger (1980) defined the coefficient for nominal
    data only; the unweighted constraint is enforced by my config
    (NoWeightIndividualConfig) and by the agreement function (chance_conger
    rejects non-identity weights).

    There is no parametric null variance. Confidence intervals and hypothesis
    tests use bootstrap and permutation methods respectively.

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer, always required).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    config : CongerKappaConfig
        The analysis configuration for Conger's kappa.
    estimator_input : EstimatorInput
        Built from ratings, carrying category_index, n_categories, and
        weights=None (from IndividualMarginalAnalyzer).

    References
    ----------
    Conger, A. J. (1980). Integration and generalization of kappas for
        multiple raters. Psychological Bulletin, 88(2), 322--328.
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    # Populated at parent class IndividualMarginalAnalyzer
    estimator_input: EstimatorInput
    # Checked at this level
    config: CongerKappaConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(self, ratings: Ratings, config: CongerKappaConfig) -> None:
        """
        Initialize Conger's kappa analyzer: validate config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze. Must contain at least two raters.
        config : CongerKappaConfig
            The analysis configuration for Conger's kappa.

        Raises
        ------
        TypeError
            If config is not CongerKappaConfig.
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
        self._validate_conger_config()

    # Public API ------------------------------------------------------------- #

    def estimate(self) -> CongerKappaResult:
        """
        Compute Conger's kappa and return the result.

        Performs point estimation and bootstrap confidence interval computation.
        Only "percentile" and "bca" CI methods are supported -- there is no
        closed-form null variance for Conger's kappa.

        Returns
        -------
        CongerKappaResult
        """
        # Guard convergence_test (not yet implemented)
        if self.config.convergence_test is not None:
            raise NotImplementedError(
                f"{_ERROR_PREFIX}: convergence_test not yet implemented"
            )

        # Compute point estimate from original data
        point = conger_kappa(self.estimator_input)

        bootstrap_distribution = self._bootstrap_distribution()
        ci = self._compute_ci(
            bootstrap_distribution, point, self._jackknife_estimates
        )

        # Build and return result
        return CongerKappaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            ci_method=self.config.ci_method,
        )

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: kappa = 0 using permutation.

        Only "permutation" is supported -- there is no closed-form null
        variance for Conger's kappa.

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

        ht_config = self.config.hypothesis_test_config
        n_categories = self.ratings.semantics.n_categories

        # Compute point estimate
        point = conger_kappa(self.estimator_input)

        if ht_config.method == "permutation":
            # Permutation test: shuffle ratings to build null distribution
            from .resampling_functions.permutation_resampling import (
                within_rater_shuffle,
                within_subject_shuffle,
                pool_and_redistribute,
            )

            scheme_map = {
                "within_rater_shuffle": within_rater_shuffle,
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
                )
                kappa_perm = conger_kappa(ei)
                null_values.append(kappa_perm)

            null_distribution = np.array(null_values)

            counts = self._permutation_p_value(
                null_distribution, point, ht_config.alternative
            )

            return HypothesisTestResult(
                alpha=ht_config.alpha,
                reject_null=self._reject_null(counts["p_value"], ht_config.alpha),
                method="permutation",
                coefficient="conger_kappa",
                point=point,
                null_distribution=null_distribution,
                n_permutations=ht_config.n_permutations,
                **counts,
            )

        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown hypothesis test method "
                f"'{ht_config.method}'. Conger's kappa has no closed-form null "
                f"variance, so only 'permutation' is available. "
                f"Got '{ht_config.method}'."
            )

    # Private helpers -------------------------------------------------------- #

    def _validate_conger_config(self) -> None:
        """
        Validate that self.config is a CongerKappaConfig.

        Raises
        ------
        TypeError
            If config is not CongerKappaConfig.
        """
        if not isinstance(self.config, CongerKappaConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be CongerKappaConfig, "
                f"got {type(self.config).__name__}"
            )

    def _estimate_from_array(self, category_index: np.ndarray) -> float:
        """
        Conger's kappa from a ratings array.

        The only coefficient-specific step in my resampling. Everything that
        draws samples is written once on Analyzer in terms of this.

        I build my own EstimatorInput because the array I am handed is a
        resample, not the original data: Conger's chance term is computed
        from each rater's own marginals, so the column positions matter.
        """
        ei = EstimatorInput(
            category_index=category_index,
            n_categories=self.ratings.semantics.n_categories,
        )
        return conger_kappa(ei)

    def _compute_ci(
        self,
        bootstrap_distribution: np.ndarray,
        point: float,
        jackknife_fn: Callable[[], np.ndarray],
    ) -> BootstrapCIResult:
        """
        Compute confidence interval from bootstrap distribution.

        Only "percentile" and "bca" are supported. There is no closed-form
        null variance for Conger's kappa, so "parametric_null" is not
        available.

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
                f"'{self.config.ci_method}'. Conger's kappa has no "
                f"closed-form null variance, so only 'percentile' and "
                f"'bca' are available."
            )
