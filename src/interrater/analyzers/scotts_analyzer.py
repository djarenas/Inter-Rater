"""Analyzer for Scott's pi coefficient."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .group_analyzer import GroupAnalyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.results import ScottsPiResult, HypothesisTestResult
    from .configs import ScottsConfig
    from ..intermediates.data_profiles.balance_profile import BalanceProfile
    from ..intermediates.data_profiles.rater_profile import RaterProfile

from ..intermediates.estimator_parts import EstimatorInput, NijMatrix
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.scotts_pi_result import (
    KAPPA_MAX_NOT_COMPUTED,
)
from .agreement_functions.group_agreement import (
    chance_fleiss,
    observed_agreement,
    scotts_pi,
)
from ..intermediates.results.ci_results import ParametricCIResult

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
    - Ratings have exactly two raters
    - Data structure is appropriate for Scott's pi

    I compute:
    - Point estimate of Scott's pi
    - Bootstrap or parametric null confidence interval
    - Hypothesis test (permutation or parametric null)
    - Return ScottsPiResult with all metadata

    Attributes
    ----------
    ratings : Ratings
        The validated ratings data (from Analyzer).
    balance_profile : BalanceProfile
        Describes the balance structure of the data (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    config : ScottsConfig
        The analysis configuration for Scott's pi.
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
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
        self._refuse_parametric_null_if_invalid(allow_weights=False)

    # Public API ------------------------------------------------------------- #

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
        from ..intermediates.results import HypothesisTestResult

        # Guard: hypothesis_test_config must be set
        if self.config.hypothesis_test_config is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis_test_config is None. "
                f"Set hypothesis_test_config in your config to use hypothesis testing."
            )

        ht_config = self.config.hypothesis_test_config

        # Compute point estimate from whichever data source is available
        nij = self.estimator_input.nij
        point = scotts_pi(self.estimator_input)
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
            p_value = self._parametric_p_value(z, ht_config.alternative)

            return HypothesisTestResult(
                p_value=p_value,
                alpha=ht_config.alpha,
                reject_null=self._reject_null(p_value, ht_config.alpha),
                method="parametric_null",
                coefficient="scotts_pi",
                point=point,
                z_statistic=z,
                null_variance=var_null,
                alternative=ht_config.alternative,
            )

        elif ht_config.method == "permutation":
            from .resampling_functions.permutation_resampling import (
                within_rater_shuffle,
                pool_and_redistribute,
            )

            scheme_map = {
                "within_rater_shuffle": within_rater_shuffle,
                "pool_and_redistribute": pool_and_redistribute,
            }

            scheme_fn = scheme_map.get(ht_config.scheme)
            if scheme_fn is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: unknown permutation scheme '{ht_config.scheme}'. "
                    f"Expected one of: {set(scheme_map.keys())}."
                )

            # Build null distribution
            null_values = []

            for permuted_array in scheme_fn(
                self.estimator_input.category_index,
                n_permutations=ht_config.n_permutations,
                random_state=ht_config.random_state,
            ):
                null_values.append(self._estimate_from_input(
                    self._with_category_index(permuted_array)
                ))

            null_distribution = np.array(null_values)
            counts = self._permutation_p_value(
                null_distribution, point, ht_config.alternative
            )

            return HypothesisTestResult(
                alpha=ht_config.alpha,
                reject_null=self._reject_null(counts["p_value"], ht_config.alpha),
                method="permutation",
                coefficient="scotts_pi",
                point=point,
                null_distribution=null_distribution,
                n_permutations=ht_config.n_permutations,
                **counts,
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

    def _compute_pj(self, nij: NijMatrix) -> np.ndarray:
        """Compute category proportions (pj) from NijMatrix."""
        n_per_subject = nij.n_per_subject
        proportions_per_subject = nij.matrix / n_per_subject[:, np.newaxis]
        return np.mean(proportions_per_subject, axis=0)

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Scott's pi from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples, bootstrap or jackknife, is written once on my base
        classes in terms of this.
        """
        return scotts_pi(data)

    def _agreement_parts(self, nij: NijMatrix) -> AgreementParts:
        """
        Po and Pe on the counts the point estimate was computed from.

        The same two functions scotts_pi combines, so the result can check
        that they rebuild the point estimate.
        """
        return AgreementParts(
            observed=observed_agreement(nij, None),
            chance=chance_fleiss(nij, None),
            kind="agreement",
            weighted=False,
        )

    def _parametric_null_ci(
        self, point: float, pj: np.ndarray
    ) -> tuple[float, float]:
        """Compute CI using Fleiss (1979) null variance."""
        from scipy import stats
        from .variance_functions import fleiss_null_variance

        nij = self.estimator_input.nij
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

    def _parametric_null_interval(self, point: float) -> ParametricCIResult:
        """
        The interval from Fleiss (1979) null variance, for
        ``ci_method="parametric_null"``: no bootstrap, so no draws.
        """
        return self._parametric_null_ci(
            point, self._compute_pj(self.estimator_input.nij)
        )

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> ScottsPiResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.
        """
        from ..intermediates.results import ScottsPiResult

        return ScottsPiResult(
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
