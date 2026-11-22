"""Analyzer for Fleiss' kappa coefficient."""

from __future__ import annotations

import numpy as np

from .group_analyzer import GroupAnalyzer

from ..data_objects.ratings import Ratings
from ..intermediates.results import FleissKappaResult, HypothesisTestResult
from ..intermediates.results.ci_results import ParametricCIResult
from .configs import FleissKappaConfig

from .agreement_functions.group_agreement import (
    chance_fleiss,
    fleiss_kappa,
    observed_agreement,
)
from ..intermediates.estimator_parts import EstimatorInput, NijMatrix, Weights
from ..intermediates.data_profiles import BalanceProfile, RaterProfile
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.fleiss_kappa_result import (
    KAPPA_MAX_NOT_COMPUTED,
)

_ERROR_PREFIX = "[FleissKappaAnalyzer] Error"


class FleissKappaAnalyzer(GroupAnalyzer):
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
    ratings : Ratings
        The ratings data (from Analyzer).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    config : FleissKappaConfig
        The analysis configuration for Fleiss' kappa.
    weights : Weights | None
        The weights object (from Analyzer).
    is_exact : bool
        Whether the exact formula applies (data is perfectly balanced).
        Computed during initialization based on balance_profile.
    """

    # Has a weighted form: Analyzer reads the weights from the semantics
    _ACCEPTS_WEIGHTS = True

    # Label in hypothesis test results
    _COEFFICIENT_KEY = "fleiss_kappa"

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
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

        # Call parent __init__ (validates group config; Analyzer builds the Weights)
        # Wrap to propagate errors with FleissKappaAnalyzer context
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        # Validate config type and compute is_exact flag
        self._validate_fleiss_config()
        self.is_exact = self.balance_profile.is_balanced
        self._refuse_parametric_null_if_invalid(allow_weights=False)

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H₀: κ = 0 by the method in config.hypothesis_test_config:
        "parametric_null" (``_parametric_null_test``) or "permutation"
        (``_permutation_test``).
        """
        if self._hypothesis_test_config().method == "parametric_null":
            return self._parametric_null_test()
        return self._permutation_test()

    def _parametric_null_test(self) -> HypothesisTestResult:
        """
        Test H₀: κ = 0 with the Fleiss (1979) null variance (Eq. 12) and the
        bias correction of Eq. 17.
        """
        from ..intermediates.results import HypothesisTestResult
        from .variance_functions import fleiss_null_variance

        ht_config = self._hypothesis_test_config()
        nij = self.estimator_input.nij
        point = fleiss_kappa(self.estimator_input)
        pj = self._compute_pj(nij)

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

        p_value = self._parametric_p_value(z, ht_config.alternative)

        return HypothesisTestResult(
            p_value=p_value,
            alpha=ht_config.alpha,
            reject_null=self._reject_null(p_value, ht_config.alpha),
            method="parametric_null",
            coefficient="fleiss_kappa",
            point=point,
            z_statistic=z,
            null_variance=var_null,
            alternative=ht_config.alternative,
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

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Fleiss' kappa from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples, bootstrap or jackknife, is written once on my base
        classes in terms of this.
        """
        return fleiss_kappa(data)

    def _agreement_parts(self, nij: NijMatrix) -> AgreementParts:
        """
        Po and Pe on the counts the point estimate was computed from.

        The same two functions fleiss_kappa combines, so the result can check
        that they rebuild the point estimate.
        """
        return AgreementParts(
            observed=observed_agreement(nij, self.weights),
            chance=chance_fleiss(nij, self.weights),
            kind="agreement",
            weighted=self.weights is not None and not self.weights.is_identity,
        )

    def _parametric_null_ci(
        self, point: float, pj: np.ndarray
    ) -> ParametricCIResult:
        """
        Compute confidence interval using Fleiss (1979) null variance.

        The half-width comes from the variance of kappa under H0: kappa = 0,
        so the interval is centred on the observed estimate but scaled as
        though the estimate were zero. The returned ParametricCIResult says
        so whenever it is printed -- see its docstring for why that matters.

        Parameters
        ----------
        point : float
            The point estimate of kappa.
        pj : np.ndarray
            Category proportions.

        Returns
        -------
        ParametricCIResult
            Bounds, the level, and the variance they were built from.
        """
        from scipy import stats
        from .variance_functions import fleiss_null_variance

        nij = self.estimator_input.nij
        n_subjects = nij.n_subjects
        n_raters = int(nij.n_per_subject[0])  # Assumes balanced

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
    ) -> FleissKappaResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.
        """
        return FleissKappaResult(
            point=point,
            confidence_interval=ci,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            bootstrap_distribution=bootstrap_distribution,
            parts=self._agreement_parts(self.estimator_input.nij),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
            weights_type=self._weights_kind(),
            weight_matrix=self._weight_matrix(),
        )
