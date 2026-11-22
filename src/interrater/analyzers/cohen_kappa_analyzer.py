"""Analyzer for Cohen's kappa coefficient."""

from __future__ import annotations

import numpy as np

from .individual_marginal_analyzer import IndividualMarginalAnalyzer

from ..data_objects.ratings import Ratings
from ..intermediates.results import CohenKappaResult, HypothesisTestResult
from ..intermediates.data_profiles.balance_profile import BalanceProfile
from ..intermediates.data_profiles.rater_profile import RaterProfile
from .configs import CohenKappaConfig

from ..intermediates.estimator_parts import Weights
from ..intermediates.estimator_parts.estimator_input import EstimatorInput
from .agreement_functions.individual_marginals import (
    cohen_kappa,
    cohen_kappa_max,
    cohen_parts,
    rater_marginals,
    shared_rows,
)
from ..intermediates.results import AgreementParts
from ..intermediates.results.ci_results import ParametricCIResult

_ERROR_PREFIX = "[CohenKappaAnalyzer] Error"


class CohenKappaAnalyzer(IndividualMarginalAnalyzer):
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

    Cohen's kappa reads rater identity -- who said what -- and counts discard
    that. If you have only counts, you want Fleiss'
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
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped. Always present here:
        this coefficient needs rater identity.
    config : CohenKappaConfig
        The analysis configuration for Cohen's kappa.
    estimator_input : EstimatorInput
        Built from ratings, carrying category_index, n_categories, and weights.
    weights : Weights | None
        The weights declared by the semantics (from Analyzer).
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

    # Has a weighted form: Analyzer reads the weights from the semantics
    _ACCEPTS_WEIGHTS = True

    # Label in hypothesis test results
    _COEFFICIENT_KEY = "cohen_kappa"

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    estimator_input: EstimatorInput
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

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: kappa_w = 0 by the method in config.hypothesis_test_config:
        "parametric_null" (``_parametric_null_test``) or "permutation"
        (``_permutation_test``).
        """
        if self._hypothesis_test_config().method == "parametric_null":
            return self._parametric_null_test()
        return self._permutation_test()

    def _parametric_null_test(self) -> HypothesisTestResult:
        """
        Test H_0: kappa_w = 0 with the Fleiss, Cohen & Everitt (1969) null
        variance.

        No bias correction: the Fleiss (1979) correction +1/(N(n-1)) was
        derived for the multi-rater case and does not apply here.
        """
        from ..intermediates.results import HypothesisTestResult
        from .variance_functions import cohen_null_variance

        ht_config = self._hypothesis_test_config()
        point = cohen_kappa(self.estimator_input)

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

        p_value = self._parametric_p_value(z, ht_config.alternative)

        return HypothesisTestResult(
            p_value=p_value,
            alpha=ht_config.alpha,
            reject_null=self._reject_null(p_value, ht_config.alpha),
            method="parametric_null",
            coefficient="cohen_kappa",
            point=point,
            z_statistic=z,
            null_variance=var_null,
            null_variance_source=self._null_variance_source(),
            alternative=ht_config.alternative,
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

    def _agreement_parts(self) -> AgreementParts:
        """
        Po and Pe from the same function cohen_kappa is built from, so the
        result can check that they rebuild the point estimate.
        """
        po, pe = cohen_parts(self.estimator_input, self.weights)
        return AgreementParts(
            observed=po,
            chance=pe,
            kind="agreement",
            weighted=self.weights is not None and not self.weights.is_identity,
        )

    def _kappa_max(self) -> float:
        """
        The largest kappa the two raters' marginals allow (Cohen 1960).

        NaN when fewer than two subjects are shared, where kappa itself is
        NaN: there is no coefficient for it to bound.
        """
        if self._get_n_subjects() < 2:
            return float("nan")
        p_a, p_b = self._compute_marginals()
        return cohen_kappa_max(p_a, p_b, self.weights)

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Cohen's kappa from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples is written once on Analyzer in terms of this. Rater
        identity matters to me, and the column positions of data carry it.
        """
        return cohen_kappa(data)

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
            variance_reference=self._null_variance_source(),
        )

    def _null_variance_source(self) -> str:
        """Fleiss, Cohen & Everitt (1969): Eq. 9 weighted, Eq. 14 unweighted."""
        from .variance_functions import (
            COHEN_NULL_VARIANCE_SOURCE,
            COHEN_UNWEIGHTED_NULL_VARIANCE_SOURCE,
        )

        weighted = self.weights is not None and not self.weights.is_identity
        return (
            COHEN_NULL_VARIANCE_SOURCE if weighted
            else COHEN_UNWEIGHTED_NULL_VARIANCE_SOURCE
        )

    def _parametric_null_interval(self, point: float) -> ParametricCIResult:
        """
        The interval from Fleiss, Cohen & Everitt (1969) null variance, for
        ``ci_method="parametric_null"``: no bootstrap, so no draws.
        """
        p_a, p_b = self._compute_marginals()
        return self._parametric_null_ci(point, p_a, p_b, self._get_n_subjects())

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> CohenKappaResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.
        """
        return CohenKappaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=self._agreement_parts(),
            kappa_max=self._kappa_max(),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
            weights_type=self._weights_kind(),
            weight_matrix=self._weight_matrix(),
        )
