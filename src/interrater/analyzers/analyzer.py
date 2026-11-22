"""Base analyzer for inter-rater agreement coefficients."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.estimator_parts import EstimatorInput, Weights
    from ..intermediates import AnalysisResult
    from ..intermediates.results import HypothesisTestResult
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
    from .configs.analyzer_config import AnalyzerConfig

_ERROR_PREFIX = "[Analyzer] Error"


class Analyzer(ABC):
    """
    I am the abstract base for all inter-rater agreement analyzers.

    Every analyzer takes Ratings and a Config, validates both, computes a
    BalanceProfile that describes the data structure, and provides an
    estimate() method that returns an AnalysisResult.

    I am responsible for:
    - Storing valid Ratings and Config objects
    - Applying the injection pipeline (if config.injection is set),
      which adds synthetic raters to Ratings before any analysis
    - Computing BalanceProfile (relevant for most coefficients)
    - Building the EstimatorInput (array, K, weights) that every estimate
      and every resample is computed from
    - Providing the public estimate() interface
    - Delegating to subclasses for coefficient-specific logic

    Subclasses override __init__ to add their own validation
    (e.g., GroupAnalyzer checks balanced_data_only, FleissKappaAnalyzer
    validates ratings type). But their constructor can utilize the parent
    constructors using super __init__.

    Attributes
    ----------
    ratings : Ratings
        The validated ratings data, after any injection.
    config : AnalyzerConfig
        The analysis configuration (FleissKappaConfig, etc.).
        Defines bootstrap resamples, confidence level, CI method, hypothesis
        test and injection. Weights come from the semantics, not from here.
    balance_profile : BalanceProfile
        Computed from ratings. Describes how subjects were rated
        (constant n per subject, distribution of subjects by rater count, etc.).
        Stored here for access by all analyzers and included in results.
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped.
    weights : Weights | None
        From the semantics, for a coefficient with a weighted form. None
        for unordered categories, and always None for a coefficient without
        a weighted form.
    estimator_input : EstimatorInput
        The ratings array, K from the semantics, and the weights. What the
        coefficient is computed from, and what every resample is drawn from.
    """

    # Attribute declaration

    # Populated at this level
    ratings: Ratings
    config: AnalyzerConfig
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    weights: Weights | None
    estimator_input: EstimatorInput

    def __init__(self, ratings: Ratings, config: AnalyzerConfig) -> None:
        """
        Validate the inputs, apply the injection, and build the profiles,
        the weights and the EstimatorInput.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : AnalyzerConfig
            The analysis configuration.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object or config is not AnalyzerConfig.
        ValueError
            If config is None.
        """
        from ..data_objects.ratings import Ratings as RatingsBase
        from .configs import AnalyzerConfig as AnalyzerConfigClass

        # Validate ratings type
        if not isinstance(ratings, RatingsBase):
            raise TypeError(
                f"{_ERROR_PREFIX}: ratings must be a Ratings object, "
                f"got {type(ratings).__name__}"
            )

        # Validate config is not None
        if config is None:
            raise ValueError(f"{_ERROR_PREFIX}: config cannot be None")

        # Validate config type
        if not isinstance(config, AnalyzerConfigClass):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be an AnalyzerConfig, "
                f"got {type(config).__name__}"
            )

        # Store inputs
        self.ratings = ratings
        self.config = config

        # Apply injection pipeline (before BalanceProfile)
        self.apply_injection()

        # Compute the profiles from the (possibly injected) ratings
        from ..intermediates.data_profiles import BalanceProfile, RaterProfile

        self.balance_profile = BalanceProfile.build_from_ratings(self.ratings)
        self.rater_profile = RaterProfile.build_from_ratings(self.ratings)

        # K and weights travel with the array from here on
        from ..intermediates.estimator_parts import EstimatorInput, Weights

        self.weights = (
            Weights.build_from_categories(self.ratings.semantics.categories)
            if self._ACCEPTS_WEIGHTS else None
        )
        self.estimator_input = EstimatorInput(
            category_index=self.ratings.category_index,
            n_categories=self.ratings.semantics.n_categories,
            weights=self.weights,
        )

    # Public API ------------------------------------------------------------- #

    def estimate(self) -> AnalysisResult:
        """
        Compute the agreement coefficient and return the result.

        The same three steps for every coefficient, so they are written once
        here: the point from the ratings, the interval (bootstrap, or the
        null variance for ``parametric_null``), then the coefficient's own
        result. Each analyzer supplies only what differs: the coefficient
        (``_estimate_from_input``), its result (``_result``), and, where it
        has one, its null variance (``_parametric_null_interval``).

        Returns
        -------
        AnalysisResult
            Subclass-specific result (FleissKappaResult, CohenKappaResult, etc.)
            holding point estimate, confidence interval, and relevant metadata.
        """
        point = self._estimate_from_input(self.estimator_input)

        if self.config.ci_method == "parametric_null":
            ci = self._parametric_null_interval(point)
            bootstrap_distribution = None
        else:
            bootstrap_distribution = self._bootstrap_distribution()
            ci = self._compute_ci(
                bootstrap_distribution, point, self._jackknife_estimates
            )

        return self._result(point, ci, bootstrap_distribution)

    @abstractmethod
    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test the null hypothesis that the agreement coefficient equals zero.

        The method, scheme and alternative come from
        ``config.hypothesis_test_config``, as everything else about a run
        comes from the config. Each leaf class implements this with the
        methods available for its coefficient:
        - "permutation": every coefficient, through ``_permutation_test``
        - "parametric_null": only Fleiss, Scott and Cohen, which have a
          published null variance, each in its own ``_parametric_null_test``

        Returns
        -------
        HypothesisTestResult
            Test statistic, p-value, and metadata.

        Raises
        ------
        NotImplementedError
            If called on the abstract base class.
        ValueError
            If config.hypothesis_test_config is None.
        """
        raise NotImplementedError(
            f"[Analyzer] Error: hypothesis_test() must be implemented by subclass, "
            f"got {type(self).__name__}"
        )

    def apply_injection(self) -> None:
        """
        Apply the injection pipeline from config to the ratings.

        Iterate through each injection step in ``self.config.injection``,
        instantiate the corresponding injector, and call ``inject()`` to
        produce a new Ratings with synthetic raters appended. Each step
        chains -- the output Ratings of one step becomes the input to the
        next. The new Ratings validates itself through its own constructor.

        Updates ``self.ratings`` in place. Does nothing if
        ``self.config.injection`` is None.

        Called automatically in ``__init__`` before BalanceProfile is built,
        so all downstream computation sees the injected data.
        """
        if self.config.injection is None:
            return

        from ..injectors.registry import get_injector_class

        for position, step in enumerate(self.config.injection):
            try:
                injector_cls = get_injector_class(step)
                injector = injector_cls(self.ratings, step, step_position=position)
                self.ratings = injector.inject()
            except Exception as e:
                raise type(e)(
                    f"[Analyzer] Error in apply_injection: "
                    f"injection step {type(step).__name__} failed: {e}"
                ) from e

    # Resampling ------------------------------------------------------------- #

    @abstractmethod
    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Compute my coefficient from an EstimatorInput.

        It takes an EstimatorInput rather than a Ratings because its callers
        may pass resamples. K and the weights come with an EstimatorInput
        object, taken from ``self.estimator_input``, so that the methods do
        not have to infer K from a resampled array.

        Parameters
        ----------
        data : EstimatorInput
            ``self.estimator_input``, a resample of it, or a leave-one-out
            subset of it.

        Returns
        -------
        float
            The coefficient, or NaN when it is undefined on this array (e.g.
            Fleiss' kappa when every rating is in one category, so chance
            agreement is 1).
        """
        # Just an extra precaution to ensure proper inheritance of the abstract method.
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: _estimate_from_input() must be implemented "
            f"by subclass"
        )

    def _with_category_index(self, category_index: np.ndarray) -> EstimatorInput:
        """
        Use the same K and same weights from my EstimatorInput attribute, but
        with a new data array. Used by the permutation tests.
        """
        import dataclasses

        return dataclasses.replace(
            self.estimator_input, category_index=category_index
        )

    def _bootstrap_distribution(self) -> np.ndarray:
        """
        Return the bootstrap distribution as a numpy array of coefficients.
        One coefficient for each of the ``n_resamples`` resamples of the
        subjects, drawn with replacement.

        Returns
        -------
        np.ndarray
            Shape (n_resamples,). NaN entries mark resamples whose
            coefficient was undefined; they are kept rather than dropped, so
            that whatever computes an interval can count them.
        """
        from .resampling_functions import subject_bootstrap

        return np.array([
            self._estimate_from_input(resampled)
            for resampled in subject_bootstrap(
                self.estimator_input,
                n_resamples=self.config.n_resamples,
                random_state=self.config.random_state,
            )
        ], dtype=np.float64)

    def _jackknife_estimates(self) -> np.ndarray:
        """
        My coefficient with each subject left out in turn.

        BCa's acceleration is computed from these values. The "a" needs
        to know how much each subject moves the coefficient.

        Returns
        -------
        np.ndarray
            Shape (n_subjects,). NaN entries mark leave-one-out subsets whose
            coefficient was undefined, for the same reason a resample can be.

        References
        ----------
        Efron, B., & Tibshirani, R. 1986. *Statistical Science*, 1(1), 54-77.
        Equation (8.5), p. 71.
        """
        from .resampling_functions import subject_jackknife

        return np.array([
            self._estimate_from_input(subset)
            for subset in subject_jackknife(self.estimator_input)
        ], dtype=np.float64)

    def _compute_ci(
        self,
        bootstrap_distribution: np.ndarray,
        point: float,
        jackknife_fn,
    ):
        """
        Build my interval from my draws: percentile or BCa.

        ``jackknife_fn`` is called only for BCa, whose acceleration needs the
        leave-one-out estimates; percentile never pays for them.
        ``parametric_null`` does not come here -- it has no draws.
        """
        from .ci_methods import bca_ci, percentile_ci

        if self.config.ci_method == "percentile":
            return percentile_ci(bootstrap_distribution, self.config.alpha)
        if self.config.ci_method == "bca":
            return bca_ci(
                bootstrap_distribution,
                point,
                self.config.alpha,
                jackknife_estimates=jackknife_fn(),
            )
        raise ValueError(
            f"{_ERROR_PREFIX}: unknown ci_method '{self.config.ci_method}' "
            f"for a bootstrap interval. Expected 'percentile' or 'bca'."
        )

    @abstractmethod
    def _result(self, point: float, ci, bootstrap_distribution):
        """
        Package my point, interval and draws into my result class.

        The last step of ``estimate()``, kept separate so that whatever
        computed the point and interval -- ``estimate()`` itself, or a
        permutation analysis that bootstrapped a pair or subset with its own
        seeds -- gets the same result, built and checked the same way.
        """
        # Just an extra precaution to ensure proper inheritance of the abstract method.
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: _result() must be implemented by subclass"
        )

    def _parametric_null_interval(self, point: float):
        """
        The interval from the null variance, for ``parametric_null``.

        Only coefficients with a published null variance override this
        (Cohen, Scott, Fleiss). Their configs and analyzers refuse
        ``parametric_null`` everywhere else, so reaching here is a bug.
        """
        raise ValueError(
            f"{_ERROR_PREFIX}: {type(self).__name__} has no null variance, "
            f"so ci_method='parametric_null' cannot be used. This should have "
            f"been refused when the config was built; this is a bug."
        )

    # Private helpers -------------------------------------------------------- #

    _COEFFICIENT_KEY: str
    """How hypothesis test results name the coefficient (e.g. "fleiss_kappa"). Set by each leaf."""

    def _hypothesis_test_config(self):
        """
        The hypothesis test settings in my config; ValueError if there are none.
        """
        ht_config = self.config.hypothesis_test_config
        if ht_config is None:
            raise ValueError(
                f"[{type(self).__name__}] Error: hypothesis_test_config is None. "
                f"Set hypothesis_test_config in your config to use hypothesis testing."
            )
        return ht_config

    def _permutation_test(self) -> HypothesisTestResult:
        """
        Test H0: coefficient = 0 by permutation. Any coefficient can use it.

        The ratings are shuffled ``n_permutations`` times by the scheme in my
        hypothesis test config, my coefficient is computed on each shuffle,
        and the p-value counts the shuffles at or beyond the observed value
        (``_permutation_p_value``).
        """
        from ..intermediates.results import HypothesisTestResult
        from .resampling_functions.permutation_resampling import (
            pool_and_redistribute,
            within_rater_shuffle,
        )

        # Check that the config object arrives with valid hypothesis test settings
        # and that you are not running a permutation test when something else was chosen.
        ht_config = self._hypothesis_test_config()
        prefix = f"[{type(self).__name__}] Error"
        if ht_config.method != "permutation":
            raise ValueError(
                f"{prefix}: hypothesis test method '{ht_config.method}' is not "
                f"available for this coefficient. Use method='permutation'."
            )

        scheme_map = {
            "within_rater_shuffle": within_rater_shuffle,
            "pool_and_redistribute": pool_and_redistribute,
        }
        scheme_fn = scheme_map.get(ht_config.scheme)
        if scheme_fn is None:
            raise ValueError(
                f"{prefix}: unknown permutation scheme '{ht_config.scheme}'. "
                f"Expected one of: {set(scheme_map.keys())}."
            )

        point = self._estimate_from_input(self.estimator_input)
        null_distribution = np.array([
            self._estimate_from_input(self._with_category_index(permuted_array))
            for permuted_array in scheme_fn(
                self.estimator_input.category_index,
                n_permutations=ht_config.n_permutations,
                random_state=ht_config.random_state,
            )
        ])
        counts = self._permutation_p_value(
            null_distribution, point, ht_config.alternative
        )

        return HypothesisTestResult(
            alpha=ht_config.alpha,
            reject_null=self._reject_null(counts["p_value"], ht_config.alpha),
            method="permutation",
            coefficient=self._COEFFICIENT_KEY,
            point=point,
            null_distribution=null_distribution,
            n_permutations=ht_config.n_permutations,
            scheme=ht_config.scheme,
            random_state=ht_config.random_state,
            **counts,
        )

    def _agreement_profile(self):
        """
        My agreement profile, from my ratings. Every result carries one.
        """
        from ..intermediates.data_profiles import AgreementProfile

        return AgreementProfile.build_from_ratings(self.ratings)

    # Weights ---------------------------------------------------------------- #

    _ACCEPTS_WEIGHTS: bool = False
    """Whether this coefficient has a weighted form. Fleiss, Gwet's AC1 (as AC2), Brennan-Prediger, Krippendorff and Cohen say yes."""

    def _refuse_weights_if_unsupported(self) -> None:
        """
        Refuse ordered or weighted categories for a coefficient with no weighted form.

        Weights are declared in the semantics, so they arrive at the Analyzer
        object with the ratings. The analyzer must look at the ratings object,
        and if its coefficient cannot use weights, then it must throw a clear
        error instead of computing an unweighted number silently.
        """
        if self._ACCEPTS_WEIGHTS:
            return
        from ..semantics.categories import DistanceCategories

        categories = self.ratings.semantics.categories
        if isinstance(categories, DistanceCategories):
            # No prefix here: the leaf constructor that called me adds its own.
            raise ValueError(
                f"this coefficient has no weighted form, but the semantics "
                f"declares categories of kind {categories._KIND!r}. To "
                f"compute it, use a semantics with kind: unordered for these "
                f"ratings."
            )

    def _refuse_parametric_null_if_invalid(self, *, allow_weights: bool) -> None:
        """
        Refuse parametric_null where its formula does not apply.

        The null variances behind it (Fleiss, Nee & Landis 1979 for Fleiss
        and Scott; Fleiss, Cohen & Everitt 1969 for Cohen) assume every
        subject has the same number of raters. The Fleiss one also has no
        weighted form, while Cohen's does, hence ``allow_weights``. Checked
        at construction, the first point where the method asked for, the
        balance of the data and the weights are all known; both the
        interval and the hypothesis test are covered.

        The analyzer object looks at the config and if its coefficient does not
        allow a parametric null, it shows a clear error.
        """
        asked = []
        if self.config.ci_method == "parametric_null":
            asked.append("ci_method='parametric_null'")
        ht = self.config.hypothesis_test_config
        if ht is not None and ht.method == "parametric_null":
            asked.append("hypothesis test method='parametric_null'")
        if not asked:
            return
        what = " and ".join(asked)
        prefix = f"[{type(self).__name__}] Error"

        profile = self.balance_profile
        if not profile.is_balanced:
            counts = profile.n_per_subject
            raise ValueError(
                f"{prefix}: {what} assumes every subject has the same number of "
                f"raters; these data have {int(counts.min())} to "
                f"{int(counts.max())}. Use ci_method='percentile' or 'bca', "
                f"or hypothesis test method='permutation'."
            )
        if not allow_weights and getattr(self, "weights", None) is not None:
            raise ValueError(
                f"{prefix}: {what} uses an unweighted null variance, but the semantics "
                f"declares a weighted kind. Use ci_method='percentile' or "
                f"'bca', or hypothesis test method='permutation'."
            )

    def _weights_kind(self) -> str:
        """How the result should name its weighting: identity, linear, quadratic or custom."""
        weights = getattr(self, "weights", None)
        return "identity" if weights is None else weights.kind

    def _weight_matrix(self) -> list[list[float]] | None:
        """The weight matrix for the result, as plain nested lists; None when unweighted."""
        weights = getattr(self, "weights", None)
        return None if weights is None else weights.matrix.tolist()

    def _reject_null(self, p_value: float, alpha: float) -> bool:
        """Decide the test as a plain Python bool, with context if it fails."""
        try:
            return bool(p_value < alpha)
        except Exception as exc:
            raise type(exc)(
                f"[{type(self).__name__}] Error: could not compare "
                f"p_value={p_value!r} with alpha={alpha!r}: {exc}"
            ) from exc

    def _permutation_p_value(
        self, null_distribution: np.ndarray, point: float, alternative: str
    ) -> dict:
        """
        Count the permutations at or beyond the observed point, and the p-value.

        Draws that came out NaN are skipped and counted, as the bootstrap
        does. This prevents NaNs from being compared as False with
        everything, which would push the p-value down.

        Returns the keyword arguments HypothesisTestResult records:
        p_value, n_usable, n_at_or_above, n_at_or_below and alternative.
        The formula is ``compute_permutation_p_value``, which the
        result uses to check them.
        """
        from ..intermediates.results.hypothesis_test_result import (
            compute_permutation_p_value,
        )

        null = np.asarray(null_distribution, dtype=float)
        usable = null[~np.isnan(null)]
        if usable.size == 0:
            raise ValueError(
                f"[{type(self).__name__}] Error: every permutation gave an "
                f"undefined coefficient (NaN), so there is no null "
                f"distribution to compare with."
            )

        # A permutation can give exactly the observed coefficient, for example
        # by reproducing the same counts, and still differ from the point in
        # the last digits, from rounding. Values within a small
        # tolerance of the point count as ties: 100 machine epsilons times
        # |point|, as in scipy.stats.permutation_test, but never less than
        # 100 epsilons.
        tol = 100 * np.finfo(float).eps * max(abs(point), 1.0)
        n_at_or_above = int(np.sum(usable >= point - tol))
        n_at_or_below = int(np.sum(usable <= point + tol))
        n_usable = int(usable.size)

        return {
            "p_value": compute_permutation_p_value(
                n_at_or_above, n_at_or_below, n_usable, alternative
            ),
            "alternative": alternative,
            "n_usable": n_usable,
            "n_at_or_above": n_at_or_above,
            "n_at_or_below": n_at_or_below,
        }

    @staticmethod
    def _parametric_p_value(z: float, alternative: str) -> float:
        """The normal-theory p-value of z for the given alternative."""
        from scipy import stats

        if alternative == "greater":
            return float(stats.norm.sf(z))
        if alternative == "less":
            return float(stats.norm.cdf(z))
        return float(min(1.0, 2.0 * stats.norm.sf(abs(z))))
