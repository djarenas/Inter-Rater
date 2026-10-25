"""Base analyzer for inter-rater agreement coefficients."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates import AnalysisResult
    from ..intermediates.results import HypothesisTestResult
    from ..intermediates.data_profiles import BalanceProfile
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
    - Providing the public estimate() interface
    - Delegating to subclasses for coefficient-specific logic

    Subclasses override __init__ to add their own validation
    (e.g., GroupAnalyzer checks balanced_data_only, WeightableGroupAnalyzer
    builds Weights, FleissKappaAnalyzer validates ratings type).

    Some subclasses may provide alternate constructors (classmethods) that accept
    NijMatrix directly instead of Ratings.

    Attributes
    ----------
    ratings : Ratings | None
        The validated ratings data. Type varies by analyzer:
        FleissKappaAnalyzer may require BalancedRatings if config.balanced_data_only=True.
        None if initialized from NijMatrix via alternate constructor.
    config : AnalyzerConfig
        The analysis configuration (FleissKappaConfig, etc.).
        Defines bootstrap resamples, confidence level, CI method, weights, etc.
    balance_profile : BalanceProfile
        Computed from ratings or nij_matrix. Describes how subjects were rated 
        (constant n per subject, distribution of subjects by rater count, fill rate, etc.). 
        Stored here for access by all analyzers and included in results.
    """

    # Attribute declaration

    # Populated at this level
    ratings: Ratings | None
    config: AnalyzerConfig
    balance_profile: BalanceProfile

    def __init__(
        self, ratings: Ratings | None, config: AnalyzerConfig
    ) -> None:
        """
        Initialize analyzer: validate inputs and compute BalanceProfile.

        Parameters
        ----------
        ratings : Ratings | None
            The ratings data to analyze. Can be None if subclass provides
            alternate constructor that uses NijMatrix directly.
        config : AnalyzerConfig
            The analysis configuration.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object (when provided) or config is not AnalyzerConfig.
        ValueError
            If config is None.
        """
        from ..data_objects.ratings import Ratings as RatingsBase
        from .configs import AnalyzerConfig as AnalyzerConfigClass

        # Validate ratings type if provided
        if ratings is not None and not isinstance(ratings, RatingsBase):
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

        # Compute BalanceProfile from ratings if available
        from ..intermediates.data_profiles import BalanceProfile

        if self.ratings is not None:
            self.balance_profile = BalanceProfile.from_ratings(self.ratings)
        else:
            # BalanceProfile will be computed by subclass from nij_matrix
            self.balance_profile = None  # type: ignore

    # Public API ------------------------------------------------------------- #

    @abstractmethod
    def estimate(self) -> AnalysisResult:
        """
        Compute the agreement coefficient and return the result.

        Returns
        -------
        AnalysisResult
            Subclass-specific result (FleissKappaResult, CohenKappaResult, etc.)
            holding point estimate, confidence interval, and relevant metadata.

        Raises
        ------
        NotImplementedError
            If called on the abstract base class.
        ValueError
            If ratings or config are invalid for this analyzer.
        """
        raise NotImplementedError(
            f"[Analyzer] Error: estimate() must be implemented by subclass, "
            f"got {type(self).__name__}"
        )

    @abstractmethod
    def hypothesis_test(self, method: str = "permutation") -> HypothesisTestResult:
        """
        Test the null hypothesis that the agreement coefficient equals zero.

        Each leaf class implements this with the methods available for that
        coefficient. Common methods:
        - "permutation": permute ratings to build null distribution (requires Ratings)
        - "parametric_null": use analytical null variance formula (coefficient-specific)

        Parameters
        ----------
        method : str, default "permutation"
            The method to use for the hypothesis test.

        Returns
        -------
        HypothesisTestResult
            Test statistic, p-value, and metadata.

        Raises
        ------
        NotImplementedError
            If called on the abstract base class, or if the requested method
            is not implemented for this coefficient.
        ValueError
            If method="permutation" but ratings is None.
        """
        raise NotImplementedError(
            f"[Analyzer] Error: hypothesis_test() must be implemented by subclass, "
            f"got {type(self).__name__}"
        )

    def apply_injection(self) -> None:
        """
        I apply the injection pipeline from config to the ratings.

        I iterate through each injection step in ``self.config.injection``,
        instantiate the corresponding injector, and call ``inject()`` to
        produce a new Ratings with synthetic raters appended. Each step
        chains -- the output Ratings of one step becomes the input to the
        next. The new Ratings validates itself through its own constructor.

        I update ``self.ratings`` in place. If ``self.config.injection``
        is None I do nothing; if it is set and ``self.ratings`` is None
        (NijMatrix path), I raise, since injection needs individual ratings.

        Called automatically in ``__init__`` before BalanceProfile is built,
        so all downstream computation sees the injected data.
        """
        if self.config.injection is None:
            return
        if self.ratings is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: injection needs individual ratings: it "
                f"builds a new rater from the real raters' answers, and a nij "
                f"matrix holds only counts per subject. Build the analyzer "
                f"from Ratings, or remove injection from the config."
            )

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
    def _estimate_from_array(self, category_index: np.ndarray) -> float:
        """
        Compute my coefficient from a ratings array, and nothing else.

        This is the one step that differs between coefficients. Everything
        that resamples -- the bootstrap, the jackknife, anything added later
        -- is expressed in terms of it, so each of them exists once here
        rather than once per analyzer.

        I take a bare array rather than a Ratings because my callers hand me
        resampled arrays, which are not Ratings and should not have to
        pretend to be. Semantics (how many categories, which weights) come
        from the analyzer, not from the array.

        Parameters
        ----------
        category_index : np.ndarray
            Shape (n_subjects, n_raters), dtype int8. A resample, a
            leave-one-out subset, or the original data.

        Returns
        -------
        float
            The coefficient, or NaN when this array cannot define it -- an
            array holding a single category drives chance agreement to 1,
            and the coefficient is then 0/0.
        """
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: _estimate_from_array() must be implemented "
            f"by subclass"
        )

    def _bootstrap_distribution(self) -> np.ndarray:
        """
        Resample subjects with replacement and recompute my coefficient.

        Returns
        -------
        np.ndarray
            Shape (n_resamples,). NaN entries mark resamples whose
            coefficient was undefined; they are kept rather than dropped, so
            that whatever computes an interval can count them.
        """
        from .resampling_functions import subject_bootstrap

        self._require_ratings("_bootstrap_distribution")
        return np.array([
            self._estimate_from_array(resampled)
            for resampled in subject_bootstrap(
                self.ratings.category_index,
                n_resamples=self.config.n_resamples,
                random_state=self.config.random_state,
            )
        ], dtype=np.float64)

    def _jackknife_estimates(self) -> np.ndarray:
        """
        Recompute my coefficient with each subject left out in turn.

        BCa's acceleration is the standardised skewness of these values. It
        measures how much each individual subject moves the coefficient,
        which is a property of the data and cannot be recovered from a
        bootstrap distribution.

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

        self._require_ratings("_jackknife_estimates")
        return np.array([
            self._estimate_from_array(subset)
            for subset in subject_jackknife(self.ratings.category_index)
        ], dtype=np.float64)

    # Private helpers -------------------------------------------------------- #

    def _require_ratings(self, caller: str) -> None:
        """
        Guard: resampling over subjects requires Ratings.

        My subject-level resamplers read ``self.ratings.category_index``,
        which is absent when I was built from a counts matrix. The public
        entry points already check this, so reaching here means a caller took
        the ratings path on an analyzer that has no ratings -- which is a
        mistake about which method to call, not about the data.

        Parameters
        ----------
        caller : str
            The method that needs ratings, named in the message so the error
            points at the call rather than at me.

        Raises
        ------
        ValueError
            If self.ratings is None.
        """
        if self.ratings is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: {caller}() requires Ratings data, but this "
                f"analyzer was built from an NijMatrix. Use the "
                f"_from_nij_matrix variant instead."
            )

    def _check_permutation_available(self) -> None:
        """
        Guard: permutation hypothesis test requires Ratings.

        Permutation tests shuffle individual rating assignments to build the
        null distribution. This requires the original Ratings object -- it
        cannot be done from an NijMatrix because rater identity is lost.

        Raises
        ------
        ValueError
            If self.ratings is None (initialized from NijMatrix).
        """
        if self.ratings is None:
            raise ValueError(
                f"[Analyzer] Error: method='permutation' requires Ratings data. "
                f"Use method='parametric_null' when initialized from NijMatrix."
            )

    # Weights ---------------------------------------------------------------- #

    _ACCEPTS_WEIGHTS: bool = False
    """Whether this coefficient has a weighted form. The weightable bases say yes."""

    def _refuse_weights_if_unsupported(self) -> None:
        """
        Refuse ordered or weighted categories for a coefficient with no weighted form.

        Weights are declared in the semantics, so they arrive with the
        ratings. A coefficient that cannot use them must say so rather than
        quietly compute the unweighted number under a weighted declaration.
        """
        if self._ACCEPTS_WEIGHTS or self.ratings is None:
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
        does: a NaN compares False with everything, so left in it would
        count as "not extreme" and push the p-value down.

        Returns the keyword arguments HypothesisTestResult records:
        p_value, n_usable, n_at_or_above, n_at_or_below and alternative.
        The formula is ``permutation_p_value_from_counts``, which the
        result uses to check them.
        """
        from ..intermediates.results.hypothesis_test_result import (
            permutation_p_value_from_counts,
        )

        null = np.asarray(null_distribution, dtype=float)
        usable = null[~np.isnan(null)]
        if usable.size == 0:
            raise ValueError(
                f"[{type(self).__name__}] Error: every permutation gave an "
                f"undefined coefficient (NaN), so there is no null "
                f"distribution to compare with."
            )

        # A shuffle that reproduces the observed table can give the same
        # value with different rounding; count it as a tie, as
        # scipy.stats.permutation_test does.
        tol = 100 * np.finfo(float).eps * max(abs(point), 1.0)
        n_at_or_above = int(np.sum(usable >= point - tol))
        n_at_or_below = int(np.sum(usable <= point + tol))
        n_usable = int(usable.size)

        return {
            "p_value": permutation_p_value_from_counts(
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
