"""Base class for permutation analysis across rater subsets.

I am NOT in the Analyzer tree. An Analyzer computes one coefficient on one
data slice. I enumerate subsets (pairs, leave-one-out, arbitrary size) and
call an estimator on each, then assemble the results.

See README.md in this directory for the full architecture: three
computations (A observed, B per-pair bootstrap, C mean bootstrap), their
costs, and why each uses its own resampling scheme.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from ...data_objects.ratings import Ratings
from ...intermediates.data_profiles.balance_profile import BalanceProfile
from ...intermediates.estimator_parts.estimator_input import EstimatorInput
from ...intermediates.results.ci_results import BootstrapCIResult

if TYPE_CHECKING:
    from ...intermediates.estimator_parts.weights import Weights
    from ..configs.permutation_config import PermutationConfig

_ERROR_PREFIX = "[PermutationAnalysis] Error"


# Type alias matching rater_combinations.Estimator
Estimator = Callable[[EstimatorInput], float]


# ------------------------------------------------------------------ #
# Module-level estimator adapters                                      #
#                                                                      #
# Each wraps a coefficient function into the Estimator signature       #
# (EstimatorInput -> float). Module-level so they are picklable if     #
# we ever move to multiprocessing. Lazy imports keep the module from   #
# pulling in every agreement function at import time.                  #
# ------------------------------------------------------------------ #


def _cohen_kappa_estimator(data: EstimatorInput) -> float:
    """Cohen's kappa, weights extracted from EstimatorInput."""
    from ..agreement_functions.individual_marginals import cohen_kappa

    return cohen_kappa(data, data.weights)


def _conger_kappa_estimator(data: EstimatorInput) -> float:
    """Conger's kappa (unweighted, individual marginals)."""
    from ..agreement_functions.individual_marginals import conger_kappa

    return conger_kappa(data)


def _scotts_pi_estimator(data: EstimatorInput) -> float:
    """Scott's pi via NijMatrix bridge. Weights from EstimatorInput."""
    from ..agreement_functions.group_agreement import scotts_pi

    return scotts_pi(data.nij, data.weights)


def _gwet_ac1_estimator(data: EstimatorInput) -> float:
    """Gwet's AC1 (no weights)."""
    from ..agreement_functions.group_agreement import gwet_ac1

    return gwet_ac1(data.nij)


def _fleiss_kappa_estimator(data: EstimatorInput) -> float:
    """Fleiss' kappa. Weights from EstimatorInput."""
    from ..agreement_functions.group_agreement import fleiss_kappa

    return fleiss_kappa(data.nij, data.weights)


def _brennan_prediger_estimator(data: EstimatorInput) -> float:
    """Brennan-Prediger. Weights from EstimatorInput."""
    from ..agreement_functions.group_agreement import brennan_prediger

    return brennan_prediger(data.nij, data.weights)


def _krippendorff_alpha_estimator(data: EstimatorInput) -> float:
    """Krippendorff's alpha. Weights from EstimatorInput."""
    from ..agreement_functions.krippendorff import krippendorff_alpha

    return krippendorff_alpha(data.nij, data.weights)


_WEIGHTABLE_COEFFICIENTS = frozenset(
    {"cohen_kappa", "fleiss_kappa", "brennan_prediger", "krippendorff_alpha"}
)
"""Coefficients with a weighted form, matching the weightable analyzers."""


_ESTIMATOR_MAP: dict[str, Estimator] = {
    "cohen_kappa": _cohen_kappa_estimator,
    "scotts_pi": _scotts_pi_estimator,
    "conger_kappa": _conger_kappa_estimator,
    "gwet_ac1": _gwet_ac1_estimator,
    "fleiss_kappa": _fleiss_kappa_estimator,
    "brennan_prediger": _brennan_prediger_estimator,
    "krippendorff_alpha": _krippendorff_alpha_estimator,
}


class PermutationAnalysis:
    """
    I am the base for permutation analysis across rater subsets.

    I am NOT in the Analyzer tree. An Analyzer computes one coefficient on
    one data slice. I enumerate subsets (pairs, leave-one-out, arbitrary
    size) and call an estimator on each, then assemble the results.

    I own:
    - Validation of ratings and config
    - Applying the injection pipeline (if config_of_coefficient.injection
      is set), which adds synthetic raters to Ratings before any
      subset enumeration
    - Building the BalanceProfile and EstimatorInput (with resolved weights)
    - Resolution of the estimator callable (EstimatorInput -> float),
      used by all three computations (A observed, B per-pair bootstrap,
      C mean bootstrap)
    - CI computation helpers shared by children

    I do NOT have a run() method -- children define their own, returning
    the result type they produce (PairwisePermutationResult,
    SubsetPermutationResult).

    Attributes
    ----------
    ratings : Ratings
        Always required. Permutation analysis reads rater identity.
    config : PermutationConfig
        Configuration specifying coefficient, nested config, and optional
        mean bootstrap config.
    balance_profile : BalanceProfile
        Structural metadata about the dataset.
    estimator_input : EstimatorInput
        Built from ratings, carrying category_index, n_categories, and
        resolved weights.
    estimator : Estimator
        The bare coefficient function (EstimatorInput -> float), used
        by all three computations. Children run their own bootstrap
        loops calling this function, rather than delegating to leaf
        analyzers.
    """

    # Attribute declaration

    ratings: Ratings
    config: PermutationConfig
    balance_profile: BalanceProfile
    estimator_input: EstimatorInput
    estimator: Estimator

    # Constructor -------------------------------------------------------- #

    def __init__(
        self, ratings: Ratings, config: PermutationConfig
    ) -> None:
        """
        Initialize permutation analysis: validate inputs, build
        EstimatorInput, resolve estimator.

        Parameters
        ----------
        ratings : Ratings
            The ratings data. Required -- permutation analysis reads
            rater identity to compute coefficients across rater subsets.
        config : PermutationConfig
            Configuration for the permutation analysis.

        Raises
        ------
        TypeError
            If ratings is not Ratings or config is not PermutationConfig.
        ValueError
            If ratings is None or coefficient is unknown.
        """
        # Validate ratings
        if ratings is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: ratings is required for permutation "
                f"analysis. Permutation analysis reads rater identity "
                f"to compute coefficients across rater subsets."
            )
        if not isinstance(ratings, Ratings):
            raise TypeError(
                f"{_ERROR_PREFIX}: ratings must be Ratings, "
                f"got {type(ratings).__name__}"
            )

        # Validate config
        from ..configs.permutation_config import PermutationConfig as _PC

        if not isinstance(config, _PC):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be PermutationConfig or "
                f"subclass (PairwisePermutationConfig, "
                f"SubsetPermutationConfig), "
                f"got {type(config).__name__}"
            )

        self.ratings = ratings
        self.config = config

        # Apply injection pipeline (before BalanceProfile)
        self.apply_injection()

        # Build BalanceProfile
        self.balance_profile = BalanceProfile.from_ratings(self.ratings)

        # Weights come from the semantics, as in the analyzers
        weights = self._resolve_weights()

        # Build EstimatorInput with resolved weights
        self.estimator_input = EstimatorInput(
            category_index=self.ratings.category_index,
            n_categories=self.ratings.semantics.n_categories,
            weights=weights,
        )

        # Resolve estimator callable
        self.estimator = self._resolve_estimator()

    # Public API --------------------------------------------------------- #

    def apply_injection(self) -> None:
        """
        I apply the injection pipeline from config_of_coefficient to the ratings.

        I read ``self.config.injection`` and iterate
        through each injection step, instantiating the corresponding
        injector and calling ``inject()`` to produce a new Ratings with
        synthetic raters appended. Each step chains -- the output Ratings
        of one step becomes the input to the next.

        I update ``self.ratings`` in place. If injection is None, I do
        nothing.

        Called automatically in ``__init__`` before BalanceProfile is built,
        so all downstream computation (EstimatorInput, subset enumeration)
        sees the injected data.
        """
        injection = self.config.config_of_coefficient.injection
        if injection is None:
            return

        from ...injectors.registry import get_injector_class

        for step in injection:
            try:
                injector_cls = get_injector_class(step)
                injector = injector_cls(self.ratings, step)
                self.ratings = injector.inject()
            except Exception as e:
                raise type(e)(
                    f"{_ERROR_PREFIX} in apply_injection: "
                    f"injection step {type(step).__name__} failed: {e}"
                ) from e

    # Private helpers ---------------------------------------------------- #

    def _resolve_weights(self) -> Weights | None:
        """
        Read the weights from the semantics of the ratings.

        Follows the same rules as the analyzers: weights are declared in the
        semantics, never in the config; a coefficient with a weighted form
        uses them; a coefficient without one refuses ordered or weighted
        categories rather than quietly ignoring them.

        Returns
        -------
        Weights | None
            The declared weights, or None for unordered categories.

        Raises
        ------
        ValueError
            If config_of_coefficient still sets weights_type or
            custom_weights_path, or if the coefficient has no weighted form
            and the semantics declares ordered or weighted categories.
        """
        from ...intermediates.estimator_parts.weights import Weights
        from ...semantics.categories import DistanceCategories

        coc          = self.config.config_of_coefficient
        weights_type = getattr(coc, "weights_type", "identity")
        custom_path  = getattr(coc, "custom_weights_path", None)
        if weights_type != "identity" or custom_path is not None:
            raise ValueError(
                f"{_ERROR_PREFIX}: weights are no longer set in "
                f"config_of_coefficient (got weights_type={weights_type!r}, "
                f"custom_weights_path={custom_path!r}). Declare them in the "
                f"semantics instead: categories kind 'ordered' with "
                f"'weighting', or kind 'weighted' with 'weights' or "
                f"'weights_file'."
            )

        categories  = self.ratings.semantics.categories
        coefficient = self.config.coefficient_to_compute
        if (
            isinstance(categories, DistanceCategories)
            and coefficient not in _WEIGHTABLE_COEFFICIENTS
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: {coefficient} has no weighted form, but "
                f"the semantics declares categories of kind "
                f"{categories._KIND!r}. To compute it, use a semantics with "
                f"kind: unordered for these ratings."
            )
        return Weights.from_categories(categories)

    def _resolve_estimator(self) -> Estimator:
        """
        Map coefficient name to the estimator callable.

        The estimator is a pure function (EstimatorInput -> float) used
        by all three computations: A (observed values via over_pairs),
        B (per-pair bootstrap loop in the child), and C (correlated
        global resample loop or independent sampling from B's output).

        The adapter functions are defined at module level so they are
        picklable if we ever move to multiprocessing. Each wraps a
        coefficient function (which may take extra arguments like weights
        or expect an NijMatrix) into the Estimator signature.

        Returns
        -------
        Estimator
            The coefficient function.

        Raises
        ------
        ValueError
            If coefficient name has no registered estimator.
        """
        coefficient = self.config.coefficient_to_compute
        estimator = _ESTIMATOR_MAP.get(coefficient)
        if estimator is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: no estimator registered for "
                f"'{coefficient}'. Known coefficients: "
                f"{sorted(_ESTIMATOR_MAP.keys())}"
            )
        return estimator

    def _compute_ci(
        self,
        distribution: np.ndarray,
        point: float,
        ci_method: str,
        alpha: float,
        jackknife_fn: Callable[[], np.ndarray] | None = None,
    ) -> BootstrapCIResult:
        """
        Compute confidence interval from a bootstrap distribution.

        Shared by children for both per-pair CIs and mean CIs.
        Delegates to the CI functions in ci_methods.py.

        Parameters
        ----------
        distribution : np.ndarray
            Bootstrap values (1-D).
        point : float
            The point estimate.
        ci_method : str
            "percentile" or "bca". Permutation results do not support
            "parametric_null" -- subsets sharing a rater are correlated,
            so there is no closed-form null variance.
        alpha : float
            Significance level (e.g., 0.05 for 95% CI).
        jackknife_fn : Callable[[], np.ndarray] | None
            Produces leave-one-out estimates when called. A callable rather
            than an array because only BCa needs them, and here they are
            expensive -- per pair, one coefficient evaluation per subject.

            **Which subjects to leave out depends on what was resampled.**
            A per-pair interval bootstraps the subjects that pair shares, so
            its jackknife leaves those out rather than all subjects: a
            subject the pair never rated cannot influence the pair's
            coefficient, so including it would contribute a deviation
            representing no influence at all. A correlated mean bootstraps
            subjects globally, so its jackknife is global too. The caller
            knows which it did, so the caller supplies the function.

            None means no jackknife exists for this quantity, which is the
            case for the independent mean: it resamples from per-pair
            distributions and never touches a subject, so there is nothing
            to leave out.

        Returns
        -------
        BootstrapCIResult
            Bounds, the method asked for against the one that ran, and how
            many resamples failed.

        Raises
        ------
        ValueError
            If ci_method is unknown, or if BCa is requested for a quantity
            that has no jackknife.
        """
        from ..ci_methods import bca_ci, percentile_ci

        if ci_method == "percentile":
            return percentile_ci(distribution, alpha)
        elif ci_method == "bca":
            if jackknife_fn is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: ci_method='bca' needs leave-one-out "
                    f"estimates to compute acceleration, and this quantity "
                    f"has none. The independent mean resamples from per-pair "
                    f"distributions rather than from subjects, so no subject "
                    f"can be left out of it. Use ci_method='percentile', or "
                    f"the correlated mean."
                )
            return bca_ci(
                distribution, point, alpha,
                jackknife_estimates=jackknife_fn(),
            )
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown ci_method '{ci_method}'. "
                f"Expected 'percentile' or 'bca'. Permutation results "
                f"do not support 'parametric_null' (subsets sharing a "
                f"rater are correlated, so there is no closed-form null "
                f"variance)."
            )
