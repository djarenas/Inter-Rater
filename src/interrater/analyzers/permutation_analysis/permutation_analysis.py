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
from ...intermediates.data_profiles.rater_profile import RaterProfile
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


def _analyzer_class(coefficient: str) -> type:
    """The analyzer that makes a coefficient's result, for each pair or subset."""
    from .. import (
        BrennanPredigerAnalyzer,
        CohenKappaAnalyzer,
        CongerKappaAnalyzer,
        FleissKappaAnalyzer,
        GwetAC1Analyzer,
        KrippendorffAlphaAnalyzer,
        ScottsAnalyzer,
    )

    return {
        "cohen_kappa": CohenKappaAnalyzer,
        "scotts_pi": ScottsAnalyzer,
        "conger_kappa": CongerKappaAnalyzer,
        "gwet_ac1": GwetAC1Analyzer,
        "fleiss_kappa": FleissKappaAnalyzer,
        "brennan_prediger": BrennanPredigerAnalyzer,
        "krippendorff_alpha": KrippendorffAlphaAnalyzer,
    }[coefficient]


class IndependentMeanDraws:
    """
    I am the independent mean's draws, built while the units are bootstrapped.

    Each draw of the independent mean picks one value from each defined
    unit's (pair's or subset's) draws and averages the picks. I take the
    units one at a time, as Computation B makes them: I pick n_resamples
    values from a unit's draws and add them to a running total. After that
    the unit's draws can be tossed, so the independent mean never needs
    every unit's draws at once.

    This is the only way the independent mean is built, whatever
    keep_draws says, so the same seed gives the same interval whether the
    draws are kept or not. The picks come from one stream seeded with
    config_of_coefficient's random_state, taken in unit order.

    A NaN pick makes that draw NaN, so it counts as failed rather than
    averaging the other units. With no defined unit, every draw is NaN.
    """

    def __init__(
        self, units: list, n_resamples: int, random_state: int | None
    ) -> None:
        self._units = list(units)
        self._defined = set(self._units)
        self._added = 0
        self._rng = np.random.default_rng(random_state)
        self._total = np.zeros(n_resamples, dtype=np.float64)

    def add(self, unit: tuple, draws: np.ndarray) -> None:
        """Pick from a unit's draws; units that are not defined are skipped."""
        if unit not in self._defined:
            return
        # The picks come off one stream, so the order fixes the result.
        if self._added >= len(self._units) or unit != self._units[self._added]:
            raise ValueError(
                f"{_ERROR_PREFIX}: independent mean got unit {unit} out of "
                f"order; this is a bug, not a property of the data."
            )
        picks = self._rng.integers(0, draws.size, size=self._total.size)
        self._total += draws[picks]
        self._added += 1

    def draws(self) -> np.ndarray:
        """The mean's draws, once every defined unit has been added."""
        if self._added != len(self._units):
            raise ValueError(
                f"{_ERROR_PREFIX}: independent mean has {self._added} of "
                f"{len(self._units)} units; this is a bug, not a property "
                f"of the data."
            )
        if not self._units:
            return np.full(self._total.size, np.nan)
        return self._total / len(self._units)


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
    rater_profile : RaterProfile
        The panel's raters, injected ones included: names, how much each
        rated, and how much each pair overlapped. Shared by every result.
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

        # Conger's kappa is for complete data only, checked after injection
        if config.coefficient_to_compute == "conger_kappa":
            from ..individual_marginal_analyzer import check_complete

            try:
                check_complete(self.ratings.category_index, "Conger's kappa")
            except ValueError as exc:
                raise ValueError(f"{_ERROR_PREFIX}: {exc}") from exc

        # Build BalanceProfile and RaterProfile
        self.balance_profile = BalanceProfile.build_from_ratings(self.ratings)
        self.rater_profile = RaterProfile.build_from_ratings(self.ratings)

        # The panel must be balanced only if the config says so
        if config.balanced_data_only and not self.balance_profile.is_balanced:
            counts = sorted(self.balance_profile.subject_counts)
            raise ValueError(
                f"{_ERROR_PREFIX}: balanced_data_only is true, but the panel "
                f"is unbalanced: subjects have {counts[0]} to {counts[-1]} "
                f"ratings. Set balanced_data_only: false (next to "
                f"coefficient_to_compute) to run anyway: each pair uses the "
                f"subjects both raters rated, each subset the subjects at "
                f"least two of its raters rated."
            )

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

        for position, step in enumerate(injection):
            try:
                injector_cls = get_injector_class(step)
                injector = injector_cls(self.ratings, step, step_position=position)
                self.ratings = injector.inject()
            except Exception as e:
                raise type(e)(
                    f"{_ERROR_PREFIX} in apply_injection: "
                    f"injection step {type(step).__name__} failed: {e}"
                ) from e

    # Private helpers ---------------------------------------------------- #

    def _check_cost(self, units_per_k: dict[int, int], unit_word: str) -> None:
        """
        Refuse, before anything runs, a run above config.max_evaluations.

        Counts the times the coefficient would be computed: the observed
        value of each unit, its bootstrap (n_resamples each), and the
        correlated mean (n_resamples draws over every unit). With BCa each
        interval adds a jackknife, counted as one evaluation per subject in
        the panel -- an upper bound, since a unit is jackknifed only over its
        own subjects. So the count is "at most", and the guard can only
        refuse a little early, never late. The independent mean adds
        nothing: it picks from draws already made.

        Parameters
        ----------
        units_per_k : dict[int, int]
            How many pairs or subsets the run will compute, per subset size.
        unit_word : str
            "pairs" or "subsets", for the message.

        Raises
        ------
        ValueError
            With the count, what it is made of, and how to lower or allow it.
        """
        coc = self.config.config_of_coefficient
        bmc = self.config.bootstrap_mean_config
        n_resamples = coc.n_resamples
        jackknife = self.ratings.n_subjects if coc.ci_method == "bca" else 0
        correlated = bmc is not None and bmc.method == "correlated"

        units = sum(units_per_k.values())
        per_unit = 1 + n_resamples + jackknife
        if correlated:
            per_unit += n_resamples + jackknife
        total = units * per_unit
        if total <= self.config.max_evaluations:
            return

        parts = [f"{units:,} {unit_word} x {n_resamples:,} resamples"]
        if correlated:
            parts.append("the correlated mean")
        if jackknife:
            parts.append(f"BCa jackknives over up to {jackknife:,} subjects")
        raise ValueError(
            f"{_ERROR_PREFIX}: this run would compute the coefficient up to "
            f"{total:,} times ({', plus '.join(parts)}), above max_evaluations "
            f"({self.config.max_evaluations:,}). Lower n_resamples or the "
            f"number of {unit_word}, or raise max_evaluations in the config "
            f"to run it on purpose."
        )

    def _unit_result(
        self,
        unit: tuple[int, ...],
        point: float,
        ci: BootstrapCIResult,
        draws: np.ndarray | None,
    ):
        """
        One pair's or subset's result, built by the coefficient's own analyzer.

        The unit's data is its usable subjects: those at least two of its
        raters rated (for a pair, both). The analyzer is built on that data
        with injection off -- the panel's injection already happened -- and
        given the point, interval and draws this analysis computed, through
        the same ``_result`` its ``estimate()`` uses. So a pair is exactly
        the result the analyzer would give on those subjects, built and
        checked the same way.

        Returns NotComputed, with the reason, when the unit cannot be a
        result: fewer than 2 usable subjects, or a subset rater who rated
        none of them (the subset would really be a smaller one).

        Parameters
        ----------
        unit : tuple[int, ...]
            Rater positions, sorted.
        point, ci : float, BootstrapCIResult
            What this analysis computed for the unit.
        draws : np.ndarray | None
            The unit's draws, or None when keep_draws is false.

        Raises
        ------
        ValueError
            If the analyzer's point on the unit's data differs from the
            analysis's point: the two computed different things, which is a
            bug.
        """
        import dataclasses

        from ...contracts import EMPTY
        from ...intermediates.not_computed import NotComputed
        from ...semantics.ratings_semantics import RatingsSemantics
        from ..agreement_functions.rater_combinations import unit_problem

        semantics = self.ratings.semantics
        names = [str(semantics.raters[i]) for i in unit]
        columns = self.ratings.category_index[:, list(unit)]
        present = columns != EMPTY
        rows = np.flatnonzero(present.sum(axis=1) >= 2)

        # The same rule that made this unit's value NaN
        reason = unit_problem(present, names)
        if reason is not None:
            if not np.isnan(point):
                raise ValueError(
                    f"{_ERROR_PREFIX}: raters {', '.join(names)}: {reason}, "
                    f"yet this analysis computed {point!r} for them; this is "
                    f"a bug."
                )
            return NotComputed(reason)

        unit_ratings = Ratings(
            columns[rows],
            RatingsSemantics(
                raters=names,
                categories=semantics.categories,
                unique_subjects=semantics.unique_subjects,
            ),
            self.ratings.subject_ids[rows],
        )
        config = dataclasses.replace(
            self.config.config_of_coefficient, injection=None
        )
        analyzer = _analyzer_class(self.config.coefficient_to_compute)(
            unit_ratings, config
        )

        # The analyzer and this analysis must agree on the unit's coefficient
        theirs = analyzer._estimate_from_array(unit_ratings.category_index)
        if not (
            (np.isnan(theirs) and np.isnan(point))
            or np.isclose(theirs, point, rtol=1e-12, atol=1e-12)
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: raters {', '.join(names)}: the analyzer "
                f"computes {theirs!r} on their subjects, but this analysis "
                f"computed {point!r}. They should be the same number; this "
                f"is a bug."
            )

        if draws is None:
            draws = NotComputed("draws not kept: keep_draws false")
        return analyzer._result(point, ci, draws)

    def _unit_results(
        self, values: dict, subjects: dict, cis: dict, draws: dict | None,
    ) -> dict:
        """
        Every pair's or subset's result, in unit order.

        ``subjects`` is Computation A's count of each unit's usable subjects;
        each result's balance profile must agree with it, since both apply
        the same rule (rated by at least two of the unit's raters).
        """
        from ...intermediates.not_computed import NotComputed

        results = {}
        for unit in sorted(values):
            result = self._unit_result(
                unit, values[unit], cis[unit],
                None if draws is None else draws[unit],
            )
            if (
                not isinstance(result, NotComputed)
                and result.balance_profile.n_subjects != subjects[unit]
            ):
                raise ValueError(
                    f"{_ERROR_PREFIX}: raters {unit}: the result was built on "
                    f"{result.balance_profile.n_subjects} subjects, but "
                    f"{subjects[unit]} are rated by at least two of them; "
                    f"this is a bug."
                )
            results[unit] = result
        return results

    def _resolve_weights(self) -> Weights | None:
        """
        Read the weights from the semantics of the ratings.

        Follows the same rules as the analyzers: weights are declared in the
        semantics; a coefficient with a weighted form
        uses them; a coefficient without one refuses ordered or weighted
        categories rather than quietly ignoring them.

        Returns
        -------
        Weights | None
            The declared weights, or None for unordered categories.

        Raises
        ------
        ValueError
            If the coefficient has no weighted form and the semantics
            declares ordered or weighted categories.
        """
        from ...intermediates.estimator_parts.weights import Weights
        from ...semantics.categories import DistanceCategories

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
        return Weights.build_from_categories(categories)

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
