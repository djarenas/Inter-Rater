"""Base class for permutation analysis across rater subsets.

PermutationAnalysis is not in the Analyzer tree. An Analyzer computes one
coefficient on one data slice; a permutation analysis enumerates pairs or
subsets of raters and computes the coefficient on each. Three computations:

- A (observed): the coefficient for each pair or subset, and their mean.
- B (per-unit bootstrap): resample subjects within each pair or subset,
  giving its distribution and interval.
- C (mean bootstrap, optional): correlated (resample subjects once, all
  units on the same draw) or independent (pick from B's draws).

README.md in this directory has the details: costs, which subjects each
unit uses, and why B and C resample differently.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from ...data_objects.ratings import Ratings
from ...intermediates.data_profiles.balance_profile import BalanceProfile
from ...intermediates.data_profiles.rater_profile import RaterProfile
from ...intermediates.estimator_parts.estimator_input import EstimatorInput
from ...intermediates.results.ci_results import BootstrapCIResult
from .rater_combinations import Estimator
from ..agreement_functions.group_agreement import (
    brennan_prediger, fleiss_kappa, gwet_ac1, scotts_pi,
)
from ..agreement_functions.individual_marginals import cohen_kappa, conger_kappa
from ..agreement_functions.krippendorff import krippendorff_alpha

if TYPE_CHECKING:
    from ...intermediates.estimator_parts.weights import Weights
    from ..configs.permutation_config import PermutationConfig

_ERROR_PREFIX = "[PermutationAnalysis] Error"


# Every coefficient function takes an EstimatorInput, K and weights included
_ESTIMATOR_MAP: dict[str, Estimator] = {
    "cohen_kappa": cohen_kappa,
    "scotts_pi": scotts_pi,
    "conger_kappa": conger_kappa,
    "gwet_ac1": gwet_ac1,
    "fleiss_kappa": fleiss_kappa,
    "brennan_prediger": brennan_prediger,
    "krippendorff_alpha": krippendorff_alpha,
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

    I do NOT have a run() method -- child classes define their own,
    returning the result type they produce (PairwisePermutationResult,
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
    rater_profile: RaterProfile
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
            If ratings is None or coefficient is unknown; if the
            coefficient is Conger's kappa and a rating is missing (after
            injection); if balanced_data_only is true and the panel is
            unbalanced; or if the coefficient has no weighted form and the
            semantics declares ordered or weighted categories.
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
        Apply config_of_coefficient's injection pipeline to the ratings.

        Reads ``self.config.config_of_coefficient.injection`` and iterates
        through each injection step, instantiating the corresponding
        injector and calling ``inject()`` to produce a new Ratings with
        synthetic raters appended. Each step chains -- the output Ratings
        of one step becomes the input to the next.

        Replaces my ``ratings`` attribute. Does nothing if injection is None.

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
        correlated mean (n_resamples draws over every unit). With BCa, each
        interval adds a jackknife, counted as one evaluation per subject in
        the panel (an upper bound, since a unit is jackknifed only over its
        own subjects). So the count is "at most", and the guard can only
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
        Build one pair's or subset's result with the coefficient's own analyzer.

        The analyzer is built on the unit's usable subjects (rated by at least
        two of its raters; for a pair, both) with injection off, since the
        panel's injection already ran. It is handed this analysis's point,
        interval and draws through its ``_result``, so each unit gets the same
        result type ``estimate()`` returns. The point is checked against the
        analyzer's own; the interval and draws are this analysis's.

        Parameters
        ----------
        unit : tuple[int, ...]
            Rater positions, sorted.
        point, ci : float, BootstrapCIResult
            What this analysis computed for the unit.
        draws : np.ndarray | None
            The unit's draws, or None when keep_draws is false.

        Returns
        -------
        PanelResult | NotComputed
            The analyzer's result, or NotComputed with the reason when the unit
            has no coefficient: fewer than 2 usable subjects, or a subset rater
            who rated none of them.

        Raises
        ------
        ValueError
            If the analyzer's point on the unit's data differs from this
            analysis's point, which is a bug.
        """
        import dataclasses

        from ...contracts import EMPTY
        from ...intermediates.not_computed import NotComputed
        from ...semantics.ratings_semantics import RatingsSemantics
        from .rater_combinations import subset_problem

        semantics = self.ratings.semantics
        names = [str(semantics.raters[i]) for i in unit]
        columns = self.ratings.category_index[:, list(unit)]
        present = columns != EMPTY
        rows = np.flatnonzero(present.sum(axis=1) >= 2)

        # The same rule that made this unit's value NaN
        reason = subset_problem(present, names)
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
        theirs = analyzer._estimate_from_input(analyzer.estimator_input)
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
            and not _analyzer_class(coefficient)._ACCEPTS_WEIGHTS
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: {coefficient} has no weighted form, but "
                f"the semantics declares categories of kind "
                f"{categories._KIND!r}. To compute it, use a semantics with "
                f"kind: unordered for these ratings."
            )
        return Weights.build_from_categories(categories)

    def _weights_record(self) -> dict:
        """
        The weighting my result records: weights_type and weight_matrix.

        Read from the weights every pair or subset was computed with, so the
        result cannot say one thing while its units did another.
        """
        weights = self.estimator_input.weights
        if weights is None:
            return {"weights_type": "identity", "weight_matrix": None}
        return {"weights_type": weights.kind, "weight_matrix": weights.matrix.tolist()}

    def _resolve_estimator(self) -> Estimator:
        """
        Map coefficient name to the estimator callable.

        The estimator is a pure function (EstimatorInput -> float) used
        by all three computations: A (observed values via estimate_pairs or
        estimate_subsets), B (per-unit bootstrap loop in the child), and C
        (correlated global resample loop or independent sampling from B's
        output).

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

        Shared by child classes for both per-pair CIs and mean CIs.
        Delegates to the CI functions in ci_methods.py.

        Parameters
        ----------
        distribution : np.ndarray
            Bootstrap values (1-D).
        point : float
            The point estimate.
        ci_method : str
            "percentile" or "bca". The package does not implement a parametric
            null for permutation analyses.
        alpha : float
            Significance level (e.g., 0.05 for 95% CI).
        jackknife_fn : Callable[[], np.ndarray] | None
            Produces leave-one-out estimates when called. A callable rather
            than an array because only BCa needs them, and here they are
            expensive -- per pair, one coefficient evaluation per subject.

            The jackknife must leave out the subjects that were resampled:
            a pair's shared subjects for a per-pair interval (see README.md,
            "BCa for a pair: Use only the subjects shared by the pair"), all
            subjects for the correlated mean. The caller knows which, so the
            caller supplies the function.

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
                f"Expected 'percentile' or 'bca'. The package implements "
                f"no parametric null for permutation analyses."
            )
