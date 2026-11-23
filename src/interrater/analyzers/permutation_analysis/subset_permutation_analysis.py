"""Subset permutation analysis -- runs a coefficient across every
C(n_raters, k) rater subset, for one or more values of k.

Runs the three computations described in README.md, once per k:

- A (observed): coefficient for each subset, subjects per subset, mean.
- B (per-subset bootstrap): resample subjects within each subset, produce
  per-subset distributions and CIs.
- C (mean bootstrap, optional): correlated (global resample, all subsets on
  the same draw) or independent (sample from B's per-subset distributions,
  picked while B runs; see IndependentMeanDraws).

Each k is independent -- "average agreement among pairs" and "average
agreement among triples" are different quantities, so each gets its own
SubsetPermutationResult with its own mean bootstrap.

All three computations use the same estimator callable
(EstimatorInput -> float), resolved by the parent.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ...contracts import EMPTY, GUARD_MAX_SUBSETS
from ...intermediates.results.ci_results import BootstrapCIResult
from .independent_mean_draws import IndependentMeanDraws
from .permutation_analysis import PermutationAnalysis

if TYPE_CHECKING:
    from ...data_objects.ratings import Ratings
    from ...intermediates.results.permutation_results.multi_k_subset_result import (
        MultiKSubsetResult,
    )
    from ...intermediates.results.permutation_results.subset_permutation_result import (
        SubsetPermutationResult,
    )
    from ..configs.subset_permutation_config import SubsetPermutationConfig

_ERROR_PREFIX = "[SubsetPermutationAnalysis] Error"


class SubsetPermutationAnalysis(PermutationAnalysis):
    """
    I run a coefficient across every C(n_raters, k) rater subset, for
    one or more values of k.

    The parent already validated the ratings and config, applied the
    injection, built the profiles and the EstimatorInput, and resolved the
    estimator. I own ``run()``, which loops over each k in ``self.ks``,
    executes Computations A, B, and C, and returns a ``MultiKSubsetResult``
    wrapping one ``SubsetPermutationResult`` per k.

    A and B run always (the analyzer config requires n_resamples >= 1).
    C runs if ``bootstrap_mean_config`` is not None.

    Attributes
    ----------
    ks : list[int]
        The subset sizes to evaluate, sorted. Normalized from
        ``config.subset_k``: int -> [int], list -> sorted copy,
        "all" -> range(2, n_raters + 1).
    """

    # Attribute declaration
    ks: list[int]

    # Constructor -------------------------------------------------------- #

    def __init__(
        self, ratings: Ratings, config: SubsetPermutationConfig
    ) -> None:
        """
        Initialize subset permutation analysis.

        Normalizes ``config.subset_k`` into a sorted list of ints and
        validates each k against the data (k <= n_raters). Every k is
        checked against GUARD_MAX_SUBSETS, and the whole run against
        config.max_evaluations, before any k runs.

        Parameters
        ----------
        ratings : Ratings
            The ratings data.
        config : SubsetPermutationConfig
            Configuration for subset permutation analysis.

        Raises
        ------
        TypeError
            If config is not SubsetPermutationConfig, or ratings is not
            Ratings.
        ValueError
            If any k > n_raters, any k gives more than GUARD_MAX_SUBSETS
            subsets, or the run would compute the coefficient more than
            config.max_evaluations times; if
            config_of_coefficient.balanced_data_only is true and a subset is
            unbalanced; or if the parent refuses the data (Conger's kappa on
            incomplete data, balanced_data_only on an unbalanced panel, or
            weights on a coefficient without a weighted form).
        """
        from ..configs.subset_permutation_config import (
            SubsetPermutationConfig as _SPC,
        )

        if not isinstance(config, _SPC):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be "
                f"SubsetPermutationConfig, "
                f"got {type(config).__name__}"
            )

        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(
                f"{_ERROR_PREFIX}: parent init failed: {e}"
            ) from e

        # self.ratings, not the argument: the parent has applied any
        # injection, and k may count the injected raters too.
        self.ks = self._normalize_ks(config.subset_k, self.ratings.n_raters)

        # Before anything runs: every k within the subset cap, and the whole
        # run within max_evaluations. Otherwise subset_k: "all" ran the small
        # k in full before a large one hit the cap.
        units_per_k = {
            k: math.comb(self.ratings.n_raters, k) for k in self.ks
        }
        for k, count in units_per_k.items():
            if count > GUARD_MAX_SUBSETS:
                raise ValueError(
                    f"{_ERROR_PREFIX}: {self.ratings.n_raters} raters taken "
                    f"{k} at a time is {count:,} subsets, past the cap of "
                    f"{GUARD_MAX_SUBSETS:,}. Nothing has run. A smaller "
                    f"subset_k, or leave-one-out, will usually answer the "
                    f"same question."
                )
        try:
            self._check_cost(units_per_k, "subsets")
        except ValueError as e:
            raise ValueError(f"{_ERROR_PREFIX}: {e}") from e

        # The coefficient's own flag, applied to each subset
        self._check_subsets_balanced()

    # Public API --------------------------------------------------------- #

    def run(self) -> MultiKSubsetResult:
        """
        Execute subset permutation analysis for all k values.

        Loops over ``self.ks``, running Computations A, B, and C for each
        k value. Wraps the per-k results into a ``MultiKSubsetResult``.

        Returns
        -------
        MultiKSubsetResult
            One SubsetPermutationResult per k, keyed by k.
        """
        from ...intermediates.results.permutation_results.multi_k_subset_result import (
            MultiKSubsetResult,
        )

        coc = self.config.config_of_coefficient

        per_k_results = {}
        for k in self.ks:
            per_k_results[k] = self._run_single_k(k)

        return MultiKSubsetResult(
            coefficient_name=self.config.coefficient_to_compute,
            ci_method=coc.ci_method,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            **self._weights_record(),
            per_k_results=per_k_results,
            ks=list(self.ks),
        )

    # Private helpers ---------------------------------------------------- #

    def _check_subsets_balanced(self) -> None:
        """
        Refuse unbalanced subsets when the coefficient asks for balance.

        ``config_of_coefficient.balanced_data_only`` has the same meaning in
        the coefficient's own analyzer. It is applied to each subset: the
        subjects at least two of its raters rated must all have the same
        number of ratings among them. This is relevant as a balanced panel
        can still give unbalanced subsets, when subjects went to different
        raters.

        Raises
        ------
        ValueError
            Naming the first unbalanced subset and how many there are.
        """
        coc = self.config.config_of_coefficient
        if not getattr(coc, "balanced_data_only", False):
            return
        from .rater_combinations import rater_subsets

        present = self.ratings.category_index != EMPTY
        names = [str(r) for r in self.ratings.semantics.raters]
        unbalanced = []
        for k in self.ks:
            for subset in rater_subsets(self.ratings.n_raters, k):
                counts = present[:, list(subset)].sum(axis=1)
                counts = counts[counts >= 2]
                if counts.size and counts.min() != counts.max():
                    unbalanced.append((subset, int(counts.min()), int(counts.max())))
        if not unbalanced:
            return
        subset, low, high = unbalanced[0]
        raise ValueError(
            f"{_ERROR_PREFIX}: config_of_coefficient.balanced_data_only is "
            f"true, but {len(unbalanced):,} subsets are unbalanced; the first, "
            f"raters {', '.join(names[r] for r in subset)}, has subjects with "
            f"{low} to {high} ratings among them. "
            f"{self.config.coefficient_to_compute} on such a subset is an "
            f"approximation. Set balanced_data_only: false inside "
            f"config_of_coefficient to allow it."
        )

    def _normalize_ks(
        self, subset_k: int | list[int] | str, n_raters: int
    ) -> list[int]:
        """
        Normalize subset_k into a sorted list and validate against data.

        Parameters
        ----------
        subset_k : int | list[int] | str
            Raw value from config. Config already validated structure
            (types, >= 2, no dupes); this method normalizes it and checks
            k <= n_raters.
        n_raters : int
            Number of raters in the data.

        Returns
        -------
        list[int]
            Sorted list of k values.

        Raises
        ------
        ValueError
            If any k > n_raters.
        """
        if isinstance(subset_k, int):
            ks = [subset_k]
        elif isinstance(subset_k, str):
            # "all" -- config already validated it's exactly "all"
            ks = list(range(2, n_raters + 1))
        else:
            # list[int] -- config already validated non-empty, all >= 2
            ks = sorted(subset_k)

        # Data-dependent validation: each k must be <= n_raters
        for k in ks:
            if k > n_raters:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subset_k={k} exceeds "
                    f"n_raters={n_raters}. Cannot form subsets "
                    f"larger than the number of raters."
                )

        return ks

    def _run_single_k(self, k: int) -> SubsetPermutationResult:
        """
        Run Computations A, B, and C for a single k value.

        Parameters
        ----------
        k : int
            Subset size for this run.

        Returns
        -------
        SubsetPermutationResult
            Observed values, per-subset bootstrap, and the optional mean
            bootstrap for this k.
        """
        from ...intermediates.results.permutation_results.subset_permutation_result import (
            SubsetPermutationResult,
        )

        coc = self.config.config_of_coefficient
        bmc = self.config.bootstrap_mean_config

        # Computation A -- always
        subset_values, subjects_per_subset, mean_over_subsets = (
            self._compute_observed(k)
        )

        # The subsets the point estimate averaged. Every draw of the mean
        # uses exactly these, so it estimates the same quantity.
        defined_subsets = sorted(
            s for s, v in subset_values.items() if not np.isnan(v)
        )

        # The independent mean picks from each subset's draws while B makes
        # them, so B can toss them when keep_draws is false.
        independent = (
            IndependentMeanDraws(defined_subsets, coc.n_resamples, coc.random_state)
            if bmc is not None and bmc.method == "independent" else None
        )

        # Computation B -- always
        subset_distributions, subset_cis = (
            self._compute_per_subset_bootstrap(k, subset_values, independent)
        )

        # Computation C -- if mean bootstrap requested
        mean_distribution = None
        mean_ci = None
        mean_method = None
        if bmc is not None:
            mean_method = bmc.method
            mean_distribution, mean_ci = self._compute_mean_bootstrap(
                k, independent, mean_over_subsets, defined_subsets
            )

        n_subsets = len(subset_values)

        subset_results = self._unit_results(
            subset_values, subjects_per_subset, subset_cis,
            subset_distributions,
        )

        return SubsetPermutationResult(
            coefficient_name=self.config.coefficient_to_compute,
            ci_method=coc.ci_method,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            **self._weights_record(),
            subset_results=subset_results,
            n_subsets=n_subsets,
            subset_k=k,
            mean_over_subsets=mean_over_subsets,
            bootstrap_mean_method=mean_method,
            bootstrap_ci_of_mean=mean_ci,
            bootstrap_distribution_of_mean=mean_distribution,
        )

    def _compute_observed(
        self, k: int,
    ) -> tuple[dict, dict, float]:
        """
        Computation A: observed values on actual data.

        Calls ``estimate_subsets`` for per-subset coefficients,
        ``subjects_per_subset`` for per-subset subject counts, and
        ``mean_over_subsets`` for the NaN-skipping mean across subsets.

        Parameters
        ----------
        k : int
            Subset size.

        Returns
        -------
        tuple[dict, dict, float]
            (subset_values, subjects_per_subset, mean_over_subsets)
        """
        from .rater_combinations import (
            mean_over_subsets,
            estimate_subsets,
            subjects_per_subset,
        )

        subset_values = estimate_subsets(
            self.estimator_input, self.estimator, k
        )
        sps = subjects_per_subset(
            self.estimator_input, k
        )
        mean = mean_over_subsets(subset_values)

        return subset_values, sps, mean

    def _compute_per_subset_bootstrap(
        self,
        k: int,
        subset_values: dict,
        independent: IndependentMeanDraws | None,
    ) -> tuple[dict | None, dict]:
        """
        Computation B: per-subset bootstrap.

        For each subset: select the k rater columns, find the subjects this
        subset can be measured on, resample those subjects with replacement,
        and compute the coefficient on each resample. Produces a bootstrap
        distribution and a confidence interval per subset.

        The subject pool is the subjects at least two of the subset's raters
        rated, the same pool Computation A measures (see README.md, "Which
        subjects a subset can be measured on").

        Seeding: ``np.random.SeedSequence`` spawns one child seed per
        subset from ``config_of_coefficient.random_state``, so each
        subset gets a reproducible but independent RNG stream.

        Subsets with fewer than 2 measurable subjects have nothing to
        resample. They receive a NaN distribution and an interval whose
        ``method_used`` is None -- no method ran, which is a different
        statement from an interval that could not be computed.

        Each subset's draws are used for its interval and, for the
        independent mean, picked from; then they are kept only if
        keep_draws is true. So with keep_draws false, memory holds one
        subset's draws at a time.

        Parameters
        ----------
        k : int
            Subset size.
        subset_values : dict
            Point estimates from Computation A, used for CI computation.
        independent : IndependentMeanDraws | None
            The independent mean, built as each subset's draws are made.
            None unless the independent mean was requested.

        Returns
        -------
        tuple[dict | None, dict]
            (subset_distributions, subset_cis)
            distributions: dict[subset -> np.ndarray of shape (n_resamples,)],
            or None when keep_draws is false
            cis: dict[subset -> BootstrapCIResult]
        """
        coc = self.config.config_of_coefficient
        n_resamples = coc.n_resamples
        alpha = coc.alpha
        ci_method = coc.ci_method
        keep_draws = self.config.keep_draws

        # Every draw and jackknife value obeys the one rule for whether a
        # subset has a coefficient (subset_problem)
        from .rater_combinations import estimate_subset

        # Spawn one child seed per subset for reproducibility
        ss = np.random.SeedSequence(coc.random_state)
        subsets = sorted(subset_values.keys())
        child_seeds = ss.spawn(len(subsets))

        distributions: dict[tuple, np.ndarray] = {}
        cis: dict[tuple, BootstrapCIResult] = {}

        for subset, seed in zip(subsets, child_seeds):
            rng = np.random.default_rng(seed)

            # Select the k rater columns
            subset_data = self.estimator_input.select_raters(
                np.array(subset)
            )

            # Subjects at least two of the subset's raters rated (README.md,
            # "Which subjects a subset can be measured on").
            present = subset_data.category_index != EMPTY
            measurable = np.flatnonzero(present.sum(axis=1) >= 2)

            if measurable.size < 2:
                # Fewer than two measurable subjects: no coefficient, so every
                # draw fails and no method runs. The result says so instead
                # of carrying bare NaNs.
                if keep_draws:
                    distributions[subset] = np.full(n_resamples, np.nan)
                cis[subset] = BootstrapCIResult(
                    lower=float("nan"),
                    upper=float("nan"),
                    alpha=alpha,
                    n_resamples=n_resamples,
                    n_failed=n_resamples,
                    method_intended=ci_method,
                    method_used=None,
                )
                continue

            # Bootstrap loop: resample measurable subjects with replacement
            draws = np.empty(n_resamples, dtype=np.float64)
            for b in range(n_resamples):
                resampled = rng.choice(
                    measurable, size=measurable.size, replace=True
                )
                draws[b] = estimate_subset(
                    subset_data.take_subjects(resampled), self.estimator
                )

            if independent is not None:
                independent.add(subset, draws)
            if keep_draws:
                distributions[subset] = draws
            point = subset_values[subset]

            # BCa's jackknife uses only the subset's pool, `measurable` (see
            # README.md, "BCa for a pair: Use only the subjects shared by
            # the pair").
            def _jackknife(sd=subset_data, sh=measurable) -> np.ndarray:
                return np.array([
                    estimate_subset(
                        sd.take_subjects(np.delete(sh, i)),
                        self.estimator,
                    )
                    for i in range(sh.size)
                ], dtype=np.float64)

            cis[subset] = self._compute_ci(
                draws, point, ci_method, alpha, jackknife_fn=_jackknife
            )

        return (distributions if keep_draws else None), cis

    def _compute_mean_bootstrap(
        self,
        k: int,
        independent: IndependentMeanDraws | None,
        mean_over_subsets: float,
        defined_subsets: list,
    ) -> tuple[np.ndarray, BootstrapCIResult]:
        """
        Computation C: mean bootstrap.

        Two paths based on ``bootstrap_mean_config.method``:

        **Correlated:** Resample subjects globally. For each draw,
        recompute all C(n_raters, k) subset coefficients on the same resampled
        subjects, then average. Preserves the correlation between subsets
        sharing raters. Cost: n_resamples x C(n_raters, k).

        **Independent:** Sample from B's per-subset distributions. For
        each draw, pick one value from each subset's bootstrap
        distribution, average. Ignores the covariance between subsets
        sharing raters, so when that covariance is positive the interval
        is too narrow.
        Cost: 0 estimator calls; the picks were made while B ran (see
        IndependentMeanDraws), so here they are only read off.

        Parameters
        ----------
        k : int
            Subset size.
        independent : IndependentMeanDraws | None
            The independent mean, built during B. Required for the
            independent path; None for correlated.
        mean_over_subsets : float
            Point estimate of the mean from A, used for CI computation.
        defined_subsets : list
            The subsets defined at the point estimate. A draw in which any
            of them is undefined is NaN, so the interval counts it as failed
            rather than averaging the other subsets.

        Returns
        -------
        tuple[np.ndarray, BootstrapCIResult]
            (mean_distribution, mean_ci)
            distribution: shape (n_resamples,)
            ci: the interval of the mean
        """
        # The mean uses the coefficient's draw count, interval method, level
        # and seed. The per-combination draws use child streams spawned from
        # that seed, so the mean's stream is distinct from all of them.
        bmc = self.config.bootstrap_mean_config
        coc = self.config.config_of_coefficient
        n_resamples = coc.n_resamples

        if bmc.method == "correlated":
            rng = np.random.default_rng(coc.random_state)
            mean_distribution = self._mean_bootstrap_correlated(
                k, defined_subsets, n_resamples, rng
            )
        elif bmc.method == "independent":
            mean_distribution = independent.draws()
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown bootstrap mean method "
                f"'{bmc.method}'. Expected 'correlated' or "
                f"'independent'."
            )

        # BCa needs a jackknife. The correlated mean resamples subjects, so
        # leaving out one subject at a time matches it. The independent mean
        # resamples per-subset draws, not subjects, so it has no jackknife and
        # _compute_ci refuses BCa for it.
        jackknife_fn = (
            (lambda: self._mean_jackknife_correlated(k, defined_subsets))
            if bmc.method == "correlated" else None
        )

        mean_ci = self._compute_ci(
            mean_distribution,
            mean_over_subsets,
            coc.ci_method,
            coc.alpha,
            jackknife_fn=jackknife_fn,
        )

        return mean_distribution, mean_ci

    def _mean_jackknife_correlated(
        self, k: int, defined_subsets: list
    ) -> np.ndarray:
        """
        Calculate leave-one-subject-out means, for the correlated mean's
        acceleration.

        For each subject in turn: remove it from the whole panel, compute every
        subset's coefficient, and average over defined_subsets. This matches
        ``_mean_bootstrap_correlated``, which resamples the whole panel.

        Parameters
        ----------
        k : int
            Subset size, passed through to estimate_subsets.
        defined_subsets : list
            The subsets defined at the point estimate; see
            ``mean_over_same_subsets``.

        Returns
        -------
        np.ndarray
            Shape (n_subjects,). NaN entries mark leave-one-out data on which
            a defined subset became undefined.
        """
        from .rater_combinations import (
            mean_over_same_subsets,
            estimate_subsets,
        )

        n_subjects = self.estimator_input.n_subjects
        all_rows = np.arange(n_subjects)

        return np.array([
            mean_over_same_subsets(
                estimate_subsets(
                    self.estimator_input.take_subjects(
                        np.delete(all_rows, i)
                    ),
                    self.estimator,
                    k,
                ),
                defined_subsets,
            )
            for i in range(n_subjects)
        ], dtype=np.float64)

    def _mean_bootstrap_correlated(
        self,
        k: int,
        defined_subsets: list,
        n_resamples: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """
        Correlated mean bootstrap: global subject resampling.

        For each draw: resample all subjects with replacement, compute
        every subset's coefficient on the same resampled data, average.
        The correlation between subsets sharing raters is preserved
        because they see the same resampled subjects.

        Parameters
        ----------
        k : int
            Subset size.
        defined_subsets : list
            The subsets defined at the point estimate. A draw in which any
            of them is undefined is NaN, not the mean of the others.
        n_resamples : int
            Number of bootstrap draws.
        rng : np.random.Generator
            Seeded random number generator.

        Returns
        -------
        np.ndarray
            Mean coefficient per draw, shape (n_resamples,).
        """
        from .rater_combinations import (
            mean_over_same_subsets,
            estimate_subsets,
        )

        n_subjects = self.estimator_input.n_subjects
        mean_draws = np.empty(n_resamples, dtype=np.float64)

        for b in range(n_resamples):
            resampled_rows = rng.choice(
                n_subjects, size=n_subjects, replace=True
            )
            resampled_data = self.estimator_input.take_subjects(
                resampled_rows
            )
            sub_values = estimate_subsets(
                resampled_data, self.estimator, k
            )
            mean_draws[b] = mean_over_same_subsets(sub_values, defined_subsets)

        return mean_draws
