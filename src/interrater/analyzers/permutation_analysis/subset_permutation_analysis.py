"""Subset permutation analysis -- runs a coefficient across every
C(n_raters, k) rater subset, for one or more values of k.

I execute the three computations described in README.md, once per k:

- A (observed): coefficient for each subset, subjects per subset, mean.
- B (per-subset bootstrap): resample subjects within each subset, produce
  per-subset distributions and CIs.
- C (mean bootstrap): correlated (global resample, all subsets on the
  same draw) or independent (sample from B's per-subset distributions).

Each k is independent -- "average agreement among pairs" and "average
agreement among triples" are different quantities, so each gets its own
SubsetPermutationResult with its own mean bootstrap.

All three computations use the same estimator callable
(EstimatorInput -> float), resolved by the parent.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ...contracts import EMPTY
from ...intermediates.results.permutation_results.subset_permutation_result import (
    MAX_SUBSETS_FOR_DISTRIBUTIONS,
)
from ...intermediates.results.ci_results import BootstrapCIResult
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

    I inherit from PermutationAnalysis -- validation, EstimatorInput,
    estimator resolution, and ``_compute_ci`` are already done. I own
    ``run()``, which loops over each k in ``self.ks``, executes
    Computations A, B, and C, and returns a ``MultiKSubsetResult``
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
        validates each k against the data (k <= n_raters).

        Parameters
        ----------
        ratings : Ratings
            The ratings data.
        config : SubsetPermutationConfig
            Configuration for subset permutation analysis.

        Raises
        ------
        TypeError
            If config is not SubsetPermutationConfig.
        ValueError
            If any k > n_raters.
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
            per_k_results=per_k_results,
            ks=list(self.ks),
        )

    # Private helpers ---------------------------------------------------- #

    def _normalize_ks(
        self, subset_k: int | list[int] | str, n_raters: int
    ) -> list[int]:
        """
        Normalize subset_k into a sorted list and validate against data.

        Parameters
        ----------
        subset_k : int | list[int] | str
            Raw value from config. Config already validated structure
            (types, >= 2, no dupes). I normalize and check k <= n_raters.
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
            Observed values, optional per-subset bootstrap, optional mean
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

        # Computation B -- always
        subset_distributions, subset_cis = (
            self._compute_per_subset_bootstrap(k, subset_values)
        )

        # Computation C -- if mean bootstrap requested
        mean_distribution = None
        mean_ci = None
        mean_method = None
        if bmc is not None:
            mean_method = bmc.method
            # The subsets the point estimate averaged. Every draw uses
            # exactly these, so it estimates the same quantity.
            defined_subsets = sorted(
                s for s, v in subset_values.items() if not np.isnan(v)
            )
            mean_distribution, mean_ci = self._compute_mean_bootstrap(
                k, subset_distributions, mean_over_subsets, defined_subsets
            )

        # Decide whether to store per-subset distributions
        n_subsets = len(subset_values)
        store_distributions = n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS

        return SubsetPermutationResult(
            coefficient_name=self.config.coefficient_to_compute,
            ci_method=coc.ci_method,
            balance_profile=self.balance_profile,
            subset_values=subset_values,
            subjects_per_subset=subjects_per_subset,
            n_subsets=n_subsets,
            subset_k=k,
            mean_over_subsets=mean_over_subsets,
            bootstrap_mean_method=mean_method,
            bootstrap_ci_of_mean=mean_ci,
            bootstrap_distribution_of_mean=mean_distribution,
            bootstrap_cis_per_subset=subset_cis,
            bootstrap_distributions_per_subset=(
                subset_distributions if store_distributions else None
            ),
        )

    def _compute_observed(
        self, k: int,
    ) -> tuple[dict, dict, float]:
        """
        Computation A: observed values on actual data.

        Calls ``over_subsets`` for per-subset coefficients,
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
        from ..agreement_functions.rater_combinations import (
            mean_over_subsets,
            over_subsets,
            subjects_per_subset,
        )

        subset_values = over_subsets(
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
    ) -> tuple[dict, dict]:
        """
        Computation B: per-subset bootstrap.

        For each subset: select the k rater columns, find the subjects this
        subset can be measured on, resample those with replacement, and
        compute the coefficient on each resample. Produces a bootstrap
        distribution and a confidence interval per subset.

        **Which subjects a subset can be measured on.** A subject
        contributes to agreement only if at least two of that subset's
        raters rated it. Fewer than two ratings is not disagreement, it is
        absence -- Fleiss makes this explicit, since the per-subject term
        ``(sum_j n_ij^2 - n_i) / (n_i (n_i - 1))`` is 0/0 at ``n_i = 1``.

        So the pool is "at least two of the k". At k=2 that means both
        raters, which is what the pairwise analysis does; above k=2 it
        admits unbalanced subsets rather than demanding a complete
        rectangle. One rule, applied at every k.

        That pool is also what Computation A measures: ``over_subsets``
        hands the whole array to the estimator and lets it handle the gaps.
        The bootstrap must resample that same pool, because an interval has
        to be around the thing the point estimate measured. Draw from the
        complete cases instead and you have quietly switched estimands; on
        partially crossed data the mistake surfaces as an interval that
        does not contain its own point estimate.

        Seeding: ``np.random.SeedSequence`` spawns one child seed per
        subset from ``config_of_coefficient.random_state``, so each
        subset gets a reproducible but independent RNG stream.

        Subsets with fewer than 2 measurable subjects have nothing to
        resample. They receive a NaN distribution and an interval whose
        ``method_used`` is None -- no method ran, which is a different
        statement from an interval that could not be computed.

        Parameters
        ----------
        k : int
            Subset size.
        subset_values : dict
            Point estimates from Computation A, used for CI computation.

        Returns
        -------
        tuple[dict, dict]
            (subset_distributions, subset_cis)
            distributions: dict[subset -> np.ndarray of shape (n_resamples,)]
            cis: dict[subset -> BootstrapCIResult]
        """
        coc = self.config.config_of_coefficient
        n_resamples = coc.n_resamples
        alpha = coc.alpha
        ci_method = coc.ci_method

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

            # At least two of these k raters rated it -- see the docstring
            # for why two, and why not all k.
            present = subset_data.category_index != EMPTY
            measurable = np.flatnonzero(present.sum(axis=1) >= 2)

            if measurable.size < 2:
                # Fewer than two measurable subjects leaves no coefficient
                # to resample. Every draw fails, so no method runs -- which
                # the result states rather than implying with a pair of
                # NaNs.
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
                draws[b] = self.estimator(
                    subset_data.take_subjects(resampled)
                )

            distributions[subset] = draws
            point = subset_values[subset]

            # The jackknife leaves out the subjects this subset was
            # computed from -- `measurable`, the same pool the bootstrap
            # draws from. A subject that contributes nothing to the
            # coefficient cannot influence it, so including it would add a
            # deviation representing no influence at all.
            def _jackknife(sd=subset_data, sh=measurable) -> np.ndarray:
                return np.array([
                    self.estimator(sd.take_subjects(np.delete(sh, i)))
                    for i in range(sh.size)
                ], dtype=np.float64)

            cis[subset] = self._compute_ci(
                draws, point, ci_method, alpha, jackknife_fn=_jackknife
            )

        return distributions, cis

    def _compute_mean_bootstrap(
        self,
        k: int,
        subset_distributions: dict,
        mean_over_subsets: float,
        defined_subsets: list,
    ) -> tuple[np.ndarray, tuple[float, float]]:
        """
        Computation C: mean bootstrap.

        Two paths based on ``bootstrap_mean_config.method``:

        **Correlated:** Resample subjects globally. For each draw,
        recompute all C(n,k) subset coefficients on the same resampled
        subjects, then average. Preserves the correlation between subsets
        sharing raters. Cost: n_resamples x C(n,k).

        **Independent:** Sample from B's per-subset distributions. For
        each draw, pick one value from each subset's bootstrap
        distribution, average. Treats subsets as independent -- the CI
        will be narrower than correlated and may underestimate variance.
        Cost: 0 estimator calls (numpy only).

        Parameters
        ----------
        k : int
            Subset size.
        subset_distributions : dict
            Per-subset bootstrap distributions from B. Required for the
            independent path; ignored for correlated.
        mean_over_subsets : float
            Point estimate of the mean from A, used for CI computation.
        defined_subsets : list
            The subsets defined at the point estimate. A draw in which any
            of them is undefined is NaN, so the interval counts it as failed
            rather than averaging the other subsets.

        Returns
        -------
        tuple[np.ndarray, tuple[float, float]]
            (mean_distribution, mean_ci)
            distribution: shape (n_resamples,)
            ci: (lower, upper)
        """
        # The mean uses the coefficient's draw count, interval method, level
        # and seed. The per-combination draws use child streams spawned from
        # that seed, so the mean's stream is distinct from all of them.
        bmc = self.config.bootstrap_mean_config
        coc = self.config.config_of_coefficient
        n_resamples = coc.n_resamples
        rng = np.random.default_rng(coc.random_state)

        if bmc.method == "correlated":
            mean_distribution = self._mean_bootstrap_correlated(
                k, defined_subsets, n_resamples, rng
            )
        elif bmc.method == "independent":
            mean_distribution = self._mean_bootstrap_independent(
                subset_distributions, defined_subsets, n_resamples, rng
            )
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown bootstrap mean method "
                f"'{bmc.method}'. Expected 'correlated' or "
                f"'independent'."
            )

        # Only the correlated mean has a jackknife. It resamples subjects, so
        # leaving one subject out is meaningful. The independent mean
        # resamples from the per-subset distributions and never touches a
        # subject, so there is nothing to leave out.
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
        Leave-one-subject-out means, for the correlated mean's acceleration.

        The counterpart of ``_mean_bootstrap_correlated``: that one resamples
        all subjects with replacement and averages every k-subset's
        coefficient on each draw, so this one drops a single subject and does
        the same. The unit left out matches the unit resampled, which is what
        makes the acceleration describe the right quantity.

        Parameters
        ----------
        k : int
            Subset size, passed through to over_subsets.
        defined_subsets : list
            The subsets defined at the point estimate; see
            ``mean_over_same_subsets``.

        Returns
        -------
        np.ndarray
            Shape (n_subjects,). NaN entries mark leave-one-out data on which
            a defined subset became undefined.
        """
        from ..agreement_functions.rater_combinations import (
            mean_over_same_subsets,
            over_subsets,
        )

        n_subjects = self.estimator_input.n_subjects
        all_rows = np.arange(n_subjects)

        return np.array([
            mean_over_same_subsets(
                over_subsets(
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
        from ..agreement_functions.rater_combinations import (
            mean_over_same_subsets,
            over_subsets,
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
            sub_values = over_subsets(
                resampled_data, self.estimator, k
            )
            mean_draws[b] = mean_over_same_subsets(sub_values, defined_subsets)

        return mean_draws

    def _mean_bootstrap_independent(
        self,
        subset_distributions: dict,
        defined_subsets: list,
        n_resamples: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """
        Independent mean bootstrap: sample from B's per-subset distributions.

        For each draw: pick one value from each defined subset's bootstrap
        distribution and average them. A draw that picks a failed (NaN)
        value for any subset is NaN, so it counts as failed rather than
        averaging the other subsets. Treats subsets as independent, so the
        CI will typically be narrower than correlated -- on real data, too
        narrow.

        Parameters
        ----------
        subset_distributions : dict
            Per-subset distributions from B.
        defined_subsets : list
            The subsets defined at the point estimate.
        n_resamples : int
            Number of bootstrap draws.
        rng : np.random.Generator
            Seeded random number generator.

        Returns
        -------
        np.ndarray
            Mean coefficient per draw, shape (n_resamples,).
        """
        subsets = defined_subsets
        if not subsets:
            return np.full(n_resamples, np.nan)
        # Stack per-subset distributions: (n_subsets, n_resamples_per_subset)
        stacked = np.array(
            [subset_distributions[s] for s in subsets]
        )
        n_subsets = stacked.shape[0]
        n_per_subset = stacked.shape[1]

        # Vectorized: draw one index per subset per bootstrap draw
        indices = rng.integers(
            0, n_per_subset, size=(n_resamples, n_subsets)
        )
        # sampled[b, s] = stacked[s, indices[b, s]]
        sampled = stacked[
            np.arange(n_subsets)[np.newaxis, :], indices
        ]
        # Mean across subsets; any NaN makes the draw NaN
        mean_draws = sampled.mean(axis=1).astype(np.float64)

        return mean_draws
