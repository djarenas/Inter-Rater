"""Pairwise permutation analysis -- runs a coefficient across every
C(n_raters, 2) rater pair.

I execute the three computations described in README.md:

- A (observed): coefficient for each pair, subjects per pair, mean.
- B (per-pair bootstrap): resample subjects within each pair, produce
  per-pair distributions and CIs.
- C (mean bootstrap): correlated (global resample, all pairs on the
  same draw) or independent (sample from B's per-pair distributions).

All three use the same estimator callable (EstimatorInput -> float),
resolved by the parent.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ...contracts import EMPTY
from ...intermediates.results.permutation_results.pairwise_permutation_result import (
    MAX_PAIRS_FOR_DISTRIBUTIONS,
)
from ...intermediates.results.ci_results import BootstrapCIResult
from .permutation_analysis import PermutationAnalysis

if TYPE_CHECKING:
    from ...data_objects.ratings import Ratings
    from ...intermediates.results.permutation_results.pairwise_permutation_result import (
        PairwisePermutationResult,
    )
    from ..configs.pairwise_permutation_config import PairwisePermutationConfig

_ERROR_PREFIX = "[PairwisePermutationAnalysis] Error"


class PairwisePermutationAnalysis(PermutationAnalysis):
    """
    I run a pairwise coefficient across every C(n_raters, 2) rater pair.

    I inherit from PermutationAnalysis -- validation, EstimatorInput,
    estimator resolution, and ``_compute_ci`` are already done. I own
    ``run()``, which executes Computations A, B, and C and returns a
    ``PairwisePermutationResult``.

    A and B run always (the analyzer config requires n_resamples >= 1).
    C runs if ``bootstrap_mean_config`` is not None.
    """

    # Constructor -------------------------------------------------------- #

    def __init__(
        self, ratings: Ratings, config: PairwisePermutationConfig
    ) -> None:
        """
        Initialize pairwise permutation analysis.

        Parameters
        ----------
        ratings : Ratings
            The ratings data.
        config : PairwisePermutationConfig
            Configuration for pairwise permutation analysis.

        Raises
        ------
        TypeError
            If config is not PairwisePermutationConfig.
        """
        from ..configs.pairwise_permutation_config import (
            PairwisePermutationConfig as _PPC,
        )

        if not isinstance(config, _PPC):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be "
                f"PairwisePermutationConfig, "
                f"got {type(config).__name__}"
            )

        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as e:
            raise type(e)(
                f"{_ERROR_PREFIX}: parent init failed: {e}"
            ) from e

    # Public API --------------------------------------------------------- #

    def run(self) -> PairwisePermutationResult:
        """
        Execute pairwise permutation analysis.

        Runs Computation A (always), B (if per-pair bootstrap requested),
        and C (if mean bootstrap requested). Assembles all results into
        a ``PairwisePermutationResult``.

        Returns
        -------
        PairwisePermutationResult
            Observed values, optional per-pair bootstrap, optional mean
            bootstrap.
        """
        from ...intermediates.results.permutation_results.pairwise_permutation_result import (
            PairwisePermutationResult,
        )

        coc = self.config.config_of_coefficient
        bmc = self.config.bootstrap_mean_config

        # Computation A -- always
        pairwise_values, subjects_per_pair, mean_over_pairs = (
            self._compute_observed()
        )

        # Computation B -- always
        pairwise_distributions, pairwise_cis = (
            self._compute_per_pair_bootstrap(pairwise_values)
        )

        # Computation C -- if mean bootstrap requested
        mean_distribution = None
        mean_ci = None
        mean_method = None
        if bmc is not None:
            mean_method = bmc.method
            # The pairs the point estimate averaged. Every draw uses exactly
            # these, so it estimates the same quantity.
            defined_pairs = sorted(
                p for p, v in pairwise_values.items() if not np.isnan(v)
            )
            mean_distribution, mean_ci = self._compute_mean_bootstrap(
                pairwise_distributions, mean_over_pairs, defined_pairs
            )

        # Decide whether to store per-pair distributions
        n_pairs = len(pairwise_values)
        store_distributions = n_pairs <= MAX_PAIRS_FOR_DISTRIBUTIONS

        return PairwisePermutationResult(
            coefficient_name=self.config.coefficient_to_compute,
            ci_method=coc.ci_method,
            balance_profile=self.balance_profile,
            pairwise_values=pairwise_values,
            subjects_per_pair=subjects_per_pair,
            n_pairs=n_pairs,
            mean_over_pairs=mean_over_pairs,
            bootstrap_mean_method=mean_method,
            bootstrap_ci_of_mean=mean_ci,
            bootstrap_distribution_of_mean=mean_distribution,
            bootstrap_cis_per_pair=pairwise_cis,
            bootstrap_distributions_per_pair=(
                pairwise_distributions if store_distributions else None
            ),
        )

    # Private helpers ---------------------------------------------------- #

    def _compute_observed(
        self,
    ) -> tuple[dict, dict, float]:
        """
        Computation A: observed values on actual data.

        Calls ``over_pairs`` for per-pair coefficients,
        ``subjects_per_subset`` for per-pair subject counts, and
        ``mean_over_subsets`` for the NaN-skipping mean across pairs.

        Returns
        -------
        tuple[dict, dict, float]
            (pairwise_values, subjects_per_pair, mean_over_pairs)
        """
        from ..agreement_functions.rater_combinations import (
            mean_over_subsets,
            over_pairs,
            subjects_per_subset,
        )

        pairwise_values = over_pairs(
            self.estimator_input, self.estimator
        )
        subjects_per_pair = subjects_per_subset(
            self.estimator_input, 2
        )
        mean_over_pairs = mean_over_subsets(pairwise_values)

        return pairwise_values, subjects_per_pair, mean_over_pairs

    def _compute_per_pair_bootstrap(
        self,
        pairwise_values: dict,
    ) -> tuple[dict, dict]:
        """
        Computation B: per-pair bootstrap.

        For each pair: select the two rater columns, find the subjects
        both raters rated, resample those subjects with replacement,
        compute the coefficient on each resample. Produces a bootstrap
        distribution and CI per pair.

        Seeding: ``np.random.SeedSequence`` spawns one child seed per
        pair from ``config_of_coefficient.random_state``, so each pair
        gets a reproducible but independent RNG stream.

        Pairs with fewer than 2 shared subjects cannot be bootstrapped
        and receive NaN distributions and CIs.

        Parameters
        ----------
        pairwise_values : dict
            Point estimates from Computation A, used for CI computation.

        Returns
        -------
        tuple[dict, dict]
            (pairwise_distributions, pairwise_cis)
            distributions: dict[pair -> np.ndarray of shape (n_resamples,)]
            cis: dict[pair -> (lower, upper)]
        """
        coc = self.config.config_of_coefficient
        n_resamples = coc.n_resamples
        alpha = coc.alpha
        ci_method = coc.ci_method

        # Spawn one child seed per pair for reproducibility
        ss = np.random.SeedSequence(coc.random_state)
        pairs = sorted(pairwise_values.keys())
        child_seeds = ss.spawn(len(pairs))

        distributions: dict[tuple, np.ndarray] = {}
        cis: dict[tuple, tuple[float, float]] = {}

        for pair, seed in zip(pairs, child_seeds):
            rng = np.random.default_rng(seed)

            # Select the two rater columns
            pair_data = self.estimator_input.select_raters(
                np.array(pair)
            )

            # Find subjects both raters rated
            present = pair_data.category_index != EMPTY
            shared = np.flatnonzero(present[:, 0] & present[:, 1])

            if shared.size < 2:
                # Two raters with fewer than two subjects in common have no
                # coefficient to resample. Every draw fails, so no method
                # runs -- which the result states rather than implying with
                # a pair of NaNs.
                distributions[pair] = np.full(n_resamples, np.nan)
                cis[pair] = BootstrapCIResult(
                    lower=float("nan"),
                    upper=float("nan"),
                    alpha=alpha,
                    n_resamples=n_resamples,
                    n_failed=n_resamples,
                    method_intended=ci_method,
                    method_used=None,
                )
                continue

            # Bootstrap loop: resample shared subjects with replacement
            draws = np.empty(n_resamples, dtype=np.float64)
            for b in range(n_resamples):
                resampled = rng.choice(
                    shared, size=shared.size, replace=True
                )
                draws[b] = self.estimator(
                    pair_data.take_subjects(resampled)
                )

            distributions[pair] = draws
            point = pairwise_values[pair]

            # The jackknife leaves out the subjects this pair was computed
            # from -- `shared`, not every subject. A subject the pair never
            # rated cannot influence the pair's coefficient, so including it
            # would contribute a deviation representing no influence at all.
            def _jackknife(pd=pair_data, sh=shared) -> np.ndarray:
                return np.array([
                    self.estimator(pd.take_subjects(np.delete(sh, i)))
                    for i in range(sh.size)
                ], dtype=np.float64)

            cis[pair] = self._compute_ci(
                draws, point, ci_method, alpha, jackknife_fn=_jackknife
            )

        return distributions, cis

    def _compute_mean_bootstrap(
        self,
        pairwise_distributions: dict,
        mean_over_pairs: float,
        defined_pairs: list,
    ) -> tuple[np.ndarray, tuple[float, float]]:
        """
        Computation C: mean bootstrap.

        Two paths based on ``bootstrap_mean_config.method``:

        **Correlated:** Resample subjects globally. For each draw,
        recompute all C(n,2) pair coefficients on the same resampled
        subjects, then average. Preserves the correlation between pairs
        sharing a rater. Cost: n_resamples x C(n,2).

        **Independent:** Sample from B's per-pair distributions. For
        each draw, pick one value from each pair's bootstrap
        distribution, average. Treats pairs as independent -- the CI
        will be narrower than correlated and may underestimate variance.
        Cost: 0 estimator calls (numpy only).

        Parameters
        ----------
        pairwise_distributions : dict
            Per-pair bootstrap distributions from B. Required for the
            independent path; ignored for correlated.
        mean_over_pairs : float
            Point estimate of the mean from A, used for CI computation.
        defined_pairs : list
            The pairs defined at the point estimate. A draw in which any of
            them is undefined is NaN, so the interval counts it as failed
            rather than averaging the other pairs.

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
                defined_pairs, n_resamples, rng
            )
        elif bmc.method == "independent":
            mean_distribution = self._mean_bootstrap_independent(
                pairwise_distributions, defined_pairs, n_resamples, rng
            )
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown bootstrap mean method "
                f"'{bmc.method}'. Expected 'correlated' or "
                f"'independent'."
            )

        # Only the correlated mean has a jackknife. It resamples subjects, so
        # leaving one subject out is meaningful. The independent mean
        # resamples from the per-pair distributions and never touches a
        # subject, so there is nothing to leave out.
        jackknife_fn = (
            (lambda: self._mean_jackknife_correlated(defined_pairs))
            if bmc.method == "correlated" else None
        )

        mean_ci = self._compute_ci(
            mean_distribution,
            mean_over_pairs,
            coc.ci_method,
            coc.alpha,
            jackknife_fn=jackknife_fn,
        )

        return mean_distribution, mean_ci

    def _mean_jackknife_correlated(self, defined_pairs: list) -> np.ndarray:
        """
        Leave-one-subject-out means, for the correlated mean's acceleration.

        The counterpart of ``_mean_bootstrap_correlated``: that one resamples
        all subjects with replacement and averages every pair's coefficient
        on each draw, so this one drops a single subject and does the same.
        The unit left out matches the unit resampled, which is what makes the
        acceleration describe the right quantity.

        Parameters
        ----------
        defined_pairs : list
            The pairs defined at the point estimate; see
            ``mean_over_same_subsets``.

        Returns
        -------
        np.ndarray
            Shape (n_subjects,). NaN entries mark leave-one-out data on which
            a defined pair became undefined.
        """
        from ..agreement_functions.rater_combinations import (
            mean_over_same_subsets,
            over_pairs,
        )

        n_subjects = self.estimator_input.n_subjects
        all_rows = np.arange(n_subjects)

        return np.array([
            mean_over_same_subsets(
                over_pairs(
                    self.estimator_input.take_subjects(
                        np.delete(all_rows, i)
                    ),
                    self.estimator,
                ),
                defined_pairs,
            )
            for i in range(n_subjects)
        ], dtype=np.float64)

    def _mean_bootstrap_correlated(
        self,
        defined_pairs: list,
        n_resamples: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """
        Correlated mean bootstrap: global subject resampling.

        For each draw: resample all subjects with replacement, compute
        every pair's coefficient on the same resampled data, average.
        The correlation between pairs sharing a rater is preserved
        because they see the same resampled subjects.

        Parameters
        ----------
        defined_pairs : list
            The pairs defined at the point estimate. A draw in which any of
            them is undefined is NaN, not the mean of the others.
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
            over_pairs,
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
            pair_values = over_pairs(resampled_data, self.estimator)
            mean_draws[b] = mean_over_same_subsets(pair_values, defined_pairs)

        return mean_draws

    def _mean_bootstrap_independent(
        self,
        pairwise_distributions: dict,
        defined_pairs: list,
        n_resamples: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """
        Independent mean bootstrap: sample from B's per-pair distributions.

        For each draw: pick one value from each defined pair's bootstrap
        distribution and average them. A draw that picks a failed (NaN)
        value for any pair is NaN, so it counts as failed rather than
        averaging the other pairs. Treats pairs as independent, so the CI
        will typically be narrower than correlated -- on real data, too
        narrow.

        Parameters
        ----------
        pairwise_distributions : dict
            Per-pair distributions from B.
        defined_pairs : list
            The pairs defined at the point estimate.
        n_resamples : int
            Number of bootstrap draws.
        rng : np.random.Generator
            Seeded random number generator.

        Returns
        -------
        np.ndarray
            Mean coefficient per draw, shape (n_resamples,).
        """
        pairs = defined_pairs
        if not pairs:
            return np.full(n_resamples, np.nan)
        # Stack per-pair distributions: (n_pairs, n_resamples_per_pair)
        stacked = np.array(
            [pairwise_distributions[p] for p in pairs]
        )
        n_pairs = stacked.shape[0]
        n_per_pair = stacked.shape[1]

        # Vectorized: draw one index per pair per bootstrap draw
        indices = rng.integers(
            0, n_per_pair, size=(n_resamples, n_pairs)
        )
        # sampled[b, p] = stacked[p, indices[b, p]]
        sampled = stacked[
            np.arange(n_pairs)[np.newaxis, :], indices
        ]
        # Mean across pairs; any NaN makes the draw NaN
        mean_draws = sampled.mean(axis=1).astype(np.float64)

        return mean_draws
