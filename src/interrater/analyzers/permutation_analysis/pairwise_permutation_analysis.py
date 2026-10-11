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
from ...intermediates.analysis_results.pairwise_permutation_result import (
    MAX_PAIRS_FOR_DISTRIBUTIONS,
)
from .permutation_analysis import PermutationAnalysis

if TYPE_CHECKING:
    from ...data_objects.ratings import Ratings
    from ...intermediates.analysis_results.pairwise_permutation_result import (
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

    A runs always. B runs if ``config_of_coefficient.n_resamples > 0``.
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
        from ...intermediates.analysis_results.pairwise_permutation_result import (
            PairwisePermutationResult,
        )

        coc = self.config.config_of_coefficient
        bmc = self.config.bootstrap_mean_config

        # Computation A -- always
        pairwise_values, subjects_per_pair, mean_over_pairs = (
            self._compute_observed()
        )

        # Computation B -- if per-pair bootstrap requested
        pairwise_distributions = None
        pairwise_cis = None
        if coc.n_resamples > 0:
            pairwise_distributions, pairwise_cis = (
                self._compute_per_pair_bootstrap(pairwise_values)
            )

        # Computation C -- if mean bootstrap requested
        mean_distribution = None
        mean_ci = None
        mean_method = None
        if bmc is not None:
            mean_method = bmc.method
            mean_distribution, mean_ci = self._compute_mean_bootstrap(
                pairwise_distributions, mean_over_pairs
            )

        # Decide whether to store per-pair distributions
        n_pairs = len(pairwise_values)
        store_distributions = (
            pairwise_distributions is not None
            and n_pairs <= MAX_PAIRS_FOR_DISTRIBUTIONS
        )

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
                distributions[pair] = np.full(n_resamples, np.nan)
                cis[pair] = (float("nan"), float("nan"))
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
            cis[pair] = self._compute_ci(draws, point, ci_method, alpha)

        return distributions, cis

    def _compute_mean_bootstrap(
        self,
        pairwise_distributions: dict | None,
        mean_over_pairs: float,
    ) -> tuple[np.ndarray, tuple[float, float]]:
        """
        Computation C: mean bootstrap.

        Two paths based on ``bootstrap_mean_config.method``:

        **Correlated:** Resample subjects globally. For each draw,
        recompute all C(n,2) pair coefficients on the same resampled
        subjects, then average. Preserves the correlation between pairs
        sharing a rater. Cost: n_bootstraps x C(n,2).

        **Independent:** Sample from B's per-pair distributions. For
        each draw, pick one value from each pair's bootstrap
        distribution, average. Treats pairs as independent -- the CI
        will be narrower than correlated and may underestimate variance.
        Cost: 0 estimator calls (numpy only). Requires B to have run.

        Parameters
        ----------
        pairwise_distributions : dict | None
            Per-pair bootstrap distributions from B. Required for the
            independent path; ignored for correlated.
        mean_over_pairs : float
            Point estimate of the mean from A, used for CI computation.

        Returns
        -------
        tuple[np.ndarray, tuple[float, float]]
            (mean_distribution, mean_ci)
            distribution: shape (n_bootstraps,)
            ci: (lower, upper)
        """
        bmc = self.config.bootstrap_mean_config
        n_bootstraps = bmc.n_bootstraps
        rng = np.random.default_rng(bmc.random_state)

        if bmc.method == "correlated":
            mean_distribution = self._mean_bootstrap_correlated(
                n_bootstraps, rng
            )
        elif bmc.method == "independent":
            mean_distribution = self._mean_bootstrap_independent(
                pairwise_distributions, n_bootstraps, rng
            )
        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown bootstrap mean method "
                f"'{bmc.method}'. Expected 'correlated' or "
                f"'independent'."
            )

        mean_ci = self._compute_ci(
            mean_distribution,
            mean_over_pairs,
            bmc.ci_method,
            bmc.alpha,
        )

        return mean_distribution, mean_ci

    def _mean_bootstrap_correlated(
        self,
        n_bootstraps: int,
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
        n_bootstraps : int
            Number of bootstrap draws.
        rng : np.random.Generator
            Seeded random number generator.

        Returns
        -------
        np.ndarray
            Mean coefficient per draw, shape (n_bootstraps,).
        """
        from ..agreement_functions.rater_combinations import (
            mean_over_subsets,
            over_pairs,
        )

        n_subjects = self.estimator_input.n_subjects
        mean_draws = np.empty(n_bootstraps, dtype=np.float64)

        for b in range(n_bootstraps):
            resampled_rows = rng.choice(
                n_subjects, size=n_subjects, replace=True
            )
            resampled_data = self.estimator_input.take_subjects(
                resampled_rows
            )
            pair_values = over_pairs(resampled_data, self.estimator)
            mean_draws[b] = mean_over_subsets(pair_values)

        return mean_draws

    def _mean_bootstrap_independent(
        self,
        pairwise_distributions: dict | None,
        n_bootstraps: int,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """
        Independent mean bootstrap: sample from B's per-pair distributions.

        For each draw: pick one value from each pair's bootstrap
        distribution, average across pairs (NaN-skipping). Treats pairs
        as independent, so the CI will typically be narrower than
        correlated -- on real data, too narrow.

        Parameters
        ----------
        pairwise_distributions : dict | None
            Per-pair distributions from B. Must not be None.
        n_bootstraps : int
            Number of bootstrap draws.
        rng : np.random.Generator
            Seeded random number generator.

        Returns
        -------
        np.ndarray
            Mean coefficient per draw, shape (n_bootstraps,).

        Raises
        ------
        ValueError
            If pairwise_distributions is None (B did not run).
        """
        if pairwise_distributions is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: 'independent' mean bootstrap "
                f"requires per-pair bootstrap (n_resamples >= 1 in "
                f"config_of_coefficient). The config should have "
                f"caught this -- if you see this error, there is a "
                f"bug in PairwisePermutationConfig validation."
            )

        pairs = sorted(pairwise_distributions.keys())
        # Stack per-pair distributions: (n_pairs, n_resamples_per_pair)
        stacked = np.array(
            [pairwise_distributions[p] for p in pairs]
        )
        n_pairs = stacked.shape[0]
        n_per_pair = stacked.shape[1]

        # Vectorized: draw one index per pair per bootstrap draw
        indices = rng.integers(
            0, n_per_pair, size=(n_bootstraps, n_pairs)
        )
        # sampled[b, p] = stacked[p, indices[b, p]]
        sampled = stacked[
            np.arange(n_pairs)[np.newaxis, :], indices
        ]
        # NaN-skipping mean across pairs for each draw
        mean_draws = np.nanmean(sampled, axis=1).astype(np.float64)

        return mean_draws
