
"""Result for pairwise permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .._bootstrap_draws import check_draws_match_interval, draw_moments
from ..ci_results import BootstrapCIResult
from ..analyzer_results import result_class_for
from .permutation_result import PermutationResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile
    from ...data_profiles.rater_profile import RaterProfile

_ERROR_PREFIX = "[PairwisePermutationResult] Error"


@dataclass(frozen=True)
class PairwisePermutationResult(PermutationResult):
    """
    I am the result of computing a coefficient on every C(n_raters, 2) pair.

    I hold two kinds of statistics, and the field names make them unmistakable:

    **Each pair (rater heterogeneity):**
    pair_results holds, for every pair, the coefficient's own result -- a
    pair of Cohen's kappa is a CohenKappaResult, with its point, interval,
    draws, parts and profiles, printed and checked like any single run --
    or NotComputed saying why the pair has none (raters 0 and 4 share 1
    subject). They answer "how do the raters compare to each other?"

    - pair_results: each pair's result, or NotComputed
    - n_pairs: C(n_raters, 2)
    - mean_over_pairs: the mean of the pairs' defined points

    **The mean's bootstrap (over draws, sampling uncertainty):**
    Fields prefixed with ``bootstrap_``. Each pair's own interval and draws
    are in its result; these are the mean's. They answer "how stable is the
    mean?"

    Two bootstrap procedures are supported, recorded in
    ``bootstrap_mean_method``:

    - "correlated": resample subjects once per draw, compute all C(n,2)
      pairwise coefficients on the same resampled subjects, then average.
      Preserves the correlation structure between pairs sharing a rater.
    - "independent": bootstrap each pair separately. Per-pair CIs and
      distributions are the same as correlated (resampling subjects within
      a pair), but the mean is computed across independently-drawn pair
      values, breaking the within-draw correlations.

    The per-pair bootstrap is the same regardless of method -- resampling
    subjects within a pair. The method only affects how per-pair draws are
    combined into the mean. The analyzer can produce two result objects
    (one per method) when both are wanted.

    Fields:

    - bootstrap_mean_method: "correlated" or "independent"
    - bootstrap_ci_of_mean: CI on the mean
    - bootstrap_distribution_of_mean: raw draws of the mean

    Observed spread and bootstrap width answer different questions with
    different remedies. High spread across pairs means raters disagree
    differently -- maybe train one of them. Wide bootstrap CIs means not
    enough subjects to pin the numbers down -- maybe collect more data.

    Pairs are tuples of rater indices (0-based, sorted). The raters' names
    are in rater_profile.rater_names.

    Attributes
    ----------
    coefficient_name : str
        Inherited. Which coefficient was computed.
    ci_method : str
        Inherited. How CIs were computed.
    balance_profile : BalanceProfile
        Inherited. Structural metadata of the full dataset.
    rater_profile : RaterProfile
        Inherited. The panel's raters; keys are positions into its
        rater_names.
    pair_results : dict[tuple[int, int], PointEstimateResult | NotComputed]
        Each pair's result, computed on the subjects both raters rated: its
        balance profile says how many, since a pair resting on 756 subjects
        and one resting on 966 are not measuring the same study. Its draws
        are NotComputed when the analysis ran with keep_draws false. A pair
        with fewer than 2 shared subjects holds NotComputed instead,
        naming the count. A pair whose coefficient is undefined (all
        ratings in one category) has a result with a NaN point.
    n_pairs : int
        Number of pairs: C(n_raters, 2).
    mean_over_pairs : float
        Mean of the pairs' defined points. Checked against them.
    bootstrap_mean_method : str | None
        Which procedure produced the mean-level bootstrap:
        "correlated" (same resampled subjects across all pairs in each
        draw, preserving pair correlations) or "independent" (each pair
        bootstrapped separately). None if no bootstrap of the mean was
        requested.
    bootstrap_ci_of_mean : BootstrapCIResult | None
        Bootstrap CI on the mean, carrying which method ran and how many
        draws failed. None means no bootstrap of the mean was requested,
        which is a different statement from an interval whose
        ``method_used`` is None -- that one was attempted and every draw
        failed.
    bootstrap_distribution_of_mean : np.ndarray | None
        Raw bootstrap draws of the mean (shape: n_resamples,). None if
        no bootstrap of the mean was requested.
    """

    # Attribute declaration

    # Inherited from parent class PermutationResult
    coefficient_name: str
    ci_method: str
    balance_profile: BalanceProfile
    rater_profile: RaterProfile

    # Each pair (rater heterogeneity)
    pair_results: dict   # values are a coefficient's result, or NotComputed
    n_pairs: int
    mean_over_pairs: float

    # Bootstrap (sampling uncertainty)
    bootstrap_mean_method: str | None
    bootstrap_ci_of_mean: BootstrapCIResult | None
    bootstrap_distribution_of_mean: np.ndarray | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate pairwise permutation result attributes."""

        # --- Each pair ------------------------------------------------------- #

        # Validate n_pairs
        if not isinstance(self.n_pairs, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_pairs must be int, "
                f"got {type(self.n_pairs).__name__}"
            )
        if self.n_pairs < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_pairs must be >= 1, "
                f"got {self.n_pairs}"
            )

        # Validate mean_over_pairs
        if not isinstance(self.mean_over_pairs, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: mean_over_pairs must be float, "
                f"got {type(self.mean_over_pairs).__name__}"
            )

        # Each pair is its coefficient's result or NotComputed, and the mean
        # is the mean of their defined points
        try:
            self._check_unit_results(
                self.pair_results, 2, self.n_pairs, self.mean_over_pairs, "pair",
            )
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # --- Bootstrap ------------------------------------------------------- #

        # Validate bootstrap_mean_method
        if self.bootstrap_mean_method is not None:
            if not isinstance(self.bootstrap_mean_method, str):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_mean_method must be str "
                    f"or None, got "
                    f"{type(self.bootstrap_mean_method).__name__}"
                )
            if self.bootstrap_mean_method not in {"correlated", "independent"}:
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_mean_method must be "
                    f"'correlated' or 'independent', "
                    f"got '{self.bootstrap_mean_method}'"
                )

        # Cross-validate: mean bootstrap fields require bootstrap_mean_method
        has_mean_bootstrap = (
            self.bootstrap_ci_of_mean is not None
            or self.bootstrap_distribution_of_mean is not None
        )
        if has_mean_bootstrap and self.bootstrap_mean_method is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: bootstrap_mean_method must be set "
                f"when mean bootstrap fields are populated"
            )
        if not has_mean_bootstrap and self.bootstrap_mean_method is not None:
            raise ValueError(
                f"{_ERROR_PREFIX}: bootstrap_mean_method is set but "
                f"no mean bootstrap fields are populated"
            )

        # Validate bootstrap_ci_of_mean. An interval validates its own
        # bounds, level, counts, and method at construction, so the only
        # thing left to check here is that it is one.
        if self.bootstrap_ci_of_mean is not None:
            if not isinstance(self.bootstrap_ci_of_mean, BootstrapCIResult):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_ci_of_mean must be "
                    f"BootstrapCIResult or None, got "
                    f"{type(self.bootstrap_ci_of_mean).__name__}"
                )

        # Validate bootstrap_distribution_of_mean
        if self.bootstrap_distribution_of_mean is not None:
            if not isinstance(self.bootstrap_distribution_of_mean, np.ndarray):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution_of_mean "
                    f"must be np.ndarray or None, got "
                    f"{type(self.bootstrap_distribution_of_mean).__name__}"
                )
            if self.bootstrap_distribution_of_mean.ndim != 1:
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution_of_mean "
                    f"must be 1-dimensional, got shape "
                    f"{self.bootstrap_distribution_of_mean.shape}"
                )
            if len(self.bootstrap_distribution_of_mean) == 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution_of_mean "
                    f"must have at least 1 value"
                )

        # The mean's draws and interval describe one bootstrap. Each pair's
        # result checked its own when it was built.
        try:
            check_draws_match_interval(
                self.bootstrap_distribution_of_mean,
                self.bootstrap_ci_of_mean,
                draws_name="bootstrap_distribution_of_mean",
                interval_name="bootstrap_ci_of_mean",
            )
        except ValueError as exc:
            raise ValueError(f"{_ERROR_PREFIX}: {exc}") from exc

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Return my printed form: the mean, then every pair.

            Pairwise Cohen's κ, 6 pairs: mean = 0.496 [0.377, 0.605] (95% CI, percentile, n=500), correlated bootstrap
              199 subjects; 4 raters; each pair on the subjects both raters rated
                raters      κ  95% CI          subjects
                0, 1    0.539  [0.383, 0.682]       199
                ...

        An independent mean adds a "!" line under the first. Above
        MAX_ROWS_PRINTED pairs, a line points to to_frame() instead of the
        list. See PermutationResult._printed.
        """
        name = result_class_for(self.coefficient_name).NAME
        return self._printed(
            title=f"Pairwise {name}, {self.n_pairs:,} pairs",
            mean=self.mean_over_pairs,
            mean_ci=self.bootstrap_ci_of_mean,
            mean_method=self.bootstrap_mean_method,
            unit="pair",
            units="pairs",
            pool="the subjects both raters rated",
        )

    def _unit_results(self) -> dict:
        return self.pair_results

    def summary(self) -> dict:
        """
        Return observed and bootstrap statistics in plain Python types.

        Returns
        -------
        dict
            Parent keys (coefficient_name, ci_method, balance, raters);
            n_pairs, mean_over_pairs; pair_results, a list with one entry per pair:
            ``{"raters": [names], "positions": [i, j], "result": its own
            summary or {"not_computed": reason}}``; and, when the mean was
            bootstrapped, bootstrap_mean_method,
            mean_ci (with its draw counts), mean_ci_lower, mean_ci_upper,
            bootstrap_draws_mean and bootstrap_se (the mean and standard
            error of the usable draws of the mean, None when too few are
            usable).
        """
        result = super().summary()

        # Observed
        result["n_pairs"] = self.n_pairs
        result["mean_over_pairs"] = self.mean_over_pairs
        result["pair_results"] = self._unit_summaries(self.pair_results)

        # Bootstrap
        if self.bootstrap_mean_method is not None:
            result["bootstrap_mean_method"] = self.bootstrap_mean_method

        if self.bootstrap_ci_of_mean is not None:
            result["mean_ci"] = self.bootstrap_ci_of_mean.summary()
            result["mean_ci_lower"] = self.bootstrap_ci_of_mean.lower
            result["mean_ci_upper"] = self.bootstrap_ci_of_mean.upper

        if self.bootstrap_distribution_of_mean is not None:
            mean, se = draw_moments(self.bootstrap_distribution_of_mean)
            result["bootstrap_draws_mean"] = mean
            result["bootstrap_se"] = se

        return result
