"""Result for subset permutation analysis."""

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

_ERROR_PREFIX = "[SubsetPermutationResult] Error"


@dataclass(frozen=True)
class SubsetPermutationResult(PermutationResult):
    """
    I am the result of computing a coefficient on every C(n_raters, k) subset.

    I hold two kinds of statistics, and the field names make them unmistakable:

    **Each subset (rater heterogeneity):**
    subset_results holds, for every subset, the coefficient's own result --
    a subset of Fleiss' kappa is a FleissKappaResult, with its point,
    interval, draws, parts and profiles, printed and checked like any single
    run -- or NotComputed saying why the subset has none. They answer "how
    does agreement scale with panel size, and which subgroups agree more or
    less?"

    - subset_results: each subset's result, or NotComputed
    - n_subsets: C(n_raters, k)
    - subset_k: the size of each subset
    - mean_over_subsets: the mean of the subsets' defined points

    **The mean's bootstrap (over draws, sampling uncertainty):**
    Fields prefixed with ``bootstrap_``. Each subset's own interval and
    draws are in its result; these are the mean's. They answer "how stable
    is the mean?"

    Two bootstrap procedures are supported, recorded in
    ``bootstrap_mean_method``:

    - "correlated": resample subjects once per draw, compute all C(n,k)
      subset coefficients on the same resampled subjects, then average.
      Preserves the correlation structure between subsets sharing raters.
    - "independent": bootstrap each subset separately. Per-subset CIs and
      distributions are the same as correlated (resampling subjects within
      a subset), but the mean is computed across independently-drawn subset
      values, breaking the within-draw correlations.

    The per-subset bootstrap is the same regardless of method -- resampling
    subjects within a subset. The method only affects how per-subset draws
    are combined into the mean.

    Fields:

    - bootstrap_mean_method: "correlated" or "independent"
    - bootstrap_ci_of_mean: CI on the mean
    - bootstrap_distribution_of_mean: raw draws of the mean

    Observed spread and bootstrap width answer different questions with
    different remedies. High spread across subsets means agreement depends
    on which raters are in the panel -- maybe some raters need calibration.
    Wide bootstrap CIs means not enough subjects to pin the numbers down --
    maybe collect more data.

    Subsets are tuples of rater indices (0-based, sorted). The raters' names
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
    weights_type : str
        Inherited. The weighting every subset was computed with; the mean is
        weighted the same way.
    weight_matrix : list[list[float]] | None
        Inherited. The matrix used; None when unweighted.
    subset_results : dict[tuple[int, ...], PointEstimateResult | NotComputed]
        Each subset's result, computed on the subjects at least two of its
        raters rated: its balance profile says how many, since a subset
        resting on 756 subjects and one resting on 966 are not measuring
        the same study. Its draws are NotComputed when the analysis ran with
        keep_draws false. A subset with fewer than 2 such subjects, or with
        a rater who rated none of them, holds NotComputed instead, saying
        which. A subset whose coefficient is undefined (all ratings in one
        category) has a result with a NaN point.
    n_subsets : int
        Number of subsets: C(n_raters, k).
    subset_k : int
        The size of each subset. Must be >= 2.
    mean_over_subsets : float
        Mean of the subsets' defined points. Checked against them.
    bootstrap_mean_method : str | None
        Which procedure produced the mean-level bootstrap:
        "correlated" (same resampled subjects across all subsets in each
        draw, preserving subset correlations) or "independent" (each subset
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
    weights_type: str
    weight_matrix: list[list[float]] | None

    # Each subset (rater heterogeneity)
    subset_results: dict   # values are a coefficient's result, or NotComputed
    n_subsets: int
    subset_k: int
    mean_over_subsets: float

    # Bootstrap (sampling uncertainty)
    bootstrap_mean_method: str | None
    bootstrap_ci_of_mean: BootstrapCIResult | None
    bootstrap_distribution_of_mean: np.ndarray | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate subset permutation result attributes."""

        # --- Each subset ----------------------------------------------------- #

        # Validate subset_k first (needed by key validation below)
        if not isinstance(self.subset_k, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: subset_k must be int, "
                f"got {type(self.subset_k).__name__}"
            )
        if self.subset_k < 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: subset_k must be >= 2, "
                f"got {self.subset_k}"
            )

        # Validate n_subsets
        if not isinstance(self.n_subsets, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_subsets must be int, "
                f"got {type(self.n_subsets).__name__}"
            )
        if self.n_subsets < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_subsets must be >= 1, "
                f"got {self.n_subsets}"
            )

        # Validate mean_over_subsets
        if not isinstance(self.mean_over_subsets, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: mean_over_subsets must be float, "
                f"got {type(self.mean_over_subsets).__name__}"
            )

        # Each subset is its coefficient's result or NotComputed, and the
        # mean is the mean of their defined points
        try:
            self._check_unit_results(
                self.subset_results, self.subset_k, self.n_subsets,
                self.mean_over_subsets, "subset",
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
            if self.bootstrap_mean_method not in {
                "correlated",
                "independent",
            }:
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
            if not isinstance(
                self.bootstrap_distribution_of_mean, np.ndarray
            ):
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

        # The mean's draws and interval describe one bootstrap. Each subset's
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
        Return my printed form: the mean, then every subset.

            Subsets of 3 raters, Fleiss' κ, 4 subsets: mean = 0.493 [0.370, 0.600] (...), correlated bootstrap
              199 subjects; 4 raters; each subset on the subjects at least two of its raters rated
                raters         κ  95% CI          subjects
                0, 1, 2    0.514  [0.367, 0.642]       199
                ...

        An independent mean adds a "!" line under the first. Above
        MAX_ROWS_PRINTED subsets, a line points to to_frame() instead of
        the list. See PermutationResult._printed.
        """
        name = result_class_for(self.coefficient_name).NAME
        word = "subset" if self.n_subsets == 1 else "subsets"
        return self._printed(
            title=f"Subsets of {self.subset_k} raters, {name}, {self.n_subsets:,} {word}",
            mean=self.mean_over_subsets,
            mean_ci=self.bootstrap_ci_of_mean,
            mean_method=self.bootstrap_mean_method,
            unit="subset",
            units="subsets",
            pool="the subjects at least two of its raters rated",
        )

    def _unit_results(self) -> dict:
        return self.subset_results

    def summary(self) -> dict:
        """
        Return observed and bootstrap statistics in plain Python types.

        Returns
        -------
        dict
            Parent keys (coefficient_name, ci_method, balance, raters);
            subset_k, n_subsets, mean_over_subsets; subset_results, a list
            with one entry per subset:
            ``{"raters": [names], "positions": [i, j], "result": its own
            summary or {"not_computed": reason}}``; and, when the mean
            was bootstrapped, bootstrap_mean_method,
            mean_ci (with its draw counts), mean_ci_lower, mean_ci_upper,
            bootstrap_draws_mean and bootstrap_se (the mean and standard
            error of the usable draws of the mean, None when too few are
            usable).
        """
        result = super().summary()

        # Observed
        result["subset_k"] = self.subset_k
        result["n_subsets"] = self.n_subsets
        result["mean_over_subsets"] = self.mean_over_subsets
        result["subset_results"] = self._unit_summaries(self.subset_results)

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
