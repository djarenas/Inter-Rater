"""Result for subset permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .._bootstrap_draws import check_draws_match_interval, draw_moments
from ..ci_results import BootstrapCIResult
from .permutation_result import PermutationResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[SubsetPermutationResult] Error"

# Default threshold for storing per-subset bootstrap distributions.
# Above this many subsets, raw distributions are not stored to save memory.
# Per-subset CIs are always stored regardless of this threshold.
# At 1000 resamples per subset, 120 subsets costs ~1 MB -- lightweight.
# Covers the common cases: C(10,3)=120, C(9,4)=126, C(7,5)=21.
MAX_SUBSETS_FOR_DISTRIBUTIONS: int = 120


@dataclass(frozen=True)
class SubsetPermutationResult(PermutationResult):
    """
    I am the result of computing a coefficient on every C(n_raters, k) subset.

    I hold two kinds of statistics, and the field names make them unmistakable:

    **Observed (over subsets, rater heterogeneity):**
    Fields with no prefix. Computed once on the actual data. They answer
    "how does agreement scale with panel size, and which subgroups agree
    more or less?"

    - subset_values: the coefficient for each subset
    - subjects_per_subset: how many subjects each subset shared
    - n_subsets: C(n_raters, k)
    - subset_k: the size of each subset
    - mean_over_subsets: the mean across subsets (NaN-skipping)

    **Bootstrap (over draws, sampling uncertainty):**
    Fields prefixed with ``bootstrap_``. Computed across resamples. They
    answer "how stable are these numbers?"

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
    - bootstrap_cis_per_subset: CI per subset (always computed and stored)
    - bootstrap_distributions_per_subset: raw draws per subset (only when
      n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS, to bound memory)

    Observed spread and bootstrap width answer different questions with
    different remedies. High spread across subsets means agreement depends
    on which raters are in the panel -- maybe some raters need calibration.
    Wide bootstrap CIs means not enough subjects to pin the numbers down --
    maybe collect more data.

    I carry no rater names. Subsets are tuples of rater indices (0-based,
    sorted). Whoever renders the result has the semantics and can supply
    names from the Ratings object.

    Attributes
    ----------
    coefficient_name : str
        Inherited. Which coefficient was computed.
    ci_method : str
        Inherited. How CIs were computed.
    balance_profile : BalanceProfile
        Inherited. Structural metadata of the full dataset.
    subset_values : dict[tuple[int, ...], float]
        Coefficient per subset. NaN for subsets that could not be computed
        (fewer than 2 shared subjects, all ratings in one category).
    subjects_per_subset : dict[tuple[int, ...], int]
        Surviving subjects per subset -- the denominator each coefficient
        was actually computed on. A subset whose coefficient rests on 756
        subjects and one resting on 966 are not measuring the same study.
    n_subsets : int
        Number of subsets: C(n_raters, k).
    subset_k : int
        The size of each subset. Must be >= 2.
    mean_over_subsets : float
        Mean across subsets, skipping NaN subsets.
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
    bootstrap_cis_per_subset : dict[tuple[int, ...], BootstrapCIResult]
        Bootstrap CI per subset, one interval object each. Subsets do not
        fail alike, and here that is routine rather than exceptional: a
        subset with fewer than two measurable subjects has nothing to
        resample and gets an interval whose ``method_used`` is None, while
        its siblings on the same data return ordinary intervals. Each one
        carries its own counts and its own method, so that difference
        survives into the result instead of arriving as an indistinguishable
        pair of NaNs. Always present: the analyzer config requires
        n_resamples >= 1, so the per-subset bootstrap always runs.
    bootstrap_distributions_per_subset : dict[tuple[int, ...], np.ndarray] | None
        Raw bootstrap draws per subset. Stored exactly when
        n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS, None above it, to bound
        memory. Use get_bootstrap_distributions_per_subset() for access
        with a helpful error when they were not stored.
    """

    # Attribute declaration

    # Inherited from parent class PermutationResult
    coefficient_name: str
    ci_method: str
    balance_profile: BalanceProfile

    # Observed, over subsets (rater heterogeneity)
    subset_values: dict
    subjects_per_subset: dict
    n_subsets: int
    subset_k: int
    mean_over_subsets: float

    # Bootstrap (sampling uncertainty)
    bootstrap_mean_method: str | None
    bootstrap_ci_of_mean: BootstrapCIResult | None
    bootstrap_distribution_of_mean: np.ndarray | None
    bootstrap_cis_per_subset: dict   # values are BootstrapCIResult
    bootstrap_distributions_per_subset: dict | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate subset permutation result attributes."""

        # --- Observed -------------------------------------------------------- #

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

        # Validate subset_values
        if not isinstance(self.subset_values, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: subset_values must be dict, "
                f"got {type(self.subset_values).__name__}"
            )

        # All keys must be k-tuples of ints
        for key in self.subset_values:
            if not isinstance(key, tuple) or len(key) != self.subset_k:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subset_values keys must be "
                    f"{self.subset_k}-tuples, got {key!r}"
                )

        # Validate n_subsets matches
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
        if len(self.subset_values) != self.n_subsets:
            raise ValueError(
                f"{_ERROR_PREFIX}: subset_values has "
                f"{len(self.subset_values)} entries but "
                f"n_subsets is {self.n_subsets}"
            )

        # Validate subjects_per_subset has same keys
        if not isinstance(self.subjects_per_subset, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: subjects_per_subset must be dict, "
                f"got {type(self.subjects_per_subset).__name__}"
            )
        if set(self.subjects_per_subset.keys()) != set(
            self.subset_values.keys()
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: subjects_per_subset must have the same "
                f"keys as subset_values"
            )

        # Validate mean_over_subsets
        if not isinstance(self.mean_over_subsets, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: mean_over_subsets must be float, "
                f"got {type(self.mean_over_subsets).__name__}"
            )

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

        # Validate bootstrap_cis_per_subset. Always present: the per-subset
        # bootstrap cannot be turned off (n_resamples >= 1).
        if not isinstance(self.bootstrap_cis_per_subset, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: bootstrap_cis_per_subset must be "
                f"a dict (the per-subset bootstrap always runs), got "
                f"{type(self.bootstrap_cis_per_subset).__name__}"
            )
        if set(self.bootstrap_cis_per_subset.keys()) != set(
            self.subset_values.keys()
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: bootstrap_cis_per_subset must have "
                f"the same keys as subset_values"
            )
        for key, ci in self.bootstrap_cis_per_subset.items():
            if not isinstance(ci, BootstrapCIResult):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap CI for subset {key} "
                    f"must be BootstrapCIResult, got {type(ci).__name__}"
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

        # Validate bootstrap_distributions_per_subset
        if self.bootstrap_distributions_per_subset is not None:
            if not isinstance(
                self.bootstrap_distributions_per_subset, dict
            ):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_distributions_per_subset "
                    f"must be dict, got "
                    f"{type(self.bootstrap_distributions_per_subset).__name__}"
                )
            if set(
                self.bootstrap_distributions_per_subset.keys()
            ) != set(self.subset_values.keys()):
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_distributions_per_subset "
                    f"must have the same keys as subset_values"
                )
            for key, dist in (
                self.bootstrap_distributions_per_subset.items()
            ):
                if not isinstance(dist, np.ndarray) or dist.ndim != 1:
                    raise ValueError(
                        f"{_ERROR_PREFIX}: bootstrap distribution for "
                        f"subset {key} must be a 1-D np.ndarray"
                    )

        # The mean's draws and interval describe one bootstrap, and so do
        # each subset's, when its draws are stored.
        try:
            check_draws_match_interval(
                self.bootstrap_distribution_of_mean,
                self.bootstrap_ci_of_mean,
                draws_name="bootstrap_distribution_of_mean",
                interval_name="bootstrap_ci_of_mean",
            )
            if self.bootstrap_distributions_per_subset is not None:
                for key, draws in self.bootstrap_distributions_per_subset.items():
                    check_draws_match_interval(
                        draws,
                        self.bootstrap_cis_per_subset[key],
                        draws_name=f"subset {key}",
                        interval_name="its interval",
                    )
        except ValueError as exc:
            raise ValueError(f"{_ERROR_PREFIX}: {exc}") from exc

        # Distributions are stored exactly when n_subsets is at or below
        # the limit, so a missing one has one explanation.
        should_store = self.n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS
        stored = self.bootstrap_distributions_per_subset is not None
        if stored != should_store:
            raise ValueError(
                f"{_ERROR_PREFIX}: bootstrap_distributions_per_subset must "
                f"be stored when n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS "
                f"({MAX_SUBSETS_FOR_DISTRIBUTIONS}) and None above it; got "
                f"n_subsets={self.n_subsets} with distributions "
                f"{'stored' if stored else 'missing'}"
            )

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Return a compact string representation.

        Format:
            Subset (k=3, fleiss_kappa): mean = 0.650
            [0.580, 0.720] (95% CI, bca, n=4000), correlated
            n_subsets=10; n_subjects=1617; n_raters=5; ci_method=bca

        The middle line is the interval printing itself, so a fallback to
        percentile, a count of failed draws, or a method that never ran at
        all appears there without this class having to know about any of
        them. It describes the mean only -- a subset that had nothing to
        resample says so in its own entry of bootstrap_cis_per_subset.
        """
        ci_line = (
            f"\n{self.bootstrap_ci_of_mean}, {self.bootstrap_mean_method}"
            if self.bootstrap_ci_of_mean is not None
            else ""
        )

        return (
            f"Subset (k={self.subset_k}, {self.coefficient_name}): "
            f"mean = {self.mean_over_subsets:.3f}"
            f"{ci_line}\n"
            f"n_subsets={self.n_subsets}; "
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"n_raters={self.balance_profile.n_raters}; "
            f"ci_method={self.ci_method}"
        )

    def get_bootstrap_distributions_per_subset(self) -> dict:
        """
        Return per-subset bootstrap distributions, or raise if not stored.

        The per-subset bootstrap always runs, but its raw draws are only
        stored when n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS, to bound
        memory. This method says so rather than returning None.

        Returns
        -------
        dict[tuple[int, ...], np.ndarray]
            Raw bootstrap draws per subset.

        Raises
        ------
        RuntimeError
            If n_subsets is above MAX_SUBSETS_FOR_DISTRIBUTIONS.
        """
        if self.bootstrap_distributions_per_subset is not None:
            return self.bootstrap_distributions_per_subset
        raise RuntimeError(
            f"{_ERROR_PREFIX}: per-subset bootstrap distributions were not "
            f"stored because n_subsets ({self.n_subsets}) > "
            f"MAX_SUBSETS_FOR_DISTRIBUTIONS ({MAX_SUBSETS_FOR_DISTRIBUTIONS}). "
            f"Per-subset CIs are in bootstrap_cis_per_subset."
        )

    def summary(self) -> dict:
        """
        Return observed and bootstrap statistics in plain Python types.

        Returns
        -------
        dict
            Parent keys (coefficient_name, ci_method, balance); the
            observed statistics (subset_values, subjects_per_subset,
            n_subsets, subset_k, mean_over_subsets); bootstrap_cis_per_subset;
            and, when the mean was bootstrapped, bootstrap_mean_method,
            mean_ci (with its draw counts), mean_ci_lower, mean_ci_upper,
            bootstrap_draws_mean and bootstrap_se (the mean and standard
            error of the usable draws of the mean, None when too few are
            usable).
            Tuple keys are converted to strings for JSON compatibility.
        """
        result = super().summary()

        # Observed
        result["subset_k"] = self.subset_k
        result["n_subsets"] = self.n_subsets
        result["mean_over_subsets"] = self.mean_over_subsets
        result["subset_values"] = {
            str(k): v for k, v in self.subset_values.items()
        }
        result["subjects_per_subset"] = {
            str(k): v for k, v in self.subjects_per_subset.items()
        }

        # Bootstrap
        if self.bootstrap_mean_method is not None:
            result["bootstrap_mean_method"] = self.bootstrap_mean_method

        if self.bootstrap_ci_of_mean is not None:
            result["mean_ci"] = self.bootstrap_ci_of_mean.summary()
            result["mean_ci_lower"] = self.bootstrap_ci_of_mean.lower
            result["mean_ci_upper"] = self.bootstrap_ci_of_mean.upper

        result["bootstrap_cis_per_subset"] = {
            str(k): v.summary()
            for k, v in self.bootstrap_cis_per_subset.items()
        }

        if self.bootstrap_distribution_of_mean is not None:
            mean, se = draw_moments(self.bootstrap_distribution_of_mean)
            result["bootstrap_draws_mean"] = mean
            result["bootstrap_se"] = se

        # Per-subset distributions are NOT included in summary (too large).
        # Access them directly via the field.
        result["subset_distributions_stored"] = (
            self.bootstrap_distributions_per_subset is not None
        )

        return result
