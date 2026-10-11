"""Result for subset permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .permutation_result import PermutationResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile

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
    - bootstrap_cis_per_subset: CI per subset (always stored if bootstrapped)
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
        bootstrapped separately). None if bootstrap was not performed.
    bootstrap_ci_of_mean : tuple[float, float] | None
        Bootstrap CI on the mean. None if bootstrap was not performed.
    bootstrap_distribution_of_mean : np.ndarray | None
        Raw bootstrap draws of the mean (shape: n_bootstraps,). None if
        bootstrap was not performed.
    bootstrap_cis_per_subset : dict[tuple[int, ...], tuple[float, float]] | None
        Bootstrap CI per subset. Always stored if bootstrap was performed,
        regardless of n_subsets. None if bootstrap was not performed.
    bootstrap_distributions_per_subset : dict[tuple[int, ...], np.ndarray] | None
        Raw bootstrap draws per subset. Only stored when
        n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS. None if not stored or
        if bootstrap was not performed. Use
        get_bootstrap_distributions_per_subset() for safe access with a
        helpful error when distributions were not stored.
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
    bootstrap_ci_of_mean: tuple[float, float] | None
    bootstrap_distribution_of_mean: np.ndarray | None
    bootstrap_cis_per_subset: dict | None
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

        # Validate bootstrap_ci_of_mean
        if self.bootstrap_ci_of_mean is not None:
            if (
                not isinstance(self.bootstrap_ci_of_mean, tuple)
                or len(self.bootstrap_ci_of_mean) != 2
            ):
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_ci_of_mean must be "
                    f"(float, float) tuple"
                )
            lower, upper = self.bootstrap_ci_of_mean
            if not isinstance(lower, (int, float)) or not isinstance(
                upper, (int, float)
            ):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_ci_of_mean bounds must "
                    f"be float, got ({type(lower).__name__}, "
                    f"{type(upper).__name__})"
                )
            if lower > upper:
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_ci_of_mean lower bound "
                    f"{lower} > upper bound {upper}"
                )

        # Validate bootstrap_cis_per_subset
        if self.bootstrap_cis_per_subset is not None:
            if not isinstance(self.bootstrap_cis_per_subset, dict):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_cis_per_subset must be "
                    f"dict, got "
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
                if not isinstance(ci, tuple) or len(ci) != 2:
                    raise ValueError(
                        f"{_ERROR_PREFIX}: bootstrap CI for subset {key} "
                        f"must be (float, float) tuple"
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
            (95% CI: [0.580, 0.720], n_bootstrap=4000, correlated)
            n_subsets=10; n_subjects=1617; n_raters=5; ci_method=bca
        """
        ci_str = ""
        if self.bootstrap_ci_of_mean is not None:
            lower, upper = self.bootstrap_ci_of_mean
            ci_str = f" (95% CI: [{lower:.3f}, {upper:.3f}]"

        bootstrap_str = ""
        if self.bootstrap_distribution_of_mean is not None:
            n_boot = len(self.bootstrap_distribution_of_mean)
            bootstrap_str = f", n_bootstrap={n_boot}"

        method_str = ""
        if self.bootstrap_mean_method is not None:
            method_str = f", {self.bootstrap_mean_method}"

        closing = ")" if ci_str else ""

        return (
            f"Subset (k={self.subset_k}, {self.coefficient_name}): "
            f"mean = {self.mean_over_subsets:.3f}"
            f"{ci_str}{bootstrap_str}{method_str}{closing}\n"
            f"n_subsets={self.n_subsets}; "
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"n_raters={self.balance_profile.n_raters}; "
            f"ci_method={self.ci_method}"
        )

    def get_bootstrap_distributions_per_subset(self) -> dict:
        """
        Return per-subset bootstrap distributions, or raise if not stored.

        Per-subset distributions are only stored when
        n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS to bound memory. This
        method provides a helpful error when distributions were not
        stored, rather than silently returning None.

        Returns
        -------
        dict[tuple[int, ...], np.ndarray]
            Raw bootstrap draws per subset.

        Raises
        ------
        RuntimeError
            If distributions were not stored (n_subsets exceeded threshold)
            or bootstrap was not performed.
        """
        if self.bootstrap_distributions_per_subset is not None:
            return self.bootstrap_distributions_per_subset
        if self.bootstrap_mean_method is None:
            raise RuntimeError(
                f"{_ERROR_PREFIX}: per-subset bootstrap distributions are "
                f"not available because bootstrap was not performed"
            )
        raise RuntimeError(
            f"{_ERROR_PREFIX}: per-subset bootstrap distributions were not "
            f"stored because n_subsets ({self.n_subsets}) > "
            f"MAX_SUBSETS_FOR_DISTRIBUTIONS ({MAX_SUBSETS_FOR_DISTRIBUTIONS}). "
            f"Per-subset CIs are still available via bootstrap_cis_per_subset. "
            f"To store raw distributions, re-run with a higher threshold "
            f"or fewer raters."
        )

    def summary(self) -> dict:
        """
        Return observed and bootstrap statistics in plain Python types.

        Returns
        -------
        dict
            Parent keys (coefficient_name, ci_method, balance) plus
            observed statistics (subset_values, subjects_per_subset,
            n_subsets, subset_k, mean_over_subsets) and bootstrap
            statistics (bootstrap_mean_method, mean_ci_lower,
            mean_ci_upper, bootstrap_n, etc.).
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
            result["mean_ci_lower"] = self.bootstrap_ci_of_mean[0]
            result["mean_ci_upper"] = self.bootstrap_ci_of_mean[1]

        if self.bootstrap_cis_per_subset is not None:
            result["bootstrap_cis_per_subset"] = {
                str(k): {"lower": v[0], "upper": v[1]}
                for k, v in self.bootstrap_cis_per_subset.items()
            }

        if self.bootstrap_distribution_of_mean is not None:
            result["bootstrap_n"] = len(
                self.bootstrap_distribution_of_mean
            )
            result["bootstrap_mean"] = float(
                np.mean(self.bootstrap_distribution_of_mean)
            )
            result["bootstrap_std"] = float(
                np.std(self.bootstrap_distribution_of_mean)
            )

        # Per-subset distributions are NOT included in summary (too large).
        # Access them directly via the field.
        result["subset_distributions_stored"] = (
            self.bootstrap_distributions_per_subset is not None
        )

        return result
