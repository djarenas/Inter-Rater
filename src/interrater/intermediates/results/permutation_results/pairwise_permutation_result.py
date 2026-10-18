
"""Result for pairwise permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..ci_results import BootstrapCIResult
from .permutation_result import PermutationResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[PairwisePermutationResult] Error"

# Default threshold for storing per-pair bootstrap distributions.
# Above this many pairs, raw distributions are not stored to save memory.
# Per-pair CIs are always stored regardless of this threshold.
MAX_PAIRS_FOR_DISTRIBUTIONS: int = 20


@dataclass(frozen=True)
class PairwisePermutationResult(PermutationResult):
    """
    I am the result of computing a coefficient on every C(n_raters, 2) pair.

    I hold two kinds of statistics, and the field names make them unmistakable:

    **Observed (over pairs, rater heterogeneity):**
    Fields with no prefix. Computed once on the actual data. They answer
    "how do the raters compare to each other?"

    - pairwise_values: the coefficient for each pair
    - subjects_per_pair: how many subjects each pair shared
    - n_pairs: C(n_raters, 2)
    - mean_over_pairs: the mean across pairs (NaN-skipping)

    **Bootstrap (over draws, sampling uncertainty):**
    Fields prefixed with ``bootstrap_``. Computed across resamples. They
    answer "how stable are these numbers?"

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
    - bootstrap_cis_per_pair: CI per pair (always stored if bootstrapped)
    - bootstrap_distributions_per_pair: raw draws per pair (only when
      n_pairs <= MAX_PAIRS_FOR_DISTRIBUTIONS, to bound memory)

    Observed spread and bootstrap width answer different questions with
    different remedies. High spread across pairs means raters disagree
    differently -- maybe train one of them. Wide bootstrap CIs means not
    enough subjects to pin the numbers down -- maybe collect more data.

    I carry no rater names. Pairs are tuples of rater indices (0-based,
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
    pairwise_values : dict[tuple[int, int], float]
        Coefficient per pair. NaN for pairs that could not be computed
        (fewer than 2 shared subjects, all ratings in one category).
    subjects_per_pair : dict[tuple[int, int], int]
        Surviving subjects per pair -- the denominator each coefficient
        was actually computed on. A pair whose coefficient rests on 756
        subjects and one resting on 966 are not measuring the same study.
    n_pairs : int
        Number of pairs: C(n_raters, 2).
    mean_over_pairs : float
        Mean across pairs, skipping NaN pairs.
    bootstrap_mean_method : str | None
        Which procedure produced the mean-level bootstrap:
        "correlated" (same resampled subjects across all pairs in each
        draw, preserving pair correlations) or "independent" (each pair
        bootstrapped separately). None if bootstrap was not performed.
    bootstrap_ci_of_mean : BootstrapCIResult | None
        Bootstrap CI on the mean, carrying which method ran and how many
        draws failed. None means no bootstrap was performed, which is a
        different statement from an interval whose ``method_used`` is None
        -- that one was attempted and every draw failed.
    bootstrap_distribution_of_mean : np.ndarray | None
        Raw bootstrap draws of the mean (shape: n_resamples,). None if
        bootstrap was not performed.
    bootstrap_cis_per_pair : dict[tuple[int, int], BootstrapCIResult] | None
        Bootstrap CI per pair, one interval object each. Pairs do not fail
        alike -- a pair with few shared subjects can lose draws that a
        well-covered pair never loses -- so each interval carries its own
        counts and its own method rather than a shared summary. Always
        stored if bootstrap was performed, regardless of n_pairs. None if
        bootstrap was not performed.
    bootstrap_distributions_per_pair : dict[tuple[int, int], np.ndarray] | None
        Raw bootstrap draws per pair. Only stored when
        n_pairs <= MAX_PAIRS_FOR_DISTRIBUTIONS. None if not stored or
        if bootstrap was not performed. Use
        get_bootstrap_distributions_per_pair() for safe access with a
        helpful error when distributions were not stored.
    """

    # Attribute declaration

    # Inherited from parent class PermutationResult
    coefficient_name: str
    ci_method: str
    balance_profile: BalanceProfile

    # Observed, over pairs (rater heterogeneity)
    pairwise_values: dict
    subjects_per_pair: dict
    n_pairs: int
    mean_over_pairs: float

    # Bootstrap (sampling uncertainty)
    bootstrap_mean_method: str | None
    bootstrap_ci_of_mean: BootstrapCIResult | None
    bootstrap_distribution_of_mean: np.ndarray | None
    bootstrap_cis_per_pair: dict | None   # values are BootstrapCIResult
    bootstrap_distributions_per_pair: dict | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate pairwise permutation result attributes."""

        # --- Observed -------------------------------------------------------- #

        # Validate pairwise_values
        if not isinstance(self.pairwise_values, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: pairwise_values must be dict, "
                f"got {type(self.pairwise_values).__name__}"
            )

        # All keys must be 2-tuples of ints
        for key in self.pairwise_values:
            if not isinstance(key, tuple) or len(key) != 2:
                raise ValueError(
                    f"{_ERROR_PREFIX}: pairwise_values keys must be "
                    f"2-tuples, got {key!r}"
                )

        # Validate n_pairs matches
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
        if len(self.pairwise_values) != self.n_pairs:
            raise ValueError(
                f"{_ERROR_PREFIX}: pairwise_values has "
                f"{len(self.pairwise_values)} entries but "
                f"n_pairs is {self.n_pairs}"
            )

        # Validate subjects_per_pair has same keys
        if not isinstance(self.subjects_per_pair, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: subjects_per_pair must be dict, "
                f"got {type(self.subjects_per_pair).__name__}"
            )
        if set(self.subjects_per_pair.keys()) != set(self.pairwise_values.keys()):
            raise ValueError(
                f"{_ERROR_PREFIX}: subjects_per_pair must have the same "
                f"keys as pairwise_values"
            )

        # Validate mean_over_pairs
        if not isinstance(self.mean_over_pairs, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: mean_over_pairs must be float, "
                f"got {type(self.mean_over_pairs).__name__}"
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

        # Validate bootstrap_cis_per_pair
        if self.bootstrap_cis_per_pair is not None:
            if not isinstance(self.bootstrap_cis_per_pair, dict):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_cis_per_pair must be "
                    f"dict, got "
                    f"{type(self.bootstrap_cis_per_pair).__name__}"
                )
            if set(self.bootstrap_cis_per_pair.keys()) != set(
                self.pairwise_values.keys()
            ):
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_cis_per_pair must have "
                    f"the same keys as pairwise_values"
                )
            for key, ci in self.bootstrap_cis_per_pair.items():
                if not isinstance(ci, BootstrapCIResult):
                    raise TypeError(
                        f"{_ERROR_PREFIX}: bootstrap CI for pair {key} must "
                        f"be BootstrapCIResult, got {type(ci).__name__}"
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

        # Validate bootstrap_distributions_per_pair
        if self.bootstrap_distributions_per_pair is not None:
            if not isinstance(self.bootstrap_distributions_per_pair, dict):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_distributions_per_pair "
                    f"must be dict, got "
                    f"{type(self.bootstrap_distributions_per_pair).__name__}"
                )
            if set(self.bootstrap_distributions_per_pair.keys()) != set(
                self.pairwise_values.keys()
            ):
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_distributions_per_pair "
                    f"must have the same keys as pairwise_values"
                )
            for key, dist in self.bootstrap_distributions_per_pair.items():
                if not isinstance(dist, np.ndarray) or dist.ndim != 1:
                    raise ValueError(
                        f"{_ERROR_PREFIX}: bootstrap distribution for "
                        f"pair {key} must be a 1-D np.ndarray"
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
            Pairwise (cohen_kappa): mean = 0.650
            [0.580, 0.720] (95% CI, bca, n=4000), correlated
            n_pairs=6; n_subjects=1617; n_raters=4; ci_method=bca

        The middle line is the interval printing itself, so a fallback to
        percentile, a count of failed draws, or a method that never ran at
        all appears there without this class having to know about any of
        them.
        """
        if self.bootstrap_ci_of_mean is not None:
            detail = str(self.bootstrap_ci_of_mean)
        elif self.bootstrap_distribution_of_mean is not None:
            n_boot = len(self.bootstrap_distribution_of_mean)
            detail = f"no CI on the mean (n={n_boot})"
        else:
            detail = ""

        if self.bootstrap_mean_method is not None:
            detail = (
                f"{detail}, {self.bootstrap_mean_method}" if detail
                else self.bootstrap_mean_method
            )

        ci_line = f"\n{detail}" if detail else ""

        return (
            f"Pairwise ({self.coefficient_name}): "
            f"mean = {self.mean_over_pairs:.3f}"
            f"{ci_line}\n"
            f"n_pairs={self.n_pairs}; "
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"n_raters={self.balance_profile.n_raters}; "
            f"ci_method={self.ci_method}"
        )

    def get_bootstrap_distributions_per_pair(self) -> dict:
        """
        Return per-pair bootstrap distributions, or raise if not stored.

        Per-pair distributions are only stored when
        n_pairs <= MAX_PAIRS_FOR_DISTRIBUTIONS to bound memory. This
        method provides a helpful error when distributions were not
        stored, rather than silently returning None.

        Returns
        -------
        dict[tuple[int, int], np.ndarray]
            Raw bootstrap draws per pair.

        Raises
        ------
        RuntimeError
            If distributions were not stored (n_pairs exceeded threshold)
            or bootstrap was not performed.
        """
        if self.bootstrap_distributions_per_pair is not None:
            return self.bootstrap_distributions_per_pair
        if self.bootstrap_mean_method is None:
            raise RuntimeError(
                f"{_ERROR_PREFIX}: per-pair bootstrap distributions are "
                f"not available because bootstrap was not performed"
            )
        raise RuntimeError(
            f"{_ERROR_PREFIX}: per-pair bootstrap distributions were not "
            f"stored because n_pairs ({self.n_pairs}) > "
            f"MAX_PAIRS_FOR_DISTRIBUTIONS ({MAX_PAIRS_FOR_DISTRIBUTIONS}). "
            f"Per-pair CIs are still available via bootstrap_cis_per_pair. "
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
            observed statistics (pairwise_values, subjects_per_pair,
            n_pairs, mean_over_pairs) and bootstrap statistics
            (bootstrap_mean_method, mean_ci_lower, mean_ci_upper,
            bootstrap_n, etc.).
            Tuple keys are converted to strings for JSON compatibility.
        """
        result = super().summary()

        # Observed
        result["n_pairs"] = self.n_pairs
        result["mean_over_pairs"] = self.mean_over_pairs
        result["pairwise_values"] = {
            str(k): v for k, v in self.pairwise_values.items()
        }
        result["subjects_per_pair"] = {
            str(k): v for k, v in self.subjects_per_pair.items()
        }

        # Bootstrap
        if self.bootstrap_mean_method is not None:
            result["bootstrap_mean_method"] = self.bootstrap_mean_method

        if self.bootstrap_ci_of_mean is not None:
            result["mean_ci"] = self.bootstrap_ci_of_mean.summary()
            result["mean_ci_lower"] = self.bootstrap_ci_of_mean.lower
            result["mean_ci_upper"] = self.bootstrap_ci_of_mean.upper

        if self.bootstrap_cis_per_pair is not None:
            result["bootstrap_cis_per_pair"] = {
                str(k): v.summary()
                for k, v in self.bootstrap_cis_per_pair.items()
            }

        if self.bootstrap_distribution_of_mean is not None:
            result["bootstrap_n"] = len(self.bootstrap_distribution_of_mean)
            result["bootstrap_mean"] = float(
                np.mean(self.bootstrap_distribution_of_mean)
            )
            result["bootstrap_std"] = float(
                np.std(self.bootstrap_distribution_of_mean)
            )

        # Per-pair distributions are NOT included in summary (too large).
        # Access them directly via the field.
        result["pair_distributions_stored"] = (
            self.bootstrap_distributions_per_pair is not None
        )

        return result
