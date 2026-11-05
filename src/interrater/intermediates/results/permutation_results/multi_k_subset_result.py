"""Result for multi-k subset permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ..analyzer_results import result_class_for
from .permutation_result import (
    INDEPENDENT_MEAN_CAVEAT,
    PermutationResult,
    interval_notes,
    raters_phrase,
)
from .subset_permutation_result import SubsetPermutationResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile
    from ...data_profiles.rater_profile import RaterProfile

_ERROR_PREFIX = "[MultiKSubsetResult] Error"


@dataclass(frozen=True)
class MultiKSubsetResult(PermutationResult):
    """
    I am the result of computing a coefficient across one or more subset sizes.

    I wrap one SubsetPermutationResult per k value. Whether the user
    requested subset_k=3, subset_k=[2, 3, 4], or subset_k="all", the
    analysis always produces me. A single k produces a dict with one key --
    the structure is identical either way, so the visualizer has one code
    path.

    Each per-k result is a complete SubsetPermutationResult with its own
    subset_results, mean_over_subsets, and (optionally) its own bootstrap
    of the mean. "Average agreement among pairs" and "average
    agreement among triples" are different quantities, so each k gets
    its own mean bootstrap.

    I carry the shared metadata (coefficient_name, ci_method,
    balance_profile) once at my level. Every per-k result shares these
    values -- I validate that on construction.

    Attributes
    ----------
    coefficient_name : str
        Inherited. Which coefficient was computed.
    ci_method : str
        Inherited. How CIs were computed.
    balance_profile : BalanceProfile
        Inherited. Structural metadata of the full dataset.
    rater_profile : RaterProfile
        Inherited. The panel's raters, shared by every k.
    per_k_results : dict[int, SubsetPermutationResult]
        One result per subset size. Keys are k values (each >= 2),
        values are complete SubsetPermutationResult objects.
    ks : list[int]
        Sorted list of k values that were computed. Must match
        per_k_results keys.
    """

    # Attribute declaration

    # Inherited from parent class PermutationResult
    coefficient_name: str
    ci_method: str
    balance_profile: BalanceProfile
    rater_profile: RaterProfile

    # Multi-k container
    per_k_results: dict  # dict[int, SubsetPermutationResult]
    ks: list  # list[int], sorted

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate multi-k subset result attributes."""

        # Validate ks
        if not isinstance(self.ks, list):
            raise TypeError(
                f"{_ERROR_PREFIX}: ks must be list, "
                f"got {type(self.ks).__name__}"
            )
        if len(self.ks) == 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: ks must not be empty"
            )
        for i, k in enumerate(self.ks):
            if not isinstance(k, int):
                raise TypeError(
                    f"{_ERROR_PREFIX}: ks elements must be int, "
                    f"got {type(k).__name__} at index {i}"
                )
            if k < 2:
                raise ValueError(
                    f"{_ERROR_PREFIX}: ks values must be >= 2, "
                    f"got {k} at index {i}"
                )
        if self.ks != sorted(self.ks):
            raise ValueError(
                f"{_ERROR_PREFIX}: ks must be sorted, got {self.ks}"
            )
        if len(self.ks) != len(set(self.ks)):
            raise ValueError(
                f"{_ERROR_PREFIX}: ks must not contain duplicates, "
                f"got {self.ks}"
            )

        # Validate per_k_results
        if not isinstance(self.per_k_results, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: per_k_results must be dict, "
                f"got {type(self.per_k_results).__name__}"
            )

        # Keys must match ks
        if set(self.per_k_results.keys()) != set(self.ks):
            raise ValueError(
                f"{_ERROR_PREFIX}: per_k_results keys "
                f"{sorted(self.per_k_results.keys())} do not match "
                f"ks {self.ks}"
            )

        # Each value must be a SubsetPermutationResult
        for k, result in self.per_k_results.items():
            if not isinstance(result, SubsetPermutationResult):
                raise TypeError(
                    f"{_ERROR_PREFIX}: per_k_results[{k}] must be "
                    f"SubsetPermutationResult, got "
                    f"{type(result).__name__}"
                )

        # Cross-validate: every per-k result must share our metadata
        for k, result in self.per_k_results.items():
            if result.coefficient_name != self.coefficient_name:
                raise ValueError(
                    f"{_ERROR_PREFIX}: per_k_results[{k}] has "
                    f"coefficient_name='{result.coefficient_name}', "
                    f"expected '{self.coefficient_name}'"
                )
            if result.ci_method != self.ci_method:
                raise ValueError(
                    f"{_ERROR_PREFIX}: per_k_results[{k}] has "
                    f"ci_method='{result.ci_method}', "
                    f"expected '{self.ci_method}'"
                )
            if result.balance_profile is not self.balance_profile:
                raise ValueError(
                    f"{_ERROR_PREFIX}: per_k_results[{k}] has a "
                    f"different balance_profile than the parent"
                )
            if result.rater_profile is not self.rater_profile:
                raise ValueError(
                    f"{_ERROR_PREFIX}: per_k_results[{k}] has a "
                    f"different rater_profile than the parent"
                )
            if result.subset_k != k:
                raise ValueError(
                    f"{_ERROR_PREFIX}: per_k_results[{k}] has "
                    f"subset_k={result.subset_k}, expected {k}"
                )

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

    def __getitem__(self, k: int) -> SubsetPermutationResult:
        """
        Access a specific k's result.

        Parameters
        ----------
        k : int
            The subset size to retrieve.

        Returns
        -------
        SubsetPermutationResult
            The result for that subset size.

        Raises
        ------
        KeyError
            If k is not in ks.
        """
        if k not in self.per_k_results:
            raise KeyError(
                f"{_ERROR_PREFIX}: k={k} not in results. "
                f"Available: {self.ks}"
            )
        return self.per_k_results[k]

    def __len__(self) -> int:
        """Return the number of k values computed."""
        return len(self.ks)

    def __iter__(self):
        """Iterate over k values in sorted order."""
        return iter(self.ks)

    def __str__(self) -> str:
        """
        Return my printed form: one line per k.

            Subsets, Fleiss' κ: 199 subjects, 4 raters
              k=2   6 subsets   mean 0.492 [0.367, 0.601]
              k=3   4 subsets   mean 0.493 [0.370, 0.600]
              k=4   1 subset    mean 0.493 [0.372, 0.600]
              intervals: 95% CI, percentile, n=300, correlated bootstrap

        Each k's full list prints from print(result.per_k_results[k]). A k
        whose bootstrap ran and produced nothing says so, rather than
        printing like a k that was never bootstrapped, and failed draws or a
        fallback are named on that k's line: an interval from a fraction of
        its draws is not the same claim as one from all of them. An
        independent mean adds the "!" line.
        """
        name = result_class_for(self.coefficient_name).NAME
        bp = self.balance_profile
        lines = [f"Subsets, {name}: {bp.n_subjects:,} subjects, {raters_phrase(self.rater_profile)}"]

        first_ci, method = None, None
        rows = []
        for k in self.ks:
            r = self.per_k_results[k]
            word = "subset" if r.n_subsets == 1 else "subsets"
            mean_ci = r.bootstrap_ci_of_mean
            if mean_ci is None:
                interval = ""
            elif not mean_ci.is_defined:
                interval = "no CI: every draw failed"
            else:
                interval = f"[{mean_ci.lower:.3f}, {mean_ci.upper:.3f}]"
                notes = interval_notes(mean_ci)
                if notes:
                    interval += f"  {'; '.join(notes)}"
            if mean_ci is not None and first_ci is None:
                first_ci, method = mean_ci, r.bootstrap_mean_method
            rows.append((f"k={k}", f"{r.n_subsets:,} {word}", f"mean {r.mean_over_subsets:.3f}", interval))

        widths = [max(len(row[i]) for row in rows) for i in range(3)]
        for row in rows:
            cells = [f"{row[i]:<{widths[i]}}" for i in range(3)] + ([row[3]] if row[3] else [])
            lines.append("  " + "   ".join(cells).rstrip())

        if first_ci is None:
            lines.append("  no interval for the mean")
        else:
            lines.append(
                f"  intervals: {first_ci.confidence_level:.0%} CI, "
                f"{first_ci.method_intended}, n={first_ci.n_resamples}, "
                f"{method} bootstrap"
            )
            if method == "independent":
                lines.append(f"  ! {INDEPENDENT_MEAN_CAVEAT}")
        return "\n".join(lines)

    def summary(self) -> dict:
        """
        Return all k results in plain Python types.

        Returns
        -------
        dict
            Parent keys (coefficient_name, ci_method, balance) plus
            ks (list of k values) and per_k (dict keyed by k, each
            value being that k's summary()).
        """
        result = super().summary()
        result["ks"] = list(self.ks)
        result["per_k"] = {
            k: self.per_k_results[k].summary()
            for k in self.ks
        }
        return result
