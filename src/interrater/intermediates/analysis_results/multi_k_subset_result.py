"""Result for multi-k subset permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .permutation_result import PermutationResult
from .subset_permutation_result import SubsetPermutationResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile

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
    subset_values, mean_over_subsets, and (optionally) its own bootstrap
    CIs and distributions. "Average agreement among pairs" and "average
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
        Return a compact multi-line string representation.

        Format (single k):
            MultiK Subset (fleiss_kappa, ks=[3]):
              k=3: mean=0.650 (10 subsets)
            n_subjects=100; n_raters=5; ci_method=percentile

        Format (multiple k):
            MultiK Subset (fleiss_kappa, ks=[2, 3, 4]):
              k=2: mean=0.720 (10 subsets)
              k=3: mean=0.650 (10 subsets, 95% CI: [0.580, 0.720])
              k=4: mean=0.610 (5 subsets)
            n_subjects=100; n_raters=5; ci_method=percentile
        """
        lines = [
            f"MultiK Subset ({self.coefficient_name}, ks={self.ks}):"
        ]

        for k in self.ks:
            r = self.per_k_results[k]
            line = f"  k={k}: mean={r.mean_over_subsets:.3f} ({r.n_subsets} subsets"
            if r.bootstrap_ci_of_mean is not None:
                lower, upper = r.bootstrap_ci_of_mean
                line += f", 95% CI: [{lower:.3f}, {upper:.3f}]"
            if r.bootstrap_mean_method is not None:
                line += f", {r.bootstrap_mean_method}"
            line += ")"
            lines.append(line)

        lines.append(
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"n_raters={self.balance_profile.n_raters}; "
            f"ci_method={self.ci_method}"
        )

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
