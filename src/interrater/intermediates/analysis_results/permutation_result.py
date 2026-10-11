"""Base result for permutation analysis (pairwise or subset)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .analysis_result import AnalysisResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[PermutationResult] Error"


@dataclass(frozen=True)
class PermutationResult(AnalysisResult):
    """
    I am the base result for permutation analysis across rater subsets.

    I am the permutation-branch counterpart of PointEstimateResult. Where
    PointEstimateResult holds one number for one coefficient, my children
    hold C(n_raters, k) numbers -- one per subset or pair -- plus summary
    statistics and bootstrap intervals. The two branches share
    AnalysisResult as their root but are otherwise independent.

    I hold only what both children genuinely share:
    - Which coefficient was computed
    - How confidence intervals were computed
    - The structural metadata of the dataset

    All statistical substance -- observed values, means, bootstrap
    distributions, per-pair or per-subset CIs -- lives on the children,
    where the vocabulary and semantics belong. Pairwise talks about
    "pairs," subset talks about "subsets." The parent does not impose
    either vocabulary.

    Subclasses
    ----------
    PairwisePermutationResult
        Coefficient computed on every C(n_raters, 2) pair of raters.
    SubsetPermutationResult
        Coefficient computed on every C(n_raters, k) subset of raters.

    Attributes
    ----------
    coefficient_name : str
        Which coefficient was computed: "cohen_kappa", "scotts_pi",
        "conger_kappa", etc.
    ci_method : str
        How confidence intervals were computed: "percentile" or "bca".
        Permutation results have no closed-form null variance (subsets
        sharing a rater are correlated), so "parametric_null" is not
        available.
    balance_profile : BalanceProfile
        Structural metadata about the full dataset (before subsetting):
        n_subjects, n_raters, balance structure.
    """

    # Attribute declaration

    coefficient_name: str
    ci_method: str
    balance_profile: BalanceProfile

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate permutation result attributes."""

        # Validate coefficient_name
        if not isinstance(self.coefficient_name, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: coefficient_name must be str, "
                f"got {type(self.coefficient_name).__name__}"
            )
        if not self.coefficient_name:
            raise ValueError(
                f"{_ERROR_PREFIX}: coefficient_name must be non-empty"
            )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or 'bca', "
                f"got '{self.ci_method}'. Permutation results have no "
                f"closed-form null variance (subsets are correlated), so "
                f"'parametric_null' is not available."
            )

        # Validate balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None"
            )

        from ..data_profiles.balance_profile import BalanceProfile

        if not isinstance(self.balance_profile, BalanceProfile):
            raise TypeError(
                f"{_ERROR_PREFIX}: balance_profile must be BalanceProfile, "
                f"got {type(self.balance_profile).__name__}"
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """
        Return shared permutation metadata in plain Python types.

        Subclasses extend this with their own statistical content.

        Returns
        -------
        dict
            Keys: coefficient_name, ci_method, balance.
        """
        return {
            "coefficient_name": self.coefficient_name,
            "ci_method": self.ci_method,
            "balance": self.balance_profile.summary(),
        }
