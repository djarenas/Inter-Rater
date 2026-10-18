"""Permutation results: a coefficient recomputed over subsets of the raters.

Every class here descends from ``PermutationResult``. Where the analyzer results
answer "what is the coefficient", these answer "how does it move when the rater
set changes" -- pairwise, leave-k-out at one k, or across several k at once.
"""

from .permutation_result import PermutationResult
from .pairwise_permutation_result import PairwisePermutationResult
from .subset_permutation_result import SubsetPermutationResult
from .multi_k_subset_result import MultiKSubsetResult

__all__ = [
    "PermutationResult",
    "PairwisePermutationResult",
    "SubsetPermutationResult",
    "MultiKSubsetResult",
]
