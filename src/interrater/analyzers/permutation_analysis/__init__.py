"""Permutation analysis across rater subsets.

An Analyzer computes one coefficient on one data slice. A permutation
analysis enumerates subsets (pairs, leave-one-out, arbitrary size),
calls an estimator on each, bootstraps sampling uncertainty, and
assembles the results.

Hierarchy:
  PermutationAnalysis (base)
    ├── PairwisePermutationAnalysis  (every C(n_raters, 2) pair)
    └── SubsetPermutationAnalysis    (every C(n_raters, k) subset)
"""

from .permutation_analysis import PermutationAnalysis
from .pairwise_permutation_analysis import PairwisePermutationAnalysis
from .subset_permutation_analysis import SubsetPermutationAnalysis

__all__ = [
    "PermutationAnalysis",
    "PairwisePermutationAnalysis",
    "SubsetPermutationAnalysis"
]
