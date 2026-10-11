"""Visualization functions for inter-rater agreement results.

v2.0 ships standalone plotting functions. A declarative config + plotter
pattern (similar to github.com/djarenas/eventus) is planned for v2.x.
"""

from .plot_pairwise_permutation import plot_pairwise_permutation
from .plot_subset_sorted import plot_subset_sorted
from .plot_subset_strip import plot_subset_strip

__all__ = [
    "plot_pairwise_permutation",
    "plot_subset_sorted",
    "plot_subset_strip"
]
