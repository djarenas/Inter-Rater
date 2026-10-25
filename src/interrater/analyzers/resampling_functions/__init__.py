"""Resampling functions for confidence intervals, null tests, and diagnostics.

Subject bootstrap and jackknife work for all coefficients. Rater jackknife and
subset are diagnostic tools. Permutation functions generate null distributions
for hypothesis testing. See README.md for rationale and developer workflow.
"""

from .subject_resampling import (
    subject_bootstrap,
    subject_jackknife,
    subject_bootstrap_from_nij_matrix,
    subject_jackknife_from_nij_matrix,
)
from .rater_resampling import rater_jackknife, rater_subset
from .permutation_resampling import within_rater_shuffle, within_subject_shuffle, pool_and_redistribute

__all__ = [
    "subject_bootstrap",
    "subject_jackknife",
    "subject_bootstrap_from_nij_matrix",
    "subject_jackknife_from_nij_matrix",
    "rater_jackknife",
    "rater_subset",
    "within_rater_shuffle",
    "within_subject_shuffle",
    "pool_and_redistribute",
]
