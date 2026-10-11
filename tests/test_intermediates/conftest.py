"""
Shared fixtures for intermediates tests.

I hold small, hand-verifiable datasets that anchor every test in this
directory. Each is small enough that every property -- counts, categories,
weights, balance -- can be checked by counting on fingers.

Two NijMatrix datasets:

- BALANCED_NIJ: 4 subjects, 3 categories, 3 raters per subject. Balanced
  because every row sums to 3.
- UNBALANCED_NIJ: 4 subjects, 2 categories, variable raters (2 to 4).
  Unbalanced because row sums differ.

One EstimatorInput dataset:

- SIMPLE_INDEX: 4 subjects, 2 raters, K=2. Fully crossed (no empties).
  Small enough that nij counts, marginals, and agreement can be verified
  on paper.

One RaterSharedMatrix dataset:

- SIMPLE_SHARED: 3 raters. Diagonal is each rater's total. Off-diagonal
  is how many subjects both raters rated. Small enough that every pair
  can be verified by inspection.

All are frozen as module-level constants so no test can mutate what another
test sees. Fixtures return fresh objects built from those constants.
"""

import numpy as np
import pytest

from interrater.contracts import INDEX_DTYPE
from interrater.intermediates.estimator_parts.nij_matrix import NijMatrix
from interrater.intermediates.estimator_parts.weights import Weights
from interrater.intermediates.estimator_parts.estimator_input import EstimatorInput
from interrater.intermediates.data_profiles.rater_shared_matrix import RaterSharedMatrix


# -- NijMatrix datasets ----------------------------------------------------- #

# BALANCED: 4 subjects, 3 categories, 3 raters per subject.
#
#   subject | cat0 | cat1 | cat2 | total
#   --------+------+------+------+------
#   s0      |  2   |  1   |  0   |  3
#   s1      |  0   |  3   |  0   |  3
#   s2      |  1   |  1   |  1   |  3
#   s3      |  0   |  0   |  3   |  3
#
# n_subjects = 4, n_categories = 3, every row sums to 3.

BALANCED_NIJ = np.array(
    [[2, 1, 0],
     [0, 3, 0],
     [1, 1, 1],
     [0, 0, 3]],
    dtype=np.int32,
)

# UNBALANCED: 4 subjects, 2 categories, variable raters per subject.
#
#   subject | cat0 | cat1 | total
#   --------+------+------+------
#   s0      |  2   |  1   |  3
#   s1      |  1   |  1   |  2
#   s2      |  3   |  1   |  4
#   s3      |  0   |  2   |  2
#
# n_subjects = 4, n_categories = 2, row sums: [3, 2, 4, 2].

UNBALANCED_NIJ = np.array(
    [[2, 1],
     [1, 1],
     [3, 1],
     [0, 2]],
    dtype=np.int32,
)


# -- EstimatorInput dataset -------------------------------------------------- #

# 4 subjects, 2 raters, K=2, fully crossed.
#
#   subject | rater0 | rater1
#   --------+--------+--------
#   s0      |   0    |   0       agree on 0
#   s1      |   1    |   1       agree on 1
#   s2      |   0    |   1       disagree
#   s3      |   1    |   0       disagree
#
# Po = 2/4 = 0.5

SIMPLE_INDEX = np.array(
    [[0, 0],
     [1, 1],
     [0, 1],
     [1, 0]],
    dtype=INDEX_DTYPE,
)

SIMPLE_K = 2


# -- RaterSharedMatrix dataset ----------------------------------------------- #

# 3 raters. Diagonal is each rater's own total (how many subjects they rated).
# Off-diagonal is how many subjects both raters rated.
#
#   rater | r0  | r1  | r2
#   ------+-----+-----+----
#   r0    | 10  |  4  |  1
#   r1    |  4  |  8  |  2
#   r2    |  1  |  2  |  9
#
# n_raters = 3, n_pairs = 3.
# Pairs (a < b): (0,1)=4, (0,2)=1, (1,2)=2.
# Thinnest pair: (0, 2, 1).

SIMPLE_SHARED = np.array(
    [[10, 4, 1],
     [ 4, 8, 2],
     [ 1, 2, 9]],
    dtype=np.int32,
)


# -- Fixtures ---------------------------------------------------------------- #


@pytest.fixture
def balanced_nij():
    """Balanced NijMatrix: 4 subjects, 3 categories, 3 raters each."""
    return NijMatrix(BALANCED_NIJ.copy())


@pytest.fixture
def unbalanced_nij():
    """Unbalanced NijMatrix: 4 subjects, 2 categories, variable raters."""
    return NijMatrix(UNBALANCED_NIJ.copy())


@pytest.fixture
def simple_input():
    """EstimatorInput: 4 subjects, 2 raters, K=2, no empties."""
    return EstimatorInput(SIMPLE_INDEX.copy(), SIMPLE_K)


@pytest.fixture
def identity_weights():
    """Identity weights for K=2."""
    return Weights.identity(2)


@pytest.fixture
def linear_weights_k3():
    """Linear weights for K=3."""
    return Weights.linear(3)


@pytest.fixture
def simple_shared():
    """RaterSharedMatrix: 3 raters, pairs (0,1)=4, (0,2)=1, (1,2)=2."""
    return RaterSharedMatrix(SIMPLE_SHARED.copy())
