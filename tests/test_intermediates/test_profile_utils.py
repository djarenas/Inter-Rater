"""
Tests for _profile_utils: prevalence, prevalence_per_rater, coverage_per_rater.

These are free functions -- numpy in, numpy out -- that compute reporting
arithmetic over what raters actually chose. They do not validate their
inputs; that is the caller's job.

Hand-verified datasets small enough that every count can be checked on
fingers.
"""

import numpy as np
import pytest

from interrater.contracts import EMPTY, INDEX_DTYPE
from interrater.intermediates._profile_utils import (
    prevalence,
    prevalence_per_rater,
    coverage_per_rater,
)

from .conftest import SIMPLE_INDEX, SIMPLE_K


# -- Test data --------------------------------------------------------------- #

# SIMPLE_INDEX from conftest (4 subjects, 2 raters, K=2, no empties):
#
#   [[0, 0],    rater0 picks: 0, 1, 0, 1 -> counts [2, 2]
#    [1, 1],    rater1 picks: 0, 1, 1, 0 -> counts [2, 2]
#    [0, 1],
#    [1, 0]]
#
# Flattened present ratings: [0, 0, 1, 1, 0, 1, 1, 0] -> counts [4, 4]
# Prevalence: [0.5, 0.5]

# With empties: 3 subjects, 2 raters, K=2.
#
#   [[0,  0],    both present
#    [1, -1],    rater1 missing
#    [0,  1]]   both present
#
# All present ratings: [0, 0, 1, 0, 1] -> counts [3, 2], prevalence [3/5, 2/5]
# Rater 0 present: [0, 1, 0] -> counts [2, 1], prevalence [2/3, 1/3]
# Rater 1 present: [0, 1]    -> counts [1, 1], prevalence [1/2, 1/2]

GAPPY = np.array(
    [[ 0,  0],
     [ 1, EMPTY],
     [ 0,  1]],
    dtype=INDEX_DTYPE,
)

# All empty: 2 subjects, 2 raters, everything is -1.
ALL_EMPTY = np.full((2, 2), EMPTY, dtype=INDEX_DTYPE)


# ============================================================================ #
#                               prevalence                                     #
# ============================================================================ #


class TestPrevalence:
    """Share of present ratings in each category."""

    def test_fully_crossed(self):
        """No empties: all 8 cells present, 4 zeros and 4 ones."""
        result = prevalence(SIMPLE_INDEX, SIMPLE_K)
        np.testing.assert_allclose(result, [0.5, 0.5], atol=1e-12)

    def test_empties_excluded(self):
        """Empties are excluded from both numerator and denominator."""
        result = prevalence(GAPPY, 2)
        # 5 present ratings: [0, 0, 1, 0, 1] -> [3, 2] / 5
        np.testing.assert_allclose(result, [3 / 5, 2 / 5], atol=1e-12)

    def test_all_empty_returns_zeros(self):
        """When every cell is empty, prevalence is zero for all categories."""
        result = prevalence(ALL_EMPTY, 2)
        np.testing.assert_array_equal(result, [0.0, 0.0])

    def test_skewed_prevalence(self):
        """One category dominates."""
        # 4 subjects, 1 rater, K=2: all pick category 0.
        idx = np.array([[0], [0], [0], [0]], dtype=INDEX_DTYPE)
        result = prevalence(idx, 2)
        np.testing.assert_allclose(result, [1.0, 0.0], atol=1e-12)

    def test_output_length_matches_k(self):
        """Output has exactly k elements, even if some categories unused."""
        idx = np.array([[0], [0]], dtype=INDEX_DTYPE)
        result = prevalence(idx, 5)
        assert len(result) == 5


# ============================================================================ #
#                          prevalence_per_rater                                #
# ============================================================================ #


class TestPrevalencePerRater:
    """Per-rater category shares."""

    def test_fully_crossed(self):
        """Each rater's prevalence computed independently."""
        result = prevalence_per_rater(SIMPLE_INDEX, SIMPLE_K)
        # Rater 0: picks [0, 1, 0, 1] -> [2, 2] / 4 = [0.5, 0.5]
        # Rater 1: picks [0, 1, 1, 0] -> [2, 2] / 4 = [0.5, 0.5]
        assert result.shape == (2, 2)
        np.testing.assert_allclose(result[0], [0.5, 0.5], atol=1e-12)
        np.testing.assert_allclose(result[1], [0.5, 0.5], atol=1e-12)

    def test_with_empties(self):
        """Empties excluded per rater."""
        result = prevalence_per_rater(GAPPY, 2)
        # Rater 0: [0, 1, 0] -> [2, 1] / 3 = [2/3, 1/3]
        # Rater 1: [0, _, 1] -> [1, 1] / 2 = [1/2, 1/2]
        assert result.shape == (2, 2)
        np.testing.assert_allclose(result[0], [2 / 3, 1 / 3], atol=1e-12)
        np.testing.assert_allclose(result[1], [1 / 2, 1 / 2], atol=1e-12)

    def test_shape(self):
        """Output shape is (n_raters, k)."""
        idx = np.array([[0, 1, 0]], dtype=INDEX_DTYPE)
        result = prevalence_per_rater(idx, 3)
        assert result.shape == (3, 3)


# ============================================================================ #
#                          coverage_per_rater                                  #
# ============================================================================ #


class TestCoveragePerRater:
    """Fraction of subjects each rater rated."""

    def test_fully_crossed(self):
        """No empties: every rater covered every subject -> 1.0."""
        result = coverage_per_rater(SIMPLE_INDEX)
        np.testing.assert_allclose(result, [1.0, 1.0], atol=1e-12)

    def test_with_empties(self):
        """Rater with an empty cell has fractional coverage."""
        result = coverage_per_rater(GAPPY)
        # 3 subjects. Rater 0 rated all 3, rater 1 rated 2.
        np.testing.assert_allclose(result, [3 / 3, 2 / 3], atol=1e-12)

    def test_output_length(self):
        """Output has one entry per rater."""
        result = coverage_per_rater(SIMPLE_INDEX)
        assert len(result) == 2
