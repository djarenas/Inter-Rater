"""
Tests for PairedRatings.

PairedRatings is exactly two raters with no gaps between them -- the shape
Cohen's kappa consumes. It inherits from BalancedRatings and runs its own
checks BEFORE the parent chain, so error messages are direct: three raters
hears "exactly 2", a gap hears "no gaps", not a complaint about balance.

Key behaviours:
- More or fewer than 2 raters is rejected.
- Gaps (empties) between 2 raters are rejected.
- A valid 2-rater, fully-crossed array is accepted.
- Error messages say PairedRatings, not BalancedRatings or Ratings.
"""

import numpy as np
import pytest

from interrater.contracts import EMPTY, INDEX_DTYPE
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.data_objects.paired_ratings import PairedRatings

from .conftest import SIMPLE_INDEX, GAPPY_INDEX


# -- Paired validation ------------------------------------------------------- #


class TestPairedValidation:
    """The two checks PairedRatings runs before the parent chain."""

    def test_accepted(self, simple_semantics):
        """SIMPLE is 2 raters, fully crossed -- accepted as a pair."""
        pr = PairedRatings(SIMPLE_INDEX.copy(), simple_semantics)
        assert pr.n_subjects == 3
        assert pr.n_raters == 2

    def test_three_raters_rejected(self, gappy_semantics):
        """More than 2 raters is rejected before any balance check."""
        with pytest.raises(ValueError, match="exactly 2"):
            PairedRatings(GAPPY_INDEX.copy(), gappy_semantics)

    def test_gaps_rejected(self, two_cat):
        """Two raters with a gap is rejected as 'no gaps'."""
        sem = RatingsSemantics(raters=2, categories=two_cat)
        arr = np.array(
            [[0,     1],
             [1, EMPTY],   # gap
             [0,     0]],
            dtype=INDEX_DTYPE,
        )
        with pytest.raises(ValueError, match="empty"):
            PairedRatings(arr, sem)


# -- Error prefix ------------------------------------------------------------ #


class TestErrorPrefix:
    """Validation errors name PairedRatings, not a parent class."""

    def test_rater_count_says_paired(self, gappy_semantics):
        """The 3-rater rejection says PairedRatings."""
        with pytest.raises(ValueError, match=r"\[PairedRatings\]"):
            PairedRatings(GAPPY_INDEX.copy(), gappy_semantics)

    def test_inherited_check_says_paired(self, two_cat):
        """Even an inherited check (shape) says PairedRatings."""
        sem = RatingsSemantics(raters=2, categories=two_cat)
        arr = np.array([0, 1, 0], dtype=INDEX_DTYPE)  # 1-D
        with pytest.raises(ValueError, match=r"\[PairedRatings\]"):
            PairedRatings(arr, sem)
