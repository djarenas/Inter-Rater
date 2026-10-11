"""
Tests for BalancedRatings.

BalancedRatings is Ratings where every subject was rated by the same *number*
of raters -- not necessarily the same raters. It inherits every check from
Ratings and adds one: validate_balanced.

Key behaviours:
- Unbalanced data (varying rater counts per subject) is rejected.
- raters_per_subject returns the constant every subject shares.
- build_from_ratings converts a valid Ratings to BalancedRatings.
- build_from_ratings rejects an unbalanced Ratings.
- Error messages say BalancedRatings, not Ratings.
"""

import numpy as np
import pytest

from interrater.contracts import INDEX_DTYPE
from interrater.data_objects.balanced_ratings import BalancedRatings
from interrater.data_objects.ratings import Ratings

from .conftest import SIMPLE_INDEX, GAPPY_INDEX, GAPPY_RATERS


# -- Extra validation -------------------------------------------------------- #


class TestBalancedValidation:
    """The one check BalancedRatings adds to Ratings."""

    def test_balanced_accepted(self, simple_semantics):
        """SIMPLE is fully crossed, therefore balanced -- accepted."""
        br = BalancedRatings(SIMPLE_INDEX.copy(), simple_semantics)
        assert br.n_subjects == 3

    def test_unbalanced_rejected(self, gappy_semantics):
        """GAPPY has varying rater counts -- rejected."""
        with pytest.raises(ValueError, match="same number"):
            BalancedRatings(GAPPY_INDEX.copy(), gappy_semantics)


# -- Public API -------------------------------------------------------------- #


class TestPublicAPI:
    """Properties and methods beyond what Ratings provides."""

    def test_raters_per_subject(self, simple_semantics):
        """The constant every subject shares -- 2 for SIMPLE."""
        br = BalancedRatings(SIMPLE_INDEX.copy(), simple_semantics)
        assert br.raters_per_subject == 2

    def test_build_from_ratings_success(self, simple_ratings):
        """A balanced Ratings converts to BalancedRatings."""
        br = BalancedRatings.build_from_ratings(simple_ratings)
        assert isinstance(br, BalancedRatings)
        assert br.n_subjects == simple_ratings.n_subjects

    def test_build_from_ratings_preserves_subject_ids(self, simple_semantics):
        """Subject ids survive the conversion -- they are not regenerated."""
        custom_ids = ["patient_A", "patient_B", "patient_C"]
        r = Ratings(SIMPLE_INDEX.copy(), simple_semantics, subject_ids=custom_ids)
        br = BalancedRatings.build_from_ratings(r)
        np.testing.assert_array_equal(br.subject_ids, custom_ids)

    def test_build_from_ratings_unbalanced_rejected(self, gappy_ratings):
        """An unbalanced Ratings cannot become BalancedRatings."""
        with pytest.raises(ValueError, match="same number"):
            BalancedRatings.build_from_ratings(gappy_ratings)


# -- Error prefix ------------------------------------------------------------ #


class TestErrorPrefix:
    """Validation errors name BalancedRatings, not a parent class."""

    def test_says_balanced_ratings(self, gappy_semantics):
        """The unbalanced rejection says BalancedRatings."""
        with pytest.raises(ValueError, match=r"\[BalancedRatings\]"):
            BalancedRatings(GAPPY_INDEX.copy(), gappy_semantics)

    def test_inherited_check_says_balanced_ratings(self, two_cat):
        """Even an inherited check (shape) says BalancedRatings."""
        sem = BalancedRatings.__bases__[0]  # just for clarity -- it's Ratings
        from interrater.semantics.ratings_semantics import RatingsSemantics

        sem = RatingsSemantics(raters=2, categories=two_cat)
        arr = np.array([0, 1, 0], dtype=INDEX_DTYPE)  # 1-D
        with pytest.raises(ValueError, match=r"\[BalancedRatings\]"):
            BalancedRatings(arr, sem)
