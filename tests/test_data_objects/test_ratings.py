"""
Tests for Ratings.

Ratings is the validated data object: a dense (n_subjects, n_raters) array of
category positions, a label for each row, and the semantics that say what all
of it means. If a Ratings exists, it is complete -- the constructor raises
rather than repair.

Key behaviours:
- Subject ids are generated (subject0, subject1, ...) when none are given.
- dtype is cast to INDEX_DTYPE; floats and non-integer types are rejected
  before the cast so the error message explains what happened.
- Shape must be 2-D with at least MIN_SUBJECTS rows.
- The array must match its semantics: rater count, n_subjects (if stated),
  values in category range.
- Subject ids must match row count and must be unique (unless the semantics
  says otherwise).
- No rater may be entirely idle; every subject must be rated at least twice.
- Ratings cannot be changed once built: it copies its input, its arrays are
  read-only, and its attributes cannot be reassigned or deleted.
- check_my_balance, check_my_crossing, check_my_pairing delegate to
  array_utils and report structural properties.
- to_nij_matrix builds the counts matrix verified against hand computation.
"""

import numpy as np
import pytest

from interrater.contracts import EMPTY, INDEX_DTYPE, MIN_SUBJECTS
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.data_objects.ratings import Ratings

from .conftest import (
    SIMPLE_INDEX,
    SIMPLE_NIJ,
)


# -- Construction ------------------------------------------------------------ #


class TestConstruction:
    """How Ratings interprets its arguments and what it stores."""

    def test_generates_subject_ids(self, simple_ratings):
        """When no subject_ids are given, positional names are generated."""
        expected = np.array(["subject0", "subject1", "subject2"], dtype=object)
        np.testing.assert_array_equal(simple_ratings.subject_ids, expected)

    def test_custom_subject_ids(self, simple_semantics):
        """Explicit subject_ids are stored as given."""
        ids = ["doc_A", "doc_B", "doc_C"]
        r = Ratings(SIMPLE_INDEX.copy(), simple_semantics, subject_ids=ids)
        np.testing.assert_array_equal(r.subject_ids, ids)

    def test_n_subjects(self, simple_ratings):
        """n_subjects is the row count."""
        assert simple_ratings.n_subjects == 3

    def test_n_raters(self, simple_ratings):
        """n_raters is the column count."""
        assert simple_ratings.n_raters == 2

    def test_n_categories(self, simple_ratings):
        """n_categories delegates to semantics."""
        assert simple_ratings.n_categories == 2

    def test_dtype_cast(self, simple_semantics):
        """An int32 array is cast to INDEX_DTYPE, not rejected."""
        arr = SIMPLE_INDEX.astype(np.int32)
        r = Ratings(arr, simple_semantics)
        assert r.category_index.dtype == INDEX_DTYPE

    def test_frozen_category_index(self, simple_ratings):
        """category_index is read-only after construction."""
        with pytest.raises(ValueError, match="read-only"):
            simple_ratings.category_index[0, 0] = 1

    def test_frozen_subject_ids(self, simple_ratings):
        """subject_ids is read-only after construction."""
        with pytest.raises(ValueError, match="read-only"):
            simple_ratings.subject_ids[0] = "tampered"


class TestCannotBeChanged:
    """
    Once built, a Ratings cannot change; a change is a new Ratings.

    Writing through it was already refused. These are the two ways that got
    through before 2.1.0: editing the array it was built from, which it
    shared, and reassigning an attribute.
    """

    def test_edit_to_callers_array_does_not_reach_it(self, simple_semantics):
        arr = SIMPLE_INDEX.copy()  # int8 already: no cast, so it used to be shared
        r = Ratings(arr, simple_semantics)
        arr[0, 0] = 1
        assert r.category_index[0, 0] == SIMPLE_INDEX[0, 0]

    def test_edit_to_callers_ids_does_not_reach_it(self, simple_semantics):
        ids = np.array(["a", "b", "c"], dtype=object)
        r = Ratings(SIMPLE_INDEX.copy(), simple_semantics, subject_ids=ids)
        ids[0] = "tampered"
        assert r.subject_ids[0] == "a"

    @pytest.mark.parametrize("attr", ["category_index", "subject_ids", "semantics", "anything_new"])
    def test_attributes_cannot_be_set(self, simple_ratings, attr):
        with pytest.raises(
            AttributeError,
            match=rf"Ratings is read-only, so {attr} cannot be set. Build a new Ratings",
        ):
            setattr(simple_ratings, attr, None)

    def test_attributes_cannot_be_deleted(self, simple_ratings):
        with pytest.raises(AttributeError, match="category_index cannot be deleted"):
            del simple_ratings.category_index

    def test_subclass_names_itself(self, simple_semantics):
        from interrater.data_objects.paired_ratings import PairedRatings

        r = PairedRatings(SIMPLE_INDEX.copy(), simple_semantics)
        with pytest.raises(AttributeError, match=r"\[PairedRatings\] Error: PairedRatings is read-only"):
            r.category_index = None


# -- Dtype validation -------------------------------------------------------- #


class TestValidationDtype:
    """dtype checks that run before any value inspection."""

    def test_float_rejected(self, simple_semantics):
        """Float arrays usually mean NaN was used for absent ratings."""
        arr = SIMPLE_INDEX.astype(np.float64)
        with pytest.raises(TypeError, match="integer array"):
            Ratings(arr, simple_semantics)

    def test_string_rejected(self, simple_semantics):
        """String arrays are raw category names, not positions."""
        arr = np.array([["include", "exclude"], ["exclude", "include"]], dtype=object)
        with pytest.raises(TypeError, match="integer array"):
            Ratings(arr, simple_semantics)


# -- Shape validation -------------------------------------------------------- #


class TestValidationShape:
    """Shape checks on the array alone."""

    def test_1d_rejected(self, simple_semantics):
        """A 1-D array has no rater dimension."""
        arr = np.array([0, 1, 0], dtype=INDEX_DTYPE)
        with pytest.raises(ValueError, match="2-dimensional"):
            Ratings(arr, simple_semantics)

    def test_too_few_subjects(self, simple_semantics):
        """Fewer than MIN_SUBJECTS rows is rejected."""
        arr = np.array([[0, 1]], dtype=INDEX_DTYPE)
        assert arr.shape[0] < MIN_SUBJECTS
        with pytest.raises(ValueError, match="at least"):
            Ratings(arr, simple_semantics)


# -- Array vs semantics ------------------------------------------------------ #


class TestValidationAgainstSemantics:
    """Checks that hold the array accountable to the declaration."""

    def test_wrong_rater_count(self, two_cat):
        """More columns than declared raters is rejected."""
        sem = RatingsSemantics(raters=2, categories=two_cat)
        arr = np.array([[0, 1, 0], [1, 0, 1]], dtype=INDEX_DTYPE)
        with pytest.raises(ValueError, match="rater column"):
            Ratings(arr, sem)

    def test_wrong_n_subjects(self, two_cat):
        """When semantics declares n_subjects, the row count must match."""
        sem = RatingsSemantics(raters=2, categories=two_cat, n_subjects=10)
        with pytest.raises(ValueError, match="subject row"):
            Ratings(SIMPLE_INDEX.copy(), sem)  # 3 rows, not 10

    def test_value_out_of_range(self, simple_semantics):
        """A category index >= K is out of range."""
        arr = SIMPLE_INDEX.copy()
        arr[0, 0] = 5  # K=2, so 5 is invalid
        with pytest.raises(ValueError, match="outside"):
            Ratings(arr, simple_semantics)

    def test_values_checked_before_the_int8_cast(self, simple_semantics):
        """int8 wraps: 256 and 257 would become 0 and 1, valid categories."""
        arr = np.array([[256, 0], [1, 257], [0, 1]])
        with pytest.raises(ValueError, match=r"outside its categories \(e.g. \[256, 257\]\)"):
            Ratings(arr, simple_semantics)

    def test_uint8_255_is_not_empty(self, two_cat):
        """A uint8 255 would become -1, EMPTY: a rating silently missing.

        Three raters, so the subject still has two ratings after the wrap
        and the twice-rated rule cannot refuse it by luck.
        """
        sem = RatingsSemantics(raters=3, categories=two_cat)
        arr = np.array([[255, 0, 1], [1, 1, 0], [0, 1, 1]], dtype=np.uint8)
        with pytest.raises(ValueError, match=r"outside its categories \(e.g. \[255\]\)"):
            Ratings(arr, sem)


# -- Subject id validation -------------------------------------------------- #


class TestValidationSubjectIds:
    """Subject id checks."""

    def test_length_mismatch(self, simple_semantics):
        """subject_ids must have one entry per row."""
        with pytest.raises(ValueError, match="subject id"):
            Ratings(SIMPLE_INDEX.copy(), simple_semantics, subject_ids=["a", "b"])

    def test_duplicate_ids_rejected(self, simple_semantics):
        """Duplicate ids with unique_subjects=True is rejected."""
        ids = ["doc_A", "doc_B", "doc_A"]
        with pytest.raises(ValueError, match="more than once"):
            Ratings(SIMPLE_INDEX.copy(), simple_semantics, subject_ids=ids)

    def test_duplicate_ids_allowed_when_not_unique(self, two_cat):
        """Duplicate ids are fine when unique_subjects=False."""
        sem = RatingsSemantics(raters=2, categories=two_cat, unique_subjects=False)
        ids = ["doc_A", "doc_B", "doc_A"]
        r = Ratings(SIMPLE_INDEX.copy(), sem, subject_ids=ids)
        assert r.n_subjects == 3


# -- Presence validation ----------------------------------------------------- #


class TestValidationPresence:
    """No rater may be entirely idle; every subject rated at least twice."""

    def test_idle_rater(self, two_cat):
        """A rater who rated nothing at all is rejected."""
        sem = RatingsSemantics(raters=3, categories=two_cat)
        # Third column all EMPTY -- idle rater
        arr = np.array(
            [[0, 1, EMPTY],
             [1, 0, EMPTY],
             [0, 1, EMPTY]],
            dtype=INDEX_DTYPE,
        )
        with pytest.raises(ValueError, match="rated nothing"):
            Ratings(arr, sem)

    def test_unrated_subject(self, simple_semantics):
        """A subject rated by nobody is rejected."""
        arr = SIMPLE_INDEX.copy()
        arr[1, :] = EMPTY  # subject1 has no ratings
        with pytest.raises(ValueError, match="rated by nobody"):
            Ratings(arr, simple_semantics)

    def test_subject_rated_once(self, simple_semantics):
        """A subject rated by one rater joins no comparison either."""
        arr = SIMPLE_INDEX.copy()
        arr[2, 1] = EMPTY  # subject2 has one rating
        with pytest.raises(ValueError, match="1 rated by only one rater"):
            Ratings(arr, simple_semantics)

    def test_nobody_and_once_in_one_message(self, simple_semantics):
        arr = SIMPLE_INDEX.copy()
        arr[1, :] = EMPTY
        arr[2, 1] = EMPTY
        with pytest.raises(
            ValueError, match="1 rated by nobody and 1 rated by only one rater",
        ):
            Ratings(arr, simple_semantics)

    def test_two_ratings_is_enough(self, gappy_semantics):
        """Two of three raters is enough; the third may be absent."""
        arr = np.array(
            [[0, 1, EMPTY],
             [1, EMPTY, 1],
             [EMPTY, 0, 0]],
            dtype=INDEX_DTYPE,
        )
        assert Ratings(arr, gappy_semantics).n_subjects == 3


# -- Structural checks ------------------------------------------------------ #


class TestCheckMethods:
    """check_my_balance, check_my_crossing, check_my_pairing."""

    def test_balance_simple(self, simple_ratings):
        """SIMPLE is fully crossed, therefore balanced."""
        assert simple_ratings.check_my_balance() is True

    def test_balance_gappy(self, gappy_ratings):
        """GAPPY has unequal rater counts per subject."""
        assert gappy_ratings.check_my_balance() is False

    def test_crossing_simple(self, simple_ratings):
        """SIMPLE has no empties -- every rater rated every subject."""
        assert simple_ratings.check_my_crossing() is True

    def test_crossing_gappy(self, gappy_ratings):
        """GAPPY has empties."""
        assert gappy_ratings.check_my_crossing() is False

    def test_pairing_simple(self, simple_ratings):
        """SIMPLE is 2 raters, fully crossed -- a pair."""
        assert simple_ratings.check_my_pairing() is True

    def test_pairing_gappy(self, gappy_ratings):
        """GAPPY has 3 raters -- not a pair."""
        assert gappy_ratings.check_my_pairing() is False


# -- to_nij_matrix ---------------------------------------------------------- #


class TestToNijMatrix:
    """The counts matrix: how many raters chose each category per subject."""

    def test_simple_nij(self, simple_ratings):
        """SIMPLE nij matches hand computation from conftest."""
        nij = simple_ratings.to_nij_matrix()
        np.testing.assert_array_equal(nij.matrix, SIMPLE_NIJ)

    def test_gappy_nij(self, gappy_ratings):
        """GAPPY nij: empties contribute nothing to counts."""
        nij = gappy_ratings.to_nij_matrix()
        # subject0: ann=include(0), bob=exclude(1), cara=include(0) -> [2,1]
        # subject1: ann=exclude(1), bob=include(0), cara=EMPTY       -> [1,1]
        # subject2: ann=include(0), bob=exclude(1), cara=exclude(1)  -> [1,2]
        # subject3: ann=exclude(1), bob=EMPTY,      cara=include(0)  -> [1,1]
        expected = np.array(
            [[2, 1],
             [1, 1],
             [1, 2],
             [1, 1]],
            dtype=np.int64,
        )
        np.testing.assert_array_equal(nij.matrix, expected)

    def test_nij_row_sums_equal_ratings_per_subject(self, gappy_ratings):
        """Each row of nij sums to the number of raters who rated that subject."""
        nij = gappy_ratings.to_nij_matrix()
        # subject0: 3 raters, subject1: 2, subject2: 3, subject3: 2
        np.testing.assert_array_equal(nij.matrix.sum(axis=1), [3, 2, 3, 2])


# -- Equality ---------------------------------------------------------------- #


class TestEquality:
    """Two Ratings are equal when they hold the same ratings and mean the same thing."""

    def test_equal_instances(self, simple_semantics):
        """Same data, same semantics."""
        a = Ratings(SIMPLE_INDEX.copy(), simple_semantics)
        b = Ratings(SIMPLE_INDEX.copy(), simple_semantics)
        assert a == b

    def test_different_data(self, simple_semantics):
        """Different category_index is a different Ratings."""
        arr = SIMPLE_INDEX.copy()
        arr[2, 1] = 0  # change a disagreement to agreement
        a = Ratings(SIMPLE_INDEX.copy(), simple_semantics)
        b = Ratings(arr, simple_semantics)
        assert a != b

    def test_not_comparable(self, simple_ratings):
        """Comparing with a non-Ratings returns NotImplemented."""
        assert simple_ratings.__eq__("not ratings") is NotImplemented


# -- Error prefix ------------------------------------------------------------ #


class TestErrorPrefix:
    """Validation errors name the concrete class, not a parent."""

    def test_says_ratings(self, simple_semantics):
        """A Ratings failing validation says Ratings in the message."""
        arr = np.array([0, 1, 0], dtype=INDEX_DTYPE)  # 1-D
        with pytest.raises(ValueError, match=r"\[Ratings\]"):
            Ratings(arr, simple_semantics)
