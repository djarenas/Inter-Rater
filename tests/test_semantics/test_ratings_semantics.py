"""
Tests for RatingsSemantics.

RatingsSemantics declares what a rating task looks like before the data is
seen: who the raters are, what the categories mean, and what the data was
expected to look like.

It sits between Categories (which it holds) and Ratings (which validates
against it). If a semantics can be built, a Ratings that matches it can be
built too -- and if a semantics rejects something, the data never gets that
far.

Key behaviours:
- Raters can be given as names or as an integer count (which generates
  rater0, rater1, ...).
- Type guards catch the traps Python's type system leaves open: a string is
  a Sequence, a bool is an int.
- Validation rejects < 2 raters, duplicate names, blank names, and bad
  n_subjects.
- from_dict / to_dict round-trips support YAML-based declarations.
"""

import pytest

from interrater.contracts import MIN_RATERS, MIN_SUBJECTS
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics

from .conftest import TWO_CAT_NAMES


# -- Construction ------------------------------------------------------------ #


class TestRaterConstruction:
    """How the raters argument is interpreted."""

    def test_from_names(self, two_cat):
        """A list of names is taken as given."""
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        assert sem.raters == ("ann", "bob")

    def test_from_count(self, two_cat):
        """An integer generates positional names."""
        sem = RatingsSemantics(raters=3, categories=two_cat)
        assert sem.raters == ("rater0", "rater1", "rater2")

    def test_n_raters_property(self, simple_semantics):
        """n_raters reflects the declared count."""
        assert simple_semantics.n_raters == 2

    def test_n_categories_property(self, simple_semantics):
        """n_categories delegates to the categories object."""
        assert simple_semantics.n_categories == 2

    def test_n_subjects_default_is_none(self, simple_semantics):
        """No subject expectation by default."""
        assert simple_semantics.n_subjects is None

    def test_n_subjects_stated(self, two_cat):
        """A stated n_subjects is stored."""
        sem = RatingsSemantics(raters=2, categories=two_cat, n_subjects=100)
        assert sem.n_subjects == 100

    def test_unique_subjects_default(self, simple_semantics):
        """Unique subjects is True by default."""
        assert simple_semantics.unique_subjects is True


# -- Type guards ------------------------------------------------------------- #


class TestTypeGuards:
    """Traps that Python's type system leaves open."""

    def test_string_raters_rejected(self, two_cat):
        """A string is a Sequence of characters -- never intended as raters."""
        with pytest.raises(TypeError, match="string"):
            RatingsSemantics(raters="ann", categories=two_cat)

    def test_bool_raters_rejected(self, two_cat):
        """A bool is an int in Python -- True would mean 1 rater."""
        with pytest.raises(TypeError, match="bool"):
            RatingsSemantics(raters=True, categories=two_cat)

    def test_non_string_rater_names(self, two_cat):
        """Rater names must be strings, not integers from a headerless file."""
        with pytest.raises(TypeError, match="strings"):
            RatingsSemantics(raters=[1, 2], categories=two_cat)

    def test_non_categories_object(self):
        """categories must be a Categories, not a list of names."""
        with pytest.raises(TypeError, match="Categories"):
            RatingsSemantics(raters=2, categories=["yes", "no"])

    def test_bool_n_subjects_rejected(self, two_cat):
        """A bool is an int -- True would mean 1 subject."""
        with pytest.raises(TypeError, match="integer or None"):
            RatingsSemantics(raters=2, categories=two_cat, n_subjects=True)

    def test_non_bool_unique_subjects(self, two_cat):
        """unique_subjects must be a boolean, not an int."""
        with pytest.raises(TypeError, match="boolean"):
            RatingsSemantics(raters=2, categories=two_cat, unique_subjects=1)


# -- Validation -------------------------------------------------------------- #


class TestValidation:
    """Value-level checks that run after types are confirmed."""

    def test_fewer_than_min_raters_by_count(self, two_cat):
        """An integer below MIN_RATERS is rejected."""
        with pytest.raises(ValueError, match="at least"):
            RatingsSemantics(raters=1, categories=two_cat)

    def test_fewer_than_min_raters_by_names(self, two_cat):
        """A single-element list is also rejected."""
        with pytest.raises(ValueError, match="at least"):
            RatingsSemantics(raters=["lonely"], categories=two_cat)

    def test_duplicate_rater_names(self, two_cat):
        """Two columns naming the same rater would be counted as independent."""
        with pytest.raises(ValueError, match="duplicate"):
            RatingsSemantics(raters=["ann", "bob", "ann"], categories=two_cat)

    def test_blank_rater_name(self, two_cat):
        """A blank name is usually a header row read as data."""
        with pytest.raises(ValueError, match="blank"):
            RatingsSemantics(raters=["ann", "  "], categories=two_cat)

    def test_n_subjects_below_minimum(self, two_cat):
        """n_subjects below MIN_SUBJECTS is rejected."""
        with pytest.raises(ValueError, match="at least"):
            RatingsSemantics(raters=2, categories=two_cat, n_subjects=1)


# -- Public API -------------------------------------------------------------- #


class TestPublicAPI:
    """Properties and methods beyond construction."""

    def test_index_of_rater(self, two_cat):
        """index_of_rater returns the column position."""
        sem = RatingsSemantics(raters=["ann", "bob", "cara"], categories=two_cat)
        assert sem.index_of_rater("bob") == 1

    def test_index_of_rater_unknown(self, simple_semantics):
        """Looking up a rater that does not exist raises KeyError."""
        with pytest.raises(KeyError, match="no rater"):
            simple_semantics.index_of_rater("nobody")

    def test_to_dict_round_trip(self, two_cat):
        """to_dict produces a mapping from_dict can rebuild from."""
        sem = RatingsSemantics(
            raters=["ann", "bob"],
            categories=two_cat,
            n_subjects=50,
            unique_subjects=False,
        )
        spec = sem.to_dict()
        rebuilt = RatingsSemantics.from_dict(spec)
        assert rebuilt == sem

    def test_from_dict_missing_raters(self):
        """A declaration without raters is incomplete."""
        with pytest.raises(KeyError, match="missing"):
            RatingsSemantics.from_dict({
                "categories": {"kind": "unordered", "names": list(TWO_CAT_NAMES)},
            })

    def test_from_dict_missing_categories(self):
        """A declaration without categories is incomplete."""
        with pytest.raises(KeyError, match="missing"):
            RatingsSemantics.from_dict({"raters": 2})


# -- Equality ---------------------------------------------------------------- #


class TestEquality:
    """Two semantics are equal when they declare the same rating task."""

    def test_equal_instances(self, two_cat):
        """Same raters, same categories, same expectations."""
        a = RatingsSemantics(raters=2, categories=two_cat)
        b = RatingsSemantics(raters=2, categories=two_cat)
        assert a == b

    def test_different_raters(self, two_cat):
        """Different raters is a different rating task."""
        a = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        b = RatingsSemantics(raters=["cara", "dev"], categories=two_cat)
        assert a != b

    def test_different_n_subjects(self, two_cat):
        """Different expectations is a different declaration."""
        a = RatingsSemantics(raters=2, categories=two_cat, n_subjects=10)
        b = RatingsSemantics(raters=2, categories=two_cat, n_subjects=20)
        assert a != b

    def test_not_comparable(self, simple_semantics):
        """Comparing with a non-semantics returns NotImplemented."""
        assert simple_semantics.__eq__("not semantics") is NotImplemented
