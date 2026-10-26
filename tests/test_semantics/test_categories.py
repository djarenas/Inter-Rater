"""
Tests for the Categories hierarchy.

Categories own K -- the number of responses a rater was allowed to give.
They are the foundation of the data-object chain: every RatingsSemantics
holds one, and every Ratings validates against it.

Three subclasses:
- UnorderedCategories: nominal, no distances, every disagreement equally wrong.
- OrderedCategories: ranked scale, distances derived from position.
- WeightedCategories: distances supplied as an explicit matrix.

The base class (Categories) is abstract; its shared validation is tested
through UnorderedCategories, the simplest concrete subclass.

The build_categories factory dispatches a plain dict to the right subclass,
which is the path YAML declarations take.
"""

import numpy as np
import pytest

from interrater.semantics.categories import (
    Categories,
    UnorderedCategories,
    OrderedCategories,
    WeightedCategories,
    build_categories,
)

from .conftest import TWO_CAT_NAMES, THREE_CAT_NAMES


# -- Base validation (tested through UnorderedCategories) -------------------- #


class TestCategoriesBaseValidation:
    """Shared validation inherited by every subclass."""

    def test_fewer_than_two_names(self):
        """A single category cannot produce agreement or disagreement."""
        with pytest.raises(ValueError, match="at least 2"):
            UnorderedCategories(("only_one",))

    def test_empty_names(self):
        """Zero categories is not a rating task."""
        with pytest.raises(ValueError, match="at least 2"):
            UnorderedCategories(())

    def test_127_categories_accepted(self):
        """The most a category position stored as int8 can hold."""
        cats = UnorderedCategories([f"c{i}" for i in range(127)])
        assert cats.n_categories == 127

    def test_more_than_127_categories_refused(self):
        """Position 128 and up would not fit the int8 the ratings use."""
        with pytest.raises(ValueError, match="at most 127 categories, got 128"):
            UnorderedCategories([f"c{i}" for i in range(128)])

    def test_cap_applies_to_every_kind(self):
        with pytest.raises(ValueError, match="at most 127"):
            OrderedCategories([f"c{i}" for i in range(300)], weighting="linear")

    def test_duplicate_names(self):
        """Two codes meaning the same thing is never intentional."""
        with pytest.raises(ValueError, match="duplicate"):
            UnorderedCategories(("yes", "no", "yes"))

    def test_labels_with_unknown_category(self):
        """A label for a category that does not exist is a typo."""
        with pytest.raises(ValueError, match="not in this set"):
            UnorderedCategories(
                TWO_CAT_NAMES, labels={"nonexistent": "display name"}
            )


# -- UnorderedCategories ----------------------------------------------------- #


class TestUnorderedCategories:
    """Nominal categories with no distances."""

    def test_construction(self, two_cat):
        """Valid construction stores names as a tuple."""
        assert two_cat.names == TWO_CAT_NAMES

    def test_n_categories(self, two_cat):
        """K is the length of names."""
        assert two_cat.n_categories == 2

    def test_index_of(self, two_cat):
        """Position in names is the integer code."""
        assert two_cat.index_of("include") == 0
        assert two_cat.index_of("exclude") == 1

    def test_index_of_unknown(self, two_cat):
        """Looking up a category that does not exist raises KeyError."""
        with pytest.raises(KeyError, match="no category"):
            two_cat.index_of("maybe")

    def test_label_for_with_label(self):
        """Labels override the category name for display."""
        cats = UnorderedCategories(TWO_CAT_NAMES, labels={"include": "Inc"})
        assert cats.label_for("include") == "Inc"

    def test_label_for_without_label(self, two_cat):
        """Without a label, the category name is its own display name."""
        assert two_cat.label_for("include") == "include"

    def test_to_dict_round_trip(self, two_cat):
        """to_dict produces a mapping build_categories can rebuild from."""
        spec = two_cat.to_dict()
        rebuilt = build_categories(spec)
        assert rebuilt == two_cat

    def test_equality_same(self):
        """Two instances with the same names are equal."""
        a = UnorderedCategories(TWO_CAT_NAMES)
        b = UnorderedCategories(TWO_CAT_NAMES)
        assert a == b

    def test_equality_different_names(self):
        """Different names means different categories."""
        a = UnorderedCategories(("yes", "no"))
        b = UnorderedCategories(("include", "exclude"))
        assert a != b

    def test_equality_different_type(self):
        """Unordered != Ordered, even with the same names."""
        u = UnorderedCategories(THREE_CAT_NAMES)
        o = OrderedCategories(THREE_CAT_NAMES)
        assert u != o

    def test_equality_not_comparable(self):
        """Comparing with a non-Categories returns NotImplemented."""
        assert UnorderedCategories(TWO_CAT_NAMES).__eq__("not categories") is NotImplemented


# -- OrderedCategories ------------------------------------------------------- #


class TestOrderedCategories:
    """Ranked scale with distances derived from position."""

    def test_construction_quadratic(self, three_cat_ordered):
        """Quadratic weighting is the default in the medical literature."""
        assert three_cat_ordered.names == THREE_CAT_NAMES
        assert three_cat_ordered.weighting == "quadratic"

    def test_construction_linear(self):
        """Linear weighting penalises distance evenly."""
        cats = OrderedCategories(THREE_CAT_NAMES, weighting="linear")
        assert cats.weighting == "linear"

    def test_invalid_weighting(self):
        """Only linear and quadratic are allowed."""
        with pytest.raises(ValueError, match="must be one of"):
            OrderedCategories(THREE_CAT_NAMES, weighting="cubic")

    def test_weight_matrix_quadratic(self, three_cat_ordered):
        """Quadratic weights for 3 categories: 1 - d^2/d_max^2."""
        w = three_cat_ordered.weight_matrix()
        # K=3, distances: 0,1,2. d_max=2. d^2: 0,1,4. w = 1 - d^2/4
        expected = np.array([
            [1.0,  0.75, 0.0],
            [0.75, 1.0,  0.75],
            [0.0,  0.75, 1.0],
        ])
        np.testing.assert_allclose(w, expected)

    def test_weight_matrix_linear(self):
        """Linear weights for 3 categories: 1 - d/d_max."""
        cats = OrderedCategories(THREE_CAT_NAMES, weighting="linear")
        w = cats.weight_matrix()
        expected = np.array([
            [1.0, 0.5, 0.0],
            [0.5, 1.0, 0.5],
            [0.0, 0.5, 1.0],
        ])
        np.testing.assert_allclose(w, expected)

    def test_weight_matrix_diagonal_is_one(self, three_cat_ordered):
        """A category agrees perfectly with itself."""
        w = three_cat_ordered.weight_matrix()
        np.testing.assert_array_equal(np.diag(w), [1.0, 1.0, 1.0])

    def test_weight_matrix_symmetric(self, three_cat_ordered):
        """How wrong a disagreement is does not depend on who said which."""
        w = three_cat_ordered.weight_matrix()
        np.testing.assert_array_equal(w, w.T)

    def test_to_dict_includes_weighting(self, three_cat_ordered):
        """The weighting survives a round trip."""
        spec = three_cat_ordered.to_dict()
        assert spec["weighting"] == "quadratic"
        rebuilt = build_categories(spec)
        assert rebuilt == three_cat_ordered

    def test_equality_different_weighting(self):
        """Same names, different weighting -- not equal."""
        lin = OrderedCategories(THREE_CAT_NAMES, weighting="linear")
        quad = OrderedCategories(THREE_CAT_NAMES, weighting="quadratic")
        assert lin != quad


# -- WeightedCategories ------------------------------------------------------ #


class TestWeightedCategories:
    """Explicit distance matrix, no ranking required."""

    @pytest.fixture
    def valid_weights(self):
        """A 2x2 weight matrix for testing."""
        return np.array([[1.0, 0.3], [0.3, 1.0]])

    def test_construction(self, valid_weights):
        """Valid construction stores the weight matrix."""
        cats = WeightedCategories(TWO_CAT_NAMES, valid_weights)
        assert cats.names == TWO_CAT_NAMES
        np.testing.assert_array_equal(cats.weight_matrix(), valid_weights)

    def test_wrong_shape(self):
        """Weight matrix must be (K, K)."""
        bad = np.array([[1.0, 0.3, 0.0], [0.3, 1.0, 0.5]])
        with pytest.raises(ValueError, match="must be"):
            WeightedCategories(TWO_CAT_NAMES, bad)

    def test_not_symmetric(self):
        """Disagreement cost cannot depend on who said which."""
        bad = np.array([[1.0, 0.3], [0.5, 1.0]])
        with pytest.raises(ValueError, match="symmetric"):
            WeightedCategories(TWO_CAT_NAMES, bad)

    def test_diagonal_not_one(self):
        """A category must agree perfectly with itself."""
        bad = np.array([[0.9, 0.3], [0.3, 1.0]])
        with pytest.raises(ValueError, match="diagonal"):
            WeightedCategories(TWO_CAT_NAMES, bad)

    def test_values_out_of_range(self):
        """Weights must lie in [0, 1]."""
        bad = np.array([[1.0, 1.5], [1.5, 1.0]])
        with pytest.raises(ValueError, match=r"\[0, 1\]"):
            WeightedCategories(TWO_CAT_NAMES, bad)

    def test_to_dict_round_trip(self, valid_weights):
        """The weight matrix survives a round trip through to_dict."""
        cats = WeightedCategories(TWO_CAT_NAMES, valid_weights)
        spec = cats.to_dict()
        rebuilt = build_categories(spec)
        assert rebuilt == cats

    def test_equality_different_weights(self):
        """Same names, different matrix -- not equal."""
        w1 = np.array([[1.0, 0.3], [0.3, 1.0]])
        w2 = np.array([[1.0, 0.7], [0.7, 1.0]])
        assert WeightedCategories(TWO_CAT_NAMES, w1) != WeightedCategories(TWO_CAT_NAMES, w2)


# -- build_categories factory ------------------------------------------------ #


class TestBuildCategories:
    """The YAML dispatch: a plain dict becomes the right subclass."""

    def test_dispatches_unordered(self):
        """kind='unordered' builds UnorderedCategories."""
        result = build_categories({"kind": "unordered", "names": list(TWO_CAT_NAMES)})
        assert isinstance(result, UnorderedCategories)

    def test_dispatches_ordered(self):
        """kind='ordered' builds OrderedCategories."""
        result = build_categories({
            "kind": "ordered", "names": list(THREE_CAT_NAMES), "weighting": "linear",
        })
        assert isinstance(result, OrderedCategories)

    def test_dispatches_weighted(self):
        """kind='weighted' builds WeightedCategories."""
        result = build_categories({
            "kind": "weighted",
            "names": list(TWO_CAT_NAMES),
            "weights": [[1.0, 0.3], [0.3, 1.0]],
        })
        assert isinstance(result, WeightedCategories)

    def test_invalid_kind(self):
        """An unknown kind is rejected."""
        with pytest.raises(ValueError, match="must be one of"):
            build_categories({"kind": "fuzzy", "names": ["a", "b"]})

    def test_missing_names(self):
        """Every kind needs names."""
        with pytest.raises(KeyError, match="names"):
            build_categories({"kind": "unordered"})
