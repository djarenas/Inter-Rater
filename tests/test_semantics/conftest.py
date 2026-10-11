"""
Shared fixtures for semantics tests.

Categories and RatingsSemantics are the foundation of the data-object chain.
These fixtures provide the building blocks that test_categories.py and
test_ratings_semantics.py need -- nothing more.
"""

import pytest

from interrater.semantics.categories import UnorderedCategories, OrderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics


# -- Category names ---------------------------------------------------------- #

TWO_CAT_NAMES:   tuple[str, ...] = ("include", "exclude")
THREE_CAT_NAMES: tuple[str, ...] = ("mild", "moderate", "severe")


# -- Fixtures ---------------------------------------------------------------- #


@pytest.fixture
def two_cat():
    """UnorderedCategories with K=2: include, exclude."""
    return UnorderedCategories(TWO_CAT_NAMES)


@pytest.fixture
def three_cat_ordered():
    """OrderedCategories with K=3: mild, moderate, severe (quadratic)."""
    return OrderedCategories(THREE_CAT_NAMES, weighting="quadratic")


@pytest.fixture
def simple_semantics(two_cat):
    """RatingsSemantics for the SIMPLE dataset: 2 raters, K=2."""
    return RatingsSemantics(raters=2, categories=two_cat)
