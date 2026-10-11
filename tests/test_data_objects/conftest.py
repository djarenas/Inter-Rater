"""
Shared fixtures for data-object tests.

Two small datasets, built bottom-up from the dependency chain:

    Categories -> RatingsSemantics -> Ratings

The SIMPLE dataset is the smallest valid case: 3 subjects, 2 raters, K=2,
fully crossed (no empties). Small enough that every property -- agreement,
balance, nij counts -- can be verified by counting on fingers.

The GAPPY dataset adds empties: 4 subjects, 3 raters, K=2, with two cells
missing. It is unbalanced and not crossed, which is the ordinary shape of a
systematic review where each subject goes to a subset of the rater pool.

Both are frozen as module-level constants so that no test can mutate what
another test sees. Fixtures return fresh objects built from those constants.
"""

import numpy as np
import pytest

from interrater.contracts import EMPTY, INDEX_DTYPE
from interrater.semantics.categories import UnorderedCategories, OrderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.data_objects.ratings import Ratings
from interrater.data_objects.balanced_ratings import BalancedRatings
from interrater.data_objects.paired_ratings import PairedRatings


# -- Category names ---------------------------------------------------------- #

TWO_CAT_NAMES:   tuple[str, ...] = ("include", "exclude")
THREE_CAT_NAMES: tuple[str, ...] = ("mild", "moderate", "severe")


# -- SIMPLE dataset: 3 subjects, 2 raters, K=2, fully crossed --------------- #
#
#   subject  | rater0 | rater1
#   ---------+--------+--------
#   subject0 |   0    |   0       both "include"
#   subject1 |   1    |   1       both "exclude"
#   subject2 |   0    |   1       disagree

SIMPLE_INDEX = np.array(
    [[0, 0],
     [1, 1],
     [0, 1]],
    dtype=INDEX_DTYPE,
)

# Expected nij matrix (3 subjects x 2 categories):
#   subject0: include=2, exclude=0
#   subject1: include=0, exclude=2
#   subject2: include=1, exclude=1
SIMPLE_NIJ = np.array(
    [[2, 0],
     [0, 2],
     [1, 1]],
    dtype=np.int64,
)


# -- GAPPY dataset: 4 subjects, 3 raters, K=2, with empties ----------------- #
#
#   subject  | ann | bob | cara
#   ---------+-----+-----+------
#   subject0 |  0  |  1  |  0      all three rated
#   subject1 |  1  |  0  | -1      cara absent
#   subject2 |  0  |  1  |  1      all three rated
#   subject3 |  1  | -1  |  0      bob absent
#
#   Ratings per subject: [3, 2, 3, 2] -- unbalanced
#   Ratings per rater:   [4, 3, 3]    -- no idle raters

GAPPY_INDEX = np.array(
    [[ 0,  1,  0],
     [ 1,  0, EMPTY],
     [ 0,  1,  1],
     [ 1, EMPTY, 0]],
    dtype=INDEX_DTYPE,
)

GAPPY_RATERS: tuple[str, ...] = ("ann", "bob", "cara")


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


@pytest.fixture
def gappy_semantics(two_cat):
    """RatingsSemantics for the GAPPY dataset: 3 named raters, K=2."""
    return RatingsSemantics(raters=list(GAPPY_RATERS), categories=two_cat)


@pytest.fixture
def simple_ratings(simple_semantics):
    """Ratings built from the SIMPLE dataset: 3x2, fully crossed."""
    return Ratings(SIMPLE_INDEX.copy(), simple_semantics)


@pytest.fixture
def gappy_ratings(gappy_semantics):
    """Ratings built from the GAPPY dataset: 4x3, with empties."""
    return Ratings(GAPPY_INDEX.copy(), gappy_semantics)
