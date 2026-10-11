"""
Shared fixtures for cleaner tests.

Three things a cleaner test always needs: categories to declare, semantics
to describe the rating task, and a DataFrame to clean. The fixtures here
are the smallest instances of each that let a test focus on its own concern
rather than on setup.

CLEAN_DF is 4 subjects x 2 raters ("ann", "bob"), K=2 ("include", "exclude"),
no nulls, no repairs, no unknowns -- the happy path in four rows.

MESSY_DF is the same shape but with a null, a missing marker ("NA"), a
repairable typo ("inc"), and an unknown ("???"), so encoding tests can see
each outcome in one frame.
"""

import numpy as np
import pandas as pd
import pytest

from interrater.semantics.categories import (
    UnorderedCategories,
    OrderedCategories,
    WeightedCategories,
)
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.cleaners.cleaner_config import CleanerConfig
from interrater.cleaners.frame_layout import FrameLayout
from interrater.cleaners.insufficient_data_policy import InsufficientDataPolicy


# -- Category names ---------------------------------------------------------- #

CAT_NAMES: tuple[str, ...] = ("include", "exclude")


# -- Clean DataFrame: 4 subjects, 2 raters, no problems --------------------- #
#
#   subject  | ann     | bob
#   ---------+---------+---------
#   doc_A    | include | exclude
#   doc_B    | exclude | exclude
#   doc_C    | include | include
#   doc_D    | exclude | include

CLEAN_DF = pd.DataFrame(
    {
        "ann": ["include", "exclude", "include", "exclude"],
        "bob": ["exclude", "exclude", "include", "include"],
    },
    index=["doc_A", "doc_B", "doc_C", "doc_D"],
)


# -- Messy DataFrame: same shape, one of each problem ----------------------- #
#
#   subject  | ann     | bob
#   ---------+---------+----------
#   doc_A    | include | exclude     clean
#   doc_B    | NA      | exclude     missing marker in ann
#   doc_C    | inc     | include     repairable typo in ann
#   doc_D    | exclude | NaN         null in bob

MESSY_DF = pd.DataFrame(
    {
        "ann": ["include",  "NA",  "inc",     "exclude"],
        "bob": ["exclude",  "exclude", "include", None],
    },
)


# -- Fixtures ---------------------------------------------------------------- #


@pytest.fixture
def two_cat():
    """UnorderedCategories with K=2: include, exclude."""
    return UnorderedCategories(CAT_NAMES)


@pytest.fixture
def semantics(two_cat):
    """RatingsSemantics for 2 named raters (ann, bob), K=2."""
    return RatingsSemantics(raters=["ann", "bob"], categories=two_cat)


@pytest.fixture
def clean_df():
    """A 4x2 DataFrame with no problems at all."""
    return CLEAN_DF.copy()


@pytest.fixture
def messy_df():
    """A 4x2 DataFrame with a null, a missing marker, a typo, and an unknown."""
    return MESSY_DF.copy()


@pytest.fixture
def default_config():
    """CleanerConfig with all defaults."""
    return CleanerConfig()
