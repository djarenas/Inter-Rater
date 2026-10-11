"""
Shared fixtures for report tests.

Reports are frozen dataclasses that record what a cleaner did. They hold
counts, not data, so every fixture here is a small hand-built instance
whose properties can be verified by reading it.

Two MatchingReport instances:
- MATCHING_CLEAN: 8 cells, nothing wrong.
- MATCHING_MESSY: 8 cells with a null, a missing marker, a repair, and
  an unknown -- one of each outcome.

Two DropReport instances:
- DROPS_NOTHING: raters_first, nothing dropped.
- DROPS_SOME: raters_first, one rater and two subjects dropped in one pass.

Two CategoryReport instances:
- CATEGORIES_ALL_USED: K=2, both used.
- CATEGORIES_WITH_UNUSED: K=3, one unused and kept.

Two CleaningReport instances (built from the above):
- REPORT_STABLE: clean matching, no drops, all used -- semantics unchanged.
- REPORT_CHANGED: messy matching, drops, unused -- semantics changed
  because a rater was dropped.
"""

import pytest

from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.cleaners.cleaner_config import CleanerConfig
from interrater.reports.matching_report import MatchingReport
from interrater.reports.drop_report import DropPass, DropReport
from interrater.reports.category_report import CategoryReport
from interrater.reports.cleaning_report import CleaningReport


# -- MatchingReport ---------------------------------------------------------- #

MATCHING_CLEAN = MatchingReport(n_cells=8, n_null=0)

MATCHING_MESSY = MatchingReport(
    n_cells=8,
    n_null=1,
    missing_matched={"NA": 1},
    repairs={"inc": 1},
    unknown_values={"???": 1},
    per_rater_null={"bob": 1},
    per_rater_missing={"ann": 1, "bob": 1},
    per_rater_unknown={"ann": 1},
    subject_missing_distribution={0: 2, 1: 1, 2: 1},
)


# -- DropReport -------------------------------------------------------------- #

DROPS_NOTHING = DropReport(method="raters_first")

DROPS_SOME = DropReport(
    method="raters_first",
    converged=True,
    passes=(DropPass(1, 1, 2),),
    dropped_raters={"carol": 1},
    dropped_subject_ids=("doc_E", "doc_F"),
    dropped_subject_distribution={1: 2},
)


# -- CategoryReport --------------------------------------------------------- #

CATEGORIES_ALL_USED = CategoryReport(
    counts={"include": 100, "exclude": 100},
    policy="error",
)

CATEGORIES_WITH_UNUSED = CategoryReport(
    counts={"include": 4, "exclude": 4, "maybe": 0},
    policy="keep",
    action="kept",
    unused=("maybe",),
)


# -- Semantics --------------------------------------------------------------- #

TWO_CAT = UnorderedCategories(("include", "exclude"))
THREE_CAT = UnorderedCategories(("include", "exclude", "maybe"))

SEMANTICS_3R = RatingsSemantics(
    raters=["ann", "bob", "carol"], categories=TWO_CAT,
)
SEMANTICS_2R = RatingsSemantics(
    raters=["ann", "bob"], categories=TWO_CAT,
)


# -- Fixtures ---------------------------------------------------------------- #


@pytest.fixture
def matching_clean():
    return MATCHING_CLEAN


@pytest.fixture
def matching_messy():
    return MATCHING_MESSY


@pytest.fixture
def drops_nothing():
    return DROPS_NOTHING


@pytest.fixture
def drops_some():
    return DROPS_SOME


@pytest.fixture
def categories_all_used():
    return CATEGORIES_ALL_USED


@pytest.fixture
def categories_with_unused():
    return CATEGORIES_WITH_UNUSED


@pytest.fixture
def report_stable():
    """CleaningReport where nothing changed: semantics in == semantics out."""
    return CleaningReport(
        config=CleanerConfig(),
        declared_semantics=SEMANTICS_2R,
        effective_semantics=SEMANTICS_2R,
        shape_in=(4, 2),
        shape_out=(4, 2),
        matching=MATCHING_CLEAN,
        drops=DROPS_NOTHING,
        categories=CATEGORIES_ALL_USED,
    )


@pytest.fixture
def report_changed():
    """CleaningReport where a rater was dropped: semantics changed."""
    return CleaningReport(
        config=CleanerConfig(),
        declared_semantics=SEMANTICS_3R,
        effective_semantics=SEMANTICS_2R,
        shape_in=(6, 3),
        shape_out=(4, 2),
        matching=MATCHING_MESSY,
        drops=DROPS_SOME,
        categories=CATEGORIES_ALL_USED,
    )
