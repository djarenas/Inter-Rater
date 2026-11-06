"""
The cleaner on deliberately damaged ratings, checked against the damage.

``fixtures/scruffy_ratings_generator.py`` takes clean ratings and breaks them
the ways real spreadsheets break -- wrong case, stray spaces, shorthand,
undeclared answers, six spellings of nothing, a rater who quit, subjects
with one rating, a category nobody chose, repeated ids -- and records every
cell it broke. ``fixtures/scruffy_cleaner_config.yaml`` is a config
permissive enough to clean it.

So the answer is known before the cleaner runs. Each test compares what the
cleaner did with what the manifest says it should have done, then the last
one takes the cleaned ratings all the way to a coefficient: the path a user
takes from a messy file to a number.
"""

from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from interrater.analyzers import KrippendorffAlphaAnalyzer
from interrater.analyzers.configs import KrippendorffAlphaConfig
from interrater.cleaners import CleanerConfig, RatingsCleaner
from interrater.contracts import EMPTY
from interrater.data_objects.ratings import Ratings
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics

from fixtures.scruffy_ratings_generator import (
    CASE,
    CATEGORIES,
    MISSING,
    NULL,
    RATERS,
    REPAIRABLE,
    SUBJECT_COL,
    UNKNOWN,
    WHITESPACE,
    make_scruffy_ratings,
)

CONFIG_FILE = Path(__file__).resolve().parents[1] / "fixtures" / "scruffy_cleaner_config.yaml"

# Damage the cleaner should undo, and damage that leaves no rating
RECOVERED = (CASE, WHITESPACE, REPAIRABLE)
ABSENT = (UNKNOWN, MISSING, NULL)


def _semantics(unique_subjects: bool = False) -> RatingsSemantics:
    """The rating task. The generator repeats two ids, so not unique by default."""
    return RatingsSemantics(
        raters=list(RATERS),
        categories=UnorderedCategories(CATEGORIES),
        unique_subjects=unique_subjects,
    )


@pytest.fixture(scope="module")
def scruffy():
    return make_scruffy_ratings(seed=0)


@pytest.fixture(scope="module")
def cleaned(scruffy):
    config = CleanerConfig.build_from_yaml(CONFIG_FILE)
    return RatingsCleaner(scruffy.frame, _semantics(), config)


@pytest.fixture(scope="module")
def expected(scruffy):
    """
    The cleaned data, worked out from the manifest alone.

    Every cell starts as its clean answer. Damage the cleaner can undo leaves
    it; damage that leaves no rating empties it. Then the rater who quit and
    the subjects with one rating go, as the config's minimums say.
    """
    answers = scruffy.clean_frame[list(RATERS)].to_numpy()
    codes = np.array([[CATEGORIES.index(a) for a in row] for row in answers])
    for cell in scruffy.damage:
        if cell.kind in ABSENT:
            codes[cell.row, RATERS.index(cell.rater)] = EMPTY

    kept_rows = [r for r in range(len(codes)) if r not in scruffy.thin_subjects]
    kept_raters = [i for i, rater in enumerate(RATERS) if rater != scruffy.quit_rater]
    return {
        "codes": codes[np.ix_(kept_rows, kept_raters)],
        "raters": [RATERS[i] for i in kept_raters],
        "subject_ids": [scruffy.frame.at[r, SUBJECT_COL] for r in kept_rows],
    }


class TestMatching:
    """Every damaged cell is counted where the manifest says it belongs."""

    def test_counts_match_the_manifest(self, scruffy, cleaned):
        matching = cleaned.report.matching
        assert matching.n_null == scruffy.count(NULL)
        assert matching.n_matched_missing == scruffy.count(MISSING)
        assert matching.n_unknown == scruffy.count(UNKNOWN)
        assert matching.n_repaired == scruffy.count(REPAIRABLE)

    def test_unknown_values_match_the_manifest(self, scruffy, cleaned):
        written = Counter(cell.written for cell in scruffy.cells(UNKNOWN))
        assert cleaned.report.matching.unknown_values == dict(written)

    def test_repairs_match_the_manifest(self, scruffy, cleaned):
        written = Counter(cell.written for cell in scruffy.cells(REPAIRABLE))
        repairs = {code: n for code, n in cleaned.report.matching.repairs.items() if n}
        assert repairs == dict(written)


class TestDrops:
    """The rater who quit and the subjects with one rating are what goes."""

    def test_the_rater_who_quit_is_dropped(self, scruffy, cleaned):
        assert list(cleaned.report.drops.dropped_raters) == [scruffy.quit_rater]

    def test_the_thin_subjects_are_dropped(self, scruffy, cleaned):
        thin_ids = tuple(scruffy.frame.at[r, SUBJECT_COL] for r in scruffy.thin_subjects)
        assert cleaned.report.drops.dropped_subject_ids == thin_ids


class TestCleanedData:
    """The ratings that come out are the original answers, cell for cell."""

    def test_every_cell_is_the_original_answer(self, cleaned, expected):
        np.testing.assert_array_equal(cleaned.ratings.category_index, expected["codes"])

    def test_raters_and_ids_survive(self, cleaned, expected):
        assert [str(r) for r in cleaned.ratings.semantics.raters] == expected["raters"]
        assert list(cleaned.ratings.subject_ids) == expected["subject_ids"]

    def test_unused_category_is_kept(self, scruffy, cleaned):
        """Nobody chose it, but it was offered, so K stays at 3."""
        categories = cleaned.report.categories
        assert categories.unused == (scruffy.unused_category,)
        assert categories.action == "kept"
        assert cleaned.ratings.n_categories == len(CATEGORIES)


class TestWhatIsFatalIsTheConfig:
    """The same file, cleaned strictly, stops instead of guessing."""

    def test_strict_defaults_stop_on_unknown_answers(self, scruffy):
        config = CleanerConfig.build_from_yaml(CONFIG_FILE)
        strict = CleanerConfig(frame_layout=config.frame_layout)
        with pytest.raises(ValueError, match="neither a declared category nor a marker for absence"):
            RatingsCleaner(scruffy.frame, _semantics(), strict)

    def test_repeated_ids_stop_a_unique_rating_task(self, scruffy):
        config = CleanerConfig.build_from_yaml(CONFIG_FILE)
        with pytest.raises(ValueError, match="appearing more than once"):
            RatingsCleaner(scruffy.frame, _semantics(unique_subjects=True), config)


def test_messy_file_to_a_coefficient(cleaned, expected):
    """
    The whole user path: messy table, cleaner, analyzer, number.

    The reference Ratings is built by hand from the manifest's answer, so the
    coefficient is checked against data the cleaner never touched.
    """
    config = KrippendorffAlphaConfig(n_resamples=50, random_state=42)
    result = KrippendorffAlphaAnalyzer(cleaned.ratings, config).estimate()

    reference = Ratings(
        expected["codes"],
        RatingsSemantics(
            raters=expected["raters"],
            categories=UnorderedCategories(CATEGORIES),
            unique_subjects=False,
        ),
        subject_ids=expected["subject_ids"],
    )
    assert np.isfinite(result.point)
    assert result.point == KrippendorffAlphaAnalyzer(reference, config).estimate().point
