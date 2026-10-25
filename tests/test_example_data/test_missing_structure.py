"""
The shuffles keep the missing-data structure of a real partially crossed study.

SR1 main is the case this matters for: 1617 abstracts, 4 reviewers, half of
all cells empty, each abstract read by 2 of the 4 (3 by three). Who was shown
what is the design of the study, so every shuffle must leave the empty cells
exactly where they are: no answer may land on an abstract its reviewer was
never shown, and every reviewer keeps their workload.

Each test also checks what that shuffle promises to keep, and that it does
move answers. The small-fixture versions of these checks live in
tests/test_analyzers/test_resampling_functions.py; these run the same
promises on the real design.
"""

import numpy as np
import pytest

from interrater.analyzers.resampling_functions import (
    pool_and_redistribute,
    within_rater_shuffle,
    within_subject_shuffle,
)

N_DRAWS = 200


@pytest.fixture(scope="module")
def sr1_main_codes(sr1_main_ratings):
    """SR1 main as a raw code array: -1 where a reviewer was not shown the abstract."""
    return sr1_main_ratings.category_index


def _per_rater_counts(codes: np.ndarray, n_categories: int) -> list[np.ndarray]:
    return [
        np.bincount(codes[codes[:, j] >= 0, j], minlength=n_categories)
        for j in range(codes.shape[1])
    ]


class TestSR1MainDesign:
    """The fixture is the partially crossed design the tests below rely on."""

    def test_half_the_cells_are_empty(self, sr1_main_codes):
        assert sr1_main_codes.shape == (1617, 4)
        assert (sr1_main_codes < 0).sum() == 3231

    def test_each_abstract_read_by_two_or_three(self, sr1_main_codes):
        per_abstract = (sr1_main_codes >= 0).sum(axis=1)
        assert set(per_abstract.tolist()) == {2, 3}


class TestShufflesKeepTheDesign:
    """Every draw keeps the empty cells, and each shuffle keeps what it promises."""

    def test_within_rater_shuffle(self, sr1_main_codes):
        """Empty cells stay put; each reviewer keeps their own answer mix."""
        empty = sr1_main_codes < 0
        k = int(sr1_main_codes.max()) + 1
        expected = _per_rater_counts(sr1_main_codes, k)
        moved = False
        for draw in within_rater_shuffle(sr1_main_codes, n_permutations=N_DRAWS):
            np.testing.assert_array_equal(draw < 0, empty)
            for got, want in zip(_per_rater_counts(draw, k), expected):
                np.testing.assert_array_equal(got, want)
            moved |= not np.array_equal(draw, sr1_main_codes)
        assert moved

    def test_within_subject_shuffle(self, sr1_main_codes):
        """Empty cells stay put; each abstract keeps the answers it got."""
        empty = sr1_main_codes < 0
        expected = np.sort(sr1_main_codes, axis=1)
        moved = False
        for draw in within_subject_shuffle(sr1_main_codes, n_permutations=N_DRAWS):
            np.testing.assert_array_equal(draw < 0, empty)
            np.testing.assert_array_equal(np.sort(draw, axis=1), expected)
            moved |= not np.array_equal(draw, sr1_main_codes)
        assert moved

    def test_pool_and_redistribute(self, sr1_main_codes):
        """Empty cells stay put; the overall mix of answers is unchanged."""
        empty = sr1_main_codes < 0
        k = int(sr1_main_codes.max()) + 1
        expected = np.bincount(sr1_main_codes[~empty], minlength=k)
        moved = False
        for draw in pool_and_redistribute(sr1_main_codes, n_permutations=N_DRAWS):
            np.testing.assert_array_equal(draw < 0, empty)
            np.testing.assert_array_equal(np.bincount(draw[~empty], minlength=k), expected)
            moved |= not np.array_equal(draw, sr1_main_codes)
        assert moved
