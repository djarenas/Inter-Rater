"""Tests for the resampling_functions sub-package.

One test file covering all 8 exported functions across the three modules:
  - subject_resampling: subject_bootstrap, subject_jackknife,
    subject_bootstrap_from_nij_matrix
  - rater_resampling: rater_jackknife, rater_subset
  - permutation_resampling: column_shuffle, within_subject_shuffle,
    pool_and_redistribute

All functions are generators. Validation runs on first next() call,
not at function-call time.
"""

from __future__ import annotations

from math import comb

import numpy as np
import pytest

from interrater.analyzers.resampling_functions.subject_resampling import (
    subject_bootstrap,
    subject_bootstrap_from_nij_matrix,
    subject_jackknife,
)
from interrater.analyzers.resampling_functions.rater_resampling import (
    rater_jackknife,
    rater_subset,
)
from interrater.analyzers.resampling_functions.permutation_resampling import (
    column_shuffle,
    pool_and_redistribute,
    within_subject_shuffle,
)


# ------------------------------------------------------------------ #
# Fixtures                                                             #
# ------------------------------------------------------------------ #


@pytest.fixture
def ratings_3x2():
    """3 subjects, 2 raters, no missing data."""
    return np.array([[0, 1], [1, 0], [1, 1]], dtype=np.int8)


@pytest.fixture
def ratings_4x3():
    """4 subjects, 3 raters, no missing data."""
    return np.array(
        [[0, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 0]], dtype=np.int8
    )


@pytest.fixture
def ratings_with_missing():
    """4 subjects, 3 raters, some missing (-1)."""
    return np.array(
        [[0, 1, -1], [1, -1, 0], [-1, 0, 1], [1, 1, 1]], dtype=np.int8
    )


@pytest.fixture
def nij_matrix_4x3():
    """4 subjects, 3 categories."""
    return np.array(
        [[2, 1, 0], [0, 2, 1], [1, 0, 2], [1, 1, 1]], dtype=np.int64
    )


# ================================================================== #
#                                                                      #
#   subject_bootstrap                                                  #
#                                                                      #
# ================================================================== #


class TestSubjectBootstrapValidation:
    """Type and value checks for subject_bootstrap."""

    def test_category_index_not_ndarray(self):
        with pytest.raises(TypeError, match="category_index must be np.ndarray"):
            list(subject_bootstrap([[0, 1]], n_resamples=1))

    def test_n_resamples_not_int(self, ratings_3x2):
        with pytest.raises(TypeError, match="n_resamples must be int"):
            list(subject_bootstrap(ratings_3x2, n_resamples=5.0))

    def test_random_state_not_int_or_none(self, ratings_3x2):
        with pytest.raises(TypeError, match="random_state must be int or None"):
            list(subject_bootstrap(ratings_3x2, random_state="bad"))

    def test_category_index_not_2d(self):
        with pytest.raises(ValueError, match="category_index must be 2D"):
            list(subject_bootstrap(np.array([0, 1, 2]), n_resamples=1))

    def test_n_resamples_zero(self, ratings_3x2):
        with pytest.raises(ValueError, match="n_resamples must be >= 1"):
            list(subject_bootstrap(ratings_3x2, n_resamples=0))

    def test_n_resamples_negative(self, ratings_3x2):
        with pytest.raises(ValueError, match="n_resamples must be >= 1"):
            list(subject_bootstrap(ratings_3x2, n_resamples=-3))


class TestSubjectBootstrapBehavior:
    """Functional tests for subject_bootstrap."""

    def test_yields_correct_count(self, ratings_3x2):
        draws = list(subject_bootstrap(ratings_3x2, n_resamples=7))
        assert len(draws) == 7

    def test_each_draw_has_correct_shape(self, ratings_3x2):
        for draw in subject_bootstrap(ratings_3x2, n_resamples=5):
            assert draw.shape == (3, 2)

    def test_draws_contain_original_rows(self, ratings_3x2):
        """Each row in a draw must be a row from the original."""
        original_set = {tuple(row) for row in ratings_3x2}
        for draw in subject_bootstrap(ratings_3x2, n_resamples=10):
            for row in draw:
                assert tuple(row) in original_set

    def test_reproducibility_same_seed(self, ratings_3x2):
        draws_a = list(subject_bootstrap(ratings_3x2, n_resamples=5, random_state=42))
        draws_b = list(subject_bootstrap(ratings_3x2, n_resamples=5, random_state=42))
        for a, b in zip(draws_a, draws_b):
            np.testing.assert_array_equal(a, b)

    def test_different_seed_different_results(self, ratings_3x2):
        draws_a = list(subject_bootstrap(ratings_3x2, n_resamples=20, random_state=1))
        draws_b = list(subject_bootstrap(ratings_3x2, n_resamples=20, random_state=2))
        any_different = any(
            not np.array_equal(a, b) for a, b in zip(draws_a, draws_b)
        )
        assert any_different

    def test_none_seed_runs(self, ratings_3x2):
        draws = list(subject_bootstrap(ratings_3x2, n_resamples=3, random_state=None))
        assert len(draws) == 3

    def test_default_random_state_is_42(self, ratings_3x2):
        """Default seed must be 42."""
        draws_default = list(subject_bootstrap(ratings_3x2, n_resamples=5))
        draws_42 = list(subject_bootstrap(ratings_3x2, n_resamples=5, random_state=42))
        for a, b in zip(draws_default, draws_42):
            np.testing.assert_array_equal(a, b)

    def test_single_resample(self, ratings_3x2):
        draws = list(subject_bootstrap(ratings_3x2, n_resamples=1))
        assert len(draws) == 1
        assert draws[0].shape == (3, 2)


# ================================================================== #
#                                                                      #
#   subject_jackknife                                                  #
#                                                                      #
# ================================================================== #


class TestSubjectJackknifeValidation:
    """Type and value checks for subject_jackknife."""

    def test_category_index_not_ndarray(self):
        with pytest.raises(TypeError, match="category_index must be np.ndarray"):
            list(subject_jackknife([[0, 1], [1, 0]]))

    def test_category_index_not_2d(self):
        with pytest.raises(ValueError, match="category_index must be 2D"):
            list(subject_jackknife(np.array([0, 1, 2])))

    def test_fewer_than_2_subjects(self):
        with pytest.raises(ValueError, match="requires >= 2 subjects"):
            list(subject_jackknife(np.array([[0, 1]], dtype=np.int8)))


class TestSubjectJackknifeBehavior:
    """Functional tests for subject_jackknife."""

    def test_yields_n_subjects_draws(self, ratings_3x2):
        draws = list(subject_jackknife(ratings_3x2))
        assert len(draws) == 3

    def test_each_draw_has_correct_shape(self, ratings_3x2):
        for draw in subject_jackknife(ratings_3x2):
            assert draw.shape == (2, 2)

    def test_each_draw_omits_one_subject(self, ratings_3x2):
        """Draw i must contain all rows except row i."""
        for i, draw in enumerate(subject_jackknife(ratings_3x2)):
            expected = np.delete(ratings_3x2, i, axis=0)
            np.testing.assert_array_equal(draw, expected)

    def test_exactly_2_subjects(self):
        data = np.array([[0, 1], [1, 0]], dtype=np.int8)
        draws = list(subject_jackknife(data))
        assert len(draws) == 2
        assert draws[0].shape == (1, 2)


# ================================================================== #
#                                                                      #
#   subject_bootstrap_from_nij_matrix                                  #
#                                                                      #
# ================================================================== #


class TestSubjectBootstrapNijValidation:
    """Type and value checks for subject_bootstrap_from_nij_matrix."""

    def test_nij_matrix_not_ndarray(self):
        with pytest.raises(TypeError, match="nij_matrix must be np.ndarray"):
            list(subject_bootstrap_from_nij_matrix([[1, 2]], n_resamples=1))

    def test_n_resamples_not_int(self, nij_matrix_4x3):
        with pytest.raises(TypeError, match="n_resamples must be int"):
            list(subject_bootstrap_from_nij_matrix(nij_matrix_4x3, n_resamples=5.0))

    def test_random_state_not_int_or_none(self, nij_matrix_4x3):
        with pytest.raises(TypeError, match="random_state must be int or None"):
            list(subject_bootstrap_from_nij_matrix(
                nij_matrix_4x3, random_state="bad"
            ))

    def test_nij_matrix_not_2d(self):
        with pytest.raises(ValueError, match="nij_matrix must be 2D"):
            list(subject_bootstrap_from_nij_matrix(
                np.array([1, 2, 3]), n_resamples=1
            ))

    def test_n_resamples_zero(self, nij_matrix_4x3):
        with pytest.raises(ValueError, match="n_resamples must be >= 1"):
            list(subject_bootstrap_from_nij_matrix(
                nij_matrix_4x3, n_resamples=0
            ))


class TestSubjectBootstrapNijBehavior:
    """Functional tests for subject_bootstrap_from_nij_matrix."""

    def test_yields_correct_count(self, nij_matrix_4x3):
        draws = list(subject_bootstrap_from_nij_matrix(
            nij_matrix_4x3, n_resamples=7
        ))
        assert len(draws) == 7

    def test_each_draw_has_correct_shape(self, nij_matrix_4x3):
        for draw in subject_bootstrap_from_nij_matrix(
            nij_matrix_4x3, n_resamples=5
        ):
            assert draw.shape == (4, 3)

    def test_draws_contain_original_rows(self, nij_matrix_4x3):
        original_set = {tuple(row) for row in nij_matrix_4x3}
        for draw in subject_bootstrap_from_nij_matrix(
            nij_matrix_4x3, n_resamples=10
        ):
            for row in draw:
                assert tuple(row) in original_set

    def test_reproducibility_same_seed(self, nij_matrix_4x3):
        draws_a = list(subject_bootstrap_from_nij_matrix(
            nij_matrix_4x3, n_resamples=5, random_state=42
        ))
        draws_b = list(subject_bootstrap_from_nij_matrix(
            nij_matrix_4x3, n_resamples=5, random_state=42
        ))
        for a, b in zip(draws_a, draws_b):
            np.testing.assert_array_equal(a, b)

    def test_default_random_state_is_42(self, nij_matrix_4x3):
        draws_default = list(subject_bootstrap_from_nij_matrix(
            nij_matrix_4x3, n_resamples=5
        ))
        draws_42 = list(subject_bootstrap_from_nij_matrix(
            nij_matrix_4x3, n_resamples=5, random_state=42
        ))
        for a, b in zip(draws_default, draws_42):
            np.testing.assert_array_equal(a, b)


# ================================================================== #
#                                                                      #
#   rater_jackknife                                                    #
#                                                                      #
# ================================================================== #


class TestRaterJackknifeValidation:
    """Type and value checks for rater_jackknife."""

    def test_category_index_not_ndarray(self):
        with pytest.raises(TypeError, match="category_index must be np.ndarray"):
            list(rater_jackknife([[0, 1], [1, 0]]))

    def test_category_index_not_2d(self):
        with pytest.raises(ValueError, match="category_index must be 2D"):
            list(rater_jackknife(np.array([0, 1, 2])))

    def test_fewer_than_2_raters(self):
        with pytest.raises(ValueError, match="requires >= 2 raters"):
            list(rater_jackknife(np.array([[0], [1]], dtype=np.int8)))


class TestRaterJackknifeBehavior:
    """Functional tests for rater_jackknife."""

    def test_yields_n_raters_draws(self, ratings_4x3):
        draws = list(rater_jackknife(ratings_4x3))
        assert len(draws) == 3

    def test_each_draw_has_correct_shape(self, ratings_4x3):
        for draw in rater_jackknife(ratings_4x3):
            assert draw.shape == (4, 2)

    def test_each_draw_omits_one_rater(self, ratings_4x3):
        """Draw i must contain all columns except column i."""
        for i, draw in enumerate(rater_jackknife(ratings_4x3)):
            expected = np.delete(ratings_4x3, i, axis=1)
            np.testing.assert_array_equal(draw, expected)

    def test_exactly_2_raters(self, ratings_3x2):
        draws = list(rater_jackknife(ratings_3x2))
        assert len(draws) == 2
        assert draws[0].shape == (3, 1)


# ================================================================== #
#                                                                      #
#   rater_subset                                                       #
#                                                                      #
# ================================================================== #


class TestRaterSubsetValidation:
    """Type and value checks for rater_subset."""

    def test_category_index_not_ndarray(self):
        with pytest.raises(TypeError, match="category_index must be np.ndarray"):
            list(rater_subset([[0, 1, 0]], k=2))

    def test_k_not_int(self, ratings_4x3):
        with pytest.raises(TypeError, match="k must be int"):
            list(rater_subset(ratings_4x3, k=2.0))

    def test_category_index_not_2d(self):
        with pytest.raises(ValueError, match="category_index must be 2D"):
            list(rater_subset(np.array([0, 1, 2]), k=2))

    def test_k_zero(self, ratings_4x3):
        with pytest.raises(ValueError, match="k must be between 1 and"):
            list(rater_subset(ratings_4x3, k=0))

    def test_k_exceeds_n_raters(self, ratings_4x3):
        with pytest.raises(ValueError, match="k must be between 1 and"):
            list(rater_subset(ratings_4x3, k=4))


class TestRaterSubsetBehavior:
    """Functional tests for rater_subset."""

    def test_yields_c_n_k_draws(self, ratings_4x3):
        # C(3, 2) = 3
        draws = list(rater_subset(ratings_4x3, k=2))
        assert len(draws) == comb(3, 2)

    def test_each_draw_has_correct_shape(self, ratings_4x3):
        for draw in rater_subset(ratings_4x3, k=2):
            assert draw.shape == (4, 2)

    def test_k_equals_1(self, ratings_4x3):
        # C(3, 1) = 3 -- each draw is a single column
        draws = list(rater_subset(ratings_4x3, k=1))
        assert len(draws) == 3
        for draw in draws:
            assert draw.shape == (4, 1)

    def test_k_equals_n_raters(self, ratings_4x3):
        # C(3, 3) = 1 -- the full data
        draws = list(rater_subset(ratings_4x3, k=3))
        assert len(draws) == 1
        np.testing.assert_array_equal(draws[0], ratings_4x3)

    def test_draws_contain_original_columns(self, ratings_4x3):
        """Each column in a draw must be a column from the original."""
        for draw in rater_subset(ratings_4x3, k=2):
            for col_idx in range(draw.shape[1]):
                col = draw[:, col_idx]
                found = any(
                    np.array_equal(col, ratings_4x3[:, j])
                    for j in range(ratings_4x3.shape[1])
                )
                assert found


# ================================================================== #
#                                                                      #
#   column_shuffle                                                     #
#                                                                      #
# ================================================================== #


class TestColumnShuffleValidation:
    """Type and value checks for column_shuffle."""

    def test_category_index_not_ndarray(self):
        with pytest.raises(TypeError, match="category_index must be np.ndarray"):
            list(column_shuffle([[0, 1]], n_permutations=1))

    def test_n_permutations_not_int(self, ratings_3x2):
        with pytest.raises(TypeError, match="n_permutations must be int"):
            list(column_shuffle(ratings_3x2, n_permutations=5.0))

    def test_random_state_not_int_or_none(self, ratings_3x2):
        with pytest.raises(TypeError, match="random_state must be int or None"):
            list(column_shuffle(ratings_3x2, random_state="bad"))

    def test_category_index_not_2d(self):
        with pytest.raises(ValueError, match="category_index must be 2D"):
            list(column_shuffle(np.array([0, 1, 2]), n_permutations=1))

    def test_n_permutations_zero(self, ratings_3x2):
        with pytest.raises(ValueError, match="n_permutations must be >= 1"):
            list(column_shuffle(ratings_3x2, n_permutations=0))


class TestColumnShuffleBehavior:
    """Functional tests for column_shuffle."""

    def test_yields_correct_count(self, ratings_3x2):
        draws = list(column_shuffle(ratings_3x2, n_permutations=7))
        assert len(draws) == 7

    def test_each_draw_has_correct_shape(self, ratings_3x2):
        for draw in column_shuffle(ratings_3x2, n_permutations=5):
            assert draw.shape == (3, 2)

    def test_preserves_per_column_values(self, ratings_4x3):
        """Each column must contain the same multiset of values."""
        for draw in column_shuffle(ratings_4x3, n_permutations=10):
            for col in range(ratings_4x3.shape[1]):
                np.testing.assert_array_equal(
                    np.sort(draw[:, col]),
                    np.sort(ratings_4x3[:, col]),
                )

    def test_preserves_missing_positions(self, ratings_with_missing):
        """Missing values (-1) should stay in their positions after shuffle.
        Actually, column_shuffle shuffles only valid (>= 0) entries within
        each column, so -1 cells keep their positions."""
        missing_mask = ratings_with_missing < 0
        for draw in column_shuffle(
            ratings_with_missing, n_permutations=10
        ):
            draw_missing = draw < 0
            np.testing.assert_array_equal(draw_missing, missing_mask)

    def test_reproducibility_same_seed(self, ratings_3x2):
        draws_a = list(column_shuffle(ratings_3x2, n_permutations=5, random_state=42))
        draws_b = list(column_shuffle(ratings_3x2, n_permutations=5, random_state=42))
        for a, b in zip(draws_a, draws_b):
            np.testing.assert_array_equal(a, b)

    def test_default_random_state_is_42(self, ratings_3x2):
        draws_default = list(column_shuffle(ratings_3x2, n_permutations=5))
        draws_42 = list(column_shuffle(
            ratings_3x2, n_permutations=5, random_state=42
        ))
        for a, b in zip(draws_default, draws_42):
            np.testing.assert_array_equal(a, b)


# ================================================================== #
#                                                                      #
#   within_subject_shuffle                                             #
#                                                                      #
# ================================================================== #


class TestWithinSubjectShuffleValidation:
    """Type and value checks for within_subject_shuffle."""

    def test_category_index_not_ndarray(self):
        with pytest.raises(TypeError, match="category_index must be np.ndarray"):
            list(within_subject_shuffle([[0, 1]], n_permutations=1))

    def test_n_permutations_not_int(self, ratings_3x2):
        with pytest.raises(TypeError, match="n_permutations must be int"):
            list(within_subject_shuffle(ratings_3x2, n_permutations=5.0))

    def test_random_state_not_int_or_none(self, ratings_3x2):
        with pytest.raises(TypeError, match="random_state must be int or None"):
            list(within_subject_shuffle(ratings_3x2, random_state="bad"))

    def test_category_index_not_2d(self):
        with pytest.raises(ValueError, match="category_index must be 2D"):
            list(within_subject_shuffle(np.array([0, 1, 2]), n_permutations=1))

    def test_n_permutations_zero(self, ratings_3x2):
        with pytest.raises(ValueError, match="n_permutations must be >= 1"):
            list(within_subject_shuffle(ratings_3x2, n_permutations=0))


class TestWithinSubjectShuffleBehavior:
    """Functional tests for within_subject_shuffle."""

    def test_yields_correct_count(self, ratings_3x2):
        draws = list(within_subject_shuffle(ratings_3x2, n_permutations=7))
        assert len(draws) == 7

    def test_each_draw_has_correct_shape(self, ratings_3x2):
        for draw in within_subject_shuffle(ratings_3x2, n_permutations=5):
            assert draw.shape == (3, 2)

    def test_preserves_per_row_values(self, ratings_4x3):
        """Each row must contain the same multiset of values."""
        for draw in within_subject_shuffle(ratings_4x3, n_permutations=10):
            for row in range(ratings_4x3.shape[0]):
                np.testing.assert_array_equal(
                    np.sort(draw[row, :]),
                    np.sort(ratings_4x3[row, :]),
                )

    def test_preserves_missing_positions(self, ratings_with_missing):
        """Missing values (-1) should stay in their positions."""
        missing_mask = ratings_with_missing < 0
        for draw in within_subject_shuffle(
            ratings_with_missing, n_permutations=10
        ):
            draw_missing = draw < 0
            np.testing.assert_array_equal(draw_missing, missing_mask)

    def test_reproducibility_same_seed(self, ratings_3x2):
        draws_a = list(within_subject_shuffle(
            ratings_3x2, n_permutations=5, random_state=42
        ))
        draws_b = list(within_subject_shuffle(
            ratings_3x2, n_permutations=5, random_state=42
        ))
        for a, b in zip(draws_a, draws_b):
            np.testing.assert_array_equal(a, b)

    def test_default_random_state_is_42(self, ratings_3x2):
        draws_default = list(within_subject_shuffle(
            ratings_3x2, n_permutations=5
        ))
        draws_42 = list(within_subject_shuffle(
            ratings_3x2, n_permutations=5, random_state=42
        ))
        for a, b in zip(draws_default, draws_42):
            np.testing.assert_array_equal(a, b)


# ================================================================== #
#                                                                      #
#   pool_and_redistribute                                              #
#                                                                      #
# ================================================================== #


class TestPoolAndRedistributeValidation:
    """Type and value checks for pool_and_redistribute."""

    def test_category_index_not_ndarray(self):
        with pytest.raises(TypeError, match="category_index must be np.ndarray"):
            list(pool_and_redistribute([[0, 1]], n_permutations=1))

    def test_n_permutations_not_int(self, ratings_3x2):
        with pytest.raises(TypeError, match="n_permutations must be int"):
            list(pool_and_redistribute(ratings_3x2, n_permutations=5.0))

    def test_random_state_not_int_or_none(self, ratings_3x2):
        with pytest.raises(TypeError, match="random_state must be int or None"):
            list(pool_and_redistribute(ratings_3x2, random_state="bad"))

    def test_category_index_not_2d(self):
        with pytest.raises(ValueError, match="category_index must be 2D"):
            list(pool_and_redistribute(np.array([0, 1, 2]), n_permutations=1))

    def test_n_permutations_zero(self, ratings_3x2):
        with pytest.raises(ValueError, match="n_permutations must be >= 1"):
            list(pool_and_redistribute(ratings_3x2, n_permutations=0))


class TestPoolAndRedistributeBehavior:
    """Functional tests for pool_and_redistribute."""

    def test_yields_correct_count(self, ratings_3x2):
        draws = list(pool_and_redistribute(ratings_3x2, n_permutations=7))
        assert len(draws) == 7

    def test_each_draw_has_correct_shape(self, ratings_3x2):
        for draw in pool_and_redistribute(ratings_3x2, n_permutations=5):
            assert draw.shape == (3, 2)

    def test_preserves_global_value_counts(self, ratings_4x3):
        """The multiset of all non-missing values must be preserved."""
        original_vals = np.sort(ratings_4x3[ratings_4x3 >= 0])
        for draw in pool_and_redistribute(ratings_4x3, n_permutations=10):
            draw_vals = np.sort(draw[draw >= 0])
            np.testing.assert_array_equal(draw_vals, original_vals)

    def test_preserves_missing_positions(self, ratings_with_missing):
        """Missing values (-1) should stay in their positions."""
        missing_mask = ratings_with_missing < 0
        for draw in pool_and_redistribute(
            ratings_with_missing, n_permutations=10
        ):
            draw_missing = draw < 0
            np.testing.assert_array_equal(draw_missing, missing_mask)

    def test_preserves_global_value_counts_with_missing(
        self, ratings_with_missing
    ):
        """Global value counts preserved even with missing data."""
        valid_mask = ratings_with_missing >= 0
        original_vals = np.sort(ratings_with_missing[valid_mask])
        for draw in pool_and_redistribute(
            ratings_with_missing, n_permutations=10
        ):
            draw_vals = np.sort(draw[valid_mask])
            np.testing.assert_array_equal(draw_vals, original_vals)

    def test_reproducibility_same_seed(self, ratings_3x2):
        draws_a = list(pool_and_redistribute(
            ratings_3x2, n_permutations=5, random_state=42
        ))
        draws_b = list(pool_and_redistribute(
            ratings_3x2, n_permutations=5, random_state=42
        ))
        for a, b in zip(draws_a, draws_b):
            np.testing.assert_array_equal(a, b)

    def test_default_random_state_is_42(self, ratings_3x2):
        draws_default = list(pool_and_redistribute(
            ratings_3x2, n_permutations=5
        ))
        draws_42 = list(pool_and_redistribute(
            ratings_3x2, n_permutations=5, random_state=42
        ))
        for a, b in zip(draws_default, draws_42):
            np.testing.assert_array_equal(a, b)
