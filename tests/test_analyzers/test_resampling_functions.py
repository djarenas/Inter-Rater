"""Tests for the resampling_functions sub-package.

One test file covering the exported functions across the two modules:
  - subject_resampling: subject_bootstrap, subject_jackknife
  - permutation_resampling: within_rater_shuffle, within_subject_shuffle,
    pool_and_redistribute

All functions are generators. Validation runs on first next() call,
not at function-call time.
"""

from __future__ import annotations


import numpy as np
import pytest

from interrater.intermediates.estimator_parts import EstimatorInput, Weights
from interrater.analyzers.resampling_functions.subject_resampling import (
    subject_bootstrap,
    subject_jackknife,
)
from interrater.analyzers.resampling_functions.permutation_resampling import (
    within_rater_shuffle,
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
def data_3x2(ratings_3x2):
    """ratings_3x2 with K = 2, unweighted."""
    return EstimatorInput(ratings_3x2, 2)


@pytest.fixture
def data_k3_weighted(ratings_3x2):
    """ratings_3x2 declared over 3 categories, one never used, with weights."""
    return EstimatorInput(ratings_3x2, 3, weights=Weights.linear(3))


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


# ================================================================== #
#                                                                      #
#   subject_bootstrap                                                  #
#                                                                      #
# ================================================================== #


class TestSubjectBootstrapValidation:
    """Type and value checks for subject_bootstrap."""

    def test_data_not_estimator_input(self, ratings_3x2):
        with pytest.raises(TypeError, match="data must be EstimatorInput"):
            list(subject_bootstrap(ratings_3x2, n_resamples=1))

    def test_n_resamples_not_int(self, data_3x2):
        with pytest.raises(TypeError, match="n_resamples must be int"):
            list(subject_bootstrap(data_3x2, n_resamples=5.0))

    def test_random_state_not_int_or_none(self, data_3x2):
        with pytest.raises(TypeError, match="random_state must be int or None"):
            list(subject_bootstrap(data_3x2, random_state="bad"))

    def test_n_resamples_zero(self, data_3x2):
        with pytest.raises(ValueError, match="n_resamples must be >= 1"):
            list(subject_bootstrap(data_3x2, n_resamples=0))

    def test_n_resamples_negative(self, data_3x2):
        with pytest.raises(ValueError, match="n_resamples must be >= 1"):
            list(subject_bootstrap(data_3x2, n_resamples=-3))


class TestSubjectBootstrapBehavior:
    """Functional tests for subject_bootstrap."""

    def test_yields_correct_count(self, data_3x2):
        draws = list(subject_bootstrap(data_3x2, n_resamples=7))
        assert len(draws) == 7

    def test_each_draw_has_correct_shape(self, data_3x2):
        for draw in subject_bootstrap(data_3x2, n_resamples=5):
            assert draw.category_index.shape == (3, 2)

    def test_draws_contain_original_rows(self, data_3x2):
        """Each row in a draw must be a row from the original."""
        original_set = {tuple(row) for row in data_3x2.category_index}
        for draw in subject_bootstrap(data_3x2, n_resamples=10):
            for row in draw.category_index:
                assert tuple(row) in original_set

    def test_draws_keep_k_and_weights(self, data_k3_weighted):
        """A draw is computed against the rating task's K, never the draw's."""
        for draw in subject_bootstrap(data_k3_weighted, n_resamples=10):
            assert draw.n_categories == 3
            assert draw.weights is data_k3_weighted.weights

    def test_reproducibility_same_seed(self, data_3x2):
        draws_a = list(subject_bootstrap(data_3x2, n_resamples=5, random_state=42))
        draws_b = list(subject_bootstrap(data_3x2, n_resamples=5, random_state=42))
        for a, b in zip(draws_a, draws_b):
            np.testing.assert_array_equal(a.category_index, b.category_index)

    def test_different_seed_different_results(self, data_3x2):
        draws_a = list(subject_bootstrap(data_3x2, n_resamples=20, random_state=1))
        draws_b = list(subject_bootstrap(data_3x2, n_resamples=20, random_state=2))
        any_different = any(
            not np.array_equal(a.category_index, b.category_index)
            for a, b in zip(draws_a, draws_b)
        )
        assert any_different

    def test_none_seed_runs(self, data_3x2):
        draws = list(subject_bootstrap(data_3x2, n_resamples=3, random_state=None))
        assert len(draws) == 3

    def test_default_random_state_is_42(self, data_3x2):
        """Default seed must be 42."""
        draws_default = list(subject_bootstrap(data_3x2, n_resamples=5))
        draws_42 = list(subject_bootstrap(data_3x2, n_resamples=5, random_state=42))
        for a, b in zip(draws_default, draws_42):
            np.testing.assert_array_equal(a.category_index, b.category_index)

    def test_single_resample(self, data_3x2):
        draws = list(subject_bootstrap(data_3x2, n_resamples=1))
        assert len(draws) == 1
        assert draws[0].category_index.shape == (3, 2)


# ================================================================== #
#                                                                      #
#   subject_jackknife                                                  #
#                                                                      #
# ================================================================== #


class TestSubjectJackknifeValidation:
    """Type and value checks for subject_jackknife."""

    def test_data_not_estimator_input(self, ratings_3x2):
        with pytest.raises(TypeError, match="data must be EstimatorInput"):
            list(subject_jackknife(ratings_3x2))

    def test_fewer_than_2_subjects(self):
        data = EstimatorInput(np.array([[0, 1]], dtype=np.int8), 2)
        with pytest.raises(ValueError, match="requires >= 2 subjects"):
            list(subject_jackknife(data))


class TestSubjectJackknifeBehavior:
    """Functional tests for subject_jackknife."""

    def test_yields_n_subjects_draws(self, data_3x2):
        draws = list(subject_jackknife(data_3x2))
        assert len(draws) == 3

    def test_each_draw_has_correct_shape(self, data_3x2):
        for draw in subject_jackknife(data_3x2):
            assert draw.category_index.shape == (2, 2)

    def test_each_draw_omits_one_subject(self, data_3x2):
        """Draw i must contain all rows except row i."""
        for i, draw in enumerate(subject_jackknife(data_3x2)):
            expected = np.delete(data_3x2.category_index, i, axis=0)
            np.testing.assert_array_equal(draw.category_index, expected)

    def test_draws_keep_k_and_weights(self, data_k3_weighted):
        for draw in subject_jackknife(data_k3_weighted):
            assert draw.n_categories == 3
            assert draw.weights is data_k3_weighted.weights

    def test_exactly_2_subjects(self):
        data = EstimatorInput(np.array([[0, 1], [1, 0]], dtype=np.int8), 2)
        draws = list(subject_jackknife(data))
        assert len(draws) == 2
        assert draws[0].category_index.shape == (1, 2)


# ================================================================== #
#                                                                      #
#   within_rater_shuffle                                               #
#                                                                      #
# ================================================================== #


class TestWithinRaterShuffleValidation:
    """Type and value checks for within_rater_shuffle."""

    def test_category_index_not_ndarray(self):
        with pytest.raises(TypeError, match="category_index must be np.ndarray"):
            list(within_rater_shuffle([[0, 1]], n_permutations=1))

    def test_n_permutations_not_int(self, ratings_3x2):
        with pytest.raises(TypeError, match="n_permutations must be int"):
            list(within_rater_shuffle(ratings_3x2, n_permutations=5.0))

    def test_random_state_not_int_or_none(self, ratings_3x2):
        with pytest.raises(TypeError, match="random_state must be int or None"):
            list(within_rater_shuffle(ratings_3x2, random_state="bad"))

    def test_category_index_not_2d(self):
        with pytest.raises(ValueError, match="category_index must be 2D"):
            list(within_rater_shuffle(np.array([0, 1, 2]), n_permutations=1))

    def test_n_permutations_zero(self, ratings_3x2):
        with pytest.raises(ValueError, match="n_permutations must be >= 1"):
            list(within_rater_shuffle(ratings_3x2, n_permutations=0))


class TestWithinRaterShuffleBehavior:
    """Functional tests for within_rater_shuffle."""

    def test_yields_correct_count(self, ratings_3x2):
        draws = list(within_rater_shuffle(ratings_3x2, n_permutations=7))
        assert len(draws) == 7

    def test_each_draw_has_correct_shape(self, ratings_3x2):
        for draw in within_rater_shuffle(ratings_3x2, n_permutations=5):
            assert draw.shape == (3, 2)

    def test_preserves_per_column_values(self, ratings_4x3):
        """Each column must contain the same multiset of values."""
        for draw in within_rater_shuffle(ratings_4x3, n_permutations=10):
            for col in range(ratings_4x3.shape[1]):
                np.testing.assert_array_equal(
                    np.sort(draw[:, col]),
                    np.sort(ratings_4x3[:, col]),
                )

    def test_preserves_per_rater_answers_with_missing(self, ratings_with_missing):
        """With missing cells, each rater keeps exactly the answers they gave.

        Only the answers move, never the empty cells, so each rater's answer
        mix and workload are unchanged. Across draws the answers do move: a
        function that returned its input unchanged would fail here.
        """
        moved = False
        for draw in within_rater_shuffle(ratings_with_missing, n_permutations=50):
            for col in range(ratings_with_missing.shape[1]):
                given = ratings_with_missing[:, col]
                np.testing.assert_array_equal(
                    np.sort(draw[draw[:, col] >= 0, col]),
                    np.sort(given[given >= 0]),
                )
            moved |= not np.array_equal(draw, ratings_with_missing)
        assert moved

    def test_preserves_missing_positions(self, ratings_with_missing):
        """Missing values (-1) should stay in their positions after shuffle.
        Actually, within_rater_shuffle shuffles only valid (>= 0) entries within
        each column, so -1 cells keep their positions."""
        missing_mask = ratings_with_missing < 0
        for draw in within_rater_shuffle(
            ratings_with_missing, n_permutations=10
        ):
            draw_missing = draw < 0
            np.testing.assert_array_equal(draw_missing, missing_mask)

    def test_reproducibility_same_seed(self, ratings_3x2):
        draws_a = list(within_rater_shuffle(ratings_3x2, n_permutations=5, random_state=42))
        draws_b = list(within_rater_shuffle(ratings_3x2, n_permutations=5, random_state=42))
        for a, b in zip(draws_a, draws_b):
            np.testing.assert_array_equal(a, b)

    def test_default_random_state_is_42(self, ratings_3x2):
        draws_default = list(within_rater_shuffle(ratings_3x2, n_permutations=5))
        draws_42 = list(within_rater_shuffle(
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

    def test_preserves_per_subject_answers_with_missing(self, ratings_with_missing):
        """With missing cells, each subject keeps exactly the answers it got.

        Only the answers move between the raters who rated that subject; no
        answer reaches a rater who was not shown it. Across draws the answers
        do move: a function that returned its input unchanged would fail here.
        """
        moved = False
        for draw in within_subject_shuffle(ratings_with_missing, n_permutations=50):
            for row in range(ratings_with_missing.shape[0]):
                given = ratings_with_missing[row, :]
                np.testing.assert_array_equal(
                    np.sort(draw[row, draw[row, :] >= 0]),
                    np.sort(given[given >= 0]),
                )
            moved |= not np.array_equal(draw, ratings_with_missing)
        assert moved

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
