"""Tests for variance_functions.py -- fleiss_null_variance, cohen_null_variance.

Deterministic tests only. Tolerance is 1e-12 for exact computations.
Hand-computed values are derived from the original paper formulas.
"""

from __future__ import annotations

import numpy as np
import pytest

from interrater.analyzers.variance_functions import (
    fleiss_null_variance,
    cohen_null_variance,
)


# ------------------------------------------------------------------ #
# Module-level _ERROR_PREFIX                                           #
# ------------------------------------------------------------------ #


class TestErrorPrefix:
    """_ERROR_PREFIX lives at module level, not inside functions."""

    def test_error_prefix_at_module_level(self):
        import interrater.analyzers.variance_functions as mod

        assert hasattr(mod, "_ERROR_PREFIX")
        assert mod._ERROR_PREFIX == "[variance_functions] Error"

    def test_error_prefix_appears_in_errors(self):
        with pytest.raises(ValueError, match=r"\[variance_functions\] Error"):
            fleiss_null_variance(np.array([0.5, 0.5]), n_subjects=0, n_raters=3)


# ------------------------------------------------------------------ #
# fleiss_null_variance                                                 #
# ------------------------------------------------------------------ #


class TestFleissNullVariance:
    """Tests for fleiss_null_variance() -- Fleiss (1979) Eq. 12."""

    def test_two_equal_categories_hand_computed(self):
        """
        pj = [0.5, 0.5], N=10, n=3.

        By hand (Eq. 12):
            pq = [0.5*0.5, 0.5*0.5] = [0.25, 0.25]
            sum_pq = 0.50
            sum_pq_sq = 0.25
            pq*(qj-pj) = [0.25*0, 0.25*0] = [0, 0]
            sum_pq_diff = 0

            num = 2 * (0.25 - 0) = 0.50
            den = 10 * 3 * 2 * 0.25 = 15.0
            Var = 0.50 / 15.0 = 1/30
        """
        pj = np.array([0.5, 0.5])
        var = fleiss_null_variance(pj, n_subjects=10, n_raters=3)
        expected = 1.0 / 30.0
        assert var == pytest.approx(expected, abs=1e-12)

    def test_three_equal_categories_hand_computed(self):
        """
        pj = [1/3, 1/3, 1/3], N=20, n=4.

        By hand:
            pq = [1/3*2/3, 1/3*2/3, 1/3*2/3] = [2/9, 2/9, 2/9]
            sum_pq = 6/9 = 2/3
            sum_pq_sq = 4/9
            qj-pj = [2/3 - 1/3, ...] = [1/3, 1/3, 1/3]
            pq*(qj-pj) = [2/9 * 1/3, ...] = [2/27, 2/27, 2/27]
            sum_pq_diff = 6/27 = 2/9

            num = 2 * (4/9 - 2/9) = 2 * 2/9 = 4/9
            den = 20 * 4 * 3 * 4/9 = 240 * 4/9 = 960/9
            Var = (4/9) / (960/9) = 4/960 = 1/240
        """
        pj = np.array([1 / 3, 1 / 3, 1 / 3])
        var = fleiss_null_variance(pj, n_subjects=20, n_raters=4)
        expected = 1.0 / 240.0
        assert var == pytest.approx(expected, abs=1e-12)

    def test_skewed_two_categories(self):
        """
        pj = [0.8, 0.2], N=50, n=5.

        By hand:
            pq = [0.8*0.2, 0.2*0.8] = [0.16, 0.16]
            sum_pq = 0.32
            sum_pq_sq = 0.1024
            qj-pj: [0.2-0.8, 0.8-0.2] = [-0.6, 0.6]
            pq*(qj-pj) = [0.16*(-0.6), 0.16*0.6] = [-0.096, 0.096]
            sum_pq_diff = 0

            num = 2 * (0.1024 - 0) = 0.2048
            den = 50 * 5 * 4 * 0.1024 = 1000 * 0.1024 = 102.4
            Var = 0.2048 / 102.4 = 0.002
        """
        pj = np.array([0.8, 0.2])
        var = fleiss_null_variance(pj, n_subjects=50, n_raters=5)
        expected = 0.002
        assert var == pytest.approx(expected, abs=1e-12)

    def test_variance_decreases_with_more_subjects(self):
        """More subjects -> smaller variance (law of large numbers)."""
        pj = np.array([0.5, 0.5])
        var_10 = fleiss_null_variance(pj, n_subjects=10, n_raters=3)
        var_100 = fleiss_null_variance(pj, n_subjects=100, n_raters=3)
        assert var_100 < var_10

    def test_variance_decreases_with_more_raters(self):
        """More raters per subject -> smaller variance."""
        pj = np.array([0.5, 0.5])
        var_2 = fleiss_null_variance(pj, n_subjects=20, n_raters=2)
        var_5 = fleiss_null_variance(pj, n_subjects=20, n_raters=5)
        assert var_5 < var_2

    def test_returns_float(self):
        pj = np.array([0.5, 0.5])
        result = fleiss_null_variance(pj, n_subjects=10, n_raters=3)
        assert isinstance(result, float)

    # --- Guards ---

    def test_n_subjects_zero_raises(self):
        with pytest.raises(ValueError, match="n_subjects must be positive"):
            fleiss_null_variance(np.array([0.5, 0.5]), n_subjects=0, n_raters=3)

    def test_n_subjects_negative_raises(self):
        with pytest.raises(ValueError, match="n_subjects must be positive"):
            fleiss_null_variance(np.array([0.5, 0.5]), n_subjects=-1, n_raters=3)

    def test_n_raters_zero_raises(self):
        with pytest.raises(ValueError, match="n_raters must be positive"):
            fleiss_null_variance(np.array([0.5, 0.5]), n_subjects=10, n_raters=0)

    def test_n_raters_negative_raises(self):
        with pytest.raises(ValueError, match="n_raters must be positive"):
            fleiss_null_variance(np.array([0.5, 0.5]), n_subjects=10, n_raters=-1)

    def test_pj_not_summing_to_one_raises(self):
        with pytest.raises(ValueError, match="pj must sum to 1"):
            fleiss_null_variance(np.array([0.5, 0.3]), n_subjects=10, n_raters=3)

    def test_single_category_raises(self):
        """pj=[1.0] -> all ratings in one category -> denominator zero."""
        with pytest.raises(ValueError, match="denominator is zero"):
            fleiss_null_variance(np.array([1.0]), n_subjects=10, n_raters=3)


# ------------------------------------------------------------------ #
# cohen_null_variance                                                  #
# ------------------------------------------------------------------ #


class TestCohenNullVariance:
    """Tests for cohen_null_variance() -- Fleiss, Cohen & Everitt (1969)."""

    def test_unweighted_equal_marginals_hand_computed(self):
        """
        p_a = p_b = [0.5, 0.5], N=20, unweighted.

        Eq. 14 (unweighted special case):
            Pe = sum(p_a_i * p_b_i) = 0.5*0.5 + 0.5*0.5 = 0.50
            Pe^2 = 0.25
            sum_i p_a_i * p_b_i * (p_a_i + p_b_i)
                = 0.25*(0.5+0.5) + 0.25*(0.5+0.5)
                = 0.25 + 0.25 = 0.50

            num = Pe + Pe^2 - sum_term = 0.50 + 0.25 - 0.50 = 0.25
            den = N * (1-Pe)^2 = 20 * 0.25 = 5.0
            Var = 0.25 / 5.0 = 0.05
        """
        p_a = np.array([0.5, 0.5])
        p_b = np.array([0.5, 0.5])
        var = cohen_null_variance(p_a, p_b, n_subjects=20)
        expected = 0.05
        assert var == pytest.approx(expected, abs=1e-12)

    def test_unweighted_asymmetric_marginals_hand_computed(self):
        """
        p_a = [0.7, 0.3], p_b = [0.6, 0.4], N=30, unweighted.

        Eq. 14:
            Pe = 0.7*0.6 + 0.3*0.4 = 0.42 + 0.12 = 0.54
            Pe^2 = 0.2916
            sum_i p_a_i * p_b_i * (p_a_i + p_b_i)
                = 0.42*(0.7+0.6) + 0.12*(0.3+0.4)
                = 0.42*1.3 + 0.12*0.7
                = 0.546 + 0.084 = 0.630

            num = 0.54 + 0.2916 - 0.63 = 0.2016
            den = 30 * (1-0.54)^2 = 30 * 0.2116 = 6.348
            Var = 0.2016 / 6.348
        """
        p_a = np.array([0.7, 0.3])
        p_b = np.array([0.6, 0.4])
        var = cohen_null_variance(p_a, p_b, n_subjects=30)
        expected = 0.2016 / 6.348
        assert var == pytest.approx(expected, abs=1e-12)

    def test_unweighted_matches_general_formula(self):
        """Unweighted (W=None) should equal passing W=identity explicitly."""
        p_a = np.array([0.4, 0.3, 0.3])
        p_b = np.array([0.5, 0.2, 0.3])
        W = np.eye(3)

        var_none = cohen_null_variance(p_a, p_b, n_subjects=50)
        var_eye = cohen_null_variance(p_a, p_b, n_subjects=50, weight_matrix=W)

        assert var_none == pytest.approx(var_eye, abs=1e-12)

    def test_weighted_three_categories_hand_computed(self):
        """
        p_a = [0.5, 0.3, 0.2], p_b = [0.4, 0.4, 0.2], N=40.
        Linear weights for K=3:
            W = [[1, 0.5, 0],
                 [0.5, 1, 0.5],
                 [0, 0.5, 1]]

        Computed step by step via the general formula (Eq. 9):
            w_bar_row = W @ p_b  (row-weighted means)
            w_bar_col = W.T @ p_a  (col-weighted means)
            Pe_w = p_a @ W @ p_b  (weighted chance)
            diff = W - w_bar_row[:,None] - w_bar_col[None,:]
            joint = p_a[:,None] * p_b[None,:]
            sum_term = sum(joint * diff^2)
            Var = (sum_term - Pe_w^2) / (N * (1 - Pe_w)^2)
        """
        p_a = np.array([0.5, 0.3, 0.2])
        p_b = np.array([0.4, 0.4, 0.2])
        W = np.array([[1.0, 0.5, 0.0], [0.5, 1.0, 0.5], [0.0, 0.5, 1.0]])
        N = 40

        # Compute expected by the same formula (but independently)
        w_bar_row = W @ p_b
        w_bar_col = W.T @ p_a
        Pe_w = float(p_a @ W @ p_b)
        diff = W - w_bar_row[:, np.newaxis] - w_bar_col[np.newaxis, :]
        joint = p_a[:, np.newaxis] * p_b[np.newaxis, :]
        sum_term = float(np.sum(joint * diff**2))
        expected = (sum_term - Pe_w**2) / (N * (1 - Pe_w) ** 2)

        var = cohen_null_variance(p_a, p_b, n_subjects=N, weight_matrix=W)
        assert var == pytest.approx(expected, abs=1e-12)

    def test_variance_decreases_with_more_subjects(self):
        """More subjects -> smaller variance."""
        p_a = np.array([0.5, 0.5])
        p_b = np.array([0.5, 0.5])
        var_20 = cohen_null_variance(p_a, p_b, n_subjects=20)
        var_200 = cohen_null_variance(p_a, p_b, n_subjects=200)
        assert var_200 < var_20

    def test_returns_float(self):
        p_a = np.array([0.5, 0.5])
        p_b = np.array([0.5, 0.5])
        result = cohen_null_variance(p_a, p_b, n_subjects=20)
        assert isinstance(result, float)

    # --- Guards ---

    def test_n_subjects_zero_raises(self):
        with pytest.raises(ValueError, match="n_subjects must be positive"):
            cohen_null_variance(
                np.array([0.5, 0.5]), np.array([0.5, 0.5]), n_subjects=0
            )

    def test_n_subjects_negative_raises(self):
        with pytest.raises(ValueError, match="n_subjects must be positive"):
            cohen_null_variance(
                np.array([0.5, 0.5]), np.array([0.5, 0.5]), n_subjects=-1
            )

    def test_p_a_not_summing_to_one_raises(self):
        with pytest.raises(ValueError, match="p_a must sum to 1"):
            cohen_null_variance(
                np.array([0.5, 0.3]), np.array([0.5, 0.5]), n_subjects=20
            )

    def test_p_b_not_summing_to_one_raises(self):
        with pytest.raises(ValueError, match="p_b must sum to 1"):
            cohen_null_variance(
                np.array([0.5, 0.5]), np.array([0.4, 0.4]), n_subjects=20
            )

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="same shape"):
            cohen_null_variance(
                np.array([0.5, 0.5]),
                np.array([0.4, 0.3, 0.3]),
                n_subjects=20,
            )

    def test_wrong_weight_matrix_shape_raises(self):
        with pytest.raises(ValueError, match="weight_matrix must be"):
            cohen_null_variance(
                np.array([0.5, 0.5]),
                np.array([0.5, 0.5]),
                n_subjects=20,
                weight_matrix=np.eye(3),
            )

    def test_single_category_raises(self):
        """Pe_w=1 -> denominator zero."""
        with pytest.raises(ValueError, match="denominator is zero"):
            cohen_null_variance(
                np.array([1.0]), np.array([1.0]), n_subjects=20
            )
