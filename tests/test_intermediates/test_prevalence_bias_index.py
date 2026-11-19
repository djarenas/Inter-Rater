"""
PrevalenceBiasIndex: Byrt, Bishop and Carlin (1993), two raters, two categories.

The table used throughout (rater_a = Ann, rater_b = Bob, positive = "include"):

                 Bob include   Bob exclude
    Ann include      a = 40        c = 5
    Ann exclude      b = 10        d = 45        N = 100

epiR's formulas give PI = a/N - d/N = -0.05 and
BI = (a + b)/N - (a + c)/N = +0.05.
"""

import numpy as np
import pytest

from interrater import EMPTY
from interrater.intermediates.data_profiles import PrevalenceBiasIndex

CATS = ("include", "exclude")
RATERS = ("Ann", "Bob")
A, B, C, D = 40, 10, 5, 45
N = A + B + C + D


def _pair_index():
    """Rows in the four cells above; 0 = include, 1 = exclude."""
    rows = [[0, 0]] * A + [[1, 0]] * B + [[0, 1]] * C + [[1, 1]] * D
    return np.array(rows, dtype=np.int8)


def _from_pair():
    return PrevalenceBiasIndex.build_from_pair(_pair_index(), CATS, RATERS)


class TestFromPair:

    def test_cells(self):
        p = _from_pair()
        assert (p.n_both_positive, p.n_both_negative) == (A, D)
        assert (p.n_only_b_positive, p.n_only_a_positive) == (B, C)
        assert p.n_disagree == B + C
        assert p.n_subjects == N

    def test_indices_match_epir(self):
        p = _from_pair()
        assert p.prevalence_index == pytest.approx(A / N - D / N)
        assert p.bias_index == pytest.approx((A + B) / N - (A + C) / N)

    def test_orientation_in_words(self):
        text = str(_from_pair())
        assert "positive = 'include'" in text
        assert "Ann vs Bob" in text
        assert "PI = -0.050: more agreement on 'exclude' than on 'include'" in text
        assert "BI = +0.050: Bob says 'include' more often than Ann" in text

    def test_swapping_raters_flips_the_bias_sign(self):
        swapped = PrevalenceBiasIndex.build_from_pair(
            _pair_index()[:, ::-1], CATS, RATERS[::-1]
        )
        assert swapped.bias_index == pytest.approx(-_from_pair().bias_index)
        assert swapped.prevalence_index == pytest.approx(_from_pair().prevalence_index)

    def test_subjects_not_rated_by_both_are_left_out(self):
        index = np.vstack([_pair_index(), [[0, EMPTY], [EMPTY, 1]]])
        assert PrevalenceBiasIndex.build_from_pair(index, CATS, RATERS).n_subjects == N

    def test_needs_two_columns(self):
        with pytest.raises(ValueError, match="exactly 2 rater columns"):
            PrevalenceBiasIndex.build_from_pair(np.zeros((5, 3), dtype=np.int8), CATS, RATERS)


class TestSelfChecks:

    def _kwargs(self, **overrides):
        base = dict(
            positive_category="include", negative_category="exclude",
            rater_a="Ann", rater_b="Bob",
            n_both_positive=A, n_both_negative=D,
            n_only_b_positive=B, n_only_a_positive=C,
        )
        base.update(overrides)
        return base

    def test_same_category_twice_refused(self):
        with pytest.raises(ValueError, match="both 'include'"):
            PrevalenceBiasIndex(**self._kwargs(negative_category="include"))

    def test_same_rater_twice_refused(self):
        with pytest.raises(ValueError, match="both 'Ann'"):
            PrevalenceBiasIndex(**self._kwargs(rater_b="Ann"))

    def test_negative_count_refused(self):
        with pytest.raises(ValueError, match=">= 0"):
            PrevalenceBiasIndex(**self._kwargs(n_both_positive=-1))

    def test_no_subjects_refused(self):
        with pytest.raises(ValueError, match="no subjects"):
            PrevalenceBiasIndex(**self._kwargs(
                n_both_positive=0, n_both_negative=0,
                n_only_b_positive=0, n_only_a_positive=0,
            ))

    def test_indices_always_between_minus_one_and_one(self):
        """They are worked out from the cells, so they cannot leave the range."""
        extreme = PrevalenceBiasIndex(**self._kwargs(
            n_both_positive=0, n_both_negative=0,
            n_only_b_positive=7, n_only_a_positive=0,
        ))
        assert extreme.bias_index == 1.0
        assert -1.0 <= extreme.prevalence_index <= 1.0
