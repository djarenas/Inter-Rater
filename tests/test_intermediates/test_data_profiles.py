"""
Tests for the data_profiles sub-package: BalanceProfile and RaterSharedMatrix.

BalanceProfile describes the balance structure of ratings -- whether every
subject was rated the same number of times, and how raters and subjects
are distributed. RaterSharedMatrix counts, for every pair of raters, how
many subjects both of them rated.

Organised bottom-up: BalanceProfile first (it has no dependencies on
RaterSharedMatrix), then RaterSharedMatrix.

Hand-verified datasets from conftest anchor every test. Each is small
enough that every property can be checked by counting on fingers.
"""

import numpy as np
import pytest

from interrater.intermediates.data_profiles.balance_profile import BalanceProfile
from interrater.intermediates.data_profiles.rater_shared_matrix import RaterSharedMatrix
from interrater.intermediates.estimator_parts.nij_matrix import NijMatrix

from .conftest import BALANCED_NIJ, UNBALANCED_NIJ, SIMPLE_SHARED


# -- Helpers ----------------------------------------------------------------- #


def _balanced_profile(**overrides):
    """Build a valid balanced BalanceProfile, overriding specific fields.

    Default: 4 subjects, 3 raters, every subject rated 3 times.
    n_per_rater = [4, 4, 4] because 12 total ratings / 3 raters.
    """
    defaults = dict(
        is_balanced=True,
        n_raters=3,
        n_raters_exact=True,
        n_per_subject=np.array([3, 3, 3, 3]),
        n_per_rater=np.array([4, 4, 4]),
        subject_counts={3: 4},
        fill_rate=1.0,
        raters_per_subject=3,
    )
    defaults.update(overrides)
    return BalanceProfile(**defaults)


def _unbalanced_profile(**overrides):
    """Build a valid unbalanced BalanceProfile, overriding specific fields.

    Default: 4 subjects, 4 raters, n_per_subject = [3, 2, 4, 2].
    subject_counts = {2: 2, 3: 1, 4: 1}.
    fill_rate = 11 / 16.
    """
    defaults = dict(
        is_balanced=False,
        n_raters=4,
        n_raters_exact=False,
        n_per_subject=np.array([3, 2, 4, 2]),
        n_per_rater=np.array([3, 3, 3, 2]),
        subject_counts={2: 2, 3: 1, 4: 1},
        fill_rate=11.0 / 16.0,
        raters_per_subject=None,
    )
    defaults.update(overrides)
    return BalanceProfile(**defaults)


# ============================================================================ #
#                              BalanceProfile                                  #
# ============================================================================ #


class TestBalanceProfileConstruction:
    """How BalanceProfile validates and stores its fields."""

    def test_balanced_valid(self):
        """A balanced profile with consistent fields is accepted."""
        bp = _balanced_profile()
        assert bp.is_balanced is True
        assert bp.n_raters == 3
        assert bp.raters_per_subject == 3

    def test_unbalanced_valid(self):
        """An unbalanced profile with consistent fields is accepted."""
        bp = _unbalanced_profile()
        assert bp.is_balanced is False
        assert bp.n_raters == 4
        assert bp.raters_per_subject is None

    def test_frozen(self):
        """BalanceProfile is immutable."""
        bp = _balanced_profile()
        with pytest.raises(AttributeError):
            bp.is_balanced = False


class TestBalanceProfileValidationRejects:
    """BalanceProfile rejects invalid field values."""

    def test_n_raters_not_positive(self):
        """n_raters must be a positive integer."""
        with pytest.raises(ValueError, match="n_raters must be positive int"):
            _balanced_profile(n_raters=0)

    def test_n_raters_not_int(self):
        """n_raters rejects float even if whole number."""
        with pytest.raises(ValueError, match="n_raters must be positive int"):
            _balanced_profile(n_raters=3.0)

    def test_n_raters_exact_not_bool(self):
        """n_raters_exact must be bool, not int."""
        with pytest.raises(TypeError, match="n_raters_exact must be bool"):
            _balanced_profile(n_raters_exact=1)

    def test_n_per_subject_not_array(self):
        """n_per_subject must be a numpy array, not a list."""
        with pytest.raises(TypeError, match="n_per_subject must be numpy"):
            _balanced_profile(n_per_subject=[3, 3, 3, 3])

    def test_n_per_subject_empty(self):
        """n_per_subject must not be empty."""
        with pytest.raises(ValueError, match="n_per_subject is empty"):
            _balanced_profile(n_per_subject=np.array([], dtype=int))

    def test_n_per_subject_has_zero(self):
        """Every subject must have at least one rater."""
        with pytest.raises(ValueError, match="must be > 0"):
            _balanced_profile(n_per_subject=np.array([3, 0, 3, 3]))

    def test_n_per_rater_not_array(self):
        """n_per_rater must be a numpy array, not a list."""
        with pytest.raises(TypeError, match="n_per_rater must be numpy"):
            _balanced_profile(n_per_rater=[4, 4, 4])

    def test_fill_rate_out_of_range(self):
        """fill_rate must be in [0, 1]."""
        with pytest.raises(ValueError, match="fill_rate must be in"):
            _balanced_profile(fill_rate=1.5)

    def test_fill_rate_not_numeric(self):
        """fill_rate must be a number, not a string."""
        with pytest.raises(TypeError, match="fill_rate must be float"):
            _balanced_profile(fill_rate="high")

    def test_n_raters_exact_without_n_raters(self):
        """n_raters_exact=True requires n_raters to be set."""
        with pytest.raises(ValueError, match="n_raters_exact=True but n_raters is None"):
            _balanced_profile(n_raters=None, n_raters_exact=True)

    def test_error_prefix(self):
        """Every error starts with [BalanceProfile] Error."""
        with pytest.raises(ValueError, match=r"\[BalanceProfile\] Error"):
            _balanced_profile(n_raters=0)


class TestBalanceProfileConsistencyRejects:
    """Cross-field consistency checks between balanced/unbalanced state."""

    def test_balanced_missing_raters_per_subject(self):
        """Balanced profile must have raters_per_subject."""
        with pytest.raises(ValueError, match="raters_per_subject must not be None"):
            _balanced_profile(raters_per_subject=None)

    def test_balanced_n_per_subject_not_constant(self):
        """Balanced but n_per_subject varies -- contradiction."""
        with pytest.raises(ValueError, match="n_per_subject is not constant"):
            _balanced_profile(n_per_subject=np.array([3, 3, 2, 3]))

    def test_unbalanced_raters_per_subject_not_none(self):
        """Unbalanced must have raters_per_subject=None."""
        with pytest.raises(ValueError, match="raters_per_subject must be None"):
            _unbalanced_profile(raters_per_subject=3)

    def test_unbalanced_subject_counts_mismatch(self):
        """subject_counts must match n_per_subject distribution."""
        with pytest.raises(ValueError, match="subject_counts does not match"):
            _unbalanced_profile(subject_counts={2: 1, 3: 2, 4: 1})


class TestBalanceProfileProperties:
    """Computed properties of BalanceProfile."""

    def test_n_subjects(self):
        """n_subjects is the length of n_per_subject."""
        bp = _balanced_profile()
        assert bp.n_subjects == 4

    def test_n_raters_lower_bound(self):
        """n_raters_lower_bound is max of n_per_subject."""
        bp = _unbalanced_profile()
        assert bp.n_raters_lower_bound == 4  # max of [3, 2, 4, 2]

    def test_n_raters_lower_bound_balanced(self):
        """For balanced data, lower bound equals raters_per_subject."""
        bp = _balanced_profile()
        assert bp.n_raters_lower_bound == 3

    def test_n_raters_display_exact(self):
        """Exact n_raters displays as 'n_raters=X'."""
        bp = _balanced_profile(n_raters_exact=True)
        assert bp.n_raters_display == "n_raters=3"

    def test_n_raters_display_lower_bound(self):
        """Approximate n_raters displays as 'n_raters>=X'."""
        bp = _unbalanced_profile(n_raters_exact=False)
        # The >= symbol is a unicode ≥
        assert "≥" in bp.n_raters_display
        assert "4" in bp.n_raters_display


class TestBalanceProfileFromNijMatrix:
    """BalanceProfile.from_nij_matrix factory method."""

    def test_from_balanced(self, balanced_nij):
        """Balanced nij matrix produces balanced profile with None rater fields."""
        bp = BalanceProfile.from_nij_matrix(balanced_nij)
        assert bp.is_balanced is True
        assert bp.n_raters is None
        assert bp.n_raters_exact is False
        assert bp.raters_per_subject == 3
        np.testing.assert_array_equal(bp.n_per_subject, [3, 3, 3, 3])
        assert bp.subject_counts == {3: 4}
        # Rater-level fields are None -- rater identity is gone in nij
        assert bp.n_per_rater is None
        assert bp.fill_rate is None
        # Lower bound is still available
        assert bp.n_raters_lower_bound == 3

    def test_from_unbalanced(self, unbalanced_nij):
        """Unbalanced nij matrix produces unbalanced profile with None rater fields."""
        bp = BalanceProfile.from_nij_matrix(unbalanced_nij)
        assert bp.is_balanced is False
        assert bp.n_raters is None
        assert bp.n_raters_exact is False
        assert bp.raters_per_subject is None
        np.testing.assert_array_equal(bp.n_per_subject, [3, 2, 4, 2])
        assert bp.subject_counts == {2: 2, 3: 1, 4: 1}
        # Rater-level fields are None
        assert bp.n_per_rater is None
        assert bp.fill_rate is None
        # Lower bound = max of n_per_subject
        assert bp.n_raters_lower_bound == 4

    def test_rater_fields_none_from_nij(self, balanced_nij):
        """n_raters, n_per_rater, and fill_rate are all None from nij."""
        bp = BalanceProfile.from_nij_matrix(balanced_nij)
        assert bp.n_raters is None
        assert bp.n_per_rater is None
        assert bp.fill_rate is None

    def test_all_zeros_raises(self):
        """Nij matrix where every row sums to zero raises."""
        nij = NijMatrix(np.zeros((3, 2), dtype=np.int32))
        with pytest.raises(ValueError, match="all subjects have zero raters"):
            BalanceProfile.from_nij_matrix(nij)


class TestBalanceProfileSummary:
    """BalanceProfile.summary() produces plain-type dicts."""

    def test_balanced_summary(self):
        """Balanced summary includes raters_per_subject, not subject_counts."""
        bp = _balanced_profile()
        s = bp.summary()
        assert s["is_balanced"] is True
        assert s["n_subjects"] == 4
        assert s["n_raters"] == 3
        assert s["n_raters_exact"] is True
        assert s["n_raters_lower_bound"] == 3
        assert s["raters_per_subject"] == 3
        assert "subject_counts" not in s

    def test_unbalanced_summary(self):
        """Unbalanced summary includes subject_counts, not raters_per_subject."""
        bp = _unbalanced_profile()
        s = bp.summary()
        assert s["is_balanced"] is False
        assert s["n_raters_lower_bound"] == 4
        assert s["subject_counts"] == {2: 2, 3: 1, 4: 1}
        assert "raters_per_subject" not in s

    def test_summary_from_nij_omits_unknown_fields(self, balanced_nij):
        """Summary from nij omits n_raters and fill_rate (both None)."""
        bp = BalanceProfile.from_nij_matrix(balanced_nij)
        s = bp.summary()
        assert "n_raters" not in s
        assert "fill_rate" not in s
        assert s["n_raters_exact"] is False
        assert s["n_raters_lower_bound"] == 3


# ============================================================================ #
#                           RaterSharedMatrix                                  #
# ============================================================================ #


class TestRaterSharedMatrixConstruction:
    """How RaterSharedMatrix validates and stores its matrix."""

    def test_valid(self):
        """A valid symmetric integer matrix is accepted."""
        m = RaterSharedMatrix(SIMPLE_SHARED.copy())
        assert m.n_raters == 3

    def test_coerces_list_to_ndarray(self):
        """A nested list is coerced to ndarray via np.asarray."""
        m = RaterSharedMatrix([[10, 4], [4, 8]])
        assert isinstance(m.shared, np.ndarray)

    def test_frozen(self):
        """RaterSharedMatrix is immutable."""
        m = RaterSharedMatrix(SIMPLE_SHARED.copy())
        with pytest.raises(AttributeError):
            m.shared = np.eye(3, dtype=int)


class TestRaterSharedMatrixValidationRejects:
    """RaterSharedMatrix rejects invalid matrices."""

    def test_not_square(self):
        """Matrix must be square."""
        with pytest.raises(ValueError, match="square matrix"):
            RaterSharedMatrix(np.array([[1, 2, 3], [3, 2, 1]], dtype=int))

    def test_too_small(self):
        """At least 2 raters required -- one rater has no pairs."""
        with pytest.raises(ValueError, match="at least 2 raters"):
            RaterSharedMatrix(np.array([[5]], dtype=int))

    def test_float_dtype(self):
        """Counts must be integers, not floats."""
        with pytest.raises(TypeError, match="integers"):
            RaterSharedMatrix(np.array([[10.0, 4.0], [4.0, 8.0]]))

    def test_not_symmetric(self):
        """How many subjects two raters share cannot depend on order."""
        with pytest.raises(ValueError, match="symmetric"):
            RaterSharedMatrix(np.array([[10, 4], [5, 8]], dtype=int))

    def test_negative_count(self):
        """Counts of subjects cannot be negative."""
        with pytest.raises(ValueError, match="negative count"):
            RaterSharedMatrix(np.array([[10, -1], [-1, 8]], dtype=int))

    def test_off_diag_exceeds_diagonal(self):
        """A pair cannot share more subjects than either member rated."""
        # rater 0 rated 5 subjects, but allegedly shares 6 with rater 1.
        with pytest.raises(ValueError, match="sharing"):
            RaterSharedMatrix(np.array([[5, 6], [6, 8]], dtype=int))

    def test_error_prefix(self):
        """Every error starts with [RaterSharedMatrix] Error."""
        with pytest.raises(ValueError, match=r"\[RaterSharedMatrix\] Error"):
            RaterSharedMatrix(np.array([[5]], dtype=int))


class TestRaterSharedMatrixProperties:
    """Computed properties of RaterSharedMatrix."""

    def test_n_raters(self, simple_shared):
        """n_raters is the side length of the matrix."""
        assert simple_shared.n_raters == 3

    def test_n_pairs(self, simple_shared):
        """n_pairs = n(n-1)/2 = 3."""
        assert simple_shared.n_pairs == 3

    def test_per_rater(self, simple_shared):
        """per_rater is the diagonal: each rater's own total."""
        np.testing.assert_array_equal(simple_shared.per_rater, [10, 8, 9])

    def test_per_rater_is_copy(self, simple_shared):
        """per_rater returns a copy, not a view into shared."""
        pr = simple_shared.per_rater
        pr[0] = 999
        np.testing.assert_array_equal(simple_shared.per_rater, [10, 8, 9])


class TestRaterSharedMatrixPublicAPI:
    """Public methods: value_for, sorted_pairs, thinnest_pair, below, summary."""

    def test_value_for(self, simple_shared):
        """value_for returns the shared count for a pair."""
        assert simple_shared.value_for(0, 1) == 4
        assert simple_shared.value_for(0, 2) == 1
        assert simple_shared.value_for(1, 2) == 2

    def test_value_for_self_raises(self, simple_shared):
        """Asking for a rater's overlap with itself is not a pair."""
        with pytest.raises(ValueError, match="not a pair with itself"):
            simple_shared.value_for(0, 0)

    def test_sorted_pairs(self, simple_shared):
        """Pairs sorted thinnest first, a < b within each pair."""
        pairs = simple_shared.sorted_pairs()
        assert len(pairs) == 3
        assert pairs[0] == (0, 2, 1)   # thinnest
        assert pairs[1] == (1, 2, 2)
        assert pairs[2] == (0, 1, 4)   # thickest

    def test_thinnest_pair(self, simple_shared):
        """Thinnest pair is (0, 2) with 1 shared subject."""
        assert simple_shared.thinnest_pair() == (0, 2, 1)

    def test_below_default_threshold(self, simple_shared):
        """below() with default threshold (MIN_OVERLAP_HARD=2)."""
        # Only pair with < 2 shared subjects: (0, 2, 1).
        result = simple_shared.below()
        assert result == [(0, 2, 1)]

    def test_below_nothing_below(self, simple_shared):
        """below(threshold=1) returns empty when no pair is below."""
        result = simple_shared.below(threshold=1)
        assert result == []

    def test_summary_keys(self, simple_shared):
        """summary() returns all expected plain-type keys."""
        s = simple_shared.summary()
        assert s["n_raters"] == 3
        assert s["n_pairs"] == 3
        assert s["per_rater"] == [10, 8, 9]
        assert len(s["pairs"]) == 3


class TestRaterSharedMatrixRender:
    """RaterSharedMatrix.render() display methods."""

    def test_render_grid(self, simple_shared):
        """Grid render contains header and rater count."""
        text = simple_shared.render(style="grid")
        assert "RaterSharedMatrix" in text
        assert "3 rater(s)" in text

    def test_render_with_names(self, simple_shared):
        """Rater names appear in rendered output."""
        text = simple_shared.render(rater_names=["ann", "bob", "cara"])
        assert "ann" in text
        assert "bob" in text

    def test_render_invalid_style_raises(self, simple_shared):
        """Unknown style is rejected."""
        with pytest.raises(ValueError, match="style must be one of"):
            simple_shared.render(style="sparkline")

    def test_render_wrong_name_count_raises(self, simple_shared):
        """Name count must match rater count."""
        with pytest.raises(ValueError, match="rater name"):
            simple_shared.render(rater_names=["ann", "bob"])


class TestRaterSharedMatrixEquality:
    """Equality and hashing for RaterSharedMatrix."""

    def test_equal(self):
        """Two matrices with the same data are equal."""
        a = RaterSharedMatrix(SIMPLE_SHARED.copy())
        b = RaterSharedMatrix(SIMPLE_SHARED.copy())
        assert a == b

    def test_not_equal(self):
        """Different data produces not-equal."""
        a = RaterSharedMatrix(SIMPLE_SHARED.copy())
        other = SIMPLE_SHARED.copy()
        other[0, 1] = 3
        other[1, 0] = 3
        b = RaterSharedMatrix(other)
        assert a != b

    def test_not_hashable(self):
        """__hash__ is None, so hash() raises TypeError."""
        m = RaterSharedMatrix(SIMPLE_SHARED.copy())
        with pytest.raises(TypeError):
            hash(m)

    def test_not_equal_to_other_type(self):
        """Comparison with a non-RaterSharedMatrix returns NotImplemented."""
        m = RaterSharedMatrix(SIMPLE_SHARED.copy())
        assert m != "not a matrix"
