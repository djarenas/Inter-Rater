"""
Tests for the estimator_parts sub-package: NijMatrix, Weights, EstimatorInput.

Ordered bottom-up by dependency. NijMatrix and Weights have no internal
dependencies; EstimatorInput depends on both.

NijMatrix is the counts matrix: (n_subjects, n_categories) of non-negative
integers. Key behaviours: construction validation (2-D, integer, non-negative),
properties, DataFrame round-trip, frozen.

Weights is a K x K symmetric matrix of credits in [0, 1] with 1.0 on the
diagonal. Key behaviours: construction validation, factory methods (identity,
linear, quadratic, from_type), is_identity property, equality.

EstimatorInput carries category_index + K + optional weights + optional
subject_ids. Key behaviours: O(1) validation, properties, nij computation,
take_subjects, select_raters, equality.
"""

import numpy as np
import pandas as pd
import pytest

from interrater.contracts import EMPTY, INDEX_DTYPE
from interrater.intermediates.estimator_parts.nij_matrix import NijMatrix
from interrater.intermediates.estimator_parts.weights import Weights
from interrater.intermediates.estimator_parts.estimator_input import EstimatorInput

from .conftest import BALANCED_NIJ, UNBALANCED_NIJ, SIMPLE_INDEX, SIMPLE_K


# ============================================================================ #
#                                 NijMatrix                                    #
# ============================================================================ #


class TestNijMatrixConstruction:
    """How NijMatrix validates and stores its matrix."""

    def test_balanced_valid(self, balanced_nij):
        """A balanced counts matrix is accepted and stored."""
        np.testing.assert_array_equal(balanced_nij.matrix, BALANCED_NIJ)

    def test_unbalanced_valid(self, unbalanced_nij):
        """An unbalanced counts matrix is accepted and stored."""
        np.testing.assert_array_equal(unbalanced_nij.matrix, UNBALANCED_NIJ)

    def test_accepts_int8(self):
        """int8 dtype is accepted."""
        m = np.array([[1, 0], [0, 1]], dtype=np.int8)
        nij = NijMatrix(m)
        np.testing.assert_array_equal(nij.matrix, m)

    def test_accepts_int16(self):
        """int16 dtype is accepted."""
        m = np.array([[1, 0], [0, 1]], dtype=np.int16)
        nij = NijMatrix(m)
        np.testing.assert_array_equal(nij.matrix, m)

    def test_accepts_int64(self):
        """int64 dtype is accepted."""
        m = np.array([[1, 0], [0, 1]], dtype=np.int64)
        nij = NijMatrix(m)
        np.testing.assert_array_equal(nij.matrix, m)

    def test_all_zeros(self):
        """A matrix of all zeros is valid -- nobody rated anything."""
        m = np.zeros((2, 2), dtype=np.int32)
        nij = NijMatrix(m)
        np.testing.assert_array_equal(nij.matrix, m)

    def test_frozen(self, balanced_nij):
        """The matrix cannot be reassigned after construction (frozen dataclass)."""
        with pytest.raises(AttributeError):
            balanced_nij.matrix = np.zeros((4, 3), dtype=np.int32)


class TestNijMatrixValidationRejects:
    """Inputs that must not produce a valid NijMatrix."""

    def test_1d_array(self):
        """A 1-D array is not a matrix."""
        with pytest.raises(ValueError, match=r"\[NijMatrix\] Error.*2-dimensional"):
            NijMatrix(np.array([1, 2, 3], dtype=np.int32))

    def test_3d_array(self):
        """A 3-D array is not a matrix."""
        with pytest.raises(ValueError, match=r"\[NijMatrix\] Error.*2-dimensional"):
            NijMatrix(np.zeros((2, 3, 4), dtype=np.int32))

    def test_float_dtype(self):
        """Float dtype is rejected -- counts must be integers."""
        with pytest.raises(ValueError, match=r"\[NijMatrix\] Error.*integers"):
            NijMatrix(np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float64))

    def test_negative_counts(self):
        """Negative counts are rejected."""
        with pytest.raises(ValueError, match=r"\[NijMatrix\] Error.*negative"):
            NijMatrix(np.array([[1, -1], [0, 1]], dtype=np.int32))

    def test_error_prefix(self):
        """All errors carry the [NijMatrix] prefix."""
        with pytest.raises(ValueError, match=r"^\[NijMatrix\] Error"):
            NijMatrix(np.array([1, 2], dtype=np.int32))


class TestNijMatrixProperties:
    """Computed properties that describe the shape of the data."""

    def test_n_subjects_balanced(self, balanced_nij):
        assert balanced_nij.n_subjects == 4

    def test_n_subjects_unbalanced(self, unbalanced_nij):
        assert unbalanced_nij.n_subjects == 4

    def test_n_categories_balanced(self, balanced_nij):
        assert balanced_nij.n_categories == 3

    def test_n_categories_unbalanced(self, unbalanced_nij):
        assert unbalanced_nij.n_categories == 2

    def test_n_per_subject_balanced(self, balanced_nij):
        """Balanced: every subject has the same number of raters."""
        expected = np.array([3, 3, 3, 3])
        np.testing.assert_array_equal(balanced_nij.n_per_subject, expected)

    def test_n_per_subject_unbalanced(self, unbalanced_nij):
        """Unbalanced: subjects have different numbers of raters."""
        expected = np.array([3, 2, 4, 2])
        np.testing.assert_array_equal(unbalanced_nij.n_per_subject, expected)


class TestNijMatrixDataFrame:
    """build_from_dataframe and to_dataframe."""

    def test_build_from_dataframe(self):
        """A DataFrame of integer counts produces a valid NijMatrix."""
        df = pd.DataFrame([[2, 1, 0], [0, 3, 0]], columns=["a", "b", "c"])
        nij = NijMatrix.build_from_dataframe(df)
        expected = np.array([[2, 1, 0], [0, 3, 0]], dtype=np.int64)
        np.testing.assert_array_equal(nij.matrix, expected)

    def test_build_from_dataframe_not_a_dataframe(self):
        """Passing a non-DataFrame raises TypeError."""
        with pytest.raises(TypeError, match=r"\[NijMatrix\] Error.*DataFrame"):
            NijMatrix.build_from_dataframe(np.array([[1, 0], [0, 1]]))

    def test_to_dataframe_shape(self, balanced_nij):
        """to_dataframe returns the right shape."""
        df = balanced_nij.to_dataframe()
        assert df.shape == (4, 3)

    def test_round_trip(self, balanced_nij):
        """to_dataframe -> build_from_dataframe recovers the same counts."""
        df = balanced_nij.to_dataframe()
        recovered = NijMatrix.build_from_dataframe(df)
        np.testing.assert_array_equal(recovered.matrix, balanced_nij.matrix)


class TestNijMatrixEquality:
    """NijMatrix declares eq=False, so identity semantics apply."""

    def test_same_object_is_equal(self, balanced_nij):
        """An object is equal to itself."""
        assert balanced_nij == balanced_nij

    def test_different_objects_not_equal(self):
        """Two NijMatrix objects with the same data are not equal (eq=False)."""
        a = NijMatrix(BALANCED_NIJ.copy())
        b = NijMatrix(BALANCED_NIJ.copy())
        assert a != b


# ============================================================================ #
#                                  Weights                                     #
# ============================================================================ #


class TestWeightsConstruction:
    """How Weights validates and stores its matrix."""

    def test_valid_custom(self):
        """A valid symmetric matrix with diagonal=1.0 is accepted."""
        m = np.array([[1.0, 0.5], [0.5, 1.0]])
        w = Weights(m, kind="custom")
        np.testing.assert_array_equal(w.matrix, m)
        assert w.kind == "custom"

    def test_coerces_to_float64(self):
        """Input is coerced to float64 regardless of input dtype."""
        m = np.array([[1, 0], [0, 1]], dtype=np.int32)
        w = Weights(m, kind="custom")
        assert w.matrix.dtype == np.float64

    def test_default_kind_is_custom(self):
        """kind defaults to 'custom' when not specified."""
        w = Weights(np.eye(2))
        assert w.kind == "custom"

    def test_frozen(self):
        """Attributes cannot be reassigned after construction."""
        w = Weights(np.eye(2))
        with pytest.raises(AttributeError):
            w.matrix = np.eye(3)


class TestWeightsValidationRejects:
    """Inputs that must not produce valid Weights."""

    def test_not_square(self):
        """A rectangular matrix is not a weight matrix."""
        with pytest.raises(ValueError, match=r"\[Weights\] Error.*square"):
            Weights(np.ones((2, 3)))

    def test_1d(self):
        """A 1-D array is not a matrix."""
        with pytest.raises(ValueError, match=r"\[Weights\] Error.*square"):
            Weights(np.array([1.0, 0.5]))

    def test_fewer_than_2_categories(self):
        """A 1x1 matrix is rejected -- needs at least 2 categories."""
        with pytest.raises(ValueError, match=r"\[Weights\] Error.*2 categories"):
            Weights(np.array([[1.0]]))

    def test_not_symmetric(self):
        """An asymmetric matrix is rejected."""
        m = np.array([[1.0, 0.3], [0.5, 1.0]])
        with pytest.raises(ValueError, match=r"\[Weights\] Error.*symmetric"):
            Weights(m)

    def test_diagonal_not_one(self):
        """A diagonal that is not 1.0 is rejected."""
        m = np.array([[0.9, 0.5], [0.5, 1.0]])
        with pytest.raises(ValueError, match=r"\[Weights\] Error.*diagonal"):
            Weights(m)

    def test_value_below_zero(self):
        """Values below 0.0 are rejected."""
        m = np.array([[1.0, -0.1], [-0.1, 1.0]])
        with pytest.raises(ValueError, match=r"\[Weights\] Error.*\[0, 1\]"):
            Weights(m)

    def test_value_above_one(self):
        """Values above 1.0 are rejected."""
        m = np.array([[1.0, 1.1], [1.1, 1.0]])
        with pytest.raises(ValueError, match=r"\[Weights\] Error.*\[0, 1\]"):
            Weights(m)

    def test_invalid_kind(self):
        """An unknown kind is rejected."""
        with pytest.raises(ValueError, match=r"\[Weights\] Error.*kind"):
            Weights(np.eye(2), kind="banana")

    def test_error_prefix(self):
        """All errors carry the [Weights] prefix."""
        with pytest.raises(ValueError, match=r"^\[Weights\] Error"):
            Weights(np.array([1.0]))


class TestWeightsFactories:
    """Class methods that build standard weight matrices."""

    def test_identity(self):
        """identity(k) produces an eye matrix."""
        w = Weights.identity(3)
        np.testing.assert_array_equal(w.matrix, np.eye(3))
        assert w.kind == "identity"

    def test_linear_k3(self):
        """linear(3) gives known values: 1 - distance / max_distance."""
        w = Weights.linear(3)
        # distance matrix for 3 categories: [[0,1,2],[1,0,1],[2,1,0]]
        # max distance = 2
        # weights: 1 - d/2 = [[1.0, 0.5, 0.0], [0.5, 1.0, 0.5], [0.0, 0.5, 1.0]]
        expected = np.array([[1.0, 0.5, 0.0],
                             [0.5, 1.0, 0.5],
                             [0.0, 0.5, 1.0]])
        np.testing.assert_allclose(w.matrix, expected, atol=1e-12)
        assert w.kind == "linear"

    def test_quadratic_k3(self):
        """quadratic(3) gives known values: 1 - distance^2 / max_distance^2."""
        w = Weights.quadratic(3)
        # distance^2: [[0,1,4],[1,0,1],[4,1,0]]
        # max = 4
        # weights: 1 - d^2/4 = [[1.0, 0.75, 0.0], [0.75, 1.0, 0.75], [0.0, 0.75, 1.0]]
        expected = np.array([[1.0, 0.75, 0.0],
                             [0.75, 1.0, 0.75],
                             [0.0, 0.75, 1.0]])
        np.testing.assert_allclose(w.matrix, expected, atol=1e-12)
        assert w.kind == "quadratic"

    def test_linear_k2(self):
        """linear(2) with only two categories: neighbours get 0, same gets 1."""
        w = Weights.linear(2)
        expected = np.array([[1.0, 0.0], [0.0, 1.0]])
        np.testing.assert_allclose(w.matrix, expected, atol=1e-12)

    def test_quadratic_k2(self):
        """quadratic(2) with only two categories: same as identity."""
        w = Weights.quadratic(2)
        expected = np.array([[1.0, 0.0], [0.0, 1.0]])
        np.testing.assert_allclose(w.matrix, expected, atol=1e-12)

    def test_from_type_identity_returns_none(self):
        """from_type('identity', k) returns None -- identity is the default."""
        assert Weights.from_type("identity", 3) is None

    def test_from_type_linear(self):
        """from_type('linear', k) returns the same as linear(k)."""
        w = Weights.from_type("linear", 3)
        expected = Weights.linear(3)
        np.testing.assert_allclose(w.matrix, expected.matrix, atol=1e-12)

    def test_from_type_quadratic(self):
        """from_type('quadratic', k) returns the same as quadratic(k)."""
        w = Weights.from_type("quadratic", 3)
        expected = Weights.quadratic(3)
        np.testing.assert_allclose(w.matrix, expected.matrix, atol=1e-12)

    def test_from_type_custom_raises(self):
        """from_type('custom', k) raises -- pass a matrix to Weights() instead."""
        with pytest.raises(ValueError, match=r"does not support 'custom'"):
            Weights.from_type("custom", 3)

    def test_from_type_unknown_raises(self):
        """from_type with an unknown string raises."""
        with pytest.raises(ValueError, match=r"unknown weights_type"):
            Weights.from_type("banana", 3)


class TestWeightsProperties:
    """Computed properties."""

    def test_n_categories(self):
        """n_categories is the side length of the matrix."""
        w = Weights.identity(4)
        assert w.n_categories == 4

    def test_is_identity_true(self):
        """Identity matrix returns True."""
        assert Weights.identity(3).is_identity is True

    def test_is_identity_false_linear(self):
        """Linear weights are not identity."""
        assert Weights.linear(3).is_identity is False

    def test_is_identity_false_quadratic(self):
        """Quadratic weights are not identity."""
        assert Weights.quadratic(3).is_identity is False

    def test_is_identity_custom_eye(self):
        """A custom matrix that happens to be the eye is still identity."""
        w = Weights(np.eye(3), kind="custom")
        assert w.is_identity is True


class TestWeightsEquality:
    """Weights defines __eq__ via array_equal, __hash__ = None."""

    def test_equal_matrices(self):
        """Two Weights with the same matrix are equal."""
        a = Weights(np.eye(2), kind="custom")
        b = Weights(np.eye(2), kind="identity")
        assert a == b  # equality is on the matrix, not kind

    def test_different_matrices(self):
        """Two Weights with different matrices are not equal."""
        a = Weights.identity(3)
        b = Weights.linear(3)
        assert a != b

    def test_not_hashable(self):
        """Weights is unhashable (mutable equality semantics)."""
        with pytest.raises(TypeError):
            hash(Weights.identity(2))


# ============================================================================ #
#                              EstimatorInput                                  #
# ============================================================================ #


class TestEstimatorInputConstruction:
    """How EstimatorInput validates and stores its parts."""

    def test_minimal_valid(self, simple_input):
        """A valid category_index + K is accepted."""
        np.testing.assert_array_equal(simple_input.category_index, SIMPLE_INDEX)
        assert simple_input.n_categories == SIMPLE_K

    def test_with_weights(self):
        """Weights matching K are accepted."""
        w = Weights.linear(2)
        data = EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=w)
        assert data.weights is w

    def test_with_subject_ids(self):
        """Subject ids matching row count are accepted."""
        ids = np.array(["a", "b", "c", "d"], dtype=object)
        data = EstimatorInput(SIMPLE_INDEX.copy(), 2, subject_ids=ids)
        np.testing.assert_array_equal(data.subject_ids, ids)

    def test_weights_default_none(self, simple_input):
        """Weights default to None (identity, not missing)."""
        assert simple_input.weights is None

    def test_subject_ids_default_none(self, simple_input):
        """Subject ids default to None."""
        assert simple_input.subject_ids is None

    def test_with_empties(self):
        """EMPTY values (-1) are valid category_index entries."""
        index = np.array([[0, EMPTY], [1, 0]], dtype=INDEX_DTYPE)
        data = EstimatorInput(index, 2)
        assert data.n_subjects == 2

    def test_frozen(self, simple_input):
        """Attributes cannot be reassigned after construction."""
        with pytest.raises(AttributeError):
            simple_input.category_index = np.zeros((4, 2), dtype=INDEX_DTYPE)


class TestEstimatorInputValidationRejects:
    """Inputs that must not produce a valid EstimatorInput."""

    def test_not_numpy(self):
        """A plain list is rejected."""
        with pytest.raises(TypeError, match=r"\[EstimatorInput\] Error.*numpy"):
            EstimatorInput([[0, 1], [1, 0]], 2)

    def test_1d_array(self):
        """A 1-D array is rejected."""
        with pytest.raises(ValueError, match=r"\[EstimatorInput\] Error.*2-dimensional"):
            EstimatorInput(np.array([0, 1], dtype=INDEX_DTYPE), 2)

    def test_float_dtype(self):
        """Float dtype is rejected -- categories are integer positions."""
        with pytest.raises(TypeError, match=r"\[EstimatorInput\] Error.*integers"):
            EstimatorInput(np.array([[0.0, 1.0]], dtype=np.float64), 2)

    def test_single_rater(self):
        """Fewer than MIN_RATERS columns is rejected."""
        with pytest.raises(ValueError, match=r"\[EstimatorInput\] Error.*2 rater"):
            EstimatorInput(np.array([[0], [1]], dtype=INDEX_DTYPE), 2)

    def test_n_categories_bool_trap(self):
        """Passing True for n_categories is rejected (bool trap)."""
        with pytest.raises(TypeError, match=r"\[EstimatorInput\] Error.*integer"):
            EstimatorInput(SIMPLE_INDEX.copy(), True)

    def test_n_categories_too_small(self):
        """n_categories < 2 is rejected."""
        with pytest.raises(ValueError, match=r"\[EstimatorInput\] Error.*at least 2"):
            EstimatorInput(SIMPLE_INDEX.copy(), 1)

    def test_weights_not_weights_object(self):
        """Passing a raw matrix instead of Weights is rejected."""
        with pytest.raises(TypeError, match=r"\[EstimatorInput\] Error.*Weights"):
            EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=np.eye(2))

    def test_weights_k_mismatch(self):
        """Weights covering a different K are rejected."""
        w = Weights.identity(3)  # K=3
        with pytest.raises(ValueError, match=r"\[EstimatorInput\] Error.*categor"):
            EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=w)

    def test_subject_ids_length_mismatch(self):
        """Subject ids with wrong length are rejected."""
        ids = np.array(["a", "b"], dtype=object)  # 2 ids for 4 rows
        with pytest.raises(ValueError, match=r"\[EstimatorInput\] Error.*subject"):
            EstimatorInput(SIMPLE_INDEX.copy(), 2, subject_ids=ids)

    def test_error_prefix(self):
        """All errors carry the [EstimatorInput] prefix."""
        with pytest.raises(TypeError, match=r"^\[EstimatorInput\] Error"):
            EstimatorInput("not an array", 2)


class TestEstimatorInputProperties:
    """Computed properties."""

    def test_n_subjects(self, simple_input):
        assert simple_input.n_subjects == 4

    def test_n_raters(self, simple_input):
        assert simple_input.n_raters == 2

    def test_is_weighted_false_default(self, simple_input):
        """No weights means unweighted."""
        assert simple_input.is_weighted is False

    def test_is_weighted_false_identity(self):
        """Identity weights are still unweighted."""
        w = Weights.identity(2)
        data = EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=w)
        assert data.is_weighted is False

    def test_is_weighted_true_linear(self):
        """Non-identity weights make is_weighted True."""
        # K=3 linear weights are NOT identity, so is_weighted is True.
        idx = np.array([[0, 1], [2, 0], [1, 2]], dtype=INDEX_DTYPE)
        w = Weights.linear(3)
        data = EstimatorInput(idx, 3, weights=w)
        assert data.is_weighted is True


class TestEstimatorInputNij:
    """The nij property computes counts from the category_index."""

    def test_nij_shape(self, simple_input):
        """nij has shape (n_subjects, n_categories)."""
        nij = simple_input.nij
        assert nij.matrix.shape == (4, 2)

    def test_nij_values(self, simple_input):
        """Hand-computed nij for the SIMPLE dataset.

        s0: both raters chose 0 -> [2, 0]
        s1: both raters chose 1 -> [0, 2]
        s2: rater0=0, rater1=1  -> [1, 1]
        s3: rater0=1, rater1=0  -> [1, 1]
        """
        expected = np.array([[2, 0], [0, 2], [1, 1], [1, 1]], dtype=np.int32)
        np.testing.assert_array_equal(simple_input.nij.matrix, expected)

    def test_nij_returns_nij_matrix(self, simple_input):
        """The nij property returns an NijMatrix instance."""
        assert isinstance(simple_input.nij, NijMatrix)


class TestEstimatorInputTakeSubjects:
    """take_subjects returns a subset of rows, carrying K and weights."""

    def test_basic_subset(self, simple_input):
        """Selecting rows 0 and 2 gives a 2-row result."""
        subset = simple_input.take_subjects(np.array([0, 2]))
        assert subset.n_subjects == 2
        np.testing.assert_array_equal(
            subset.category_index,
            SIMPLE_INDEX[np.array([0, 2])],
        )

    def test_carries_k(self, simple_input):
        """K is preserved -- it must not be re-derived from the subset."""
        subset = simple_input.take_subjects(np.array([0]))
        assert subset.n_categories == SIMPLE_K

    def test_carries_weights(self):
        """Weights travel with the subset."""
        w = Weights.linear(2)
        data = EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=w)
        subset = data.take_subjects(np.array([0, 1]))
        assert subset.weights is w

    def test_carries_subject_ids(self):
        """Subject ids are subsetted along with rows."""
        ids = np.array(["a", "b", "c", "d"], dtype=object)
        data = EstimatorInput(SIMPLE_INDEX.copy(), 2, subject_ids=ids)
        subset = data.take_subjects(np.array([1, 3]))
        np.testing.assert_array_equal(subset.subject_ids, ["b", "d"])

    def test_repeated_rows(self, simple_input):
        """Repeated indices produce repeated rows (bootstrap draw)."""
        draw = simple_input.take_subjects(np.array([0, 0, 1]))
        assert draw.n_subjects == 3
        np.testing.assert_array_equal(
            draw.category_index[0], draw.category_index[1]
        )

    def test_no_shared_memory(self, simple_input):
        """The subset does not share memory with the original."""
        subset = simple_input.take_subjects(np.array([0, 1]))
        assert not np.shares_memory(
            subset.category_index, simple_input.category_index
        )


class TestEstimatorInputSelectRaters:
    """select_raters returns a subset of columns, carrying K and weights."""

    def test_basic_select(self, simple_input):
        """Selecting one column (of two) gives a 1-column... wait, that's
        below MIN_RATERS. Select both to verify the mechanic."""
        # Select both raters in reversed order.
        subset = simple_input.select_raters(np.array([1, 0]))
        assert subset.n_raters == 2
        # Columns are reversed.
        np.testing.assert_array_equal(
            subset.category_index[:, 0], SIMPLE_INDEX[:, 1]
        )

    def test_carries_k(self, simple_input):
        """K is preserved -- a subset that misses a category keeps the study's K."""
        subset = simple_input.select_raters(np.array([0, 1]))
        assert subset.n_categories == SIMPLE_K

    def test_carries_weights(self):
        """Weights travel with the rater subset."""
        w = Weights.linear(2)
        data = EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=w)
        subset = data.select_raters(np.array([0, 1]))
        assert subset.weights is w

    def test_carries_subject_ids(self):
        """Subject ids are preserved -- only columns change."""
        ids = np.array(["a", "b", "c", "d"], dtype=object)
        data = EstimatorInput(SIMPLE_INDEX.copy(), 2, subject_ids=ids)
        subset = data.select_raters(np.array([0, 1]))
        np.testing.assert_array_equal(subset.subject_ids, ids)


class TestEstimatorInputEquality:
    """EstimatorInput defines __eq__ by comparing all fields."""

    def test_equal_minimal(self):
        """Two inputs with same index and K are equal."""
        a = EstimatorInput(SIMPLE_INDEX.copy(), 2)
        b = EstimatorInput(SIMPLE_INDEX.copy(), 2)
        assert a == b

    def test_equal_with_weights(self):
        """Two inputs with same index, K, and weights are equal."""
        w = Weights.linear(2)
        a = EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=w)
        b = EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=w)
        assert a == b

    def test_equal_with_subject_ids(self):
        """Two inputs with same index, K, and subject_ids are equal."""
        ids = np.array(["a", "b", "c", "d"], dtype=object)
        a = EstimatorInput(SIMPLE_INDEX.copy(), 2, subject_ids=ids.copy())
        b = EstimatorInput(SIMPLE_INDEX.copy(), 2, subject_ids=ids.copy())
        assert a == b

    def test_different_k(self):
        """Different K makes them unequal."""
        a = EstimatorInput(SIMPLE_INDEX.copy(), 2)
        b = EstimatorInput(SIMPLE_INDEX.copy(), 3)
        assert a != b

    def test_different_index(self):
        """Different category_index makes them unequal."""
        other = SIMPLE_INDEX.copy()
        other[0, 0] = 1
        a = EstimatorInput(SIMPLE_INDEX.copy(), 2)
        b = EstimatorInput(other, 2)
        assert a != b

    def test_one_has_weights(self):
        """One with weights and one without are unequal."""
        a = EstimatorInput(SIMPLE_INDEX.copy(), 2)
        b = EstimatorInput(SIMPLE_INDEX.copy(), 2, weights=Weights.linear(2))
        assert a != b

    def test_one_has_subject_ids(self):
        """One with subject_ids and one without are unequal."""
        ids = np.array(["a", "b", "c", "d"], dtype=object)
        a = EstimatorInput(SIMPLE_INDEX.copy(), 2)
        b = EstimatorInput(SIMPLE_INDEX.copy(), 2, subject_ids=ids)
        assert a != b

    def test_not_hashable(self):
        """EstimatorInput is unhashable."""
        with pytest.raises(TypeError):
            hash(EstimatorInput(SIMPLE_INDEX.copy(), 2))

    def test_not_equal_to_other_type(self):
        """Comparing with a non-EstimatorInput returns NotImplemented."""
        data = EstimatorInput(SIMPLE_INDEX.copy(), 2)
        assert data != "not an EstimatorInput"
