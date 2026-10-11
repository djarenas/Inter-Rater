"""Standardized object to feed to estimators and resamplers"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from ...array_utils import category_counts
from ...contracts import MIN_RATERS
from .nij_matrix import NijMatrix
from .weights import Weights
_ERROR_PREFIX = "[EstimatorInput] Error"


@dataclass(frozen=True, eq=False)
class EstimatorInput:
    """
    I am everything a coefficient needs to be computed, in one object.

    An array of category positions, the number of categories they index into,
    and -- when a coefficient cares -- how far apart those categories are. I
    exist so that those three cannot be separated. 
    
    This abstraction serves to prevent a possible bug: if K were read off the data instead
    of carried alongside, a bootstrap draw that happens to miss a rare category would be computed
    against a smaller K than the rating task had, and AC1 and
    Brennan-Prediger, whose chance terms count categories, would move from
    bootstrap-draw to draw for no reason but the inference.

    I check almost nothing. Shape, dtype, that the weights match K -- all of it
    O(1). I never look at my contents: ``Ratings`` checked every value against
    the categories once, at the boundary, and a draw taken from a validated
    array is valid by construction. At 0.4 microseconds against a 40 microsecond row gather, 
    I am about 1% of a draw.

    I hold no names. A rater is a column index and a category is a position;
    what either is called matters only when a result is rendered, which is the
    analyser's job and happens once.

    Attributes
    ----------
    category_index:
        A numpy array of size n_subjects by n_raters, holding integers. Each
        entry is a position in the rating task's categories, or EMPTY.
    n_categories:
        K. From the semantics, never inferred from the array.
    weights:
        A ``Weights``, or None for identity.

        None is not a missing value. An unweighted coefficient *is* a weighted
        one with identity weights, so None is the honest default rather than a
        placeholder. It also keeps the fast path fast: unweighted observed
        agreement is one pass over the counts, while weighted agreement is a
        sum over every pair of categories, so a coefficient that knows it is
        unweighted can skip work rather than multiply by ones.
    subject_ids:
        A numpy array of length n_subjects, holding one label per row, or
        None. Carried only when something downstream
        reads them -- intra-rater work needs to know which rows are the same
        subject. None means the builder chose not to carry them, not that the
        data lacks them; gathering them costs about 44% on top of the codes
        gather, which is not worth paying on thousands of draws that will
        never look.

    Examples
    --------
    >>> import numpy as np
    >>> data = EstimatorInput(np.array([[0, 1], [1, 1]], dtype=np.int8), 2)
    >>> data.n_subjects, data.n_raters, data.n_categories
    (2, 2, 2)
    >>> data.take_subjects(np.array([0, 0, 1])).n_subjects
    3
    """

    # Attributes -------------------------------------------------------------  #

    category_index: np.ndarray
    n_categories:   int
    weights:        Weights | None = None
    subject_ids:    np.ndarray | None = None

    # Constructor ------------------------------------------------------------  #

    def __post_init__(self) -> None:
        self._validate()

    def _validate(self) -> None:
        """
        Shape, dtype, and that my parts belong together. Nothing more.

        Every check here is O(1) -- reading ``.ndim``, ``.shape`` and
        ``.dtype``, never the array itself. That is the whole reason I am cheap
        enough to build inside a resampling loop.
        """
        index = self.category_index

        if not isinstance(index, np.ndarray):
            raise TypeError(
                f"{_ERROR_PREFIX}: category_index must be a numpy array, got "
                f"{type(index).__name__}."
            )

        if index.ndim != 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: category_index must be 2-dimensional, "
                f"one row per subject and one column per rater, got shape "
                f"{index.shape}."
            )

        if index.dtype.kind not in "iu":
            raise TypeError(
                f"{_ERROR_PREFIX}: category_index must hold integers, got "
                f"dtype {index.dtype}. Category names are resolved to "
                f"positions by a cleaner long before anything reaches me."
            )

        if index.shape[1] < MIN_RATERS:
            raise ValueError(
                f"{_ERROR_PREFIX}: needs at least {MIN_RATERS} rater columns, "
                f"got {index.shape[1]}. Agreement is undefined for a single "
                f"rater."
            )

        if isinstance(self.n_categories, bool) or not isinstance(
            self.n_categories, (int, np.integer)
        ):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_categories must be an integer, got "
                f"{type(self.n_categories).__name__}."
            )

        if self.n_categories < 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_categories must be at least 2, got "
                f"{self.n_categories}. A single category cannot produce "
                f"agreement or disagreement."
            )

        if self.weights is not None:
            if not isinstance(self.weights, Weights):
                raise TypeError(
                    f"{_ERROR_PREFIX}: weights must be a Weights or None, got "
                    f"{type(self.weights).__name__}. Build one with "
                    f"Weights.quadratic(k), Weights.linear(k), or from a "
                    f"matrix directly."
                )
            if self.weights.n_categories != self.n_categories:
                raise ValueError(
                    f"{_ERROR_PREFIX}: weights cover "
                    f"{self.weights.n_categories} categor(ies) but "
                    f"n_categories is {self.n_categories}. Weights and K that "
                    f"disagree usually mean they came from different category "
                    f"sets."
                )

        if self.subject_ids is not None:
            if len(self.subject_ids) != index.shape[0]:
                raise ValueError(
                    f"{_ERROR_PREFIX}: given {len(self.subject_ids):,} subject "
                    f"id(s) for {index.shape[0]:,} row(s). Ids that do not "
                    f"line up with rows are worse than no ids at all."
                )

    # Public API -------------------------------------------------------------  #

    @property
    def n_subjects(self) -> int:
        return int(self.category_index.shape[0])

    @property
    def n_raters(self) -> int:
        return int(self.category_index.shape[1])

    @property
    def is_weighted(self) -> bool:
        """Whether a coefficient reading me would apply anything but identity."""
        return self.weights is not None and not self.weights.is_identity

    @property
    def nij(self) -> NijMatrix:
        """
        Me as counts: how many raters put each subject in each category.

        What every coefficient that treats raters as interchangeable actually
        reads. Computed on demand rather than held, for two reasons. A pairwise
        coefficient never wants it, so building it at construction would tax
        every draw of a bootstrap that was never going to look. And an analyser
        converts once and passes the result down, so there is no repeat access
        for a cache to save.

        **The conversion only runs one way.** Counts have no columns: which
        rater said what is discarded here and cannot be recovered. Resampling
        my rows is still a subject bootstrap, exactly, but leaving out a rater
        has to happen before this point, on me, while I still have raters.
        That is why an analyser converts only once it knows which kind of
        operation is ahead.
        """
        return NijMatrix(category_counts(self.category_index, self.n_categories))

    def take_subjects(self, rows: np.ndarray) -> EstimatorInput:
        """
        Me, restricted to the given rows, carrying K and weights along.

        The resampling primitive. ``rows`` may repeat an index, which is the
        point: drawing whole subject rows with replacement preserves the
        correlation between coefficients that share a rater without anything
        having to model it.

        Numpy's fancy indexing always copies, so the array I return shares no
        memory with mine and a draw cannot disturb what it was drawn from.

        Subject ids come along only if I am carrying them. Gathering them when
        nothing will read them is a tax on every draw.
        """
        return EstimatorInput(
            category_index = self.category_index[rows],
            n_categories   = self.n_categories,
            weights        = self.weights,
            subject_ids    = (
                None if self.subject_ids is None else self.subject_ids[rows]
            ),
        )

    def select_raters(self, keep: np.ndarray) -> EstimatorInput:
        """
        Me, restricted to the given rater columns.

        For leave-one-out and rater permutation. Note that what comes back may
        no longer be balanced, and on sparse data may hold subjects with too
        few ratings to contribute to anything -- I do not check, because the
        caller is the one who knows whether that matters and needs to report
        how much survived.
        """
        return EstimatorInput(
            category_index = self.category_index[:, keep],
            n_categories   = self.n_categories,
            weights        = self.weights,
            subject_ids    = self.subject_ids,
        )

    # Comparison and display -------------------------------------------------  #

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, EstimatorInput):
            return NotImplemented
        if self.n_categories != other.n_categories:
            return False
        if not np.array_equal(self.category_index, other.category_index):
            return False
        if (self.weights is None) != (other.weights is None):
            return False
        if self.weights is not None and not np.array_equal(
            self.weights, other.weights
        ):
            return False
        if (self.subject_ids is None) != (other.subject_ids is None):
            return False
        return self.subject_ids is None or np.array_equal(
            self.subject_ids, other.subject_ids
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        parts = [
            f"{self.n_subjects:,} x {self.n_raters}",
            f"K={self.n_categories}",
        ]
        if self.is_weighted:
            parts.append("weighted")
        if self.subject_ids is not None:
            parts.append("with ids")
        return f"EstimatorInput({', '.join(parts)})"
