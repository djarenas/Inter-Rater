"""The data object: a validated array of ratings and what it means."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .. import array_utils as _a
from ..contracts import INDEX_DTYPE, LABEL_DTYPE
from ..semantics.ratings_semantics import RatingsSemantics
from . import _ratings_validation as _v


class Ratings:
    """
    I am a validated array of ratings. If I exist, I am complete.

    I hold three things: a dense ``(n_subjects, n_raters)`` array of category
    positions, a label for each row, and the semantics that say what all of it
    means. That is the whole of me. I do not clean, analyse, transform, or
    count -- a numpy array does not know what it is, and my job is to be the
    array that does.

    My constructor raises rather than repair. Anything holding me may assume my
    values are in range, my shape matches what my semantics declare, and no
    declared rater of mine sat idle.

    I describe my contents, not my history. Nothing about what a cleaner
    discarded to make me survives in me; that lives in the cleaner report, alongside
    the config that produced it.

    I compute nothing about myself either. How skewed my categories are,
    whether every subject was rated the same number of times, which raters
    overlap -- those are facts *about* me, not part of me. Keeping them out
    means I stay cheap to build; this is a deliberate choice as it matters when a bootstrapper is making
    thousands of arrays a second.

    I allow empties freely. Raters need not have rated the same subjects, which
    is ordinary in some use cases (such as systematic review) where each subject goes to a
    subset of raters. An empty says only that no rating exists there, never
    why: 'unsure' and 'abstain' are categories the rating task offered and a
    rater chose, and they count like any other answer.

    Attributes
    ----------
    category_index:
        (n_subjects, n_raters) int8, read-only. Each cell is a position in
        ``semantics.categories.names``, or EMPTY. Column order is rater order
        as declared; I never sort.
    subject_ids:
        (n_subjects,) object. One label per row, in input order. Generated as
        ``subject0, subject1, ...`` when the data had none, since a row must be
        nameable even if the file did not name it. May repeat only if the
        semantics allow it.
    semantics:
        Who the raters are, what the categories mean, and what the data was
        expected to look like.

    Examples
    --------
    >>> import numpy as np
    >>> from interrater.semantics import RatingsSemantics, UnorderedCategories
    >>> sem = RatingsSemantics(
    ...     raters=["ann", "bob"],
    ...     categories=UnorderedCategories(("include", "exclude")),
    ... )
    >>> r = Ratings(np.array([[0, 0], [1, 1], [0, 1]], dtype=np.int8), sem)
    >>> r.n_subjects, r.n_raters
    (3, 2)
    >>> r.subject_ids[0]
    'subject0'
    """

    # Attributes -------------------------------------------------------------  #

    category_index: np.ndarray
    # (n_subjects, n_raters) int8, C-contiguous, read-only after construction.
    # A cell is a position in semantics.categories.names, or EMPTY (-1).
    # K is owned by the semantics and is never inferred from these values:
    # category_index.max() + 1 is only ever a lower bound on K.

    subject_ids:    np.ndarray
    # (n_subjects,) object, holding real Python strings so long ids survive.
    # Row order is meaningful -- it is what distinguishes a bootstrap draw from
    # the sample it came from.

    semantics:      RatingsSemantics
    # The declaration my contents were checked against.

    # Constructor ------------------------------------------------------------  #

    def __init__(
        self,
        category_index: np.ndarray,
        semantics:      RatingsSemantics,
        subject_ids:    Sequence[object] | None = None,
    ) -> None:
        name = type(self).__name__

        self.category_index = np.ascontiguousarray(category_index)
        self.semantics      = semantics

        _v.validate_dtype(self.category_index, name)
        self.category_index = self.category_index.astype(INDEX_DTYPE, copy=False)
        _v.validate_shape(self.category_index, name)

        self.subject_ids = self._as_subject_ids(
            subject_ids, self.category_index.shape[0]
        )

        self._validate()
        self._freeze()

    @staticmethod
    def _as_subject_ids(
        subject_ids: Sequence[object] | None, n_rows: int
    ) -> np.ndarray:
        """Absent ids become positional ones; given ids are taken as given."""
        if subject_ids is None:
            return np.array(
                [f"subject{i}" for i in range(n_rows)], dtype=LABEL_DTYPE
            )
        return np.asarray(subject_ids, dtype=LABEL_DTYPE)

    def _validate(self) -> None:
        """
        Every check I make, in order.

        Each is told the concrete class name, so a PairedRatings failing one
        of these says PairedRatings. The check living on me is an
        implementation detail the caller has no reason to know.
        """
        name = type(self).__name__
        _v.validate_matches_raters(self.category_index, self.semantics, name)
        _v.validate_matches_n_subjects(self.category_index, self.semantics, name)
        _v.validate_values_in_categories(self.category_index, self.semantics, name)
        _v.validate_subject_ids_match_shape(self.category_index, self.subject_ids, name)
        _v.validate_subject_ids_unique(self.subject_ids, self.semantics, name)
        _v.validate_no_idle_raters(self.category_index, self.semantics, name)
        _v.validate_no_unrated_subjects(self.category_index, self.subject_ids, name)

    def _freeze(self) -> None:
        """Read-only, so nothing computed from me can be invalidated later."""
        self.category_index = self.category_index.view()
        self.category_index.flags.writeable = False
        self.subject_ids = self.subject_ids.view()
        self.subject_ids.flags.writeable = False

    # Public API -------------------------------------------------------------  #

    @property
    def n_subjects(self) -> int:
        return int(self.category_index.shape[0])

    @property
    def n_raters(self) -> int:
        return int(self.category_index.shape[1])

    @property
    def n_categories(self) -> int:
        """K, from the semantics. Never inferred from my values."""
        return self.semantics.n_categories

    def check_my_balance(self) -> bool:
        """
        Was every subject rated the same *number* of times?

        Not by the same raters. Subject 1 may have gone to Ann and Bob while
        subject 2 went to Cara and Dev; what matters is that both were seen
        twice. This is what exact Fleiss' kappa asks for, and it is weaker
        than a complete rectangle.

        I compute this rather than hold it, because I hold only what I was
        given. Ask me once if you need to know; don't ask me repeatedly.
        """
        return _a.is_balanced(self.category_index)

    def check_my_crossing(self) -> bool:
        """
        Did every rater rate every subject? No empties anywhere.

        Implies balance. Beyond that it makes every pair fully overlapping, so
        pairwise coefficients rest on the same subjects and are comparable to
        one another rather than merely each valid on its own footing.
        """
        return _a.is_crossed(self.category_index)

    def check_my_pairing(self) -> bool:
        """
        Am I two raters with no gaps between them?

        The shape Cohen's kappa consumes. Two raters with gaps are not paired:
        a paired comparison is defined on the subjects both of them rated, so
        reducing to that overlap is part of what makes a pair rather than
        something done to one afterwards.
        """
        return self.n_raters == 2 and self.check_my_crossing()

    def to_nij_matrix(self):
        """
        I build the counts matrix from my ratings.

        For each subject, I count how many raters chose each category,
        ignoring empties. The result discards rater identity -- only
        tallies remain.

        Returns
        -------
        NijMatrix
            Shape (n_subjects, n_categories), entry (i, j) is how many
            of my raters chose category j for subject i.
        """
        from ..intermediates.estimator_parts.nij_matrix import NijMatrix

        K = self.n_categories
        N = self.n_subjects
        nij = np.zeros((N, K), dtype=np.int64)

        for i in range(N):
            row = self.category_index[i]
            valid = row[row != -1]
            if len(valid) > 0:
                nij[i] = np.bincount(valid, minlength=K)

        return NijMatrix(matrix=nij)

    # Comparison and display -------------------------------------------------  #

    def __eq__(self, other: object) -> bool:
        """
        Two of me are equal when we hold the same ratings and mean the same
        thing.

        Contents only. Where the data came from is not part of the comparison,
        because I do not carry that. Positional on both axes: rater order is
        declared order, and subject order distinguishes a bootstrap draw from
        the sample it came from.
        """
        if not isinstance(other, Ratings):
            return NotImplemented
        return (
            self.semantics == other.semantics
            and self.category_index.shape == other.category_index.shape
            and np.array_equal(self.category_index, other.category_index)
            and np.array_equal(self.subject_ids, other.subject_ids)
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        filled = int((self.category_index != -1).sum())
        total  = self.n_subjects * self.n_raters
        return (
            f"Ratings(\n"
            f"    n_subjects = {self.n_subjects:,}\n"
            f"    raters     = {list(self.semantics.raters)}\n"
            f"    categories = {list(self.semantics.categories.names)}\n"
            f"    filled     = {filled:,} / {total:,} ({filled / total:.1%})\n"
            f")"
        )
