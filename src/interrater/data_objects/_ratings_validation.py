"""
Validation for Ratings.

Everything here raises. Semantics declares what a rating task should look like but
never sees an array; Ratings holds both, so this is where a declaration meets
the data that is supposed to satisfy it.

Kept apart from the arithmetic so that a bootstrapper can import the maths
without dragging error machinery along.
"""

from __future__ import annotations

import numpy as np

from .. import array_utils as _a
from ..contracts import EMPTY, MIN_SUBJECTS
from ..semantics.ratings_semantics import RatingsSemantics

def _prefix(cls_name: str) -> str:
    """
    Name the class the caller asked for, not the one defining the check.

    A PairedRatings that fails a base check -- a bad dtype, a value outside the
    declared categories -- should say PairedRatings. It is what the user named,
    and the check happening to live on the parent is an implementation detail
    they have no reason to know.
    """
    return f"[{cls_name}] Error"


# The array on its own --------------------------------------------------------  #


def validate_dtype(category_index: np.ndarray, cls_name: str) -> None:
    if category_index.dtype.kind == "f":
        raise TypeError(
            f"{_prefix(cls_name)}: category_index must be an integer array, got dtype "
            f"{category_index.dtype}. Floats usually mean NaN was used for "
            f"absent ratings; use EMPTY ({EMPTY}) instead, which keeps the "
            f"array integral and indexable."
        )
    if category_index.dtype.kind not in "iu":
        raise TypeError(
            f"{_prefix(cls_name)}: category_index must be an integer array, got dtype "
            f"{category_index.dtype}. Category names are resolved to positions "
            f"by the cleaner; by the time an array reaches me it holds "
            f"integers."
        )


def validate_shape(category_index: np.ndarray, cls_name: str) -> None:
    if category_index.ndim != 2:
        raise ValueError(
            f"{_prefix(cls_name)}: category_index must be 2-dimensional "
            f"(n_subjects, n_raters), got shape {category_index.shape}."
        )
    if category_index.shape[0] < MIN_SUBJECTS:
        raise ValueError(
            f"{_prefix(cls_name)}: needs at least {MIN_SUBJECTS} subjects, got "
            f"{category_index.shape[0]}."
        )


# The array against the declaration -------------------------------------------  #


def validate_matches_raters(
    category_index: np.ndarray, semantics: RatingsSemantics, cls_name: str
) -> None:
    """The check that makes declaring a rating task worth the trouble."""
    found = category_index.shape[1]
    if found != semantics.n_raters:
        raise ValueError(
            f"{_prefix(cls_name)}: was given {found} rater column(s) but its semantics declare "
            f"{semantics.n_raters}: {list(semantics.raters)}. A column that "
            f"was not declared is a mistake to look into, not a rater to "
            f"invent; if the extra column is real, add it to the semantics."
        )


def validate_matches_n_subjects(
    category_index: np.ndarray, semantics: RatingsSemantics, cls_name: str
) -> None:
    if semantics.n_subjects is None:
        return
    found = category_index.shape[0]
    if found != semantics.n_subjects:
        raise ValueError(
            f"{_prefix(cls_name)}: was given {found:,} subject row(s) but its semantics expect "
            f"{semantics.n_subjects:,}. Either the file is not the one the "
            f"rating task describes, or a cleaner dropped rows it should not "
            f"have. "
            f"Set n_subjects to None in the semantics for no expectation."
        )


def validate_values_in_categories(
    category_index: np.ndarray, semantics: RatingsSemantics, cls_name: str
) -> None:
    k       = semantics.n_categories
    illegal = (category_index < EMPTY) | (category_index >= k)
    n_bad   = int(illegal.sum())
    if n_bad:
        examples = np.unique(category_index[illegal])[:3].tolist()
        raise ValueError(
            f"{_prefix(cls_name)}: category_index holds {n_bad:,} value(s) outside its "
            f"categories (e.g. {examples}). Valid values are 0..{k - 1}, "
            f"indexing {list(semantics.categories.names)}, or EMPTY ({EMPTY}). "
            f"Use a cleaner to map raw responses onto the categories first."
        )


# Subject ids -----------------------------------------------------------------  #


def validate_subject_ids_match_shape(
    category_index: np.ndarray, subject_ids: np.ndarray, cls_name: str
) -> None:
    if len(subject_ids) != category_index.shape[0]:
        raise ValueError(
            f"{_prefix(cls_name)}: has {category_index.shape[0]:,} subject row(s) but "
            f"{len(subject_ids):,} subject id(s). Every row needs a label."
        )


def validate_subject_ids_unique(
    subject_ids: np.ndarray, semantics: RatingsSemantics, cls_name: str
) -> None:
    if not semantics.unique_subjects:
        return
    uniques, counts = np.unique(subject_ids, return_counts=True)
    repeated = uniques[counts > 1]
    if repeated.size:
        raise ValueError(
            f"{_prefix(cls_name)}: has {repeated.size:,} subject id(s) appearing more than "
            f"once (e.g. {repeated[:3].tolist()}). That is usually an export "
            f"that duplicated rows. If the rating task deliberately re-shows "
            f"subjects, set unique_subjects=False in the semantics."
        )


# Nothing may be entirely absent ----------------------------------------------  #


def validate_no_idle_raters(
    category_index: np.ndarray, semantics: RatingsSemantics, cls_name: str
) -> None:
    rated = (category_index != EMPTY).sum(axis=0)
    idle  = np.flatnonzero(rated == 0)
    if idle.size:
        names = [semantics.raters[i] for i in idle[:3]]
        raise ValueError(
            f"{_prefix(cls_name)}: has {idle.size} declared rater(s) who rated nothing at all "
            f"(e.g. {names}). A rater with no ratings contributes to no "
            f"coefficient and makes per-rater rates undefined. Drop them from "
            f"the semantics, or find out why their column is empty."
        )


def validate_no_unrated_subjects(
    category_index: np.ndarray, subject_ids: np.ndarray, cls_name: str
) -> None:
    rated   = (category_index != EMPTY).sum(axis=1)
    unrated = np.flatnonzero(rated == 0)
    if unrated.size:
        raise ValueError(
            f"{_prefix(cls_name)}: has {unrated.size:,} subject(s) rated by nobody (e.g. "
            f"{subject_ids[unrated][:3].tolist()}). An all-empty row joins no "
            f"comparison; drop these in the cleaner."
        )


# What the narrower types add -------------------------------------------------  #
#
# These take the concrete class name, because a PairedRatings that fails the
# balance check should say PairedRatings rather than the class that happens to
# define the check.


def validate_balanced(category_index: np.ndarray, cls_name: str) -> None:
    """Every subject rated the same number of times."""
    counts   = _a.n_per_subject(category_index)
    distinct = np.unique(counts)
    if distinct.size != 1:
        offenders = np.flatnonzero(counts != distinct[0])[:3]
        raise ValueError(
            f"{_prefix(cls_name)}: every subject must be rated the same number "
            f"of times, but found {distinct.size} distinct counts "
            f"{distinct[:5].tolist()} (e.g. rows {offenders.tolist()}). Hold "
            f"this as Ratings and accept the approximated Fleiss' kappa, or "
            f"drop the odd subjects first."
        )


def validate_crossed(category_index: np.ndarray, cls_name: str) -> None:
    """Every rater rated every subject."""
    gaps = np.flatnonzero((category_index == EMPTY).any(axis=1))
    if gaps.size:
        raise ValueError(
            f"{_prefix(cls_name)}: every rater must have rated every subject, "
            f"but {gaps.size:,} subject(s) have at least one empty (e.g. rows "
            f"{gaps[:3].tolist()}). Hold this as BalancedRatings if the rater "
            f"count is constant, or as Ratings otherwise."
        )


def validate_exactly_two_raters(category_index: np.ndarray, cls_name: str) -> None:
    """Exactly two rater columns."""
    found = category_index.shape[1]
    if found != 2:
        raise ValueError(
            f"{_prefix(cls_name)}: exactly 2 raters are needed, got {found}. "
            f"A paired comparison is defined on two raters and the subjects "
            f"they share."
        )
