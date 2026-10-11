"""
Arithmetic over the presence of ratings.

Every function here reads one thing: whether a cell holds a rating. None of
them look at *which* category was chosen, with one exception noted below.

The rule that keeps this file from drifting is not presence-versus-values, it
is shared-versus-specific: a function belongs here when more than one layer
needs it, and beside its layer when only that layer does. Structure is asked by
a data object validating itself, by a profile reporting, and by an extractor
finding a pair. ``category_counts`` is the exception to the presence rule --
it reads values, but an estimator needs the same projection a profile does, and
one definition is better than two that can drift.

Takes arrays, returns arrays or scalars, imports nothing but numpy and
``contracts``. A function that needs to know what a rating *means* -- what a
category is called, whether a scale is ordered -- belongs with the layer that
knows. Prevalence and coverage read values but only the profile reports them,
so they live there.

Nothing here validates. Callers hold an array that already satisfies the
contract in ``interrater.contracts``.
"""

from __future__ import annotations

import numpy as np

from .contracts import EMPTY

# Presence --------------------------------------------------------------------  #


def present_mask(category_index: np.ndarray) -> np.ndarray:
    """
    True where a rating exists.

    The one primitive everything else here is built from. ``!= EMPTY`` rather
    than ``>= 0`` so that a second sentinel, if one is ever wanted, changes
    this line and nothing else.
    """
    return category_index != EMPTY


def n_per_subject(category_index: np.ndarray) -> np.ndarray:
    """
    How many raters rated each subject.

    Whether this is constant is exactly what ``is_balanced`` answers, and a
    constant is what exact Fleiss' kappa requires.
    """
    # Axis = 1 -> Sum over rows
    return present_mask(category_index).sum(axis=1).astype(np.int32)


def n_per_rater(category_index: np.ndarray) -> np.ndarray:
    """
    How many subjects each rater rated.

    The diagonal of ``pairwise_overlap``.
    """
    return present_mask(category_index).sum(axis=0).astype(np.int32)


def category_counts(category_index: np.ndarray, k: int) -> np.ndarray:
    """
    Calculates the n_ij matrix as a np.ndarray
    n_ij matrix: How many raters put each subject in each category. i: Subject index; j: Category index. 

    Rows (i) are subjects, columns are categories. 
    Empty ratings contribute nothing to n_ij; therefore, a row sums 
    to the number of raters who actually rated that subject, 
    not to the number of raters in the rating task.
    """
    n_subjects = category_index.shape[0]
    present    = present_mask(category_index)
    row_index  = np.repeat(np.arange(n_subjects), present.sum(axis=1))
    flat       = row_index * k + category_index[present].astype(np.int64)
    counts     = np.bincount(flat, minlength=n_subjects * k)
    return counts.reshape(n_subjects, k).astype(np.int32)


# Shared subjects -------------------------------------------------------------  #


def pairwise_overlap(category_index: np.ndarray) -> np.ndarray:
    """
    Co-rated subject counts, (n_raters, n_raters).

    Cell (a, b) is how many subjects both a and b rated; the diagonal is each
    rater's own total. Every pairwise coefficient rests on one of the
    off-diagonal counts, so a small entry means a coefficient that is valid
    and not worth much. Zero is allowed: disjoint assignment is legitimate,
    just uninformative.

    A matrix multiply rather than a loop over pairs -- the mask is boolean, so
    the product of a column and a row counts the rows where both are true.
    """
    present = present_mask(category_index).astype(np.int32)
    return present.T @ present


def shared_subject_mask(
    category_index: np.ndarray, rater_a: int, rater_b: int
) -> np.ndarray:
    """
    True for subjects both raters rated.

    The rows a paired comparison is defined on. Extracting a pair means taking
    these rows and those two columns, which is what leaves the extracted array
    with no empties at all.
    """
    present = present_mask(category_index)
    return present[:, rater_a] & present[:, rater_b]


# Structure -------------------------------------------------------------------  #


def is_balanced(category_index: np.ndarray) -> bool:
    """
    Whether every subject was rated by the same *number* of raters.

    Not the same raters. Subject 1 may have gone to Ann and Bob while subject
    2 went to Cara and Dev; what matters is that both were seen twice. This is
    the condition exact Fleiss' kappa actually asks for, and it is weaker than
    a complete rectangle -- which is why it is worth asking separately.
    """
    return len(np.unique(n_per_subject(category_index))) == 1


def is_crossed(category_index: np.ndarray) -> bool:
    """
    Whether every rater rated every subject. No empties anywhere.

    Implies balanced. Beyond that it makes every pair fully overlapping, so
    pairwise coefficients rest on the same subjects and are comparable to each
    other rather than merely each valid on its own footing.
    """
    return bool(present_mask(category_index).all())


def raters_per_subject(category_index: np.ndarray) -> int | None:
    """
    The constant every subject shares, or None when there is not one.

    The n that exact Fleiss' kappa needs.
    """
    counts = np.unique(n_per_subject(category_index))
    return int(counts[0]) if len(counts) == 1 else None


def subject_count_distribution(category_index: np.ndarray) -> dict[int, int]:
    """
    How many subjects were rated how many times: ``{2: 1614, 3: 3}``.

    What ``is_balanced`` reduces to a boolean. The distinction matters when
    reporting: a dataset with 1,614 subjects at two raters and 3 at three is
    unbalanced, but calling it that without the shape suggests an irregular
    design rather than three stragglers worth looking at.
    """
    values, howmany = np.unique(n_per_subject(category_index), return_counts=True)
    return {int(v): int(n) for v, n in zip(values, howmany)}


# Repeats ---------------------------------------------------------------------  #


def repeated_subject_positions(
    subject_ids: np.ndarray,
) -> dict[object, np.ndarray]:
    """
    Row indices of every subject id appearing more than once.

    Repeats arise two ways -- a rating task re-showing a subject deliberately,
    and a bootstrap draw taking the same row twice. This does not distinguish
    them; it reports the structure.
    """
    order      = np.argsort(subject_ids, kind="stable")
    sorted_ids = subject_ids[order]
    out: dict[object, np.ndarray] = {}
    start = 0
    for i in range(1, len(sorted_ids) + 1):
        if i == len(sorted_ids) or sorted_ids[i] != sorted_ids[start]:
            if i - start > 1:
                out[sorted_ids[start]] = np.sort(order[start:i])
            start = i
    return out
