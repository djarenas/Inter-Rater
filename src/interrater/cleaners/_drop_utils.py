"""
Dropping rows and columns that hold too little data.

Free functions over the code array. Nothing here raises or knows what a config
is: the caller supplies two thresholds and an order, and gets back boolean
masks saying what survived.

Dropping is not symmetric. Removing a rater lowers the rating count of every
subject that rater covered, and removing a subject lowers the count of every
rater who saw it. So one sweep can leave data below the minimums, and repeating
it can keep finding new casualties. Repeated until a sweep drops nothing, it
ends at the largest data meeting both minimums: if two sets of rows and columns
each meet them, so does their union, so that largest set is unique and no sweep
ever removes part of it. The order of the axes changes the path, not the end.
"""

from __future__ import annotations

import numpy as np

from ..contracts import EMPTY


def rating_counts(codes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Ratings held by each subject (rows) and each rater (columns)."""
    present = codes != EMPTY
    return present.sum(axis=1), present.sum(axis=0)


def sweep_once(
    codes:      np.ndarray,
    min_rater:   int,
    min_subject: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    One sweep. Returns (subject_keep, rater_keep, subject_counts, rater_counts).

    Raters are evaluated on the full array; subjects on what the raters left
    behind. Repeating sweeps (sweep_until_settled) makes the order irrelevant
    to where they end.

    The counts are the ones each subject and rater was judged on, so for the
    second axis they already reflect what the first removed. A report that
    says how much a casualty held when it went has to use these, not counts
    taken from the array before the sweep.
    """
    _, per_rater = rating_counts(codes)
    rater_keep   = per_rater >= min_rater
    per_subject, _ = rating_counts(codes[:, rater_keep])
    subject_keep = per_subject >= min_subject
    return subject_keep, rater_keep, per_subject, per_rater


def sweep_until_settled(
    codes:        np.ndarray,
    min_rater:    int,
    min_subject:  int,
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int]], np.ndarray, np.ndarray]:
    """
    Sweep until a pass drops nothing.

    Returns cumulative masks over the original array, the per-pass counts,
    and what each casualty held in the pass that removed it (subjects, then
    raters; -1 for those kept). Every pass but the last removes at least one
    row or column, so it ends within n_subjects + n_raters + 1 passes. A
    cascade can still run a long way, which is why the caller is told how many
    passes it took.
    """
    subject_keep = np.ones(codes.shape[0], dtype=bool)
    rater_keep   = np.ones(codes.shape[1], dtype=bool)
    subject_held = np.full(codes.shape[0], -1, dtype=np.int64)
    rater_held   = np.full(codes.shape[1], -1, dtype=np.int64)
    per_pass: list[tuple[int, int]] = []

    while True:
        live = codes[np.ix_(subject_keep, rater_keep)]
        s_now, r_now, s_count, r_count = sweep_once(live, min_rater, min_subject)

        dropped_subjects = int((~s_now).sum())
        dropped_raters   = int((~r_now).sum())
        per_pass.append((dropped_raters, dropped_subjects))

        if dropped_subjects == 0 and dropped_raters == 0:
            return subject_keep, rater_keep, per_pass, subject_held, rater_held

        s_index = np.flatnonzero(subject_keep)
        r_index = np.flatnonzero(rater_keep)
        subject_held[s_index[~s_now]] = s_count[~s_now]
        rater_held[r_index[~r_now]]   = r_count[~r_now]
        subject_keep[s_index[~s_now]] = False
        rater_keep[r_index[~r_now]]   = False


def distribution(counts: np.ndarray) -> dict[int, int]:
    """How many entries hold how many ratings: ``{0: 12, 1: 380}``."""
    if counts.size == 0:
        return {}
    values, howmany = np.unique(counts, return_counts=True)
    return {int(v): int(n) for v, n in zip(values, howmany)}
