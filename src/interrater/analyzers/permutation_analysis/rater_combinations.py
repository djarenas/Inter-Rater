"""
Run a coefficient over subsets of the raters.

Nothing here computes agreement. Each function picks subsets of the rater
columns and calls an estimator on each, so swapping Cohen's kappa for
Fleiss' changes what is measured, not how the subsets are picked. Pairs are
subsets of size two; leaving one rater out is a subset of size n_raters - 1.

The estimator is any coefficient function, since every one takes an
``EstimatorInput``. Analyzers pass the function itself, not a lambda,
because these run inside resampling loops.

**Each subset may be computed on fewer subjects.** Within a subset, a subject
left with fewer than two ratings in the subset drops out, and on sparse data
how many varies by subset. In SR1 main (1,617 subjects, nearly all rated by
two of the four reviewers), removing one rater leaves 650 to 861 subjects
with a single rating, depending on which rater. A difference between two
subsets can therefore reflect how much data each kept rather than how the
raters agree. This is the purpose of ``subjects_per_subset``: it gives each
subset's subject count, keyed like ``estimate_subsets``, so an analyzer can
report the two side by side.
"""

from __future__ import annotations

import itertools
import math
from typing import Callable

import numpy as np

from ...contracts import EMPTY, GUARD_MAX_SUBSETS, MIN_RATERS, MIN_SUBJECTS
from ...intermediates.estimator_parts.estimator_input import EstimatorInput

_ERROR_PREFIX = "[rater_combinations] Error"

Estimator = Callable[[EstimatorInput], float]


# Enumeration ------------------------------------------------------------  #


def rater_subsets(n_raters: int, size: int) -> list[tuple[int, ...]]:
    """
    Every combination of ``size`` rater columns, as sorted tuples of indices.

    Indices are used instead of names. Code that renders the result objects
    can use semantics to get names.

    Combinations, not permutations. Reporting (0, 1) and (1, 0) separately
    would double every result for no information. They are sorted so that a
    subset has one spelling and can be used as a dictionary key.
    """
    if size < MIN_RATERS:
        raise ValueError(
            f"{_ERROR_PREFIX}: a subset needs at least {MIN_RATERS} raters, "
            f"got {size}. Agreement is undefined for a single rater."
        )
    if size > n_raters:
        raise ValueError(
            f"{_ERROR_PREFIX}: asked for subsets of {size} rater(s) from a "
            f"panel of {n_raters}."
        )

    count = math.comb(n_raters, size)
    if count > GUARD_MAX_SUBSETS:
        raise ValueError(
            f"{_ERROR_PREFIX}: {n_raters} raters taken {size} at a time is "
            f"{count:,} subsets, past the cap of {GUARD_MAX_SUBSETS:,}. "
            f"That is almost certainly not what was intended -- a smaller "
            f"size will usually answer the same question."
        )

    return list(itertools.combinations(range(n_raters), size))


# Counting ---------------------------------------------------------------  #


def subjects_per_subset(
    data: EstimatorInput, size: int
) -> dict[tuple[int, ...], int]:
    """
    How many subjects each subset keeps: those rated at least twice by its
    raters.

    Same keys as ``estimate_subsets`` at the same size, so each count lines up
    with its coefficient.
    """
    present = data.category_index != EMPTY
    return {
        subset: int((present[:, list(subset)].sum(axis=1) >= 2).sum())
        for subset in rater_subsets(data.n_raters, size)
    }


# Whether a subset has a coefficient -------------------------------------  #


def subset_problem(present: np.ndarray, names: list[str]) -> str | None:
    """
    Why a subset (or pair) has no coefficient, as text, or None if it has one.

    The one rule, used wherever a subset is measured: its values, its
    bootstrap, the bootstrap of the mean, and its result object. Written
    once so those cannot disagree about which subsets count.

    A subset is measured on its usable subjects, those at least two of its
    raters rated (for a pair, both). It has no coefficient when:

    - fewer than MIN_SUBJECTS subjects are usable, or
    - one of its raters rated none of them. The subset would really be a
      smaller one, and its number would be that smaller subset's.

    Parameters
    ----------
    present : np.ndarray
        Boolean (n_subjects, k): which of the subset's raters rated each
        subject.
    names : list[str]
        The subset's raters, for the reason.

    Returns
    -------
    str | None
        The reason, e.g. "raters a, b share 1 subject", or None.
    """
    usable = present.sum(axis=1) >= 2
    n_usable = int(usable.sum())
    if n_usable < MIN_SUBJECTS:
        word = "subject" if n_usable == 1 else "subjects"
        if len(names) == 2:
            return f"raters {', '.join(names)} share {n_usable} {word}"
        return (
            f"raters {', '.join(names)}: {n_usable} {word} rated by at least "
            f"two of them"
        )
    rated = present[usable].any(axis=0)
    if not rated.all():
        idle = names[int(np.flatnonzero(~rated)[0])]
        return (
            f"rater {idle} rated none of the {n_usable} subjects the others "
            f"share"
        )
    return None


def estimate_subset(unit_data: EstimatorInput, estimator: Estimator) -> float:
    """
    ``estimator`` on one pair's or subset's data, or NaN when it has no
    coefficient under ``subset_problem``.

    ``unit_data`` holds only the subset's rater columns. Every place that
    measures a subset -- observed, bootstrapped, or jackknifed -- calls this,
    so a resample that leaves a rater idle fails like any other draw rather
    than returning a smaller subset's number.
    """
    present = unit_data.category_index != EMPTY
    names = [str(i) for i in range(present.shape[1])]
    if subset_problem(present, names) is not None:
        return float("nan")
    return estimator(unit_data)


# Running an estimator ---------------------------------------------------  #


def estimate_subsets(
    data: EstimatorInput, estimator: Estimator, size: int
) -> dict[tuple[int, ...], float]:
    """
    Run ``estimator`` on every subset of ``size`` raters.

    Keyed by the subset, so ``(0, 2)`` holds what the estimator returned for
    raters 0 and 2 alone. A subset whose data cannot produce a coefficient --
    fewer than two usable subjects or an idle rater (``subset_problem``), or
    the estimator itself returns NaN (for kappa-type coefficients, every
    rating in one category) -- comes back as NaN rather than being
    omitted, because a missing key is easy to overlook and a subset that
    could not be computed is a fact about the design.

    The subset is taken with ``select_raters``, which carries K and the weights
    along. K in particular must not be re-derived: a subset that happens to
    miss a rare category would otherwise be computed against a smaller K than
    the rating task offered, and any coefficient whose chance term counts
    categories would move for no reason but the subsetting.
    """
    return {
        subset: estimate_subset(data.select_raters(np.array(subset)), estimator)
        for subset in rater_subsets(data.n_raters, size)
    }


def estimate_pairs(
    data: EstimatorInput, estimator: Estimator
) -> dict[tuple[int, int], float]:
    """
    Run ``estimator`` on every pair of raters. The size-two case.

    What version 1 of this package called permutating pairs. Worth its own name
    because it is the common case and because pairs are the only subsets on
    which the two-rater coefficients -- Cohen's kappa, Scott's pi -- are
    defined at all.

    Each pair is analysed on the subjects both raters rated, which the
    estimator handles. On sparse data that is a different set of subjects for
    every pair, which is what ``subjects_per_subset`` is for.
    """
    return estimate_subsets(data, estimator, 2)


# Summarising ------------------------------------------------------------  #


def mean_over_subsets(values: dict) -> float:
    """
    The mean of a subset dictionary, skipping the subsets that could not be
    computed.

    **The values are not independent.** Two subsets sharing a rater are
    usually positively correlated, so a standard error that treats them as
    independent is too small. Version 1 of this package made that mistake.

    A bootstrap that resamples subjects and recomputes every subset within each
    draw preserves the correlation without having to model it. So this returns
    a point estimate and leaves intervals to the resampler, which is the only
    thing holding enough of the picture to get them right.
    """
    array  = np.array(list(values.values()), dtype=np.float64)
    usable = ~np.isnan(array)
    if not usable.any():
        return float("nan")
    return float(array[usable].mean())


def mean_over_same_subsets(values: dict, subsets: list) -> float:
    """
    The mean over exactly these subsets, or NaN if any of them is undefined.

    A bootstrap or jackknife value of the mean has to estimate the same
    quantity as the point estimate: the average over the subsets that were
    defined on the real data. ``mean_over_subsets`` skips whatever is
    undefined, so a resample that loses one subset would average the others
    -- a different quantity, counted as an ordinary draw. Here such a
    resample is NaN, and the interval counts it as failed, as it does for a
    single coefficient. Subsets not listed are ignored even when defined.

    Parameters
    ----------
    values : dict
        Coefficient per subset, as from ``estimate_subsets`` or
        ``estimate_pairs``.
    subsets : list
        The subsets to average: those defined at the point estimate.

    Returns
    -------
    float
        Their mean, or NaN if any is NaN or the list is empty.
    """
    if not subsets:
        return float("nan")
    array = np.array([values[s] for s in subsets], dtype=np.float64)
    if np.isnan(array).any():
        return float("nan")
    return float(array.mean())
