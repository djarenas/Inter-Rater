"""
Running a coefficient over subsets of the raters, rather than over all of them.

A panel coefficient is one number for the whole group. It cannot tell you that
one rater was the problem, because it has already averaged that rater away. On
one real screening round the six pairs span 0.58 to 0.82 with a single reviewer
accounting for the spread, and the panel coefficient is 0.71 with no hint that
anything is uneven; on another the three pairs span 0.05. Those are different
situations calling for different responses, and the single number is identical
in kind for both.

Nothing here computes agreement. Each function enumerates subsets of the rater
columns and calls an estimator on each, so swapping Cohen's kappa for Fleiss'
changes what is being measured and not how the enumeration works. Pairs are the
subsets of size two; leaving one rater out is a subset of size n-1; anything
between is available for the asking.

**The estimator is a callable taking an ``EstimatorInput``.** The coefficients
in ``group_agreement`` take an ``NijMatrix`` instead, so an analyser using one
of those supplies a small module-level adapter rather than a lambda -- these
results are often computed inside a resampling loop, and a closure cannot be
shipped to a worker process.

**Dropping raters guts sparse data, and the damage is uneven.** Every subject
left with fewer than two ratings contributes to nothing, so the coefficient for
a subset is computed on whatever survived rather than on the study. Measured on
one real dataset of 1,617 subjects, each rated by two of four reviewers,
removing a rater leaves between 650 and 861 subjects holding a single rating --
and *which* rater was removed decides how many. A comparison across subsets is
therefore confounded with how much data each rater was carrying unless the
surviving count travels alongside the coefficient.

That is why ``subjects_per_subset`` exists and why an analyser is expected to
call it. These functions return plain dictionaries keyed the same way, so the
analyser can join them into a result object that keeps the two together; at
this level they are loose numbers, and keeping them from being misused is the
analyser's job.
"""

from __future__ import annotations

import itertools
import math
from typing import Callable

import numpy as np

from ...array_utils import n_per_subject
from ...contracts import EMPTY, MIN_RATERS
from ...intermediates.estimator_parts.estimator_input import EstimatorInput

_ERROR_PREFIX = "[rater_combinations] Error"

MAX_SUBSETS = 10_000
"""
How many subsets a single call will enumerate before refusing.

Subsets grow as a binomial coefficient: ten raters give 45 pairs and 252
subsets of size five, but twenty raters give 184,756 of size ten. The cap is a
guard against asking for a number of coefficients nobody intended rather than a
statement about what is meaningful, so it is deliberately loose.
"""

Estimator = Callable[[EstimatorInput], float]


# Enumeration ------------------------------------------------------------  #


def rater_subsets(n_raters: int, size: int) -> list[tuple[int, ...]]:
    """
    Every combination of ``size`` rater columns, as sorted tuples of indices.

    Combinations, not permutations: a subset of raters has no order, and
    reporting (0, 1) and (1, 0) separately would double every result for no
    information. Sorted so that a subset has one spelling and can be used as a
    dictionary key without two entries meaning the same thing.

    Indices, not names. Whoever renders the result has the semantics and can
    supply them.
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
    if count > MAX_SUBSETS:
        raise ValueError(
            f"{_ERROR_PREFIX}: {n_raters} raters taken {size} at a time is "
            f"{count:,} subsets, past the cap of {MAX_SUBSETS:,}. That is "
            f"almost certainly not what was intended -- a smaller size, or "
            f"leave-one-out, will usually answer the same question."
        )

    return list(itertools.combinations(range(n_raters), size))


# Counting ---------------------------------------------------------------  #


def subjects_per_subset(
    data: EstimatorInput, size: int
) -> dict[tuple[int, ...], int]:
    """
    How many subjects survive each subset: those still rated at least twice.

    The denominator every coefficient in this module was actually computed on,
    and the thing that makes results comparable across subsets or reveals that
    they are not. A subset whose coefficient rests on 756 subjects and one
    resting on 966 are not measuring the same study, and nothing about the two
    numbers says so.

    Returns the same keys as ``over_subsets`` at the same size, so the two can
    be joined without matching anything up by hand.
    """
    present = data.category_index != EMPTY
    return {
        subset: int((present[:, list(subset)].sum(axis=1) >= 2).sum())
        for subset in rater_subsets(data.n_raters, size)
    }


# Running an estimator ---------------------------------------------------  #


def over_subsets(
    data: EstimatorInput, estimator: Estimator, size: int
) -> dict[tuple[int, ...], float]:
    """
    Run ``estimator`` on every subset of ``size`` raters.

    Keyed by the subset, so ``(0, 2)`` holds what the estimator returned for
    raters 0 and 2 alone. A subset whose data cannot produce a coefficient --
    too few shared subjects, every rating in one category -- comes back as NaN
    rather than being omitted, because a missing key is easy to overlook and a
    subset that could not be computed is a fact about the design.

    The subset is taken with ``select_raters``, which carries K and the weights
    along. K in particular must not be re-derived: a subset that happens to
    miss a rare category would otherwise be computed against a smaller K than
    the rating task offered, and any coefficient whose chance term counts
    categories would move for no reason but the subsetting.
    """
    return {
        subset: estimator(data.select_raters(np.array(subset)))
        for subset in rater_subsets(data.n_raters, size)
    }


def over_pairs(
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
    return over_subsets(data, estimator, 2)


def leave_one_rater_out(
    data: EstimatorInput, estimator: Estimator
) -> dict[int, float]:
    """
    Run ``estimator`` on the panel with each rater removed in turn.

    **Keyed by the rater who was dropped**, not by the subset that remained,
    because the question is about that rater. Subtracting each value from the
    full-panel coefficient gives a jackknife influence measure: what the panel
    estimate becomes without this person.

    That is a different question from how a rater agrees with each other rater
    one at a time, which is what ``over_pairs`` answers. A rater can agree
    reasonably with everybody and still be the one holding the panel
    coefficient down, and a rater can disagree sharply with one colleague
    without moving the panel at all.

    Two cautions, both from the module docstring and both sharper here. The
    surviving subject count differs by which rater was removed, so the
    differences are confounded on sparse data unless that count is reported
    beside them. And dropping a rater from balanced data usually breaks the
    balance, so a panel coefficient that required it -- exact Fleiss' kappa --
    no longer applies to what is left. An analyser promising exact Fleiss has
    to notice that; this function will not.
    """
    subsets = rater_subsets(data.n_raters, data.n_raters - 1)
    dropped = lambda subset: (set(range(data.n_raters)) - set(subset)).pop()
    return {
        dropped(subset): estimator(data.select_raters(np.array(subset)))
        for subset in subsets
    }


def subjects_per_leave_one_out(data: EstimatorInput) -> dict[int, int]:
    """
    How many subjects survive each rater's removal, keyed by the rater dropped.

    The companion to ``leave_one_rater_out``, and on sparse data the more
    important half of the pair. Removing a reviewer who covered a great deal
    destroys more of the dataset than removing one who covered little, so a
    large change in the coefficient may say more about how much work that
    person did than about how well they did it.
    """
    present = data.category_index != EMPTY
    out: dict[int, int] = {}
    for rater in range(data.n_raters):
        keep = [r for r in range(data.n_raters) if r != rater]
        out[rater] = int((present[:, keep].sum(axis=1) >= 2).sum())
    return out


# Summarising ------------------------------------------------------------  #


def mean_over_subsets(values: dict) -> float:
    """
    The mean of a subset dictionary, skipping the subsets that could not be
    computed.

    **The values are not independent.** Two subsets sharing a rater are
    positively correlated, so the variance of this mean is not the variance of
    an average of independent things. That is what version 1 of this package
    assumed when it put an analytic standard error on the average of pairwise
    kappas, and it is why those intervals were too narrow -- on real data, by
    about 28%.

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
        Coefficient per subset, as from ``over_subsets`` or ``over_pairs``.
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
