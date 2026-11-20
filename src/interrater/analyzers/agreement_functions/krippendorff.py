"""
Krippendorff's alpha, which is built the other way up.

Alpha measures *disagreement* and compares it with the disagreement chance
would produce:

    alpha = 1 - Do / De

- Do, observed disagreement: the disagreeing pairs of ratings within each
  subject rated at least twice, each subject divided by n_i - 1, summed and
  divided by n:

      Do = (1/n) Sum_i Sum_{a,b} d_ab c_ia c_ib / (n_i - 1)

- De, expected disagreement: how often two ratings drawn without
  replacement from all the pairable ratings would differ:

      De = Sum_{a,b} d_ab n_a n_b / [ n (n - 1) ]

  c_ia: ratings of subject i in category a. n_a: pairable ratings in
  category a. n: pairable ratings (``pairable_ratings``). d_ab: the cost of
  confusing a with b, 1 - W_ab, zero on the diagonal.

Two differences from the kappa family. Pairs are pooled across subjects, so
a subject with more ratings weighs more. And De draws without replacement,
which is where the n - 1 comes from.

References
- Krippendorff2004: Krippendorff, Klaus. *Content Analysis: An Introduction to
  Its Methodology.* 2nd ed., Sage, 2004, chapter 11.
- Krippendorff2011: Krippendorff, Klaus. "Computing Krippendorff's
  alpha-reliability." Departmental Papers, Annenberg School for Communication,
  University of Pennsylvania, 2011.
"""

from __future__ import annotations

import numpy as np

from ...intermediates.estimator_parts import EstimatorInput, NijMatrix, Weights

# Shared parts -----------------------------------------------------------  #


def pairable_ratings(nij: NijMatrix) -> int:
    """
    n: the number of ratings that belong to a subject rated more than once.

    Alpha's sample size, and not the same as the number of ratings. A subject
    rated by one rater contributes a rating but no pair, so it cannot disagree
    with anything and takes no part in the calculation. Krippendorff calls the
    ratings that remain *pairable*.

    This is the denominator the chance term is corrected by, which is why it
    has to be counted rather than assumed: using the full rating count would
    understate the correction on data with many singly-rated subjects.
    """
    per_subject = nij.n_per_subject
    return int(per_subject[per_subject >= 2].sum())


def observed_disagreement(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Do: the proportion of pairable rating pairs that disagree.

    Every pair of ratings of the same subject is counted, pooled across
    subjects. For a subject rated n_i times there are n_i(n_i - 1) ordered
    pairs, of which the disagreeing ones are counted with whatever credit the
    weights allow:

        Do = Sum_i Sum_{a != b} d_ab c_ia c_ib / [ n_i - 1 ]   /  n

    where d_ab is the cost of confusing category a with category b -- one minus
    the weight -- and n is the count of pairable ratings.

    Note the ``n_i - 1`` rather than ``n_i (n_i - 1)``. Each subject's pairs are
    divided by one less than its rating count and the whole sum by the total
    pairable ratings, which is what weights subjects by how many pairs they
    offered instead of equally. That is the pooling described in the module
    docstring, and it is the first of the two places alpha parts company with
    Fleiss.

    Identity weights make d_ab one for every disagreement and zero on the
    diagonal, which is the unweighted case.
    """
    counts = nij.matrix.astype(np.float64)
    per    = nij.n_per_subject.astype(np.float64)

    usable = per >= 2
    if not usable.any():
        return float("nan")

    counts, per = counts[usable], per[usable]
    n = float(pairable_ratings(nij))      # one definition of alpha's sample size

    if weights is None or weights.is_identity:
        # Disagreeing pairs are all pairs less the agreeing ones.
        agreeing    = (counts * (counts - 1)).sum(axis=1)
        disagreeing = per * (per - 1) - agreeing
    else:
        cost        = 1.0 - weights.matrix
        disagreeing = np.einsum("ia,ab,ib->i", counts, cost, counts)

    return float((disagreeing / (per - 1.0)).sum() / n)


def expected_disagreement(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    De: the disagreement that shuffling the same ratings would produce.

    Alpha's chance model is not a claim about how raters behave. It takes the
    ratings actually given, forgets which subject each belonged to, and asks
    how often two drawn at random would differ:

        De = Sum_{a != b} d_ab n_a n_b / [ n (n - 1) ]

    where n_a is how many pairable ratings went to category a. Drawing without
    replacement is what puts ``n - 1`` in the denominator, and it is the second
    place alpha parts company with the kappa family: Fleiss' chance term is a
    population quantity with no sample-size correction in it, so the two
    diverge on small datasets and converge as n grows.

    Only pairable ratings are counted, for the same reason as in the observed
    term -- a rating that could not disagree with anything should not
    contribute to what disagreement chance would produce either.
    """
    counts = nij.matrix.astype(np.float64)
    per    = nij.n_per_subject

    usable = per >= 2
    if not usable.any():
        return float("nan")

    per_category = counts[usable].sum(axis=0)
    n            = float(pairable_ratings(nij))   # the same n as the observed term
    if n < 2:
        return float("nan")

    if weights is None or weights.is_identity:
        agreeing    = (per_category * (per_category - 1.0)).sum()
        disagreeing = n * (n - 1.0) - agreeing
    else:
        cost        = 1.0 - weights.matrix
        disagreeing = float(per_category @ cost @ per_category)

    return float(disagreeing / (n * (n - 1.0)))


# Coefficient ------------------------------------------------------------  #


def krippendorff_alpha(data: EstimatorInput) -> float:
    """
    Krippendorff's alpha: one minus observed disagreement over expected.

    Krippendorff (2004), chapter 11. Reaches 1 when every pair of ratings of a
    subject agrees, 0 when disagreement is what shuffling would have produced,
    and goes negative when raters disagree more than chance.

    Counts ``data`` once (``data.nij``, K wide) and reads ``data.weights``.
    Unlike the AC1 family it has a published weighted form:
    the cost of confusing two categories is one minus their weight, entering
    both terms identically. Identity weights give the nominal case.

    Returns NaN when expected disagreement is zero, which happens when every
    pairable rating went to one category. Bootstrap
    draws reach this state on heavily skewed data, which is why the resampler
    counts NaNs rather than treating them as zeros.
    """
    nij = data.nij
    do = observed_disagreement(nij, data.weights)
    de = expected_disagreement(nij, data.weights)

    if np.isnan(do) or np.isnan(de) or np.isclose(de, 0.0):
        return float("nan")
    return float(1.0 - do / de)
