"""
Coefficients that give each rater their own distribution.

Every function here reads which rater said what. That is the whole of what
separates them from ``group_agreement``: those treat the raters as
interchangeable and read only counts, while these ask how often *this* rater
said "include" as opposed to *that* one, and subtract the difference out before
measuring agreement.

The distinction is one modelling choice appearing twice:

                    pooled marginals      individual marginals
    two raters      Scott's pi            Cohen's kappa
    many raters     Fleiss' kappa         Conger's kappa

Pooled treats the raters as one distribution, so a rater who says "include"
twice as often as another is counted as disagreeing. Individual lets each rater
keep their own rate and removes the difference first, treating a systematic
gap in threshold as a known feature of the raters rather than as disagreement.

Neither is more correct. The gap between them is informative in its own right:
zero when two raters use the categories at the same rate, and growing as their
rates diverge, so it measures how differently two raters draw the line as
distinct from how erratic they are. Those are different problems with different
remedies, and a single coefficient cannot tell them apart.

These take an ``EstimatorInput`` rather than a projection of it, because the
array is what they read and K and the weights travel with it. There is no
counts-shaped summary for this family the way ``NijMatrix`` serves the other:
counts discard rater identity, which is precisely what these need.

Nothing here validates, and that is the contract. Cohen's kappa assumes it was
handed exactly two rater columns and will return a number for any other shape;
the analyser that calls it is responsible for the extraction. That is a mistake
only a developer of this package can make, not a user, because nothing in this
module is exported at the package top level.

References
- Cohen1960: Cohen, Jacob. "A coefficient of agreement for nominal scales."
  Educational and Psychological Measurement 20.1 (1960): 37-46.
- Cohen1968: Cohen, Jacob. "Weighted kappa: nominal scale agreement with
  provision for scaled disagreement or partial credit." Psychological Bulletin
  70.4 (1968): 213-220.
- Conger1980: Conger, Anthony J. "Integration and generalization of kappas for
  multiple raters." Psychological Bulletin 88.2 (1980): 322-328.
- Fleiss1969: Fleiss, Joseph L., Jacob Cohen, and Brian S. Everitt. "Large
  sample standard errors of kappa and weighted kappa." Psychological Bulletin
  72.5 (1969): 323-327.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np

from ...contracts import EMPTY
from ...intermediates.estimator_parts import Weights
from ...intermediates.estimator_parts.estimator_input import EstimatorInput
from .group_agreement import agreement_corrected_by_chance

# Shared parts -----------------------------------------------------------  #
# Marginals are the thing this family has and the other does not. Everything
# below is built from them.


def rater_marginals(data: EstimatorInput, rows: np.ndarray | None = None) -> np.ndarray:
    """
    p_rj: how often each rater chose each category, as a proportion.

    A numpy array of size n_raters by K. Row r is rater r's own distribution
    over the categories, summing to one across the categories they actually
    used.

    ``rows`` restricts the calculation to a subset of subjects, which is what a
    pairwise coefficient needs: Cohen's kappa is defined on the subjects both
    raters saw, and a rater's marginal must be computed on those subjects
    rather than on everything they ever rated. A rater who was strict on the
    shared subjects and lenient elsewhere would otherwise bring the wrong
    distribution to the comparison.

    Empty cells contribute to neither numerator nor denominator, so this says
    how a rater used the categories, not how much they covered.

    Length K in every row, taken from the input rather than from the data, so a
    category this rater never chose is a zero rather than an absence.
    """
    index = data.category_index if rows is None else data.category_index[rows]
    k     = data.n_categories

    out = np.zeros((index.shape[1], k), dtype=np.float64)
    for r in range(index.shape[1]):
        column = index[:, r]
        column = column[column != EMPTY]
        if column.size == 0:
            continue
        out[r] = np.bincount(column.astype(np.int64), minlength=k) / column.size
    return out


def shared_rows(data: EstimatorInput, rater_a: int, rater_b: int) -> np.ndarray:
    """
    Row indices of the subjects both raters rated.

    Where a paired comparison is defined. Reducing to the overlap is part of
    what makes a pair rather than something done to one afterwards -- a subject
    only one of them saw contributes to neither agreement nor disagreement, and
    including it would mean dividing by a denominator that counts subjects
    nobody compared.
    """
    present = data.category_index != EMPTY
    return np.flatnonzero(present[:, rater_a] & present[:, rater_b])


# Chance terms -----------------------------------------------------------  #


def chance_cohen(
    marginals_a: np.ndarray, marginals_b: np.ndarray, weights: Weights | None = None
) -> float:
    """
    Pe as Cohen (1960): the product of the two raters' own proportions.

    Two raters assigning independently at their own rates agree on category j
    with probability p_aj * p_bj, so the chance term is the sum of those
    products. Scott's pi squares the *average* of the two rates instead --
    arithmetic mean against geometric, which is the whole of the difference
    between the two coefficients.

    They coincide exactly when the two raters use the categories at the same
    rate, and diverge as those rates separate. Cohen's is the one that treats a
    rater who is systematically stricter as a known feature to remove rather
    than as disagreement to measure.

    With weights the same model applies, but near misses earn partial credit
    too, so the term becomes p_a . W . p_b. Identity weights reduce it to the
    sum of products above.
    """
    if weights is None or weights.is_identity:
        return float((marginals_a * marginals_b).sum())
    return float(marginals_a @ weights.matrix @ marginals_b)


def chance_conger(marginals: np.ndarray, weights: Weights | None = None) -> float:
    """
    Pe as Conger (1980): pooled chance, less the variance across raters.

    Conger's kappa is Fleiss' with individual marginals. The chance term takes
    the mean rate at which the panel used each category, as Fleiss does, and
    then subtracts the variance of those rates across raters -- which is what
    removes the effect of raters having different thresholds.

        Pe = Sum_j p_bar_j^2  -  Sum_j var_r(p_rj) / n_raters

    When every rater uses the categories identically the variance is zero and
    this reduces to Fleiss' term exactly, which is why the two coefficients
    agree on a panel with matched marginals and separate as the raters diverge.

    The variance uses the sample form, dividing by n_raters - 1, following
    Conger. With a single rater there is nothing to vary and the correction is
    dropped.

    **Takes no weights.** The weighted form of Conger's kappa has a definition
    of its own and is not implemented here; the parameter is present only to
    fail loudly rather than silently if one is ever passed.
    """
    if weights is not None and not weights.is_identity:
        raise NotImplementedError(
            "[chance_conger] weighted Conger's kappa is not implemented. Its "
            "weighted form has a definition of its own, and applying a weight "
            "matrix to the observed agreement while this term ignored it "
            "would give a number that is neither."
        )

    n_raters = marginals.shape[0]
    pooled   = float((marginals.mean(axis=0) ** 2).sum())
    if n_raters < 2:
        return pooled
    return pooled - float(marginals.var(axis=0, ddof=1).sum()) / n_raters


# Maximum agreement -----------------------------------------------------  #


def max_observed_agreement(
    marginals_a: np.ndarray, marginals_b: np.ndarray, weights: Weights | None = None
) -> float:
    """
    The most observed agreement two raters could reach with these marginals.

    Cohen (1960): when each rater's rates are fixed, agreement cannot exceed
    the overlap of the two distributions, Sum_j min(p_aj, p_bj). A pair whose
    rates differ cannot reach kappa = 1 however carefully they rate, and this
    says how far short they must fall.

    With weights, near misses earn credit too, and the most agreement is the
    largest weighted agreement of any table whose rows sum to rater a's rates
    and whose columns sum to rater b's. That is a small transportation
    problem, solved exactly as a linear program. With identity weights it
    reduces to the sum of minima, which is used directly.

    Parameters
    ----------
    marginals_a, marginals_b : np.ndarray
        Each rater's proportions over the K categories, summing to 1.
    weights : Weights | None
        Agreement weights, 1.0 on the diagonal.

    Returns
    -------
    float
        The maximum of Po.
    """
    a = np.asarray(marginals_a, dtype=np.float64)
    b = np.asarray(marginals_b, dtype=np.float64)
    if weights is None or weights.is_identity:
        return float(np.minimum(a, b).sum())

    from scipy.optimize import linprog

    k = a.size
    # One variable per cell of the K x K table, row-major: x[j, l] at j*K + l
    rows = np.kron(np.eye(k), np.ones((1, k)))   # Sum_l x[j, l] = a_j
    cols = np.kron(np.ones((1, k)), np.eye(k))   # Sum_j x[j, l] = b_l
    solved = linprog(
        c=-weights.matrix.reshape(-1),
        A_eq=np.vstack([rows, cols]),
        b_eq=np.concatenate([a, b]),
        bounds=(0, None),
        method="highs",
    )
    if not solved.success:
        raise RuntimeError(
            f"[max_observed_agreement] the linear program did not solve: "
            f"{solved.message}"
        )
    return float(-solved.fun)


def cohen_kappa_max(
    marginals_a: np.ndarray, marginals_b: np.ndarray, weights: Weights | None = None
) -> float:
    """
    The largest Cohen's kappa these marginals allow (Cohen 1960).

    kappa_max = (Po_max - Pe) / (1 - Pe), with Pe the chance term Cohen's
    kappa itself uses. Weighted when weights are given.
    """
    return agreement_corrected_by_chance(
        max_observed_agreement(marginals_a, marginals_b, weights),
        chance_cohen(marginals_a, marginals_b, weights),
    )


# Coefficients -----------------------------------------------------------  #


def cohen_kappa(data: EstimatorInput, weights: Weights | None = None) -> float:
    """
    Cohen's kappa: agreement between two raters, each keeping their own rates.

    Cohen (1960), and Cohen (1968) for the weighted form. Defined on the
    subjects both raters rated; subjects only one of them saw are dropped
    before anything is computed.

    **Assumes exactly two rater columns and does not check.** An analyser
    extracts the pair -- ``EstimatorInput.select_raters`` does it -- and is
    responsible for handing over the right shape. Given three columns this
    reads the first two and returns a number for a question nobody asked.

    What separates this from Scott's pi is whose rates are used for the chance
    term, not how agreement is measured. Two raters who use the categories
    identically give the same answer under both.

    Returns NaN when the pair shares no subjects, which is legitimate rather
    than exceptional: disjoint assignment is a real design, just one that
    cannot produce a pairwise coefficient.
    """
    return agreement_corrected_by_chance(*cohen_parts(data, weights))


def cohen_parts(data: EstimatorInput, weights: Weights | None = None) -> tuple[float, float]:
    """
    Po and Pe for Cohen's kappa: the two numbers ``cohen_kappa`` combines.

    Both on the subjects both raters rated, Pe from each rater's marginals on
    those subjects. (NaN, NaN) when fewer than two subjects are shared, which
    is when ``cohen_kappa`` is NaN. Kept apart so a result can carry the two
    terms and check that they rebuild the coefficient; ``cohen_kappa`` is
    built from this, so the two cannot drift.
    """
    rows = shared_rows(data, 0, 1)
    if rows.size < 2:
        return float("nan"), float("nan")

    index = data.category_index[rows]
    po    = _observed_pairwise(index, data.n_categories, weights)

    marginals = rater_marginals(data, rows)
    pe        = chance_cohen(marginals[0], marginals[1], weights)
    return po, pe


def conger_kappa(data: EstimatorInput) -> float:
    """
    Conger's kappa: many raters, each keeping their own rates.

    Conger (1980). The individual-marginals counterpart to Fleiss' kappa, and
    the many-rater counterpart to Cohen's -- the fourth corner of the table in
    the module docstring.

    Observed agreement is Fleiss', unchanged; only the chance term differs.
    That makes the gap between this and Fleiss' kappa a panel-level measure of
    the same thing the Cohen-minus-Scott gap measures for a pair: how much the
    raters differ in where they draw the line, as opposed to how erratic they
    are.

    Conger motivated the coefficient as an approximation to the average of all
    pairwise Cohen's kappas. It is not that average, and the two can be
    computed separately and compared -- see ``rater_combinations`` -- so if
    they diverge the approximation is not holding for that panel, which is
    itself worth knowing.

    **Assumes every rater rated every subject, and does not check**, as
    nothing in this module does. Conger (1980) defines the coefficient for
    that design: each rater's marginals must come from the same subjects, or
    their spread measures who saw which subjects rather than how the raters
    differ. The analysers refuse incomplete data before calling this. On two
    raters it is exactly Cohen's kappa.

    **Takes no weights**, because its weighted form is not implemented. There
    is no parameter to pass one to, which is a TypeError at the call site
    rather than a number that quietly ignored half of what it was given.
    """
    return agreement_corrected_by_chance(*conger_parts(data))


def conger_parts(data: EstimatorInput) -> tuple[float, float]:
    """
    Po and Pe for Conger's kappa: the two numbers ``conger_kappa`` combines.

    Assumes complete data, as ``conger_kappa`` does. Kept apart so a result
    can carry the two terms and check that they rebuild the coefficient;
    ``conger_kappa`` is built from this, so the two cannot drift.
    """
    from .group_agreement import observed_agreement

    return observed_agreement(data.nij), chance_conger(rater_marginals(data))


def light_kappa(data: EstimatorInput) -> float:
    """
    Light's kappa: the mean of all C(n, 2) pairwise Cohen's kappas.

    Light (1971). Each pair of raters gets its own Cohen's kappa -- computed
    on the subjects they share, with their own marginals -- and the headline
    number is the unweighted mean of those kappas.

    This keeps rater identity twice: once through Cohen's kappa (each pair
    uses its own marginals) and once through the averaging (every pair counts
    equally). Conger's kappa approximates this mean but is not identical to
    it, and comparing the two tells you whether the approximation holds for
    your panel.

    **Assumes every rater rated every subject, and does not check**, as
    nothing in this module does. Light (1971) defines the coefficient for that
    design; with missing ratings, pairs would be averaged over different
    subjects. The analysers refuse incomplete data before calling this.

    A pair whose kappa is undefined (every rating in one category, or fewer
    than two subjects) is NaN from ``cohen_kappa`` and left out of the mean.
    If every pair is NaN, the mean is NaN -- returned, not warned about, since
    a bootstrap can meet it in many draws.

    **Takes no weights**, following Light (1971). There is no weighted form
    of Light's kappa.

    References
    ----------
    Light, R. J. (1971). Measures of response agreement for qualitative
        data: Some generalizations and alternatives. Psychological
        Bulletin, 76(5), 365--377.
    """
    kappas = [
        agreement_corrected_by_chance(po, pe)
        for po, pe in light_pair_parts(data).values()
    ]
    defined = [k for k in kappas if not np.isnan(k)]
    return float(np.mean(defined)) if defined else float("nan")


def light_pair_parts(data: EstimatorInput) -> dict[tuple[int, int], tuple[float, float]]:
    """
    Po and Pe of every pair's Cohen's kappa: what ``light_kappa`` averages.

    Keyed by the pair's column positions, in ``combinations`` order. Assumes
    complete data, as ``light_kappa`` does. A pair with fewer than two
    subjects has (NaN, NaN), as in ``cohen_parts``.
    ``light_kappa`` is built from this, so a result carrying these terms can
    check that they rebuild it.
    """
    return {
        (i, j): cohen_parts(data.select_raters(np.array([i, j])))
        for i, j in combinations(range(data.n_raters), 2)
    }


# Helper methods ---------------------------------------------------------  #


def _observed_pairwise(
    index: np.ndarray, n_categories: int, weights: Weights | None = None
) -> float:
    """
    Po for two raters: the proportion of subjects they agreed on.

    Separate from ``group_agreement.observed_agreement`` only because it can be
    much cheaper. With two raters every subject contributes exactly one pair,
    so agreement is a comparison rather than a sum over category pairs. The
    weighted case still has to look up a credit per subject, but it is one
    lookup rather than a double sum.

    Assumes no empty cells, which holds because the caller reduced to the
    subjects both raters rated before calling.
    """
    a, b = index[:, 0], index[:, 1]
    if weights is None or weights.is_identity:
        return float((a == b).mean())
    return float(weights.matrix[a.astype(np.int64), b.astype(np.int64)].mean())
