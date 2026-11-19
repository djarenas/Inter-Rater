"""
Group coefficients from counts: Fleiss' kappa, Scott's pi, Gwet's AC1,
Brennan-Prediger and percent agreement.

Each coefficient is (Po - Pe) / (1 - Pe). Po is the same for all of them
(``observed_agreement``); only the chance term Pe differs.

Each reads an NijMatrix: for each subject, how many raters chose each
category. Therefore, they use pooled marginals; there is no information
about which rater chose which category.

- K (the number of categories) is calculated from the width of the NijMatrix.
  This means that an analyzer must build the matrix with a K width;
  ``EstimatorInput.nij`` does, with K from the semantics.
- Nothing here validates. The analyzer checks the data once, before; these
  functions run inside bootstrap loops.
- A ``weights`` argument appears only where the coefficient has a weighted
  form: Fleiss, Brennan-Prediger, percent agreement. Gwet's AC1 has none
  (its weighted form is AC2). Scott's pi accepts weights only to refuse
  non-identity ones.

Formulas are in README.md.

References
- Brennan, R. L., & Prediger, D. J. (1981). Educational and Psychological
  Measurement, 41(3), 687-699.
- Fleiss, J. L. (1971). Psychological Bulletin, 76(5), 378-382.
- Gwet, K. L. (2014). Handbook of Inter-Rater Reliability (4th ed.).
"""

from __future__ import annotations

import numpy as np

from ...intermediates.estimator_parts import NijMatrix, Weights

_SCOTT_ERROR_PREFIX = "[scotts_pi] Error"

# Shared parts -----------------------------------------------------------  #
# These are methods shared by different kappa statistics. 
# Documentation specifies how the function deals with the case where not 
# all subjects were rated by all raters.

def observed_agreement(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Po: the proportion of agreeing rater pairs, averaged over subjects.

    For one subject rated by n raters there are n(n-1) ordered pairs of
    ratings, of which sum_j n_ij(n_ij - 1) agree. Fleiss (1971) Equation 2
    defines this per subject and Equation 3 averages it, and the average is
    what every coefficient in this module uses as its observed agreement.

    Fleiss writes n as a constant, the same for every subject. Here each
    subject is divided by its own count, which agrees with Fleiss exactly when
    the data is balanced. An analyser that promises *exact*
    Fleiss must therefore require balanced data; this function will not stop
    it, because refusing here would also refuse the coefficients that have no
    such requirement.

    Subjects rated fewer than twice are skipped.

    With weights, a disagreement still earns whatever credit the weight matrix
    gives it, so the count of agreeing pairs becomes a weighted sum over every
    pair of categories:

        P_i = ( c_i . W . c_i  -  n_i ) / [ n_i (n_i - 1) ]

    where c_i is subject i's row of counts. The subtraction removes the pairing
    of a rating with itself, which is not a pair of raters. Identity weights
    reduce this exactly to the unweighted form, which is why there is one
    formula here and not two -- but the unweighted path is kept separate anyway
    because it is one pass over the counts where the weighted one is a sum over
    every pair of categories.
    """
    nij_matrix   = nij.matrix                       # i: subject index; j: category index
    per_subj = nij.n_per_subject.astype(np.float64)   # n_i: raters per subject

    usable = per_subj >= 2  # Mask of subjects with at least two raters
    if not usable.any():
        return float("nan")

    # Remove entries without at least two raters
    nij_matrix   = nij_matrix[usable].astype(np.float64)
    per_subj = per_subj[usable]

    # Agreement per subject (P_i)
    if weights is None or weights.is_identity:
        agreeing = (nij_matrix * (nij_matrix - 1)).sum(axis=1)
    else:
        agreeing = np.einsum("ia,ab,ib->i", nij_matrix, weights.matrix, nij_matrix)
        agreeing = agreeing - nij_matrix.sum(axis=1)

    P_i = agreeing / (per_subj * (per_subj - 1))

    # Average over subjects
    return float(P_i.mean())


def category_proportions(nij: NijMatrix) -> np.ndarray:
    """
    p_j: the share of ratings that went to each category, weighted by subject.

    Fleiss (1971) Equation 1, which is `p_j = (1/Nn) Sum_i n_ij` with a
    constant n. When n is constant there is only one way to read that. When it
    is not, there are two, and they give different answers:

    - **averaged per subject**, the mean of each subject's own proportions,
      where every subject counts once however many raters saw it.
    - **pooled per rating**, `Sum_i n_ij / Sum_i n_i`, where a subject rated
      four times contributes twice what a subject rated twice does

    This package uses the first option, averaged per subject. 
    The reason is consistency rather than authority: the
    observed agreement above is already a mean over subjects, each weighted
    equally, so weighting Pe by rating instead would have the two halves of one
    formula disagree about what they are averaging over.

    Subjects rated fewer than twice are skipped.

    Has one component per column of the NijMatrix, so a category nobody chose
    gets a zero, provided the matrix was built K wide. This is relevant for
    Gwet's AC1 and Brennan-Prediger.

    Weights do not enter. These are proportions of ratings.
    """
    counts   = nij.matrix                              # i: subject; j: category
    per_subj = nij.n_per_subject.astype(np.float64)    # n_i: raters per subject
    usable   = per_subj >= 2                           # Mask for subjects
    if not usable.any():
        return np.zeros(nij.n_categories, dtype=np.float64)

    return (counts[usable] / per_subj[usable, None]).mean(axis=0)



# Chance terms -----------------------------------------------------------  #


def chance_fleiss(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Pe as Fleiss (1971) Equation 5: the sum of squared category proportions.

    The model is that raters assign at the overall rates p_j and ignore the
    subject, so two raters agree by chance with probability sum_j p_j^2. Scott's
    pi uses the same expression for two raters, which is why they are one
    function here.

    This is the term that makes kappa behave oddly under skew. At 90% in one
    category, Pe is around 0.82 and kappa is left dividing by 0.18; therefore, a
    dataset where raters agreed on 93% of subjects scores lower than one where
    they agreed on 91% but spread their ratings more evenly.

    With weights the same model applies, but two raters assigning at the
    overall rates are credited for near misses as well as matches, so the term
    becomes p . W . p. Identity weights reduce it to the sum of squares above.
    """
    pj = category_proportions(nij)
    if weights is None or weights.is_identity:
        return float((pj ** 2).sum())
    return float(pj @ weights.matrix @ pj)


def chance_gwet_ac1(nij: NijMatrix) -> float:
    """
    Pe as Gwet's AC1: the mean of p_j(1 - p_j), scaled by K - 1.

    Gwet's argument is that chance agreement happens when a rater is uncertain,
    and uncertainty is highest for categories used about half the time. A
    category used 90% of the time is being chosen on purpose, not at random, so
    it should contribute little to the chance term. The result is a coefficient
    that does not collapse under skew the way kappa does.

    This function is an example of why it is important that NijMatrix
    gets properly built from the analyzers with the right number of
    categories.

    **Takes no weights, this package has no weighted form for it.** The
    weighted version of AC1 is Gwet's AC2, a separate coefficient with a
    published definition of its own. Handing weights to this one is a
    TypeError.
    """
    pj = category_proportions(nij)
    k = nij.n_categories
    return float((pj * (1.0 - pj)).sum() / (k - 1))


def chance_brennan_prediger(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Pe as Brennan and Prediger: 1/K, regardless of what the raters did.

    Chance is a uniform guess among the K categories.

    It assumes chance means picking
    uniformly among the categories on offer, so the observed proportions do not
    enter at all. Unaffected by skew because it never looks at it.

    With weights the formula becomes sum(W) / K^2, where W is the K x K weight
    matrix and K is the number of categories. Identity weights reduce this to
    K / K^2 = 1/K, so the unweighted case is a special case of the weighted
    one. Off-diagonal entries inflate the sum, reflecting that near-misses earn
    partial credit and therefore raise the baseline that chance alone would
    produce.

    The weighted formula is from Gwet (2014).
    """
    k = nij.n_categories
    if weights is None or weights.is_identity:
        return float(1.0 / k)
    return float(weights.matrix.sum() / (k * k))


# Coefficients -----------------------------------------------------------  #

def agreement_corrected_by_chance(observed: float, chance: float) -> float:
    """
    (Po - Pe) / (1 - Pe): how much of the room above chance was taken up.

    Returns NaN when chance agreement is total. That happens when every rating
    in the data went to one category, and it is not a failure to handle but a
    question with no answer: if chance alone would have produced perfect
    agreement, there is no room above it to have taken up. 
    
    Returning NaN is especially relevant for bootstrap draws. It would 
    be misleading to return zero.
    """
    if np.isclose(chance, 1.0):
        return float("nan")
    return float((observed - chance) / (1.0 - chance))


def percent_agreement(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Po, uncorrected. The proportion of rater pairs that agreed.

    Worth reporting beside a kappa rather than instead of one: a kappa of 0.42 next to 93% agreement is
    telling you something about the prevalence, and not just about the raters.
    """
    return observed_agreement(nij, weights)


def fleiss_kappa(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Fleiss' kappa: agreement among many raters, corrected for chance.

    Fleiss (1971), Equation 7. Defined for the case where every subject is
    rated the same number of times, though not necessarily by the same raters
    -- Fleiss' own data drew six psychiatrists per patient from a pool of 43.

    **The constant-rater requirement is real and is not checked here.** Fleiss
    Equation 2 divides by n(n-1) with a fixed n, and the introduction says
    plainly that the case of varying numbers of raters is not considered. 
    
    On unbalanced data the result reports itself as an approximation.

    Examples
    --------
    Fleiss (1971) Table 1: thirty patients, six raters each, five diagnoses.
    The paper reports 0.430.
    """
    return agreement_corrected_by_chance(
        observed_agreement(nij, weights), chance_fleiss(nij, weights)
    )


def gwet_ac1(nij: NijMatrix) -> float:
    """
    Gwet's AC1: chance-corrected agreement that holds up under skew.

    Same observed agreement as Fleiss', a different answer to how much of it
    chance would explain. On heavily skewed data the two diverge sharply --
    across the six systematic-review datasets in example_data, by 0.17 to
    0.36 --
    and the gap is entirely the chance term.

    **Unweighted only.** The weighted coefficient is Gwet's AC2, defined
    separately and not implemented here, so this takes no weights argument at
    all and offering one is a TypeError at the call site.
    """
    return agreement_corrected_by_chance(
        observed_agreement(nij), chance_gwet_ac1(nij)
    )


def brennan_prediger(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Brennan and Prediger's coefficient: chance is a uniform guess among K.

    Its chance term is independent of how the raters actually used the
    categories.

    With weights, both halves change: observed agreement credits near-misses
    via the weight matrix, and the chance term becomes sum(W) / K^2 instead
    of 1/K. The unweighted coefficient is the special case where W is the
    identity.
    """
    return agreement_corrected_by_chance(
        observed_agreement(nij, weights), chance_brennan_prediger(nij, weights)
    )


def scotts_pi(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Scott's pi: agreement between raters using pooled category rates.

    Scott (1955) proposed this two-rater coefficient; it generalizes to many
    raters identically to how Fleiss' kappa does. Both use the same chance term
    (the sum of squared pooled category proportions); therefore, at two raters they are one
    function. The name honors Scott's priority; the implementation is shared
    with Fleiss' kappa.

    **Scott's pi is defined for nominal (unweighted) data only.** If weights
    are provided, they must be None or identity; any other weighting raises
    ValueError. For ordinal agreement with pooled marginals, use a weighted
    generalization of Fleiss' kappa instead.

    Parameters
    ----------
    nij : NijMatrix
        The counts matrix.
    weights : Weights | None
        Must be None or identity weights. Non-identity weights raise ValueError.

    Returns
    -------
    float
        Scott's pi coefficient.

    Raises
    ------
    ValueError
        If weights are provided and non-identity.

    References
    ----------
    Scott, William A. 1955. "Reliability of content analysis: The case of
    nominal scale coding." *Public opinion quarterly* 19.3: 321-325.
    """
    # Scott's pi is unweighted; reject non-identity weights
    if weights is not None and not weights.is_identity:
        raise ValueError(
            f"{_SCOTT_ERROR_PREFIX}: Scott's pi is defined for nominal (unweighted) data only. "
            f"weights must be None or identity, got {weights.kind}."
        )
    
    return agreement_corrected_by_chance(
        observed_agreement(nij, None),  # Always unweighted
        chance_fleiss(nij, None)
    )
