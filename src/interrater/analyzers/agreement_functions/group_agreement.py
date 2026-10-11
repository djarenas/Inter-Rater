"""
Coefficients for many raters, computed from how many raters chose each category.

Every function here reads the same thing: for each subject, how many raters
placed it in each category. None of them care *which* rater said what, which is
what separates them from the pairwise family -- Cohen's kappa cannot be
computed from these counts at all, because the counts have thrown away the
identity of the raters.

They also share a skeleton. Observed agreement is computed once and identically
for all of them; only the chance term differs:

    coefficient = (Po - Pe) / (1 - Pe)

On real data with Po = 0.8941, Fleiss (Pe = 0.784), Gwet's AC1 (0.216),
Brennan-Prediger (0.500) and percent agreement (0) are the same arithmetic with
four different answers to "how much agreement would chance alone produce". That
question is the whole of the disagreement between them, and it is why they can
differ by 0.35 on identical ratings.

These take an ``NijMatrix`` and, optionally, a ``Weights``. Not an
``EstimatorInput``: counts are what they actually read, and an analyser
converts once and passes the result down rather than having each function
recompute it. K travels with the matrix -- its width *is* K -- so a coefficient
can never be computed against a smaller K than the rating task had, which is
what would happen if a bootstrap draw that missed a rare category were allowed
to infer K from itself.

Nothing here validates, and that is the contract rather than an omission. They
assume the counts came from an array that satisfied ``Ratings`` at the
boundary, and they assume the caller knows which coefficient it is invoking and
what that coefficient supports. The analysers are that caller: each one was
written to compute a particular coefficient and already knows its
preconditions, so asking a function to describe itself at runtime would be
querying a fact the caller was authored with.

Whether a coefficient *applies* to the data -- whether Fleiss' constant-rater
requirement is met, whether there are enough raters for the question being
asked -- is therefore checked once, by the analyser, before any of this runs.

**A weights argument appears only where it is read.** Fleiss' kappa, percent
agreement, and Brennan-Prediger take one; Gwet's AC1 does not, because there is
no weighted form of it here -- the weighted AC1 is Gwet's AC2, a separate
coefficient. Offering weights to a coefficient that has none is therefore a
TypeError at the call site rather than a number that quietly ignored half of
what it was given. Signature uniformity is not needed for resampling: a
resampler receives a one-argument callable, and binding whatever else a
coefficient wants is the analyser's job.

References
- Fleiss1969: Fleiss, Joseph L., Jacob Cohen, and Brian S. Everitt.
  "Large sample standard errors of kappa and weighted kappa."
  Psychological bulletin 72.5 (1969): 323.
- Fleiss1971: Fleiss, Joseph L.
  "Measuring nominal scale agreement among many raters."
  Psychological bulletin 76.5 (1971): 378.


"""

from __future__ import annotations

import numpy as np

from ...intermediates.estimator_parts import NijMatrix, Weights

_ERROR_PREFIX = "[scotts_pi] Error"

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
    the data is balanced and generalises it otherwise -- the generalisation
    being Landis and Koch's, not Fleiss'. An analyser that promises *exact*
    Fleiss must therefore require balanced data; this function will not stop
    it, because refusing here would also refuse the coefficients that have no
    such requirement.

    Subjects rated fewer than twice are skipped. A single rating forms no pair,
    so there is nothing for it to agree or disagree with, and dividing by
    n(n-1) = 0 would be undefined rather than merely uninformative.

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
    formula disagree about what they are averaging over. It also matches Gwet's
    implementation, which matters for anyone comparing.

    The choice is invisible on balanced data. An
    analyser computing *Fleiss'* kappa should require balance and sidestep the
    question; one that accepts unbalanced data is choosing among
    generalisations Fleiss did not give, and should say so.

    Length K, taken from the input rather than the data, so a category nobody
    chose still has an entry holding zero. Gwet's AC1 and Brennan-Prediger put
    the number of categories in their chance term, so dropping an empty one
    would move them.

    Subjects rated fewer than twice are skipped, as they are in the observed
    agreement -- the two must average over the same subjects or they are not
    describing the same data.

    Guaranteed to have K components, since an NijMatrix is K columns wide by
    construction and a category nobody chose is a column of zeros rather than a
    column that is missing.

    Weights do not enter. These are proportions of ratings, and how much credit
    a disagreement earns has nothing to say about how often a category was
    chosen.
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
    category, Pe is around 0.82 and kappa is left dividing by 0.18 -- so a
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

    Note K appearing explicitly. That is why an unused category cannot be
    quietly dropped from the scheme: doing so changes this number, while
    leaving Fleiss' untouched. NijMatrix is K columns wide for that reason.

    **Takes no weights, because there is no weighted form of it here.** The
    weighted version of AC1 is Gwet's AC2, a separate coefficient with a
    published definition of its own. Handing weights to this one is a
    TypeError, which is the right answer: the alternative is a weighted
    observed agreement over an unweighted chance term, which is neither
    coefficient and looks entirely plausible.
    """
    pj = category_proportions(nij)
    k = nij.n_categories
    return float((pj * (1.0 - pj)).sum() / (k - 1))


def chance_brennan_prediger(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Pe as Brennan and Prediger: 1/K, regardless of what the raters did.

    The most severe answer available. It assumes chance means picking
    uniformly among the categories on offer, so the observed proportions do not
    enter at all. Unaffected by skew because it never looks at it, which is
    either the point or the problem depending on whether you think a category
    nobody uses should still count as a possibility.

    With weights the formula becomes sum(W) / K^2, where W is the K x K weight
    matrix and K is the number of categories. Identity weights reduce this to
    K / K^2 = 1/K, so the unweighted case is a special case of the weighted
    one. Off-diagonal entries inflate the sum, reflecting that near-misses earn
    partial credit and therefore raise the baseline that chance alone would
    produce.

    The weighted formula is from Gwet (2014), implemented in irrCAC's
    ``bp.coeff.raw`` as ``pe <- sum(weights.mat) / (q^2)``.
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
    agreement, there is no room above it to have taken up. Bootstrap draws
    reach this state on skewed data, which is why the resampler counts NaNs
    rather than treating them as zeros.
    """
    if np.isclose(chance, 1.0):
        return float("nan")
    return float((observed - chance) / (1.0 - chance))


def percent_agreement(nij: NijMatrix, weights: Weights | None = None) -> float:
    """
    Po, uncorrected. The proportion of rater pairs that agreed.

    Not chance-corrected and not pretending to be. Worth reporting beside a
    kappa rather than instead of one: a kappa of 0.42 next to 93% agreement is
    telling you something about the prevalence, not about the raters.
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
    plainly that the case of varying numbers of raters is not considered. An
    analyser promising exact Fleiss must require balanced data; given
    unbalanced data this returns the Landis and Koch generalisation, which is a
    reasonable number but not the one Fleiss defined.

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
    measured across five real screening datasets, by between 0.17 and 0.36 --
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

    Independent of how the raters actually used the categories, which makes it
    stable under skew and insensitive to it in a way some find too blunt.
    Included mainly because it is the clearest demonstration that the chance
    term is a modelling choice rather than a fact: same Po, same data, a
    different number because of a different assumption.

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
    (sum of squared category proportions), so at two raters they are one
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
            f"{_ERROR_PREFIX}: Scott's pi is defined for nominal (unweighted) data only. "
            f"weights must be None or identity, got {weights.kind}."
        )
    
    return agreement_corrected_by_chance(
        observed_agreement(nij, None),  # Always unweighted
        chance_fleiss(nij, None)
    )
