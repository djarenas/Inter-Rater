# agreement_functions/

Coefficient arithmetic for inter-rater agreement. Called by analysers, not by
users directly.

Four modules, one contract: each function takes an ``NijMatrix`` or an
``EstimatorInput`` and returns a float. None of them validate -- whether a
coefficient applies to the data is the analyser's question, asked once before
any of this runs. That division is deliberate: these functions are called
inside bootstrap loops thousands of times, and re-checking per draw would be
both slow and redundant. It also means they are not safe to call on arbitrary
input, which is why they are not exported at the package top level.

| module | what it holds | input type |
| --- | --- | --- |
| ``group_agreement`` | coefficients that treat raters as interchangeable | ``NijMatrix`` |
| ``individual_marginals`` | coefficients that give each rater their own distribution | ``EstimatorInput`` |
| ``krippendorff`` | Krippendorff's alpha, built from disagreement | ``NijMatrix`` |
| ``rater_combinations`` | running a coefficient over subsets of the raters | ``EstimatorInput`` |

The input-type split is structural, not arbitrary. ``NijMatrix`` holds counts
-- how many raters put each subject in each category -- and discards which
rater said what. ``EstimatorInput`` keeps the full array of rater-by-subject
assignments. A coefficient that reads only counts (Fleiss, AC1,
Brennan--Prediger, alpha) takes the matrix. One that reads rater identity
(Cohen, Conger, Light) takes the input. The rater-combinations module takes
the input because it needs to drop rater columns before anything is counted.


## Notation

Following Fleiss (1971):

| symbol | meaning |
| --- | --- |
| *N* | number of subjects |
| *n<sub>i</sub>* | number of raters who rated subject *i* |
| *K* | number of categories the rating task offered |
| *n<sub>ij</sub>* | number of raters who put subject *i* in category *j* |
| *p<sub>j</sub>* | share of ratings that went to category *j* |
| *P<sub>o</sub>* | observed agreement |
| *P<sub>e</sub>* | agreement expected by chance |


## The kappa skeleton

Every coefficient in ``group_agreement`` and ``individual_marginals`` has the
same shape:

```
coefficient = (Po - Pe) / (1 - Pe)
```

The numerator is how much agreement was observed beyond chance; the denominator
is how much room there was above chance to begin with. The coefficient is the
fraction of the available room that was actually taken up.

``agreement_corrected_by_chance(observed, chance)`` computes this, returning
NaN when Pe = 1 -- the case where chance alone would have produced perfect
agreement and there is no room above it to measure.

### Observed agreement

Computed once and identically for all group coefficients. For subject *i* rated
by *n<sub>i</sub>* raters there are *n<sub>i</sub>*(*n<sub>i</sub>* - 1) ordered
pairs of ratings, of which Sum_j *n<sub>ij</sub>*(*n<sub>ij</sub>* - 1) agree:

```
P_i = Sum_j n_ij (n_ij - 1) / [ n_i (n_i - 1) ]

Po  = mean of P_i over subjects
```

Subjects rated fewer than twice are skipped: a single rating forms no pair, so
dividing by n_i(n_i - 1) = 0 would be undefined rather than uninformative.

### Four chance terms

Only the chance term differs. That is the whole of the disagreement between
these coefficients:

| coefficient | *P<sub>e</sub>* | what "by chance" means |
| --- | --- | --- |
| percent agreement | 0 | nothing is chance |
| Fleiss' kappa | Sum *p<sub>j</sub>*^2 | raters assign at the observed rates and ignore the subject |
| Gwet's AC1 | Sum *p<sub>j</sub>*(1 - *p<sub>j</sub>*) / (*K* - 1) | chance comes from rater uncertainty, highest for categories used about half the time |
| Brennan--Prediger | 1 / *K* (unweighted) or Sum(*W*) / *K*^2 (weighted) | chance means picking uniformly among the categories on offer |

On one real dataset -- 881 abstracts, three raters, 12% included -- all four
share Po = 0.894 and give:

| coefficient | Pe | value |
| --- | --- | --- |
| percent agreement | 0.000 | 0.894 |
| Fleiss' kappa | 0.784 | 0.509 |
| Brennan--Prediger | 0.500 | 0.788 |
| Gwet's AC1 | 0.216 | 0.865 |

Identical ratings, a spread of 0.36. The chance term is a modelling assumption,
not a fact about the data, and the four answers are four different assumptions
rather than four measurements.

### Why kappa falls under skew

Fleiss' Pe = Sum *p<sub>j</sub>*^2 grows as the ratings concentrate in one
category. At 90% in a single category it is around 0.82, so kappa divides by
1 - 0.82 = 0.18 -- a small denominator that makes a small change in observed
agreement move the result a long way.

The consequence is that kappa can rank datasets differently from raw agreement.
Across five real screening datasets, the review with the *lowest* observed
agreement of the top three has the *highest* kappa, because its ratings were
less concentrated. Gwet's AC1 was built to avoid this and preserves the
agreement ordering.

Neither is wrong. They answer different questions, and reporting Po and Pe
alongside any kappa is what makes the difference visible.

### Scott's pi

Scott's pi (Scott 1955) uses the same chance term as Fleiss: the sum of
squared category proportions. At two raters, Fleiss' kappa equals Scott's pi
exactly -- they are the same formula, generalised and then specialised back.

Scott's pi is defined for nominal (unweighted) data only. ``scotts_pi(nij,
weights)`` with non-identity weights raises ValueError. For ordinal agreement
with pooled marginals, use weighted Fleiss' kappa instead.


## Weighting

On an ordered scale, disagreeing by one category is less wrong than disagreeing
by three. A weight matrix says how much credit a disagreement still earns:
entry (*i*, *j*) is what a rater who said *i* gets when another said *j*, with
1.0 on the diagonal.

Unweighted agreement is not a separate idea. It is weighting with an identity
matrix, where every disagreement costs everything. That is why there are not
weighted and unweighted versions of each formula -- there is one formula, and
identity is the default position on it:

```
Po = ( c_i . W . c_i  -  n_i ) / [ n_i (n_i - 1) ]     for subject i

Pe = p . W . p                                          Fleiss' chance term
```

where *c<sub>i</sub>* is subject *i*'s row of counts. The subtraction removes
the pairing of a rating with itself, which is not a pair of raters. Put
*W* = *I* and both reduce exactly to the unweighted forms above.

The two usual weightings derive from rank distance:

| | weight for categories *d* steps apart | neighbours on a 4-point scale |
| --- | --- | --- |
| linear | 1 - *d* / *d*<sub>max</sub> | 0.67 |
| quadratic | 1 - *d*^2 / *d*<sub>max</sub>^2 | 0.89 |

Quadratic is the usual choice in the medical literature and forgives near
misses much more readily. Neither is more correct, and the choice changes the
coefficient.

### Which coefficients accept weights

Not all of them, and the signatures say so rather than a docstring:

| coefficient | weights | why |
| --- | --- | --- |
| percent agreement | yes | it is Po, which weights |
| Fleiss' kappa | yes | Pe = *p* . *W* . *p* |
| Brennan--Prediger | yes | Pe = Sum(*W*) / *K*^2, from Gwet (2014) |
| Scott's pi | **no** | defined for nominal data only; ValueError on non-identity weights |
| Gwet's AC1 | **no** | the weighted version is AC2, a separate coefficient |
| Cohen's kappa | yes | Pe = *p<sub>a</sub>* . *W* . *p<sub>b</sub>* |
| Conger's kappa | **no** | its weighted form has a definition of its own; NotImplementedError |
| Light's kappa | **no** | no weighted form in Light (1971) |
| Krippendorff's alpha | yes | cost = 1 - *W*, entering both Do and De |

Brennan--Prediger's weighted chance term Sum(*W*) / *K*^2 has a clean
interpretation: with identity weights Sum(*W*) = *K*, so it reduces to 1 / *K*.
Off-diagonal entries inflate the sum, reflecting that near-misses earn partial
credit and raise the baseline chance alone would produce. Unlike Fleiss'
weighted chance term, this one is still independent of the observed proportions
-- the hallmark of B-P.


## Pooled vs individual marginals

Every group coefficient uses *pooled* marginals: the raters are treated as one
distribution, so a difference in how often each rater says "include" is counted
as disagreement.

The alternative is *individual* marginals -- letting each rater keep their own
rate and subtracting the difference out first. That produces different
coefficients, and they live in ``individual_marginals`` because they need rater
identity, which the counts have discarded:

|  | pooled | individual |
| --- | --- | --- |
| two raters | Scott's pi | Cohen's kappa |
| many raters | Fleiss' kappa | Conger's kappa |

The gap between the pooled and individual versions measures how differently the
raters draw the line, as distinct from how erratic they are. Those are different
problems with different remedies, and a single coefficient cannot tell them
apart.


## Individual marginals coefficients

These live in ``individual_marginals.py`` and take an ``EstimatorInput``
rather than an ``NijMatrix``, because they read the array -- which rater said
what -- not just the counts.

### Rater marginals

``rater_marginals(data, rows=None)`` computes p_rj: how often each rater chose
each category, as a proportion. Returns an array of size n_raters by K. The
``rows`` parameter restricts the calculation to a subset of subjects, which is
what Cohen's kappa needs -- it is defined on the subjects both raters saw, and
a rater's marginal must be computed on those subjects rather than on everything
they ever rated.

### Cohen's kappa

Agreement between two raters, each keeping their own rates. Cohen (1960), and
Cohen (1968) for the weighted form.

Assumes exactly two rater columns and does not check -- the analyser extracts
the pair with ``select_raters`` and is responsible for handing over the right
shape. Defined on the subjects both raters rated; subjects only one saw are
dropped via ``shared_rows`` before anything is computed. Returns NaN when the
pair shares fewer than two subjects.

The chance term is the product of the two raters' own proportions:
Pe = Sum_j p_aj * p_bj (unweighted), or p_a . W . p_b (weighted).

### Conger's kappa

The many-rater individual-marginals counterpart to Fleiss' kappa. Conger
(1980). Observed agreement is Fleiss', unchanged; only the chance term differs.
The chance term takes the mean category rate across raters and subtracts the
variance of those rates, which is what removes the effect of raters having
different thresholds.

Takes no weights -- its weighted form has a definition of its own and is not
implemented here. Passing one raises NotImplementedError.

### Light's kappa

The mean of all C(n, 2) pairwise Cohen's kappas. Light (1971). Each pair gets
its own Cohen's kappa, and the headline number is the unweighted mean. Conger's
kappa approximates this mean but is not identical to it; comparing the two says
whether the approximation holds for a particular panel.

Takes no weights, following Light (1971).


## Krippendorff's alpha

Alpha is built the other way up. Where the kappa family measures agreement and
corrects for chance, alpha measures *disagreement* and compares it with the
disagreement chance would produce:

```
alpha = 1 - Do / De
```

Those are not the same expression rearranged. Two differences make them
genuinely separate.

### The unit is a pair of ratings, not a subject

Alpha pools every pairable pair of ratings across the whole dataset and asks
how often those pairs disagree. Fleiss computes agreement within each subject
and then averages over subjects, weighting every subject equally. On balanced
data the two routes almost coincide; on unbalanced data they separate, because
a subject rated by six raters offers fifteen pairs and one rated by two offers
one.

### The chance term is a sampling correction, not a model of rater behaviour

Fleiss asks what agreement would look like if raters assigned at the observed
rates and ignored the subject. Alpha asks what disagreement would look like if
the same multiset of ratings were shuffled across subjects at random, which
gives a finite-sample correction with an ``n - 1`` in the denominator. On
small datasets that correction is visible; on large ones it vanishes.

``pairable_ratings(nij)`` counts the ratings that belong to a subject rated
more than once -- alpha's sample size, and not the same as the total number of
ratings. A subject rated once contributes a rating but no pair, so it cannot
disagree with anything and takes no part in the calculation.

### Weights

Alpha has a published weighted form: the cost of confusing two categories is
one minus their weight, entering both Do and De identically. Identity weights
give the nominal case.

### Relationship to kappa on balanced data

Expect alpha and Fleiss' kappa to be close but not equal, and to separate as
the data gets thinner or more unbalanced. On real balanced datasets with 770+
pairable ratings, the two match to 3+ decimal places. On a 200-rating balanced
dataset the finite-sample correction produces a gap of about 0.002. On
unbalanced data (pair-pooling vs subject-averaging) the gap can be larger. On
the Fleiss (1971) 30-patient data (180 pairable ratings), alpha is 0.433 vs
kappa 0.430 -- the correction has the same bite as on the 200-rating synthetic
data.


## Rater combinations

Running a coefficient over subsets of the raters, rather than over all of them.

A panel coefficient is one number for the whole group. It cannot tell you that
one rater was the problem, because it has already averaged that rater away. On
one real screening round the six pairs span 0.58 to 0.82 with a single
reviewer accounting for the spread, and the panel coefficient is 0.71 with no
hint that anything is uneven.

### The Estimator callable

``Estimator = Callable[[EstimatorInput], float]`` -- the type alias used
throughout this module. The estimator is whatever coefficient is being
measured: swapping Cohen's kappa for Fleiss' changes what is being measured and
not how the enumeration works. An analyser using a ``group_agreement``
coefficient (which takes an ``NijMatrix``) supplies a small module-level adapter
rather than a lambda -- these results are often computed inside a resampling
loop, and a closure cannot be shipped to a worker process.

### Subsets and pairs

``rater_subsets(n_raters, size)`` enumerates every combination of ``size``
rater columns as sorted tuples of indices. Combinations, not permutations: a
subset of raters has no order.

``over_subsets(data, estimator, size)`` runs the estimator on every subset.
``over_pairs(data, estimator)`` is the size-two case, which is the common one
and the only size on which two-rater coefficients (Cohen's kappa, Scott's pi)
are defined at all.

### Leave-one-rater-out

``leave_one_rater_out(data, estimator)`` runs the estimator on the panel with
each rater removed in turn. Keyed by the rater who was dropped, not by the
subset that remained, because the question is about that rater. Subtracting
each value from the full-panel coefficient gives a jackknife influence measure.

That is a different question from pairwise analysis. A rater can agree
reasonably with everybody and still be the one holding the panel coefficient
down, and a rater can disagree sharply with one colleague without moving the
panel at all.

### Sparse-data survival

Dropping raters guts sparse data, and the damage is uneven. Every subject left
with fewer than two ratings contributes to nothing, so the coefficient for a
subset is computed on whatever survived. Measured on one real dataset of 1,617
subjects, each rated by two of four reviewers, removing a rater leaves between
650 and 861 subjects holding a single rating -- and *which* rater was removed
decides how many.

``subjects_per_subset(data, size)`` and ``subjects_per_leave_one_out(data)``
return the count of subjects still rated at least twice after each
extraction -- the denominator every coefficient was actually computed on, and
the thing that makes results comparable across subsets or reveals that they
are not.

### The MAX_SUBSETS guard

Subsets grow as a binomial coefficient: ten raters give 45 pairs and 252
subsets of size five, but twenty raters give 184,756 of size ten.
``MAX_SUBSETS = 10_000`` is a guard against asking for a number of coefficients
nobody intended. It is deliberately loose.

### Summarising

``mean_over_subsets(values)`` averages a subset dictionary, skipping NaN
entries. The values are *not* independent -- two subsets sharing a rater are
positively correlated -- so the variance of this mean is not the variance of
an average of independent things. A bootstrap that resamples subjects and
recomputes every subset within each draw preserves the correlation without
having to model it. The function returns a point estimate and leaves intervals
to the resampler.


## Unbalanced data: the p_j choice

Fleiss (1971) defines his coefficient for the case where every subject is rated
the same number of times. The introduction says plainly that varying numbers
are not considered, and the closing paragraph of Fleiss, Nee and Landis (1979)
credits that case to Landis and Koch instead.

Applying the formula anyway means choosing among generalisations Fleiss did not
give. One choice in particular is silent and consequential: how do you compute
p_j when subjects have different numbers of raters?

```
pooled per rating      p_j = Sum_i n_ij / Sum_i n_i
averaged per subject   p_j = (1/N) Sum_i ( n_ij / n_i )
```

Pooling weights every *rating* equally, so a subject rated by four raters pulls
twice as hard as one rated by two. Averaging weights every *subject* equally,
however many raters saw it.

**This package averages per subject.** The reason is internal consistency
rather than authority: Po is already a mean over subjects with each weighted
equally, so weighting Pe by rating would have the two halves of one formula
disagree about what they are averaging over. It also matches Gwet's
implementation, which matters for anyone comparing numbers across packages.

### How much it matters

Nothing, when the data is balanced -- the two definitions are algebraically
identical when every n_i is the same, which is why every implementation agrees
there.

On unbalanced data it depends on the data, and the range is wide:

| data | gap in kappa |
| --- | --- |
| 81 subjects, one of them with fewer raters | 0.01 |
| 1,617 subjects, three of them with an extra rater | 0.0004 |
| constructed: 1,000 subjects x 2 raters, 20 x 50 raters, leaning opposite ways | **0.55** |

That last case gives 0.80 under one definition and 0.25 under the other, on
identical ratings. Two different conclusions, not a rounding difference.

**The mechanism says when to worry.** Pooling lets a handful of heavily rated
subjects dominate p_j; averaging gives every subject one vote. The gap opens
only when rater counts vary widely **and** the thinly and heavily rated
subjects were rated differently. Either condition alone is harmless, which is
why real screening data shows a gap of 0.01 and a constructed case shows 0.55.

An analyser computing *Fleiss'* kappa should require balanced data and sidestep
the question entirely. One that accepts unbalanced data is choosing a
generalisation, and should say which.


## Subjects rated once

A subject rated by a single rater forms no pair. The formula divides by
n_i(n_i - 1), which is zero -- undefined rather than merely uninformative --
so such subjects are excluded from both Po and p_j.

Both, not one. If the two averaged over different sets of subjects they would
not be describing the same data.


## Validation

Every coefficient here is checked against published worked examples and against
an independent implementation:

- **Fleiss (1971) Table 1** -- 30 patients, 6 raters, 5 diagnoses. The paper
  reports kappa = 0.430; this implementation gives 0.4302.
- **Gwet's AC1**, transcribed literally from the AHRQ methods appendix
  (NBK82266) -- agrees to machine precision on crossed data.
- **``irrCAC``**, the Python port of Gwet's R library -- Fleiss, AC1 and
  Brennan--Prediger across five real datasets, fifteen comparisons, all agreeing
  to within that library's five-digit rounding. Weighted Fleiss with quadratic
  weights agrees on the same five.

That cross-check has already earned its keep once: it found a disagreement on
the two unbalanced datasets and nowhere else, which turned out to be the
pooled-versus-subject-weighted choice described above. The implementation was
changed to match.


## References

- Brennan, R. L., & Prediger, D. J. (1981). Coefficient kappa: some uses,
  misuses, and alternatives. *Educational and Psychological Measurement*,
  41(3), 687--699.
- Cohen, J. (1960). A coefficient of agreement for nominal scales.
  *Educational and Psychological Measurement*, 20(1), 37--46.
- Cohen, J. (1968). Weighted kappa: nominal scale agreement with provision for
  scaled disagreement or partial credit. *Psychological Bulletin*, 70(4),
  213--220.
- Conger, A. J. (1980). Integration and generalization of kappas for multiple
  raters. *Psychological Bulletin*, 88(2), 322--328.
- Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters.
  *Psychological Bulletin*, 76(5), 378--382.
- Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). Large sample standard
  errors of kappa and weighted kappa. *Psychological Bulletin*, 72(5),
  323--327.
- Fleiss, J. L., Nee, J. C. M., & Landis, J. R. (1979). Large sample variance
  of kappa in the case of different sets of raters. *Psychological Bulletin*,
  86(5), 974--977.
- Gwet, K. L. *Handbook of Inter-Rater Reliability.* Advanced Analytics.
- Krippendorff, K. (2004). *Content Analysis: An Introduction to Its
  Methodology.* 2nd ed., Sage, chapter 11.
- Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
  Departmental Papers, Annenberg School for Communication, University of
  Pennsylvania.
- Landis, J. R., & Koch, G. G. (1977). The measurement of observer agreement
  for categorical data. *Biometrics*, 33(1), 159--174.
- Light, R. J. (1971). Measures of response agreement for qualitative data:
  some generalizations and alternatives. *Psychological Bulletin*, 76(5),
  365--377.
- Scott, W. A. (1955). Reliability of content analysis: the case of nominal
  scale coding. *Public Opinion Quarterly*, 19(3), 321--325.
