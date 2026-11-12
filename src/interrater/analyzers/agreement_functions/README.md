# agreement_functions/

Coefficient arithmetic for inter-rater agreement. They are meant to be called by analysers,
and not directly by users.

Four modules, one contract: each function takes an ``NijMatrix`` or an
``EstimatorInput`` and returns a float. None of them validate: they run inside
bootstrap loops thousands of times, and the analyser has already checked the
data once.

| module | what it holds | input type |
| --- | --- | --- |
| ``group_agreement`` | coefficients that treat raters as interchangeable | ``NijMatrix`` |
| ``individual_marginals`` | coefficients that give each rater their own distribution | ``EstimatorInput`` |
| ``krippendorff`` | Krippendorff's alpha, built from disagreement | ``NijMatrix`` |
| ``rater_combinations`` | running a coefficient over subsets of the raters | ``EstimatorInput`` |

``NijMatrix`` holds counts (how many raters put each subject in each
category) and drops which rater said what. ``EstimatorInput`` keeps the full
rater-by-subject array. Coefficients that read only counts (Fleiss, AC1,
Brennan-Prediger, alpha) take the matrix; those that read rater identity
(Cohen, Conger, Light) take the input. ``rater_combinations`` takes the input
because it drops rater columns before anything is counted.


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


## The kappa form

Every coefficient in ``group_agreement`` and ``individual_marginals`` is

```
coefficient = (Po - Pe) / (1 - Pe)
```

the share of the agreement possible above chance that was observed.
``agreement_corrected_by_chance(observed, chance)`` computes it and returns
NaN when Pe = 1.

### Observed agreement

The same for every group coefficient. Subject *i*, rated by *n<sub>i</sub>*
raters, has *n<sub>i</sub>*(*n<sub>i</sub>* - 1) ordered pairs of ratings, of
which Sum_j *n<sub>ij</sub>*(*n<sub>ij</sub>* - 1) agree:

```
P_i = Sum_j n_ij (n_ij - 1) / [ n_i (n_i - 1) ]

Po  = mean of P_i over subjects
```

Subjects rated fewer than twice are skipped: they form no pair, and the
denominator is zero.

### Four chance terms

Only the chance term differs:

| coefficient | *P<sub>e</sub>* | what "by chance" means |
| --- | --- | --- |
| percent agreement | 0 | nothing is chance |
| Fleiss' kappa | Sum *p<sub>j</sub>*^2 | raters assign at the observed rates and ignore the subject |
| Gwet's AC1 | Sum *p<sub>j</sub>*(1 - *p<sub>j</sub>*) / (*K* - 1) | chance comes from rater uncertainty, highest for categories used about half the time |
| Brennan-Prediger | 1 / *K* (unweighted) or Sum(*W*) / *K*^2 (weighted) | picking uniformly among the categories on offer |

On one real dataset (881 abstracts, three raters, 12% included), all four
share Po = 0.894:

| coefficient | Pe | value |
| --- | --- | --- |
| percent agreement | 0.000 | 0.894 |
| Fleiss' kappa | 0.784 | 0.509 |
| Brennan-Prediger | 0.500 | 0.788 |
| Gwet's AC1 | 0.216 | 0.865 |

The same ratings give values 0.36 apart, because each chance term is a
different assumption.

### Kappa under skew

Fleiss' Pe = Sum *p<sub>j</sub>*^2 grows as ratings concentrate in one
category. At 90% in one category it is about 0.82, so kappa divides by 0.18,
and a small change in observed agreement moves it a long way. Kappa can
therefore rank datasets differently from raw agreement: across five real
screening datasets, the one with the lowest observed agreement of the top
three has the highest kappa, because its ratings were less concentrated.
Gwet's AC1 keeps the agreement ordering. Reporting Po and Pe with any kappa
shows which effect is at work.

### Scott's pi

Scott (1955). Same chance term as Fleiss; at two raters Fleiss' kappa equals
Scott's pi. Nominal only: ``scotts_pi(nij, weights)`` with non-identity
weights raises ValueError. For ordinal data with pooled marginals, use
weighted Fleiss' kappa.


## Weighting

A weight matrix gives the credit a disagreement still earns on an ordered
scale: entry (*i*, *j*) is the credit when one rater said *i* and another
said *j*, with 1.0 on the diagonal. Unweighted is the identity matrix, so
each formula has one form:

```
Po = ( c_i . W . c_i  -  n_i ) / [ n_i (n_i - 1) ]     for subject i

Pe = p . W . p                                          Fleiss' chance term
```

where *c<sub>i</sub>* is subject *i*'s row of counts; subtracting *n_i*
removes each rating paired with itself. With *W* = *I* both reduce to the
unweighted forms above.

The two usual weightings, from rank distance *d*:

| | weight | neighbours on a 4-point scale |
| --- | --- | --- |
| linear | 1 - *d* / *d*<sub>max</sub> | 0.67 |
| quadratic | 1 - *d*^2 / *d*<sub>max</sub>^2 | 0.89 |

Quadratic is the usual choice in the medical literature and gives more credit
to near misses. The choice changes the coefficient.

### Which coefficients accept weights

| coefficient | weights | why |
| --- | --- | --- |
| percent agreement | yes | it is Po, which weights |
| Fleiss' kappa | yes | Pe = *p* . *W* . *p* |
| Brennan-Prediger | yes | Pe = Sum(*W*) / *K*^2, from Gwet (2014) |
| Scott's pi | **no** | nominal only; ValueError on non-identity weights |
| Gwet's AC1 | **no** | the weighted version is AC2, a separate coefficient |
| Cohen's kappa | yes | Pe = *p<sub>a</sub>* . *W* . *p<sub>b</sub>* |
| Conger's kappa | **no** | its weighted form has its own definition; NotImplementedError |
| Light's kappa | **no** | no weighted form in Light (1971) |
| Krippendorff's alpha | yes | cost = 1 - *W*, in both Do and De |

Brennan-Prediger's weighted Pe = Sum(*W*) / *K*^2 reduces to 1 / *K* with
identity weights (Sum(*W*) = *K*). Partial credit raises it, and it still does
not depend on the observed proportions.


## Pooled vs individual marginals

Group coefficients pool the raters into one distribution, so a rater who says
"include" more often than another adds to the disagreement. Individual
marginals give each rater their own rate and remove that difference first.
They need rater identity, so they live in ``individual_marginals``:

|  | pooled | individual |
| --- | --- | --- |
| two raters | Scott's pi | Cohen's kappa |
| many raters | Fleiss' kappa | Conger's kappa |

The gap between the pooled and individual versions measures how differently
the raters draw the line, as distinct from how erratic they are.


## Individual marginals coefficients

In ``individual_marginals.py``; they take an ``EstimatorInput``.

### Rater marginals

``rater_marginals(data, rows=None)`` gives p_rj, how often each rater chose
each category, as an n_raters x K array. ``rows`` restricts it to some
subjects: Cohen's kappa is defined on the subjects both raters saw, so each
rater's marginal is computed on those.

### Cohen's kappa

Two raters, each with their own rates. Cohen (1960); Cohen (1968) for the
weighted form. Assumes exactly two rater columns (the analyser extracts the
pair with ``select_raters``). Subjects only one rater saw are dropped via
``shared_rows``; NaN when the pair shares fewer than two subjects. Chance term:
Pe = Sum_j p_aj * p_bj, or p_a . W . p_b weighted.

### Conger's kappa

Conger (1980), the many-rater counterpart. Observed agreement is Fleiss';
the chance term is the mean category rate across raters minus the variance of
those rates, which removes the effect of raters having different thresholds.
No weights (NotImplementedError).

### Light's kappa

Light (1971): the unweighted mean of all C(n, 2) pairwise Cohen's kappas.
Conger's kappa approximates this mean; comparing the two shows whether the
approximation holds for a panel. No weights.

### The most agreement the marginals allow (kappa_max)

With each rater's rates fixed, observed agreement cannot exceed
Sum_j min(p_aj, p_bj) (Cohen 1960), so a pair whose rates differ cannot reach
kappa = 1. ``cohen_kappa_max`` puts that ceiling on the kappa scale,
(Po_max - Pe) / (1 - Pe), with Cohen's chance term.

With weights, the ceiling is the largest weighted agreement of any table whose
rows match rater a's rates and whose columns match rater b's: a small
transportation problem, solved exactly as a linear program by
``max_observed_agreement``. With identity weights it is the sum of minima.

It is defined for Cohen's kappa and each pair in Light's. Conger's kappa with
three or more raters has no standard closed form for it, and the pooled or
marginal-free coefficients (Fleiss, Scott, Gwet, Brennan-Prediger) have none.


## Krippendorff's alpha

Alpha measures disagreement and compares it with the disagreement chance would
produce:

```
alpha = 1 - Do / De
```

It differs from the kappa family in two ways.

- **The unit is a pair of ratings.** Alpha pools every pairable pair across
  the dataset; Fleiss averages agreement within each subject, so every subject
  counts equally. On balanced data the two nearly coincide; on unbalanced
  data they separate, since a subject rated by six raters gives fifteen pairs
  and one rated by two gives one.
- **The chance term is a sampling correction.** De is the disagreement if
  the same ratings were shuffled across subjects, which gives an ``n - 1`` in
  the denominator: visible on small datasets, negligible on large ones.

``pairable_ratings(nij)`` counts the ratings on subjects rated more than once,
which is alpha's sample size. A subject rated once has no pair and takes no
part. Do and De both take their n from it, so the two terms cannot disagree
about the sample size.

**Weights.** The cost of confusing two categories is one minus their weight,
in both Do and De. Identity weights give the nominal case.

**Alpha vs Fleiss' kappa.** Close but not equal, and further apart on thin or
unbalanced data. On balanced datasets with 770+ pairable ratings they match to
3+ decimals; on a 200-rating balanced dataset they differ by about 0.002. On
Fleiss (1971) (180 pairable ratings), alpha is 0.433 and kappa 0.430.


## Rater combinations

A coefficient over subsets of the raters. A panel coefficient averages every
rater in; on one real screening round the six pairs span 0.58 to 0.82, one
reviewer accounts for the spread, and the panel coefficient is 0.71.

### The Estimator callable

``Estimator = Callable[[EstimatorInput], float]``, the coefficient being
measured. Changing the estimator changes what is measured, not how subsets are
enumerated. A ``group_agreement`` coefficient (which takes an ``NijMatrix``) is
wrapped in a module-level adapter rather than a lambda, so it can be pickled
if the loops are ever run in parallel.

### Subsets and pairs

``rater_subsets(n_raters, size)`` lists every combination of ``size`` raters
as sorted index tuples. ``over_subsets(data, estimator, size)`` runs the
estimator on each; ``over_pairs(data, estimator)`` is size two, the only size
on which Cohen's kappa and Scott's pi are defined.

### When a subset has no coefficient

``unit_problem(present, names)`` is the one rule, and returns the reason when
there is no coefficient:

- fewer than two of the subset's subjects are usable (rated by at least two
  of its raters), or
- one of its raters rated none of the usable subjects; the subset is really a
  smaller one.

``estimate_unit(unit_data, estimator)`` applies it: NaN when the rule says no,
the estimator otherwise. ``over_subsets``, every bootstrap draw and every
jackknife value call it, so the values, the intervals and the mean agree on
which subsets count. The reason is what a permutation result's
``NotComputed`` says.

### Leave-one-rater-out

``leave_one_rater_out(data, estimator)`` runs the estimator with each rater
removed in turn, keyed by the rater dropped. Subtracting each value from the
full-panel coefficient gives a jackknife influence measure. It differs from
pairwise analysis: a rater can agree reasonably with everyone and still hold
the panel coefficient down, or disagree sharply with one colleague without
moving it.

### Subjects that survive

Removing raters leaves some subjects with fewer than two ratings, and they
count for nothing. On one real dataset (1,617 subjects, each rated by two of
four reviewers), removing a rater leaves 650 to 861 subjects with a single
rating, depending on which rater. ``subjects_per_subset(data, size)`` and
``subjects_per_leave_one_out(data)`` give the subjects still rated at least
twice, the count each coefficient was computed on.

### GUARD_MAX_SUBSETS

Ten raters give 45 pairs and 252 subsets of five; twenty raters give 184,756
subsets of ten. ``GUARD_MAX_SUBSETS = 10_000`` (in ``contracts``) stops a
request for more coefficients than anyone meant. It is deliberately loose.

### Summarising

``mean_over_subsets(values)`` averages a subset dictionary, skipping NaN. The
values are positively correlated when subsets share a rater, so this mean has
no simple variance; its interval comes from a bootstrap that recomputes every
subset in each draw. ``mean_over_same_subsets(values, subsets)`` is the
resampler's version: the mean over exactly the subsets defined at the point
estimate, or NaN (a failed draw) if any of them is undefined in that resample.


## Unbalanced data: the p_j choice

Fleiss (1971) assumes every subject is rated the same number of times, and
Fleiss, Nee and Landis (1979) leave the varying case to Landis and Koch. With
varying numbers of raters, p_j can be computed two ways:

```
pooled per rating      p_j = Sum_i n_ij / Sum_i n_i
averaged per subject   p_j = (1/N) Sum_i ( n_ij / n_i )
```

Pooling weights every rating equally, so a subject rated by four raters
counts twice as much as one rated by two. Averaging weights every subject
equally.

**This package averages per subject.** Po is already a mean over subjects,
so Pe should average over the same thing. It also matches Gwet's
implementation.

On balanced data the two are identical. On unbalanced data the gap depends on
the data:

| data | gap in kappa |
| --- | --- |
| 81 subjects, one of them with fewer raters | 0.01 |
| 1,617 subjects, three of them with an extra rater | 0.0004 |
| constructed: 1,000 subjects x 2 raters, 20 x 50 raters, leaning opposite ways | **0.53** |

The constructed case: of the 1,000 subjects with two raters, 950 got A from
both and 50 got one A and one B; each of the 20 subjects with fifty raters
got 5 A and 45 B. Averaging per subject, as this package does, gives 0.35;
pooling per rating gives 0.88 (pinned in
`test_agreement_functions.py::TestCategoryProportions`). The gap opens only
when rater counts vary widely and the thinly and heavily rated subjects were
rated differently; either alone is harmless.


## Subjects rated once

A subject with one rating forms no pair (n_i(n_i - 1) = 0). It is excluded
from both Po and p_j, so the two describe the same subjects.


## Validation

Checked against published examples and an independent implementation:

- **Fleiss (1971) Table 1**: 30 patients, 6 raters, 5 diagnoses. Published
  kappa 0.430; this implementation 0.4302.
- **Gwet's AC1**, transcribed from the AHRQ methods appendix (NBK82266):
  agrees to machine precision on crossed data.
- **``irrCAC``**, the Python port of Gwet's R library: Fleiss, AC1 and
  Brennan-Prediger on five real datasets, fifteen comparisons, all within that
  library's five-digit rounding. Weighted Fleiss (quadratic) agrees on the
  same five.

The ``irrCAC`` comparison first disagreed on the two unbalanced datasets only;
the cause was the pooled-versus-per-subject p_j choice above, and the
implementation was changed to match.


## References

- Brennan, R. L., & Prediger, D. J. (1981). Coefficient kappa: some uses,
  misuses, and alternatives. *Educational and Psychological Measurement*,
  41(3), 687-699.
- Cohen, J. (1960). A coefficient of agreement for nominal scales.
  *Educational and Psychological Measurement*, 20(1), 37-46.
- Cohen, J. (1968). Weighted kappa: nominal scale agreement with provision for
  scaled disagreement or partial credit. *Psychological Bulletin*, 70(4),
  213-220.
- Conger, A. J. (1980). Integration and generalization of kappas for multiple
  raters. *Psychological Bulletin*, 88(2), 322-328.
- Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters.
  *Psychological Bulletin*, 76(5), 378-382.
- Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). Large sample standard
  errors of kappa and weighted kappa. *Psychological Bulletin*, 72(5),
  323-327.
- Fleiss, J. L., Nee, J. C. M., & Landis, J. R. (1979). Large sample variance
  of kappa in the case of different sets of raters. *Psychological Bulletin*,
  86(5), 974-977.
- Gwet, K. L. *Handbook of Inter-Rater Reliability.* Advanced Analytics.
- Krippendorff, K. (2004). *Content Analysis: An Introduction to Its
  Methodology.* 2nd ed., Sage, chapter 11.
- Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
  Departmental Papers, Annenberg School for Communication, University of
  Pennsylvania.
- Landis, J. R., & Koch, G. G. (1977). The measurement of observer agreement
  for categorical data. *Biometrics*, 33(1), 159-174.
- Light, R. J. (1971). Measures of response agreement for qualitative data:
  some generalizations and alternatives. *Psychological Bulletin*, 76(5),
  365-377.
- Scott, W. A. (1955). Reliability of content analysis: the case of nominal
  scale coding. *Public Opinion Quarterly*, 19(3), 321-325.
