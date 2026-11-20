# agreement_functions/

Coefficient arithmetic. They should only be called by the analyzers, and not manually by users.

| module | holds | input |
| --- | --- | --- |
| ``group_agreement`` | Fleiss, Scott, Gwet's AC1, Brennan-Prediger, percent agreement | ``EstimatorInput`` |
| ``krippendorff`` | Krippendorff's alpha | ``EstimatorInput`` |
| ``individual_marginals`` | Cohen, Conger, Light, kappa_max | ``EstimatorInput`` |


## Contract

- Each function returns a float, NaN when the coefficient is undefined.
- They are not meant to be called manually because nothing here validates. There is no validation because these functions run inside bootstrap loops. Analyzers check data before calling these functions.
- **One input.** Every coefficient takes an ``EstimatorInput``, so all of
  them are the same ``Callable[[EstimatorInput], float]``.
- **K** comes from ``data.n_categories``. Group coefficients count the data
  once (``data.nij``, K wide) and hand that NijMatrix to their parts (Po,
  Pe, Do, De), which take K from its width.
- **Weights** come from ``data.weights``; ``None`` means identity. A
  coefficient without a weighted form refuses non-identity weights.


## Notation

Following Fleiss (1971):

| symbol | meaning |
| --- | --- |
| *N* | subjects |
| *n<sub>i</sub>* | raters who rated subject *i* |
| *K* | categories the rating task offered |
| *n<sub>ij</sub>* | raters who put subject *i* in category *j* |
| *p<sub>j</sub>* | share of ratings in category *j* |
| *W* | weight matrix, 1 on the diagonal |
| *P<sub>o</sub>*, *P<sub>e</sub>* | observed and chance agreement |


## group_agreement

Every coefficient here is ``(Po - Pe) / (1 - Pe)``
(``agreement_corrected_by_chance``, NaN when Pe = 1). Po is the same for all;
only Pe differs.

```
Po_i = ( c_i . W . c_i  -  n_i ) / [ n_i (n_i - 1) ]     c_i = row i of the counts
Po   = mean of Po_i over subjects rated at least twice
```

With *W* = *I* this is Sum_j *n<sub>ij</sub>*(*n<sub>ij</sub>* - 1) / [*n<sub>i</sub>*(*n<sub>i</sub>* - 1)].

| coefficient | Pe | weights |
| --- | --- | --- |
| percent agreement | 0 | yes |
| Fleiss' kappa | *p* . *W* . *p* | yes |
| Scott's pi | Sum *p<sub>j</sub>*^2 | no, ValueError on non-identity |
| Gwet's AC1 | Sum *p<sub>j</sub>*(1 - *p<sub>j</sub>*) / (*K* - 1) | no (the weighted form is AC2) |
| Brennan-Prediger | Sum(*W*) / *K*^2, which is 1/*K* unweighted | yes |

At two raters Scott's pi equals Fleiss' kappa.

**p<sub>j</sub> is averaged per subject**, (1/*N*) Sum_i *n<sub>ij</sub>* / *n<sub>i</sub>*,
not pooled per rating. Po is a mean over subjects, so Pe is too. On balanced
data the two are identical. On unbalanced data they can differ a lot;
``test_agreement_functions.py::TestCategoryProportions`` pins a case where
kappa is 0.35 one way and 0.88 the other.

Subjects rated once form no pair and are left out of both Po and
*p<sub>j</sub>*.


## krippendorff

```
alpha = 1 - Do / De
```

- Do and De count disagreement over every pairable pair of ratings, not per
  subject, so a subject rated by six raters weighs more than one rated by two.
- De is the disagreement if the same ratings were shuffled across subjects;
  this puts an ``n - 1`` in its denominator.
- ``pairable_ratings(nij)`` is n: the ratings on subjects rated at least
  twice. Do and De both use it.
- Weights: the cost of a disagreement is 1 - *W*.


## individual_marginals

Each rater keeps their own category rates. ``rater_marginals(data, rows)``
gives them as an n_raters x K array, optionally on some rows only.

| coefficient | definition | weights |
| --- | --- | --- |
| Cohen's kappa | two raters, subjects both rated (``shared_rows``); Pe = *p<sub>a</sub>* . *W* . *p<sub>b</sub>*; NaN on fewer than 2 shared subjects | yes |
| Conger's kappa | Fleiss' Po; Pe = Sum *p̄<sub>j</sub>*^2 - Sum var<sub>r</sub>(*p<sub>rj</sub>*) / n_raters, sample variance | no, NotImplementedError |
| Light's kappa | mean of the pairwise Cohen's kappas | no |

Conger's and Light's assume complete data; the analyzer refuses anything
else.

**kappa_max.** With each rater's rates fixed, Po cannot exceed
Sum_j min(*p<sub>aj</sub>*, *p<sub>bj</sub>*) (Cohen 1960).
``cohen_kappa_max`` puts that on the kappa scale. With weights the ceiling is
a transportation problem, solved as a linear program by
``max_observed_agreement``. Defined for Cohen's kappa and each pair in
Light's, not for Conger's at three or more raters.


## Validation

- Fleiss (1971) Table 1: published kappa 0.430; here 0.4302.


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
- Gwet, K. L. (2014). *Handbook of Inter-Rater Reliability* (4th ed.).
  Advanced Analytics.
- Krippendorff, K. (2004). *Content Analysis: An Introduction to Its
  Methodology.* 2nd ed., Sage, chapter 11.
- Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
  Departmental Papers, Annenberg School for Communication, University of
  Pennsylvania.
- Light, R. J. (1971). Measures of response agreement for qualitative data:
  some generalizations and alternatives. *Psychological Bulletin*, 76(5),
  365-377.
- Scott, W. A. (1955). Reliability of content analysis: the case of nominal
  scale coding. *Public Opinion Quarterly*, 19(3), 321-325.
