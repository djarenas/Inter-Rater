# Example Data

Validation and test datasets for interrater.

- **V1** is interrater 1.x, as published. **V2** is this package.
- V2 bootstrap results use `n_resamples=4000`, `random_state=42` and
  `ci_method="bca"`, unless a row says otherwise.
- Where V2 gives a "parametric null" row, the interval comes from Fleiss'
  variance under the null hypothesis (Fleiss, Nee & Landis 1979, Eq. 12). It
  is narrower than the bootstrap interval on every dataset here with more
  than two raters.

## Folder layout

Every dataset lives in its own folder with the YAML it needs to load: a data
file, a semantics file and a cleaner config.

```
example_data/
  Fleiss1971/          counts only, as published (turned into Ratings in memory)
  Gwet_example/        ratings
  Synthetic_2raters/   ratings
  SR1_OSF2018/         prelim + main, plus permutation configs and weighted examples
  SR2_JSRC2019/        ratings
  SR3_JGIM2019/        ratings
  SR4_FIS_CVD/         ratings
  SR5_SciRep2020/      ratings
  Statistical_background/  ratings, every table on guide/statistical_background.md
```

**File conventions**

- `.csv` files have a header row. `.txt` files are tab-separated, with no header.
- A header row and `use_column_names` (in the cleaner config) are separate things:
  - `use_column_names: true`: column labels must match the rater names in the semantics.
  - `use_column_names: false`: raters are taken by position, and any labels are ignored.
  - The subject id column is named by `subject_id_col` and found either way.
- In practice here: the `.csv` files use `true`, and the headerless `.txt`
  files use `false`, since they have no labels to match.

## Original Research Data

### SR1_OSF2018/SR1_prelim_data.csv

Abstract screening, pilot round, from a systematic review of health outcomes
associated with food insecurity status in the USA (Arenas et al. 2018). Four
reviewers rated the same abstracts as include, maybe or exclude.

| Property | Value |
|---|---|
| Subjects | 199 |
| Raters | 4 (rater0, rater1, rater2, rater3) |
| Categories | 3: -1 (exclude), 0 (maybe), 1 (include) |
| Balance | Balanced (all subjects rated by all 4 raters) |
| Prevalence | Highly skewed toward exclude (~83% of ratings are -1) |

**Note:** the paper reports 200 pilot abstracts; this file has 199. One was
probably a duplicate and removed, and only the 199 were kept.

**Validation results:**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | ~0.49 | [0.32, 0.66] | 0.34 | V1's interval was computed incorrectly |
| V2 parametric null | 0.493 | [0.446, 0.540] | 0.094 | |
| V2 bootstrap | 0.493 | [0.383, 0.606] | 0.223 | BCa bootstrap |
| Gwet AC1 | 0.828 | [0.777, 0.869] | 0.092 | BCa bootstrap |
| Brennan-Prediger | 0.779 | [0.721, 0.828] | 0.107 | BCa bootstrap |
| Krippendorff's alpha | 0.493 | [0.383, 0.606] | 0.223 | BCa bootstrap, exact formula |

- The parametric null interval is 42% the width of the bootstrap interval.
- With 83% of ratings "exclude", Fleiss' κ (0.493) is far below B-P (0.779)
  and AC1 (0.828): its chance term grows with the skew.
- Krippendorff's α equals Fleiss' κ to 3 decimals. On balanced data, α's
  expected disagreement is Fleiss' times n/(n − 1) (Artstein & Poesio 2008,
  §2.6.1; Gwet 2015, Eq. 3), so α = κ + (1 − κ)/n, with n the number of
  ratings. Here n = 796, so the two differ by 0.0006.

**Weighted versions of the same data**

- Unless a section says otherwise, the categories in these examples are
  unordered: every disagreement counts the same, with no partial credit.
- SR1 prelim's categories have a natural order (include, maybe, exclude),
  so two extra semantics files analyse the same ratings with weights, each
  with its own reading of what "maybe" meant. The data
  file and cleaner config are the same.
- In a semantics file, `kind: ordered` gives linear or quadratic weights,
  and `kind: weighted` takes a custom matrix.

| Semantics file | Weights | Fleiss' κ | Brennan-Prediger | Krippendorff's α |
|---|---|---|---|---|
| `SR1_prelim_semantics.yaml` | none (unordered) | 0.493 [0.383, 0.606] | 0.779 [0.721, 0.828] | 0.493 [0.383, 0.606] |
| `SR1_prelim_ordered_semantics.yaml` | linear | 0.554 [0.440, 0.669] | 0.748 [0.683, 0.808] | 0.555 [0.440, 0.669] |
| `SR1_prelim_weighted_semantics.yaml` | custom (`SR1_prelim_weights.yaml`) | 0.549 [0.438, 0.661] | 0.728 [0.658, 0.791] | 0.549 [0.438, 0.661] |

(95% BCa intervals, `n_resamples=4000`, `random_state=42`.)

- **Linear** reads "maybe" as a 50/50 call, halfway between: include–maybe and maybe–exclude each keep 0.5 of the credit.
- **Custom** reads it by what happened to the abstract. A maybe went on to full-text review, where it was read carefully to rule out exclusion, so it behaved much like an include: include–maybe 0.8, maybe–exclude 0.2, include–exclude 0.0. The matrix and its reasoning are in `SR1_prelim_weights.yaml`. `SR1_prelim_weighted_inline_semantics.yaml` holds the same matrix written inline and gives identical results.

Weighting raises Fleiss' κ and Krippendorff's α but lowers Brennan-Prediger. That is expected: B-P's chance term is the average weight, sum(W)/K², so partial credit raises chance agreement as well as observed agreement. In this package, Gwet's AC1, Scott's pi, Conger's kappa and Light's kappa have no weighted form and refuse the ordered and weighted files; use `SR1_prelim_semantics.yaml` for them.

---

### SR1_OSF2018/SR1_main_data.txt

Abstract screening, main round, from the same review (Arenas et al. 2018):
the remaining 1,617 abstracts, each rated by at least two of the four
reviewers.

| Property | Value |
|---|---|
| Subjects | 1617 |
| Raters | 4 (rater0, rater1, rater2, rater3) |
| Categories | 3: -1 (exclude), 0 (maybe), 1 (include) |
| Format | Tab-separated, no header |
| Missing marker | `o` (rater did not rate this subject) |
| Balance | Partially crossed: each subject rated by 2 of the 4 raters (3 subjects by 3). 50% of cells filled |
| Prevalence | Highly skewed toward exclude (~90.2% of ratings are -1) |
| Rater workloads | rater0: 652, rater1: 861, rater2: 860, rater3: 864 |

**Validation results:**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.62 | [0.55, 0.69] | 0.14 | |
| V2 parametric null | — | — | — | not computed: needs the same number of raters per subject (here 2–3) |
| V2 bootstrap | 0.619 | [0.553, 0.680] | 0.127 | BCa bootstrap |
| Gwet AC1 | 0.926 | [0.911, 0.939] | 0.028 | BCa bootstrap, Gwet (2014) formula for missing ratings |
| Brennan-Prediger | 0.898 | [0.879, 0.916] | 0.037 | BCa bootstrap, Gwet (2014) formula for missing ratings |
| Krippendorff's alpha | 0.618 | [0.553, 0.679] | 0.126 | BCa bootstrap, exact formula (unbalanced data) |

- The data is unbalanced (1,614 subjects rated twice, 3 rated three
  times), so each coefficient uses its formula for that case: Fleiss' κ
  the approximation that weights each subject equally, B-P and AC1 Gwet's
  (2014) formula for missing ratings, and α its own definition, which
  needs no balance.
- With 90% of ratings "exclude", Fleiss' κ (0.619) and α (0.618) are far
  below B-P (0.898) and AC1 (0.926).

---

### SR2_JSRC2019/SR2_data.txt

Abstract screening from a systematic review of patient-centered needs
assessments performed by free health clinics (Arenas et al. 2019a).

| Property | Value |
|---|---|
| Subjects | 81 |
| Raters | 4 (rater0, rater1, rater2, rater3) |
| Categories | 3: -1 (exclude), 0 (uncertain), 1 (include) |
| Format | Tab-separated, no header |
| Missing marker | `o` (rater did not rate this subject) |
| Balance | Nearly balanced: 80 subjects rated by all 4 raters, 1 subject by 2 |
| Prevalence | Highly skewed toward exclude (~80.4% of ratings are -1) |

**Validation results:**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.42 | [0.34, 0.52] | 0.18 | |
| V2 parametric null | — | — | — | not computed: needs the same number of raters per subject (here 2–4) |
| V2 bootstrap | 0.431 | [0.282, 0.607] | 0.325 | BCa bootstrap |
| Gwet AC1 | 0.782 | [0.694, 0.849] | 0.155 | BCa bootstrap, Gwet (2014) formula for missing ratings |
| Brennan-Prediger | 0.725 | [0.623, 0.802] | 0.179 | BCa bootstrap, Gwet (2014) formula for missing ratings |

- With 80% of ratings "exclude", Fleiss' κ (0.431) is far below B-P (0.725)
  and AC1 (0.782).

---

### SR3_JGIM2019/SR3_data.txt

Abstract screening from a systematic review and meta-analysis of depression,
anxiety and sleep disorders in US adults with food insecurity (Arenas et al.
2019b). Three reviewers screened each of the 770 abstracts.

| Property | Value |
|---|---|
| Subjects | 770 |
| Raters | 3 (rater0, rater1, rater2) |
| Categories | 2: -1 (exclude), 1 (include) |
| Format | Tab-separated, no header |
| Balance | Balanced (all subjects rated by all 3 raters) |
| Prevalence | Highly skewed toward exclude (~90.3% of ratings are -1) |
| Notes | Binary: include/exclude, no uncertain |

**Validation results:**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.68 | [0.63, 0.72] | 0.09 | |
| V2 parametric null | 0.675 | [0.634, 0.716] | 0.082 | |
| V2 bootstrap | 0.675 | [0.602, 0.740] | 0.138 | BCa bootstrap |
| Gwet AC1 | 0.931 | [0.911, 0.945] | 0.034 | BCa bootstrap |
| Brennan-Prediger | 0.886 | [0.856, 0.908] | 0.052 | BCa bootstrap |

- The parametric null interval is 59% the width of the bootstrap interval.
- With 90% of ratings "exclude", Fleiss' κ (0.675) is far below B-P (0.886)
  and AC1 (0.931).

---

### SR4_FIS_CVD/SR4_data.txt

Abstract screening from a systematic search on food insecurity and
cardiometabolic risk factors (PROSPERO CRD42020149560). The same screening
was used for three reviews: hypertension (Beltrán et al. 2020), type 2
diabetes (Beltrán et al. 2022) and dyslipidemia (Arenas et al. 2022). Four
reviewers screened each of the 787 abstracts.

| Property | Value |
|---|---|
| Subjects | 787 |
| Raters | 4 (rater0, rater1, rater2, rater3) |
| Categories | 2: -1 (exclude), 1 (include) |
| Format | Tab-separated, no header |
| Balance | Balanced (all subjects rated by all 4 raters) |
| Prevalence | Highly skewed toward exclude (~82.1% of ratings are -1) |
| Notes | Binary: include/exclude, no uncertain |

**Validation results:**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.68 | [0.66, 0.71] | 0.05 | |
| V2 parametric null | 0.708 | [0.679, 0.737] | 0.058 | |
| V2 bootstrap | 0.708 | [0.663, 0.750] | 0.087 | BCa bootstrap |
| Gwet AC1 | 0.879 | [0.853, 0.899] | 0.046 | BCa bootstrap |
| Brennan-Prediger | 0.828 | [0.798, 0.855] | 0.057 | BCa bootstrap |

- The parametric null interval is 67% the width of the bootstrap interval.
- With 82% of ratings "exclude", Fleiss' κ (0.708) is below B-P (0.828) and
  AC1 (0.879).

---

### SR5_SciRep2020/SR5_data.txt

Abstract screening from a systematic review and meta-analysis of cocaine,
cardiomyopathy and heart failure (Arenas et al. 2020). Three reviewers
screened each of the 881 abstracts.

| Property | Value |
|---|---|
| Subjects | 881 |
| Raters | 3 (rater0, rater1, rater2) |
| Categories | 2: -1 (exclude), 1 (include) |
| Format | Tab-separated, no header |
| Balance | Balanced (all subjects rated by all 3 raters) |
| Prevalence | Highly skewed toward exclude (~87.7% of ratings are -1) |
| Notes | Binary: include/exclude, no uncertain |

**Validation results:**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.51 | [0.47, 0.55] | 0.08 | |
| V2 parametric null | 0.509 | [0.471, 0.547] | 0.076 | |
| V2 bootstrap | 0.509 | [0.443, 0.574] | 0.131 | BCa bootstrap |
| Gwet AC1 | 0.865 | [0.840, 0.886] | 0.046 | BCa bootstrap |
| Brennan-Prediger | 0.788 | [0.755, 0.820] | 0.065 | BCa bootstrap |

- The parametric null interval is 58% the width of the bootstrap interval.
- With 88% of ratings "exclude", Fleiss' κ (0.509) is far below B-P (0.788)
  and AC1 (0.865).

---

### Fleiss1971/Fleiss1971_psychiatric_data_nij.csv

The classic Fleiss (1971) psychiatric diagnosis dataset, as counts.

| Property | Value |
|---|---|
| Subjects | 30 |
| Raters | 6 psychiatrists per subject, identities not given |
| Categories | 5: Depression, Personality_disorder, Schizophrenia, Neurosis, Other |
| Format | CSV with header |
| Columns | (index), Depression, Personality_disorder, Schizophrenia, Neurosis, Other |
| Row sums | All equal 6 |
| Source | Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. *Psychological Bulletin*, 76(5), 378--382. |

The paper gives counts only: for each patient, how many of the six
psychiatrists chose each diagnosis, not which one. No ratings file is kept,
because one would have to invent who rated what. To analyse it,
`example_data/made_up_ratings.py` turns the counts into a Ratings in memory,
with made-up rater slots. Use the result only with coefficients that treat
raters as interchangeable (Fleiss, Scott, AC1, Brennan-Prediger,
Krippendorff).

**Reference values:**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| Fleiss (1971) original | 0.430 | | | Its variance formula was corrected by Fleiss, Nee & Landis (1979) |
| V2 parametric null | 0.430 | [0.382, 0.478] | 0.096 | |
| V2 bootstrap | 0.430 | [0.335, 0.549] | 0.214 | BCa bootstrap |
| Gwet AC1 | 0.448 | [0.349, 0.564] | 0.215 | BCa bootstrap |
| Brennan-Prediger | 0.444 | [0.344, 0.561] | 0.217 | BCa bootstrap, Pe = 1/5 = 0.2 |

- The parametric null interval is 45% the width of the bootstrap interval.
- The ratings are spread over 5 categories with no strong skew, so the
  chance terms are similar and so are the coefficients: Fleiss' κ (0.430),
  B-P (0.444), AC1 (0.448).

---

### Synthetic_2raters/synthetic_2raters_data.txt

Synthetic dataset for testing 2-rater agreement coefficients.

| Property | Value |
|---|---|
| Subjects | 100 |
| Raters | 2 (rater0, rater1) |
| Categories | 3: -1 (exclude), 0 (uncertain), 1 (include) |
| Format | Tab-separated, no header |
| Balance | Balanced (all subjects rated by both raters) |
| Prevalence | Moderate: 1 (include) ~48.0%, -1 (exclude) ~36.5%, 0 (uncertain) ~15.5% |
| Notes | Synthetic, not from a real study |

**Validation results (V2):**

| Coefficient | Method | Value | 95% CI | CI Width |
|---|---|---|---|---|
| Cohen's kappa | bootstrap | 0.641 | [0.497, 0.757] | 0.260 |
| Fleiss' kappa | parametric null | 0.641 | [0.493, 0.788] | 0.295 |
| Fleiss' kappa | bootstrap | 0.641 | [0.499, 0.759] | 0.260 |
| Scott's pi | parametric null | 0.641 | [0.493, 0.788] | 0.295 |
| Scott's pi | bootstrap | 0.641 | [0.499, 0.759] | 0.260 |
| Gwet AC1 | bootstrap | 0.683 | [0.546, 0.786] | 0.240 |
| Brennan-Prediger | bootstrap | 0.670 | [0.535, 0.775] | 0.240 |

- Fleiss' κ and Scott's π are the same coefficient for two raters on
  balanced data.
- Cohen's κ agrees with them to 3 decimals (0.64082 vs 0.64073). The two
  raters' shares are close (rater0: 0.36, 0.15, 0.49; rater1: 0.37, 0.16,
  0.47 for exclude, uncertain, include), so their chance terms are close.
- Here the parametric null interval is wider than the bootstrap interval.
- With no strong skew, the coefficients are close: Cohen's and Fleiss' κ
  (0.641), B-P (0.670), AC1 (0.683).

---

## Gwet Handbook Data

Source: Gwet, K. L. (2014). *Handbook of Inter-Rater Reliability* (4th ed.). Advanced Analytics, LLC.

### Gwet_example/gwet_cont3x3_data.csv

Ratings format of the cont3x3abstractors dataset.

| Property | Value |
|---|---|
| Subjects | 100 |
| Raters | 2 (Rater1, Rater2) |
| Categories | 3: Ectopic, AIU, NIU |
| Balance | Balanced (all subjects rated by both raters) |
| Source table | Gwet Handbook Table 4.1 (cont3x3abstractors) |
| Domain | Emergency room pregnancy classification |

Expanded from the contingency table below: rater identity is real (rows are
Rater1, columns Rater2), subject order is not.

**Original contingency table:**

|  | Ectopic | AIU | NIU |
|---|---|---|---|
| Ectopic | 13 | 0 | 0 |
| AIU | 0 | 20 | 7 |
| NIU | 0 | 4 | 56 |

**Reference values (irrCAC, asymptotic CI):**

| Coefficient | Value | 95% CI | Pa | Pe |
|---|---|---|---|---|
| Gwet AC1 | 0.84933 | [0.764, 0.935] | 0.890 | 0.270 |
| Cohen kappa | 0.79641 | | | |
| Scott pi | 0.79624 | | | |
| Brennan-Prediger | 0.83500 | | | |

**V2 results (bootstrap CI):**

| Coefficient | Value | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| Cohen's kappa | 0.796 | [0.647, 0.893] | 0.246 | BCa bootstrap |
| Gwet AC1 | 0.849 | [0.745, 0.919] | 0.174 | BCa bootstrap |
| Fleiss' kappa | 0.796 | [0.647, 0.893] | 0.246 | BCa bootstrap |
| Scott's pi | 0.796 | [0.647, 0.893] | 0.246 | BCa bootstrap |
| Brennan-Prediger | 0.835 | [0.715, 0.910] | 0.195 | BCa bootstrap |
| Krippendorff's alpha | 0.797 | [0.649, 0.894] | 0.245 | BCa bootstrap, exact formula |

- V2 point estimates match irrCAC: Cohen's κ, AC1, Scott's π (equal to
  Fleiss' κ for two raters on balanced data) and B-P.
- Cohen's κ agrees with Fleiss' κ to 3 decimals: the two raters' shares are
  close, so their chance terms are close.
- Krippendorff's α is 0.001 above Fleiss' κ: on balanced data
  α = κ + (1 − κ)/n (see SR1 prelim), with n = 200 ratings here.
- For AC1, irrCAC's asymptotic interval (width 0.171) and the V2 bootstrap
  interval (width 0.174) are about the same width.

---

## Derived files

### SR1_OSF2018/SR1_prelim_data_with_consensus.csv

SR1 prelim plus a sixth column, `consensus`, holding the adjudicated decision
for each subject. Columns: `(internal) id, rater0..rater3, consensus`.

- Not used by any test.
- Not an agreement dataset as it stands: consensus is an outcome, not a
  fifth rater, and treating it as one would count the raters who produced
  it twice.
- It is here for work on how agreement relates to the final decision.

---

## Statistical_background/

Every table on `guide/statistical_background.md`, as ratings files that load
through the cleaner. Synthetic: built from the page's cells, or simulated
with seed 42.

| Files | Page section |
|---|---|
| `professors_data.csv` | Percent agreement; Cohen's kappa |
| `strict_grader_data.csv` | Pooled or individual rates; Scott's quote |
| `pathologists_data.csv`, `prevalence_*_data.csv` | When the condition is rare |
| `biopsies_2raters_data.csv`, `biopsy_panel_data.csv` | More than two raters |
| `abstracts_data.csv` | When not everyone rates everything |

- `generate_tables.py` builds the CSVs. `--check` fails if a committed file
  differs from what it builds.
- `reproduce_page.py` runs interrater on each table and checks the page's
  numbers at the precision the page prints. The simulated panel is seeded,
  so its observed columns are checked the same way; its theoretical columns
  are computed from the design, and each theoretical kappa must lie inside
  the observed 95% interval.
- Each prevalence table has the fewest slides at which every cell is a whole
  number, so the sweep is exact.
- `tests/test_example_data/test_statistical_background.py` runs both.

---

## References

- Arenas, D. J., Zhou, S., Thomas, A., Wang, J., Vilá-Arroyo, G., & Bash, K.
  (2018). Negative health outcomes associated with food insecurity status in
  the United States of America: A systematic review of peer-reviewed
  studies. *OSF Preprints (NutriXiv).*
  https://osf.io/preprints/nutrixiv/uk9xw_v1
- Arenas, D. J., Noche, R., Thomas, A., Hallak, D., & Raman, S. (2019a).
  Systematic review of patient-centered needs assessments performed by free
  health clinics. *Journal of Student-Run Clinics*, 5(1).
  https://doi.org/10.59586/jsrc.v5i1.79
- Arenas, D. J., Thomas, A., Wang, J., & DeLisser, H. M. (2019b). A
  systematic review and meta-analysis of depression, anxiety, and sleep
  disorders in US adults with food insecurity. *Journal of General Internal
  Medicine*, 34, 2874–2882. https://doi.org/10.1007/s11606-019-05202-4
- Arenas, D. J., Beltran, S., Zhou, S., & Goldberg, L. R. (2020). Cocaine,
  cardiomyopathy, and heart failure: A systematic review and meta-analysis.
  *Scientific Reports*, 10, 19795. https://doi.org/10.1038/s41598-020-76273-1
- Arenas, D. J., Beltrán, S., Montgomery, C., Pharel, M., Lopez-Hinojosa, I.,
  Vilá-Arroyo, G., & DeLisser, H. M. (2022). A systematic review and
  meta-analysis of food insecurity and dyslipidemia. *Journal of the American
  Board of Family Medicine*, 35(4), 656–667.
  https://doi.org/10.3122/jabfm.2022.04.210413
- Artstein, R., & Poesio, M. (2008). Inter-coder agreement for computational
  linguistics. *Computational Linguistics*, 34(4), 555–596.
  https://doi.org/10.1162/coli.07-034-R2
- Beltrán, S., Pharel, M., Montgomery, C. T., López-Hinojosa, I. J., Arenas,
  D. J., & DeLisser, H. M. (2020). Food insecurity and hypertension: A
  systematic review and meta-analysis. *PLOS ONE*, 15(11), e0241628.
  https://doi.org/10.1371/journal.pone.0241628
- Beltrán, S., Arenas, D. J., Pharel, M., Montgomery, C., Lopez-Hinojosa, I.,
  & DeLisser, H. M. (2022). Food insecurity, type 2 diabetes, and
  hyperglycaemia: A systematic review and meta-analysis. *Endocrinology,
  Diabetes & Metabolism*, 5, e00315. https://doi.org/10.1002/edm2.315
- Gwet, K. L. (2015). *On Krippendorff's alpha coefficient* (revised
  October 5, 2015). Advanced Analytics, LLC.
  https://www.agreestat.com/papers/onkrippendorffalpha_rev10052015.pdf
