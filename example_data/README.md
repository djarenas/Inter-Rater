# Example Data

Validation and test datasets for the interrater v2 package.

**All V2 bootstrap results use:** `n_resamples=4000`, `random_state=42`, `ci_method="bca"` unless otherwise noted. These are reference runs, not the settings the test suite uses -- `tests/test_example_data/` runs percentile at B=250 for speed, and asserts point estimates rather than CI bounds.

## Folder layout

Every dataset lives in its own folder with the YAML it needs to load. The
ratings pathway needs three files, the nij pathway needs one.

```
example_data/
  Fleiss1971/          nij only (the original data is counts)
  Gwet_example/        both pathways, for cross-validating them against each other
  Synthetic_2raters/   both pathways
  SR1_OSF2018/         prelim + main, plus permutation configs and weighted examples
  SR2_JSRC2019/        ratings
  SR3_JGIM2019/        ratings
  SR4_FIS_CVD/         ratings
  SR5_CocCardio/       ratings
  Unbalanced_nij/      nij only, the one unbalanced counts matrix
```

**File conventions.** `.csv` files have a header row; `.txt` files do not.
That is a property of the file, not of `CleanerConfig.frame_layout.use_column_names`
-- that flag controls whether DataFrame column labels must match declared
rater names. The Gwet CSV shows the two are independent: it has headers
(`Subject,Rater1,Rater2`) and `use_column_names: false`, because its raters
are positional.

## Original Research Data

### SR1_OSF2018/SR1_prelim_data.csv

Preliminary screening data from Systematic Review 1.

| Property | Value |
|---|---|
| Subjects | 199 |
| Raters | 4 (rater0, rater1, rater2, rater3) |
| Categories | 3: -1 (exclude), 0 (uncertain), 1 (include) |
| Balance | Balanced (all subjects rated by all 4 raters) |
| Prevalence | Highly skewed toward exclude (~83% of ratings are -1) |

**Validation results (Fleiss' kappa):**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published (buggy CI) | ~0.49 | [0.32, 0.66] | 0.34 | Wrong CI implementation |
| V2 parametric null | 0.493 | [0.446, 0.540] | 0.094 | Correct null variance (Fleiss 1979 Eq. 12), but too narrow |
| V2 bootstrap | 0.493 | [0.383, 0.606] | 0.223 | Empirical bootstrap CI, realistic width |
| Gwet AC1 | 0.828 | [0.777, 0.869] | 0.092 | BCa bootstrap |
| Brennan-Prediger | 0.779 | [0.721, 0.828] | 0.107 | BCa bootstrap |
| Krippendorff's alpha | 0.493 | [0.383, 0.606] | 0.223 | BCa bootstrap, exact formula |

The parametric null CI is ~42% the width of the bootstrap CI, illustrating that Fleiss's null variance underestimates uncertainty around the point estimate. Under heavy prevalence skew (~83% exclude), Fleiss' κ (0.493) is dragged down by inflated marginal-based chance, while B-P (0.779) and AC1 (0.828) remain high. Krippendorff's α (0.493) matches Fleiss κ exactly on this balanced dataset -- the two use different units and chance models but converge when every subject has the same number of raters.

**Weighted versions of the same data.** The categories have a natural order (include, uncertain, exclude), so the same ratings can also be analysed with weights. Two extra semantics files do this, each with its own reading of what "uncertain" meant. The data file and cleaner config are unchanged.

| Semantics file | Weights | Fleiss' κ | Brennan-Prediger | Krippendorff's α |
|---|---|---|---|---|
| `SR1_prelim_semantics.yaml` | none (unordered) | 0.493 [0.383, 0.606] | 0.779 [0.721, 0.828] | 0.493 [0.383, 0.606] |
| `SR1_prelim_ordered_semantics.yaml` | linear | 0.554 [0.440, 0.669] | 0.748 [0.683, 0.808] | 0.555 [0.440, 0.669] |
| `SR1_prelim_weighted_semantics.yaml` | custom (`SR1_prelim_weights.yaml`) | 0.549 [0.438, 0.661] | 0.728 [0.658, 0.791] | 0.549 [0.438, 0.661] |

(95% BCa intervals, `n_resamples=4000`, `random_state=42`.)

- **Linear** reads "uncertain" as a 50/50 call, halfway between: include–uncertain and uncertain–exclude each keep 0.5 of the credit.
- **Custom** reads it by what happened to the abstract. An uncertain went on to full-text review, where it was read carefully to rule out exclusion, so it behaved much like an include: include–uncertain 0.8, uncertain–exclude 0.2, include–exclude 0.0. The matrix and its reasoning are in `SR1_prelim_weights.yaml`. `SR1_prelim_weighted_inline_semantics.yaml` holds the same matrix written inline and gives identical results.

Weighting raises Fleiss' κ and Krippendorff's α but lowers Brennan-Prediger. That is expected: B-P's chance term is the average weight, sum(W)/K², so partial credit raises chance agreement as well as observed agreement. Gwet's AC1, Scott's pi, Conger's kappa and Light's kappa have no weighted form and refuse the ordered and weighted files; use `SR1_prelim_semantics.yaml` for them.

---

### SR1_OSF2018/SR1_main_data.txt

Main screening data from Systematic Review 1 (FIS + Psych combined).

| Property | Value |
|---|---|
| Subjects | 1617 |
| Raters | 4 (rater0, rater1, rater2, rater3) |
| Categories | 3: -1 (exclude), 0 (uncertain), 1 (include) |
| Format | Tab-separated, no header |
| Missing marker | `o` (rater did not rate this subject) |
| Balance | Partially crossed -- each subject rated by exactly 2 of the 4 raters (3 subjects have 3). 50% fill |
| Prevalence | Highly skewed toward exclude (~90.2% of ratings are -1) |
| Rater workloads | rater0: 652, rater1: 861, rater2: 860, rater3: 864 |

**Validation results (Fleiss' kappa):**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.62 | [0.55, 0.69] | 0.14 | |
| V2 parametric null | — | — | — | not computed: needs the same number of raters per subject (here 2–3) |
| V2 bootstrap | 0.619 | [0.553, 0.680] | 0.127 | BCa bootstrap, realistic width |
| Gwet AC1 | 0.926 | [0.911, 0.939] | 0.028 | BCa bootstrap, approximation (unbalanced) |
| Brennan-Prediger | 0.898 | [0.879, 0.916] | 0.037 | BCa bootstrap, approximation (unbalanced) |
| Krippendorff's alpha | 0.618 | [0.553, 0.679] | 0.126 | BCa bootstrap, exact formula (unbalanced data) |

**This is the only shipped dataset where a rater subset's measurable pool and
its complete cases differ**, and it is what the subset permutation tests exist
for. A subject contributes to a subset's agreement if at least two of that
subset's raters rated it; before 2.0.1 the bootstrap used only subjects *all*
k raters rated. On every other dataset here the two sets are identical, so
nothing caught it. At k=3 on this one:

| subset | measurable (>=2 of 3) | complete case (all 3) |
|---|---|---|
| (0, 1, 2) | 756 | 0 |
| (0, 1, 3) | 760 | 0 |
| (0, 2, 3) | 757 | 2 |
| (1, 2, 3) | 967 | 1 |

Three of the four had nothing to resample at all under the old rule. At k=2
the two pools coincide on all six pairs, as they must -- "at least two of the
k" is one rule, and at k=2 it means both.

Fleiss κ, B-P, and AC1 all use the unbalanced-data approximation on this dataset. Krippendorff's α (0.618) handles the unbalanced data exactly via pair-pooling -- no approximation needed. It is nearly identical to Fleiss κ (0.619) because only 3 of 1617 subjects deviate from 2 raters, so the data is barely unbalanced. Under extreme skew (~90% exclude), Fleiss' κ (0.619) and α (0.618) are much lower than B-P (0.898) and AC1 (0.926).

---

### SR2_JSRC2019/SR2_data.txt

Screening data from Systematic Review 2.

| Property | Value |
|---|---|
| Subjects | 81 |
| Raters | 4 (rater0, rater1, rater2, rater3) |
| Categories | 3: -1 (exclude), 0 (uncertain), 1 (include) |
| Format | Tab-separated, no header |
| Missing marker | `o` (rater did not rate this subject) |
| Balance | Nearly balanced -- 80 subjects rated by all 4 raters, 1 subject rated by 2 raters |
| Prevalence | Highly skewed toward exclude (~80.4% of ratings are -1) |

**Validation results (Fleiss' kappa):**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.42 | [0.34, 0.52] | 0.18 | |
| V2 parametric null | — | — | — | not computed: needs the same number of raters per subject (here 2–4) |
| V2 bootstrap | 0.431 | [0.282, 0.607] | 0.325 | BCa bootstrap, realistic width |
| Gwet AC1 | 0.782 | [0.694, 0.849] | 0.155 | BCa bootstrap, approximation (unbalanced) |
| Brennan-Prediger | 0.725 | [0.623, 0.802] | 0.179 | BCa bootstrap, approximation (unbalanced) |

Under skew (~80% exclude), Fleiss' κ (0.431) is much lower than B-P (0.725) and AC1 (0.782).

---

### SR3_JGIM2019/SR_FIS_Psych_data.txt

Full-text screening data from SR1, FIS + Psych databases combined.

| Property | Value |
|---|---|
| Subjects | 770 |
| Raters | 3 (rater0, rater1, rater2) |
| Categories | 2: -1 (exclude), 1 (include) |
| Format | Tab-separated, no header |
| Balance | Balanced (all subjects rated by all 3 raters) |
| Prevalence | Highly skewed toward exclude (~90.3% of ratings are -1) |
| Notes | No uncertain (0) ratings -- binary include/exclude only |

**Validation results (Fleiss' kappa):**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.68 | [0.63, 0.72] | 0.09 | |
| V2 parametric null | 0.675 | [0.634, 0.716] | 0.082 | Fleiss 1979 Eq. 12, but too narrow |
| V2 bootstrap | 0.675 | [0.602, 0.740] | 0.138 | BCa bootstrap, realistic width |
| Gwet AC1 | 0.931 | [0.911, 0.945] | 0.034 | BCa bootstrap |
| Brennan-Prediger | 0.886 | [0.856, 0.908] | 0.052 | BCa bootstrap |

The parametric null CI is ~59% the width of the bootstrap CI. With binary categories and ~90% exclude, Fleiss' κ (0.675) is again much lower than B-P (0.886) and AC1 (0.931).

---

### SR4_FIS_CVD/SR_FIS_CVD_data.txt

Full-text screening data from SR1, FIS CVD database.

| Property | Value |
|---|---|
| Subjects | 787 |
| Raters | 4 (rater0, rater1, rater2, rater3) |
| Categories | 2: -1 (exclude), 1 (include) |
| Format | Tab-separated, no header |
| Balance | Balanced (all subjects rated by all 4 raters) |
| Prevalence | Highly skewed toward exclude (~82.1% of ratings are -1) |
| Notes | No uncertain (0) ratings -- binary include/exclude only |

**Validation results (Fleiss' kappa):**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.68 | [0.66, 0.71] | 0.05 | |
| V2 parametric null | 0.708 | [0.679, 0.737] | 0.058 | Fleiss 1979 Eq. 12, but too narrow |
| V2 bootstrap | 0.708 | [0.663, 0.750] | 0.087 | BCa bootstrap, realistic width |
| Gwet AC1 | 0.879 | [0.853, 0.899] | 0.046 | BCa bootstrap |
| Brennan-Prediger | 0.828 | [0.798, 0.855] | 0.057 | BCa bootstrap |

The parametric null CI is ~67% the width of the bootstrap CI. Under ~82% exclude skew, Fleiss' κ (0.708) is lower than B-P (0.828) and AC1 (0.879).

---

### SR5_CocCardio/SR_cocaine_cardio_data.txt

Full-text screening data from SR2, Cocaine-Cardio database.

| Property | Value |
|---|---|
| Subjects | 881 |
| Raters | 3 (rater0, rater1, rater2) |
| Categories | 2: -1 (exclude), 1 (include) |
| Format | Tab-separated, no header |
| Balance | Balanced (all subjects rated by all 3 raters) |
| Prevalence | Highly skewed toward exclude (~87.7% of ratings are -1) |
| Notes | No uncertain (0) ratings -- binary include/exclude only |

**Validation results (Fleiss' kappa):**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| V1 published | 0.51 | [0.47, 0.55] | 0.08 | |
| V2 parametric null | 0.509 | [0.471, 0.547] | 0.076 | Fleiss 1979 Eq. 12, but too narrow |
| V2 bootstrap | 0.509 | [0.443, 0.574] | 0.131 | BCa bootstrap, realistic width |
| Gwet AC1 | 0.865 | [0.840, 0.886] | 0.046 | BCa bootstrap |
| Brennan-Prediger | 0.788 | [0.755, 0.820] | 0.065 | BCa bootstrap |

The parametric null CI is ~58% the width of the bootstrap CI. Under ~88% exclude skew, Fleiss' κ (0.509) is much lower than B-P (0.788) and AC1 (0.865).

---

### Fleiss1971/Fleiss1971_psychiatric_data_nij.csv

Nij (counts matrix) format of the classic Fleiss (1971) psychiatric diagnosis dataset.

| Property | Value |
|---|---|
| Subjects | 30 |
| Raters | 6 (psychiatrists) |
| Categories | 5: Depression, Personality_disorder, Schizophrenia, Neurosis, Other |
| Format | CSV with header |
| Columns | (index), Depression, Personality_disorder, Schizophrenia, Neurosis, Other |
| Row sums | All equal 6 (number of raters) |
| Source | Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters. *Psychological Bulletin*, 76(5), 378--382. |
| Notes | Nij format only -- no ratings format (original data is counts, not individual rater assignments) |

**Reference values:**

| Method | kappa | 95% CI | CI Width | Notes |
|---|---|---|---|---|
| Fleiss (1971) original | 0.430 | | | Variance formula in this paper is known to be wrong |
| V2 parametric null | 0.430 | [0.382, 0.478] | 0.096 | Fleiss 1979 Eq. 12, but too narrow |
| V2 bootstrap | 0.430 | [0.335, 0.549] | 0.214 | BCa bootstrap (n_bootstrap=4000), realistic width |
| Gwet AC1 | 0.448 | [0.349, 0.564] | 0.215 | BCa bootstrap (n_bootstrap=4000) |
| Brennan-Prediger | 0.444 | [0.344, 0.561] | 0.217 | BCa bootstrap (n_bootstrap=4000), Pe = 1/5 = 0.2 |

The parametric null CI is ~45% the width of the bootstrap CI. With 5 categories and no extreme skew toward any single category, all three coefficients are remarkably close: Fleiss κ (0.430), B-P (0.444), and AC1 (0.448) -- the chance terms diverge little when prevalence is spread across categories.

Note: irrCAC's built-in `distrib_6raters()` contains only the first 15 of the 30 subjects (κ = 0.414). This file contains all 30 subjects from the original paper.

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
| Notes | Synthetic data -- not from a real study |

**Validation results (V2):**

| Coefficient | Method | Value | 95% CI | CI Width |
|---|---|---|---|---|
| Cohen's kappa | bootstrap | 0.641 | [0.497, 0.757] | 0.260 |
| Fleiss' kappa | parametric null | 0.641 | [0.493, 0.788] | 0.295 |
| Fleiss' kappa | bootstrap | 0.641 | [0.499, 0.759] | 0.260 |
| Scott's pi | parametric null | 0.641 | [0.493, 0.788] | 0.295 |
| Scott's pi | bootstrap | 0.641 | [0.499, 0.759] | 0.260 |
| Gwet AC1 (ratings) | bootstrap | 0.683 | [0.546, 0.786] | 0.240 |
| Gwet AC1 (nij) | bootstrap | 0.683 | [0.546, 0.786] | 0.240 |
| Brennan-Prediger | bootstrap | 0.670 | [0.535, 0.775] | 0.240 |

Cohen's κ (0.641) matches Fleiss' κ and Scott's π to 3 decimal places on this balanced 2-rater data -- the individual marginals (p_a = [0.36, 0.15, 0.49], p_b = [0.37, 0.16, 0.47]) are close enough that Pe_Cohen ≈ Pe_Fleiss, so the coefficients converge. They differ at the 4th decimal place. Fleiss' kappa and Scott's pi are identical for balanced 2-rater data, as expected. With only 2 raters, the parametric null CI is actually *wider* than the bootstrap CI here -- the opposite of the pattern seen with more raters. Gwet AC1 from both pathways (ratings and nij matrix) gives the same point estimate and, at this seed, the same interval: both pathways resample the same subjects. With moderate prevalence (no extreme skew), all four chance-corrected coefficients are close: Cohen κ / Fleiss κ (0.641), B-P (0.670), and AC1 (0.683).

---

## Gwet Handbook Data

Source: Gwet, K. L. (2014). *Handbook of Inter-Rater Reliability* (4th ed.). Advanced Analytics, LLC.

These datasets are provided in two formats to validate both analyzer pathways:
the ratings format (subject x rater) for `estimate()`, and the nij format
(subject x category counts) for `estimate_from_nij_matrix()`.

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

### Gwet_example/gwet_cont3x3_data_nij.csv

Nij (counts matrix) format of the same dataset.

| Property | Value |
|---|---|
| Subjects | 100 |
| Columns | Subject, Ectopic, AIU, NIU |
| Row sums | All equal 2 (number of raters) |

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
| Cohen's kappa | 0.796 | [0.647, 0.893] | 0.246 | BCa bootstrap (n_bootstrap=4000) |
| Gwet AC1 | 0.849 | [0.745, 0.919] | 0.174 | BCa bootstrap (n_bootstrap=4000) |
| Fleiss' kappa | 0.796 | [0.647, 0.893] | 0.246 | BCa bootstrap (n_bootstrap=4000) |
| Scott's pi | 0.796 | [0.647, 0.893] | 0.246 | BCa bootstrap (n_bootstrap=4000) |
| Brennan-Prediger | 0.835 | [0.715, 0.910] | 0.195 | BCa bootstrap (n_bootstrap=4000) |
| Krippendorff's alpha | 0.797 | [0.649, 0.894] | 0.245 | BCa bootstrap (n_bootstrap=4000), exact formula |

V2 point estimates match irrCAC: Cohen κ 0.796 vs 0.79641, AC1 0.849 vs 0.84933, Fleiss κ 0.796 vs Scott π 0.79624 (equivalent for 2 balanced raters), B-P 0.835 vs 0.83500. Cohen's κ (0.796) matches Fleiss' κ (0.796) to 3 decimal places -- on this 2-rater data the individual marginals are close, so Pe_Cohen ≈ Pe_Fleiss. Krippendorff's α (0.797) is nearly identical to Fleiss κ (0.796) on this balanced 2-rater data, as expected -- the two use different units and chance models but converge when every subject has the same number of raters and the dataset is large enough. The irrCAC asymptotic CI [0.764, 0.935] (width 0.171) is close to the V2 bootstrap CI [0.745, 0.919] (width 0.174) for AC1.

---

## Alternate-format and derived files

Datasets above that ship in a second format, plus one derived file. These
exist so the two analyzer pathways can be validated against each other: the
same data through `estimate()` and `estimate_from_nij_matrix()` must give the
same point estimate.

### Fleiss1971/Fleiss1971_psychiatric_data_nij.txt

Tab-separated twin of the CSV, with a header row and no subject-ID column.
Same 30 x 5 counts. The CSV is what the tests load.

### Synthetic_2raters/synthetic_2raters_data_nij.csv

Nij form of the synthetic 2-rater data: 100 subjects x 3 categories, all row
sums 2. Columns are `0, 1, 2` -- category *indices*, not the `-1 / 0 / 1`
labels of the ratings file.

### SR1_OSF2018/SR1_prelim_data_with_consensus.csv

SR1 prelim plus a sixth column, `consensus`, holding the adjudicated decision
for each subject. Columns: `(internal) id, rater0..rater3, consensus`.

Not used by any test, and not an agreement dataset as it stands -- consensus
is an outcome, not a fifth rater, and treating it as one would double-count
the raters that produced it. It is here for work on how agreement relates to
the decision finally reached.

---

## Unbalanced_nij/unbalanced_data_nij.txt

Synthetic unbalanced counts matrix. The only dataset here where an nij row
sum varies.

| Property | Value |
|---|---|
| Subjects | 50 |
| Categories | 3 |
| Format | Tab-separated, no header |
| Row sums | 2, 3 or 4 -- subjects rated by different numbers of raters |
| Balance | Unbalanced |
| Notes | Synthetic -- not from a real study |

**No reference values yet.** This dataset has no published or internally
validated coefficient values, and no test uses it. The `unbalanced_nij`
fixture in `tests/test_example_data/conftest.py` loads it and nothing
consumes the fixture.

It covers a real gap -- the nij pathway on unbalanced data, where Fleiss, B-P
and AC1 all fall back to an approximation and Krippendorff does not. Closing
it means establishing reference values first, which for synthetic data means
deciding what they should be validated *against*.
