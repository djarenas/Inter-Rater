# Data Profiles

Structural metadata about the ratings data. Profiles describe the shape
of the data -- how subjects and raters are distributed -- so that
downstream code (reports, visualization, interpretation) can adapt
without re-inspecting the raw arrays.

| Class               | What It Holds                                              | Why It Exists                                                                 |
|---------------------|------------------------------------------------------------|-------------------------------------------------------------------------------|
| `AgreementProfile`  | Category prevalence, category-specific agreement, per-rater marginals, and the prevalence/bias indices, on the subjects rated at least twice | A coefficient is one number; this is what a reader needs to understand it. |
| `BalanceProfile`    | Whether ratings are balanced, per-subject counts, optional rater-level stats | Reports and visualizations need to know data shape, not just coefficient values. |
| `RaterSharedMatrix` | (n_raters, n_raters) symmetric matrix of shared subject counts | Overlap diagnostics: which rater pairs share enough subjects for a stable estimate. |
| `PrevalenceBiasIndex` | Two raters, two categories: the four cells, the prevalence index and the bias index, with the positive category and rater order | Shows why kappa is low despite high agreement: a dominant category (prevalence) or one rater saying "positive" more often (bias). |

## AgreementProfile

What the data says about agreement before any coefficient, computed on the
subjects every coefficient uses (rated at least twice):

- `prevalence` -- the package's p_j, averaged per subject: the same numbers
  Fleiss', Scott's and Gwet's chance terms use, so the profile explains the
  coefficient rather than describing different data.
- `specific_agreement` -- for each category, the share of rater pairs
  involving it that agree on it. Unweighted. On balanced data,
  Po = Sum_j p_j * s_j; on unbalanced data that identity does not hold,
  because Po averages per subject and this pools pairs.
- `rater_marginals` -- each rater's shares, as Cohen's and Conger's chance
  terms use them; `NotComputed` from counts.
- `prevalence_bias` -- a `PrevalenceBiasIndex` for two raters and two
  categories, otherwise `NotComputed` with the reason.

It prints two labelled tables: how often each category was used (an
"overall" row, then one row per rater under it), and how often raters agreed
on each category. On unbalanced data the overall row is the per-subject
average the coefficients use, so it is not exactly the mean of the rater
rows.

Built with `build_from_ratings()` or `build_from_nij_matrix()`. A nij matrix carries no
category names, so categories are labelled "column 0", "column 1", ...

## BalanceProfile

A frozen dataclass with two construction paths:

- **`build_from_ratings()`** -- all fields are exact. Rater identity is known,
  so `n_raters`, `n_per_rater`, and `fill_rate` are populated.
- **`build_from_nij_matrix()`** -- subject-level fields are exact, but
  rater-level fields (`n_raters`, `n_per_rater`, `fill_rate`) are `None`.
  Rater identity is gone in a counts matrix; honest absence is better
  than fabricated data. Use `n_raters_lower_bound` (always available)
  for what is known.

## PrevalenceBiasIndex

Byrt, Bishop and Carlin (1993), for two raters and two categories. The signs
follow R's epiR:

- prevalence index = (a - d) / N: above zero, more agreement on the positive
  category than on the negative one;
- bias index = (b - c) / N = P(rater_b says positive) - P(rater_a says
  positive): above zero, rater_b says positive more often.

"Positive" is the first declared category. Both the category and the rater
order travel with the numbers and are printed in words, because the signs
mean nothing without them. The four cells are stored and the indices worked
out from them, so they cannot disagree.

Two construction paths, as for `BalanceProfile`:

- **`build_from_pair()`** -- two raters' answers on the subjects both rated. Both
  indices are computed.
- **`build_from_counts()`** -- per-subject category counts. The prevalence index is
  computed; the bias index is `NotComputed`, because counts do not say which
  rater gave which answer.

For two categories, PABAK (prevalence- and bias-adjusted kappa) equals
2 Po - 1, which is unweighted Brennan-Prediger, so it has no field of its own.

## RaterSharedMatrix

Entry (i, j) is how many subjects raters i and j both rated. Diagonal
is each rater's total. Provides `thinnest_pair`, `below(threshold)`,
and `render()` for diagnostics.

## Dependencies

- `array_utils` (is_balanced, n_per_subject, n_per_rater, subject_count_distribution)
