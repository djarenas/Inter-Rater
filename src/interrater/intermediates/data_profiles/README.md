# Data Profiles

Structural metadata about the ratings. Profiles describe the shape
of the data (i.e. how subjects and raters are distributed). This helps 
downstream code (reports, visualization, interpretation) to adapt
without re-inspecting the raw arrays.

| Class               | What It Holds                                              | Why It Exists                                                                 |
|---------------------|------------------------------------------------------------|-------------------------------------------------------------------------------|
| `AgreementProfile`  | Category prevalence, category-specific agreement, per-rater marginals, and the prevalence/bias indices, on the subjects rated at least twice | A coefficient is one number; this is what a reader needs to understand it. |
| `BalanceProfile`    | Whether ratings are balanced, and how many ratings each subject got | Reports and visualizations need to know data shape, not just coefficient values. |
| `RaterProfile`      | The raters: names, how many subjects each rated, and their overlap (a `RaterSharedMatrix`) | Who rated what. Needs rater identity, so it exists only for data built from ratings. |
| `RaterSharedMatrix` | (n_raters, n_raters) symmetric matrix of shared subject counts | Overlap diagnostics: which rater pairs share enough subjects for a stable estimate. |
| `PrevalenceBiasIndex` | Two raters, two categories: the four cells, the prevalence index and the bias index, with the positive category and rater order | Shows why kappa is low despite high agreement: a dominant category (prevalence) or one rater saying "positive" more often (bias). |

## AgreementProfile

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

Computed on the subjects every coefficient uses (rated at least twice).

It prints two labelled tables: how often each category was used (an
"overall" row, then one row per rater under it), and how often raters agreed
on each category. On unbalanced data the overall row is the per-subject
average the coefficients use, so it is not exactly the mean of the rater
rows.

Built with `build_from_ratings()` or `build_from_nij_matrix()`. In the case of nij matrix, as it carries no category names, the categories are labelled "column 0", "column 1", ...

## BalanceProfile

The subject side only: `is_balanced`, `n_per_subject`, `subject_counts`
({raters: subjects}) and `raters_per_subject` (when balanced), with
`n_subjects` and `n_raters_lower_bound` as properties. All of that is in a
counts matrix as well as in ratings, so every result carries one. Two
construction paths, both exact:

- **`build_from_ratings()`**
- **`build_from_nij_matrix()`** -- a rater who rated nothing is invisible in
  counts, so `n_raters_lower_bound` is all a counts matrix can say about the
  raters.

Plotted by `BalanceProfilePlotter` (see `visualizers/README.md`).

## RaterProfile

The rater side: `rater_names`, `n_subjects`, `n_per_rater` and
`rater_overlap` (a `RaterSharedMatrix`), with `n_raters` and `fill_rate`
derived from them so they cannot disagree. It checks that the overlap's
diagonal is each rater's own total.

Built only with `build_from_ratings()`. A counts matrix keeps how many raters
chose each category, not who they were, so a result built from counts holds
`NotComputed("built from counts: no rater identity")` in its `rater_profile`.
Cohen's, Conger's and Light's kappas and every permutation analysis require
one.

The printed tables, `to_frame()`, summaries and plots name raters through
`rater_profile.rater_names`. Plotted by `RaterProfilePlotter` (see
`visualizers/README.md`): `draw_overlap()` and `draw_counts()`.

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
is each rater's total. Built with `build_from_ratings()`; a `RaterProfile`
holds one as `rater_overlap`. Provides `thinnest_pair`, `below(threshold)`,
and `render()` for diagnostics.

## Dependencies

- `array_utils` (is_balanced, n_per_subject, n_per_rater, subject_count_distribution, pairwise_overlap)
