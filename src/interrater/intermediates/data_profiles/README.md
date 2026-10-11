# Data Profiles

Structural metadata about the ratings data. Profiles describe the shape
of the data -- how subjects and raters are distributed -- so that
downstream code (reports, visualization, interpretation) can adapt
without re-inspecting the raw arrays.

| Class               | What It Holds                                              | Why It Exists                                                                 |
|---------------------|------------------------------------------------------------|-------------------------------------------------------------------------------|
| `BalanceProfile`    | Whether ratings are balanced, per-subject counts, optional rater-level stats | Reports and visualizations need to know data shape, not just coefficient values. |
| `RaterSharedMatrix` | (n_raters, n_raters) symmetric matrix of shared subject counts | Overlap diagnostics: which rater pairs share enough subjects for a stable estimate. |

## BalanceProfile

A frozen dataclass with two construction paths:

- **`from_ratings()`** -- all fields are exact. Rater identity is known,
  so `n_raters`, `n_per_rater`, and `fill_rate` are populated.
- **`from_nij_matrix()`** -- subject-level fields are exact, but
  rater-level fields (`n_raters`, `n_per_rater`, `fill_rate`) are `None`.
  Rater identity is gone in a counts matrix; honest absence is better
  than fabricated data. Use `n_raters_lower_bound` (always available)
  for what is known.

## RaterSharedMatrix

Entry (i, j) is how many subjects raters i and j both rated. Diagonal
is each rater's total. Provides `thinnest_pair`, `below(threshold)`,
and `render()` for diagnostics.

## Dependencies

- `array_utils` (is_balanced, n_per_subject, n_per_rater, subject_count_distribution)
