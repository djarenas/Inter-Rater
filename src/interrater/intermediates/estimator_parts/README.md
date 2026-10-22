# Estimator Parts

The three objects that sit between `Ratings` (the validated input) and the
coefficient estimators. They strip away everything a coefficient does not
need -- names, semantics, validation -- so the hot path runs on bare arrays
and integers.

| Class            | What It Holds                                                      | Why It Exists                                                                 |
|------------------|--------------------------------------------------------------------|-------------------------------------------------------------------------------|
| `Weights`        | K x K float matrix in [0, 1], symmetric, 1 on diagonal             | Says how much credit a disagreement earns. Identity = unweighted.             |
| `NijMatrix`      | (n_subjects, K) int counts matrix                                  | Rater-anonymous tallies. The sole input to group coefficients.                |
| `EstimatorInput` | category_index + n_categories + weights + optional subject_ids     | Bundles what a coefficient needs so K cannot be separated from the array. O(1) validation -- cheap enough for bootstrap loops. |

## Weights

A frozen dataclass. Built from class methods (`identity`, `linear`,
`quadratic`), or from the categories a semantics declares
(`from_categories`), which is how the analyzers get theirs. Remembers its
`kind` so a report can say *how* a coefficient was weighted.

## NijMatrix

The counts matrix every group coefficient reads. Entry (i, j) is how many
raters chose category j for subject i. Rater identity is discarded --
only tallies remain. Converts to and from pandas DataFrames.

## EstimatorInput

Everything a coefficient needs to be computed, in one object. Validates
only shape and consistency, all O(1) -- never inspects array contents.
Provides `take_subjects` (the bootstrap primitive) and `select_raters`
(for leave-one-out and permutation). Computes `nij` on demand rather
than at construction, because pairwise coefficients never want it.

## Dependencies

- `contracts` (EMPTY, MIN_RATERS)
- `array_utils` (category_counts)
