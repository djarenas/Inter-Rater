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
(`build_from_categories`), which is how the analyzers get theirs. The object remembers its
`kind` -> a report can say *how* a coefficient was weighted.

## NijMatrix

The counts matrix every group coefficient reads. Entry (i, j) is how many
raters chose category j for subject i. Rater identity is discarded, and 
only tallies remain. Converts to and from pandas DataFrames.

## EstimatorInput

Everything a coefficient needs to be computed, in one object. Validates
only shape and consistency, all O(1) operations. It never inspects array contents.

Three methods make a new `EstimatorInput` from an existing one:

- `take_subjects(rows)`: the given rows, repeats allowed. This is the
  bootstrap primitive: a draw is a copy, so it cannot disturb the original.
- `select_raters(keep)`: the given rater columns. Used for the pairs and
  subsets of a permutation analysis, and for Light's pairs. The result may be
  unbalanced or hold subjects with too few ratings; the caller checks that.
- `nij`: the counts matrix, computed when asked rather than at construction,
  because pairwise coefficients never want it. The conversion runs one way:
  counts have no rater columns, so a rater must be selected before it.

## Dependencies

- `contracts` (MIN_RATERS)
- `array_utils` (category_counts)
- `semantics`, imported inside `Weights.build_from_categories`
- `pandas`, imported inside the `NijMatrix` DataFrame conversions
