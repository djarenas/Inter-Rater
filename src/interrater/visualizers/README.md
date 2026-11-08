# Visualizers

Plots of results and of the data behind them. 

Every plot draws on an `ax` you pass, or makes a new figure, and returns the
`Figure`. None calls `plt.show()` or saves a file. matplotlib is imported
only when a plot is first used.

## What each plot shows

| Plot | Takes | Shows |
|------|-------|-------|
| `plot_bootstrap_distribution(result)` | a coefficient's result from `estimate()` | The usable bootstrap draws, the point, and the interval labelled with the method that produced it. For BCa, the plain percentile interval as an outline and the two corrections (z0, a) in words. |
| `plot_null_distribution(result)` | a permutation `hypothesis_test()` result | The null distribution, the observed value, the shuffles counted toward p, and the p-value with its arithmetic. |
| `plot_pairwise_permutation(result)` | a `PairwisePermutationResult` | Every pair under both of its raters, with per-pair intervals and the mean with its interval. |
| `plot_subset_sorted(result)` | one `SubsetPermutationResult` (`result[k]` from a multi-k run) | Every subset of one size, sorted by value, with intervals and the mean. |
| `plot_subset_strip(result)` | a `MultiKSubsetResult` | Every subset at every size, and how the mean moves with panel size. |
| `RaterProfilePlotter` | `Ratings`, or any result built from ratings | `draw_overlap()`: subjects shared by each pair of raters. `draw_counts()`: subjects each rater rated. |
| `BalanceProfilePlotter` | `Ratings`, an `NijMatrix`, or any result | `draw()`: how many subjects got 1, 2, 3, ... raters. |

In the three permutation plots, a pair, subset or mean whose bootstrap ran
and produced no interval (every draw failed) is drawn hollow, with a legend
entry saying so, rather than without a whisker like one that was never
bootstrapped.

## The profile plotters

The two profile plotters are classes, one per profile. Each holds its profile
and is built by the constructor that fits what you have:

```python
from interrater import RaterProfilePlotter, BalanceProfilePlotter

raters = RaterProfilePlotter.build_from_ratings(ratings)
raters.draw_overlap()      # heatmap of shared subjects
raters.draw_counts()       # bar per rater

BalanceProfilePlotter.build_from_result(result).draw()
```

- **`RaterProfilePlotter`** needs rater identity. `build_from_result` on a
  result built from counts refuses when built, not when drawn.
  - `draw_overlap()`: the full matrix in declared rater order. The diagonal
    is each rater's own total, white and outside the colour scale, so the
    colours compare pairs only. A pair sharing fewer than 2 subjects is
    hatched: the analysis gives it no coefficient. `annotate=False` drops
    the cell numbers on large panels.
  - `draw_counts()`: one bar per rater, in the heatmap's order, with a dashed
    line at all subjects.
- **`BalanceProfilePlotter`** works from counts too, since a balance profile
  is subject-side only. `draw()` has one bar for every count of raters from
  the fewest to the most, so a count no subject had shows as a zero bar, and
  each bar's count is written on it.

The overlap heatmap is the context a pairwise or subset table needs on
partially crossed data: a pair resting on 8 shared subjects and one resting
on 400 print the same way.

## Functions now, classes later

The five result plots are functions, as in 2.0. The two profile plotters are
classes because they take several kinds of input (ratings, counts, results)
through named constructors. Converting the five to the same pattern is on
the 2.2 list.
