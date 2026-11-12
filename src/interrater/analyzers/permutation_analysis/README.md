# Permutation Analysis

Permutation analyses run a coefficient across subsets of the rater panel.
They are standalone, not in the Analyzer tree. An Analyzer computes one
coefficient on one data slice. A PermutationAnalysis enumerates subsets
(pairs, leave-one-out, any size), calls an estimator on each, and assembles
the results.


## Hierarchy

```
PermutationAnalysis                (parent: shared validation, estimator resolution)
  ├── PairwisePermutationAnalysis  (every C(n,2) pair)
  └── SubsetPermutationAnalysis    (every C(n,k) subset)
```

Scott's pi is pairwise only. It is defined for two raters; on three or more
it is Fleiss' kappa, so subset permutation refuses it and points to
`fleiss_kappa`.

Weights come from the semantics, as in the analyzers. Every pair or subset
is computed with them, and the result records them (`weights_type`,
`weight_matrix`).


## Three computations

### A: observed values (always, no bootstrap)

The estimator on each pair, on the actual data:

- `over_pairs(estimator_input, estimator)`: each pair's value (NaN for a
  pair with no coefficient; see *When a pair or subset has no number*)
- `subjects_per_subset(estimator_input, 2)`: each pair's subject count,
  checked against the subjects its result was built on
- `mean_over_subsets(values)`: `mean_over_pairs`

Cost: C(n,2) estimator calls. Needs nothing from the config beyond the
coefficient.

### B: per-subset bootstrap (always)

For each subset: keep the subjects it can be measured on, resample within
them, compute the coefficient each draw. Gives each subset's draws and
interval: how stable is this subset's coefficient?

Cost: C(n,k) x n_resamples estimator calls. Set by `config_of_coefficient`
(n_resamples, alpha, ci_method); it always runs, since n_resamples >= 1.

#### Which subjects a subset can be measured on

**A subject counts toward a subset's agreement only if at least two of the
subset's raters rated it.** One rating is absence, not disagreement: Fleiss'
per-subject term

```
(sum_j n_ij^2 - n_i) / (n_i (n_i - 1))
```

is 0/0 at `n_i = 1`. So the pool is "at least two of the k", one rule at
every k:

| k | "at least two of the k" means | effect |
|---|-------------------------------|--------|
| 2 | both raters rated the subject | complete case, as Cohen's kappa requires |
| 3+ | two or more of the k rated it | unbalanced subsets allowed |

**The bootstrap resamples that same pool**, so the interval is around what
the point estimate measured. Computation A reaches the same subjects another
way (`over_subsets` passes the whole array and the estimator skips the gaps).
Resampling the complete cases of a 3-rater subset instead changes the
estimand; on partially crossed data it shows as an interval that does not
contain its own point.

#### When a pair or subset has no number

One rule, written once, in `unit_problem` in
`agreement_functions/rater_combinations.py`. A pair or subset has no
coefficient when:

- fewer than two of its subjects are usable (rated by at least two of its
  raters), or
- one of its raters rated none of the usable subjects. The subset is really
  a smaller one.

Every place that measures a unit goes through `estimate_unit`, which returns
NaN when the rule says no: the observed value (A), every bootstrap draw and
jackknife value (B), and every draw of the mean (C). A resample that leaves a
rater idle fails like any other draw. The unit's result is `NotComputed` with
the reason, and the mean leaves it out.

Before 2.1.0 the estimator ran on whatever it got, so a subset with one
usable subject, or an idle rater, could get a number and enter the mean.

### C: bootstrap of the mean (if `bootstrap_mean_config` is set)

Gives the mean's draws and interval: how stable is the average across pairs?

**Correlated:** resample subjects once per draw, compute every pair on the
same resample, average. Pairs sharing a rater move together, and this keeps
that correlation, so the interval is wider than if the pairs were
independent. Cost: n_resamples x C(n,2) estimator calls.

**Independent:** pick one value from each pair's draws (from B) and average.
This treats the pairs as independent, so the interval is narrower, usually
too narrow on real data. Version 1's analytic standard error made the same
assumption and was roughly 30% too narrow. Cost: no estimator calls; it is
built during B (see *Keeping the draws*).

Set by `bootstrap_mean_config` (method only). The draw count, interval
method, level and seed come from `config_of_coefficient`, so the mean and the
per-pair intervals are built the same way. The independent mean needs
`ci_method: percentile`.

#### Keeping the draws

Each unit's draws are kept by default, for every pair or subset: they can be
plotted, re-analysed, or used to see how far the bounds move with another
seed. (Before 2.1.0 they were kept only up to 20 pairs or 120 subsets, after
being held through the whole run.)

`keep_draws: false` is for panels whose draws will not fit. Each unit is
bootstrapped, its interval built, its share of the mean added, and its draws
discarded, so memory holds one unit's draws at a time. Every value and
interval is the same; each unit's `bootstrap_distribution` is
`NotComputed("draws not kept: keep_draws false")`.

The same seed gives the same numbers under either setting. The independent
mean picks its values as each unit is bootstrapped, before the draws can be
discarded, in the same way under both settings. The correlated mean
resamples subjects from its own seed stream and never reads the unit draws.
The mean's own draws are always kept.

#### Every draw averages the same pairs

The point averages the pairs defined on the real data. Each draw of the mean,
correlated or independent, averages exactly those pairs; a draw in which any
is undefined is NaN and counted as failed, as for a single coefficient. The
correlated jackknife follows the same rule (`mean_over_same_subsets` in
`agreement_functions`).

#### BCa is available for the correlated mean only

BCa's acceleration comes from leave-one-out estimates, and the unit left out
must be the unit that was resampled.

- **Per-subset intervals** resample the subjects in the subset's pool, so
  the jackknife leaves out one of those. Leaving out a subject outside the
  pool would add a zero deviation, shrink the acceleration toward zero, and
  still report BCa.
- **The correlated mean** resamples subjects globally, so its jackknife is
  global: drop one subject, recompute every subset, average.
- **The independent mean** never resamples subjects, so there is nothing to
  leave out and no acceleration. `ci_method="bca"` for it is refused.


## Each pair or subset is a result

Each unit becomes the result its coefficient's own analyzer gives on the
unit's usable subjects: the analyzer is built on those subjects with
injection off (the panel's injection already happened), and handed the
point, interval and draws computed here, through the same `_result` its
`estimate()` uses. The analysis checks that the analyzer computes the same
point; a different number is a bug.

The point matches a direct run; the interval does not. Run the analyzer
alone on a pair's subjects and the point, Po and Pe, kappa_max and profiles
are identical. The interval and draws differ because each unit draws from its
own child seed (`np.random.SeedSequence(random_state).spawn(n_units)`), so
the units' draws are independent of each other.

See `intermediates/results/permutation_results/README.md` for what the
result holds and checks.


## Why B is separate from C

The correlated mean computes every pair in every draw, so per-pair
distributions come out of it too. They are not the same as B's on sparse
data:

- **B (per pair):** pair (0,1) shares 170 subjects; each draw resamples 170
  from those 170. How would this pair's kappa vary if we re-drew its own
  subjects?
- **C (global):** each draw resamples all 1,617 subjects, and pair (0,1)
  gets however many of them both raters rated, maybe 155, maybe 185. How
  would it vary if we re-drew the whole study?

On balanced data the two are the same. On sparse data their variances
differ, and per-pair resampling is the standard choice for per-pair
intervals. The global draws are for the mean.


## Total cost

| Configuration              | Estimator calls                            |
|----------------------------|--------------------------------------------|
| A + B                      | C(n,2) x (1 + n_resamples)                |
| A + B + C independent      | C(n,2) x (1 + n_resamples)                |
| A + B + C correlated       | C(n,2) x (1 + 2 x n_resamples)            |

| n_raters | C(n,2) | A + B (n_resamples = 4,000) | A + B + C correlated |
|----------|--------|-----------------------------|----------------------|
| 2        | 1      | 4,001                       | 8,001                |
| 4        | 6      | 24,006                      | 48,006               |
| 6        | 15     | 60,015                      | 120,015              |
| 10       | 45     | 180,045                     | 360,045              |
| 20       | 190    | 760,190                     | 1,520,190            |

### The cost guard

Before anything runs, the analysis counts these calls and refuses a run above
`max_evaluations` in the config. The default is `GUARD_MAX_EVALUATIONS` in
`contracts`, 2,000,000: at the measured 50 to 260 microseconds per call,
about 2 to 9 minutes. Raise it to run a larger analysis on purpose.

With BCa, each interval also runs a jackknife, counted as one call per
subject in the panel. That is an upper bound (a unit is jackknifed only over
its own subjects), so the guard can refuse a little early, never late, and
its message says "up to".

Subset permutation checks every k before the first one runs: each against
`GUARD_MAX_SUBSETS` (10,000 subsets), and their total against
`max_evaluations`.
