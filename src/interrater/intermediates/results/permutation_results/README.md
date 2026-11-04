# Permutation Results

One coefficient, recomputed over subsets of the raters. Where the analyzer
results answer *what is the coefficient*, these answer *how much does it move
when the rater set changes* -- across all pairs, across leave-k-out subsets at
one k, or across several k at once.

```
PermutationResult               coefficient_name, ci_method, balance_profile
├── PairwisePermutationResult   all C(n,2) pairs
├── SubsetPermutationResult     leave-k-out at a single k
└── MultiKSubsetResult          several k values at once
```

| Class                       | What It Holds                                                                 | Why It Exists                                                        |
|-----------------------------|--------------------------------------------------------------------------------|----------------------------------------------------------------------|
| `PermutationResult`         | coefficient_name, ci_method, balance_profile                                   | Shared base. Deliberately carries no point estimate -- see below.    |
| `PairwisePermutationResult` | pair_results (each pair's own result, or NotComputed), n_pairs, mean_over_pairs; optionally bootstrap_mean_method, bootstrap_ci_of_mean, bootstrap_distribution_of_mean | All C(n,2) pairwise coefficients, each a full result.                |
| `SubsetPermutationResult`   | subset_results (each subset's own result, or NotComputed), n_subsets, subset_k, mean_over_subsets; the same three optional fields for the mean | Leave-one-out or leave-k-out at a single k.                          |
| `MultiKSubsetResult`        | per_k_results (dict[int, SubsetPermutationResult]), ks                          | Stability curve: how the coefficient moves as more raters are removed. |

## No single point estimate

`PermutationResult` has no `point` field, since a permutation analysis produces a *collection*
of coefficients (i.e. one per pair, one per subset). The mean over them is
reported as its own field (`mean_over_pairs`, `mean_over_subsets`) rather than
promoted to "the" estimate.

That distinction is important because **the values in these collections are not
independent**. Two pairs sharing a rater are positively correlated, so the
variance of their mean is not the variance of an average of independent
things. I was guilty of doing this wrongly in version 1; the analytic standard error produced intervals roughly 30% too narrow on real data.

## Each pair or subset is a result

`pair_results` and `subset_results` map each unit's raters, e.g. `(0, 2)`,
to the result its coefficient's own analyzer gives on that unit's subjects:
a `CohenKappaResult` for a pair under Cohen's kappa, a `FleissKappaResult`
for a subset under Fleiss'. This helps each pair carry everything a single analysis
would have produced: point, interval, draws, Po and Pe, kappa_max, its own balance and
agreement profiles, ...

**The container checks its units** when it is built, and raises if any
fails:

- each unit is the coefficient's own result class, or `NotComputed`;
- each was built with the container's `ci_method`;
- each computed unit has a bootstrap interval -- the per-unit bootstrap
  always runs, so one without is a bug;
- the draws are kept for every computed unit or for none -- `keep_draws`
  applies to the whole analysis;
- the mean is the mean of the units' defined points.

**A unit with no coefficient holds `NotComputed`**, which says why: "raters
ann, bob share 1 subject", or "rater cat rated none of the 40 subjects the
others share". It is not a result full of NaNs. The rule is in
`analyzers/permutation_analysis/README.md`.


## Every interval here is a `BootstrapCIResult`

`bootstrap_ci_of_mean` holds a `BootstrapCIResult`, and so does every
computed unit's `confidence_interval`. A tuple can only
say where the bounds are, and these intervals routinely need to say more than
that.


## The mean's bootstrap is optional

Each pair or subset is always bootstrapped; its interval is in its result.
The bootstrap of the mean is requested, not assumed: `bootstrap_mean_method`,
`bootstrap_ci_of_mean` and `bootstrap_distribution_of_mean` are `None`
without it. The spread of the units is itself a diagnostic, and absence is
recorded as `None` rather than as an empty container.

The correlated mean resamples subjects and recomputes every pair or subset
within each draw, which preserves the correlation above without having to
model it. `bootstrap_mean_method` records which scheme produced the
distribution.

## Subjects per pair or subset

Each unit's result carries the number of subjects it was computed from, in
its `balance_profile.n_subjects`, and the printed table lists it.
Under a partially crossed design two raters may overlap on very few subjects,
and a coefficient from a thin overlap should not be read beside one from a
thick overlap without that context travelling with it.

To see every pair's overlap at once, `RaterSharedMatrix` (in
`intermediates/data_profiles`) can be used: it counts the subjects each pair
of raters shares, and prints them as a grid, or as a list of the thinnest
pairs first.

**What counts as a subject here is a rule, not an obvious fact.** A subject
contributes to a k-rater subset when **at least two of those k raters rated
it** -- not when all k did. Fewer than two ratings is not disagreement, it is
absence: Fleiss' per-subject term `(sum_j n_ij^2 - n_i) / (n_i (n_i - 1))` is
0/0 at `n_i = 1`. At k = 2 the rule reduces to both raters, which is the
complete-case pairing Cohen's kappa requires, so it is one rule at every k
rather than a special case for pairs.

That pool is what the bootstrap draws from, and the reason is estimand
matching: an interval has to be around the thing the point estimate measured.
Resampling the complete cases of a 3-rater subset instead silently switches
estimands, and on partially crossed data the symptom is an interval that does
not contain its own point estimate. On this package's SR1 main dataset the two
pools differ enormously -- 756 measurable subjects against 0 complete cases
for one 3-rater subset -- so the mistake is not a rounding difference.

A subset with fewer than two measurable subjects has no coefficient, and
holds `NotComputed`.

## Getting at the pieces

`MultiKSubsetResult` is a container: `result[3]` returns the
`SubsetPermutationResult` for k = 3, `len(result)` counts the k values, and
iterating yields them in sorted order. The `ks` field lists them.

Each unit's `bootstrap_distribution` holds its raw draws, for callers who
want to plot or re-analyse them rather than read the summarised interval.
With `keep_draws: false` it is `NotComputed("draws not kept: keep_draws
false")` (see `analyzers/permutation_analysis/README.md`).

`to_frame()` gives one row per pair or subset -- raters, estimate,
ci_lower, ci_upper, point_outside_interval, method_used, n_failed,
subjects, not_computed -- to save with `to_csv`. Empty cells are NaN, so it
reads back from CSV unchanged. The printed table is made from the same rows.

A unit's point outside its own interval raises no warning; with many pairs
there would be one per pair. Its row says `POINT OUTSIDE INTERVAL`, and the
frame's column says True.

`summary()` on each class returns plain types. `pair_results` and
`subset_results` hold each unit's own `summary()`, or
`{"not_computed": reason}`. When the mean was bootstrapped,
`bootstrap_draws_mean` and `bootstrap_se` describe its usable draws, as in
`PointEstimateResult`. The mean's draws must agree with the mean's interval;
each unit's draws are checked against its interval by its own result.

## Dependencies

- `..analysis_result` (AnalysisResult)
- `..analyzer_results` (each coefficient's result class, for the units)
- `..ci_results` (BootstrapCIResult)
- `...not_computed` (NotComputed)
- `...data_profiles` (BalanceProfile)
- `numpy` (bootstrap distributions)
