# Permutation Analysis

Permutation analyses run a coefficient across subsets of the rater panel.
They are standalone and NOT in the Analyzer tree. An Analyzer
computes one coefficient on one data slice. A PermutationAnalysis
enumerates subsets (pairs, leave-one-out, arbitrary size) and calls an
estimator on each, then assembles the results.


## Hierarchy

```
PermutationAnalysis                (parent -- shared validation, estimator resolution)
  ├── PairwisePermutationAnalysis  (every C(n,2) pair)
  └── SubsetPermutationAnalysis    (every C(n,k) subset)
```

Scott's pi is pairwise only. It is defined for two raters; on a subset of
three or more it is Fleiss' kappa under another name, so subset
permutation refuses it and points to `fleiss_kappa`.


## Three Computations

### Computation A -- Observed values (always, no bootstrap)

Run the estimator on each pair using the actual data. No resampling.

- `over_pairs(estimator_input, estimator)` --> each pair's value (NaN for
  a pair with no coefficient; see *When a pair or subset has no number*)
- `subjects_per_subset(estimator_input, 2)` --> each pair's subject count,
  checked against the subjects its result was built on
- `mean_over_subsets(values)` --> mean_over_pairs

Cost: C(n,2) estimator calls.

No config imnformation needed beyond the coefficient itself.


### Computation B -- Per-subset bootstrap (always)

For each subset: slice data to the subjects that subset can be measured on,
resample within those subjects, compute the coefficient each draw.

Produces: per-subset distributions, per-subset CIs.

Answers: "how stable is *this subset's* coefficient?"

Cost: C(n,k) x n_resamples estimator calls.

Controlled by: `config_of_coefficient` (n_resamples, alpha, ci_method).
It always runs: the coefficient's config requires n_resamples >= 1.

#### Which subjects a subset can be measured on

**A subject contributes to a subset's agreement only if at least two of that
subset's raters rated it.** Fewer than two ratings is not disagreement, it is
absence. Fleiss makes this explicit -- the per-subject agreement term is

```
(sum_j n_ij^2 - n_i) / (n_i (n_i - 1))
```

which is 0/0 at `n_i = 1`. Two is where the quantity starts existing, not
where we chose to draw a line.

So the pool is **"at least two of the k"**, and that is one rule at every k:

| k | "at least two of the k" means | effect |
|---|-------------------------------|--------|
| 2 | both raters rated the subject | complete case, which is what Cohen's kappa requires |
| 3+ | two or more of the k rated it | unbalanced subsets allowed, no complete rectangle demanded |

The pairwise case is not an exception to the rule, it is the rule at k=2.

**The bootstrap must resample that same pool**, because an interval has to be
around the thing the point estimate measured. Computation A gets there by a
different route -- `over_subsets` hands the whole array to the estimator and
lets it handle the gaps -- but it lands on the same set of contributing
subjects. Draw from the complete cases of a 3-rater subset instead and you
have quietly switched estimands; on partially crossed data the mistake
surfaces as an interval that does not contain its own point estimate.

#### When a pair or subset has no number

**One rule, written once**, in `unit_problem` in
`agreement_functions/rater_combinations.py`. A pair or subset has no
coefficient when:

- fewer than two of its subjects are usable (rated by at least two of its
  raters), or
- one of its raters rated none of the usable subjects. The subset is really
  a smaller one, and its number would be that smaller subset's.

Every place that measures a unit goes through `estimate_unit`, which applies
the rule and returns NaN when it says no: the observed value (A), every
bootstrap draw and jackknife value (B), and every draw of the mean (C). So a
resample that leaves a rater idle fails like any other draw, rather than
reporting the smaller subset's number. The unit's result is `NotComputed`,
with the rule's reason, and the mean leaves it out.

Until 2.1.0 the values ran the estimator on whatever they got, so a subset
with one usable subject, or with an idle rater, could still get a number and
go into the mean.


### Computation C -- Mean bootstrap (if bootstrap_mean_config is set)

Produces: mean distribution, CI on the mean.

Answers: "how stable is the *average* across pairs?"

Two paths:

**Correlated:** Resample subjects globally, compute all pairs on the same
resampled data, average. Preserves the correlation between pairs sharing a
rater -- when one pair drops, pairs sharing a rater tend to drop with it.
This is the Pythagorean intuition: correlated uncertainties add linearly
(sum of the sides), so the CI on the mean is wider than if pairs were
independent.

Cost: n_resamples x C(n,2) estimator calls.

**Independent:** Sample from B's pre-computed per-pair distributions,
average. Treats pairs as if they were independent. Uncertainties add in
quadrature (Pythagorean -- the hypotenuse), so the CI on the mean is
narrower. On real data this is typically too narrow -- pairs sharing a
rater ARE correlated, and ignoring it underestimates the variance of the
mean. The v1 package's analytic SE made this assumption and was ~28% too
narrow.

Cost: 0 estimator calls. It is built during Computation B: as each pair's
draws are made, one value per mean draw is picked from them and added to a
running total (see *Keeping the draws*).

Controlled by: `bootstrap_mean_config` (method only). The draw count,
interval method, level and seed are the coefficient's own
(`config_of_coefficient`), so the mean and the per-pair intervals are built
the same way. The independent mean needs `ci_method: percentile` there.

#### Keeping the draws

Every pair or subset is bootstrapped, and its draws are the evidence behind
its interval: they can be plotted, re-analysed, or used to check how much the
interval's bounds would move under another seed. So by default the result
keeps them, for every pair or subset.

Until 2.1.0 they were kept only up to 20 pairs or 120 subsets. Above that
they were computed, held through the whole run, and thrown away at the end.
That freed memory after the run, never during it, and it meant a 10-rater
panel could not show a single pair's distribution.

`keep_draws: false` is for panels where the draws will not fit. Each pair or
subset is bootstrapped, its interval built, its share of the mean added, and
then its draws are tossed. Memory holds one unit's draws at a time, whatever
the number of units. Every value, interval and mean interval is the same as
with `keep_draws: true`; only the draws are missing: each unit's
`bootstrap_distribution` is `NotComputed("draws not kept: keep_draws
false")`.

**The same seed gives the same numbers either way.** The independent mean
samples from each unit's draws, so it is built as the units are bootstrapped:
as each unit's draws are made, one value per mean draw is picked from them and
added to a running total, before the draws can be tossed. This is done the
same way under both settings. If it were not, `keep_draws` would change the
independent interval, and a setting about storage must not change a result.
The correlated mean never reads the per-unit draws: it resamples subjects
itself, one draw at a time, from its own seed stream, so `keep_draws` cannot
change it either.

The mean's own draws, one array, are always kept.

#### Every draw averages the same pairs

The point estimate averages the pairs defined on the real data. Each draw
of the mean -- correlated or independent -- averages exactly those pairs,
and a draw in which any of them is undefined is NaN, counted as failed, as
for a single coefficient. Averaging the others instead would make that draw
a different quantity. The correlated jackknife follows the same rule. See
`mean_over_same_subsets` in `agreement_functions`.

#### BCa is available for the correlated mean only

BCa's acceleration is the standardised skewness of leave-one-out estimates,
and **the unit left out must be the unit that was resampled**.

- **Per-subset intervals** resample the subjects in that subset's measurable
  pool, so the jackknife leaves one of those out. Drop a subject the subset
  cannot be measured on and the coefficient does not move, so that subject
  enters the jackknife as a zero deviation. Pad the jackknife with zeros and
  the standardised skewness shrinks toward zero -- the acceleration correction
  quietly switches itself off, leaving the bias correction to do the work of
  both while the interval still calls itself BCa.
- **The correlated mean** resamples subjects globally, so its jackknife is
  global too: drop one subject, recompute every subset, average.
- **The independent mean** resamples from the per-subset distributions and
  never touches a subject. There is nothing to leave out, so there is no
  acceleration to compute. Requesting `ci_method="bca"` for it raises rather
  than approximating one.


## Each pair or subset is a result

After the three computations, each unit becomes the result its coefficient's
own analyzer gives on that unit's usable subjects: the analyzer is built on
those subjects, with injection off (the panel's injection already happened),
and handed the point, interval and draws this analysis computed, through the
same `_result` its `estimate()` uses. The analysis checks that the analyzer
computes the same point on those subjects; a different number is a bug.

**The point matches a direct run; the interval does not, on purpose.** Run
the analyzer alone on a pair's subjects and the point, Po and Pe, kappa_max
and profiles are identical. The interval and draws differ, because each unit
draws from its own child seed, spawned from `random_state`
(`np.random.SeedSequence(random_state).spawn(n_units)`), while a direct run
draws from `random_state` itself. Separate streams keep the units' draws
independent of each other.

See `intermediates/results/permutation_results/README.md` for what the
result holds and checks.


## Why Computation B is separate from Computation C

The correlated mean loop (Computation C) computes kappa for every pair in
every draw, so per-pair distributions fall out as a side product. It is
tempting to skip Computation B and use those instead.

The two resampling schemes produce different per-pair distributions on
sparse data:

- **Per-pair resampling (B):** Pair (0,1) shares 756 subjects. Resample
  756 from those 756 with replacement. Always exactly 756 subjects per
  draw. Asks: "how would this pair's kappa vary if we re-drew from the
  pair's own subjects?"

- **Global resampling (C side product):** Resample from all 1,617
  subjects. Pair (0,1) gets however many of the resampled subjects both
  raters happened to rate -- maybe 740, maybe 770. The effective sample
  size varies across draws. Asks: "how would this pair's kappa vary if we
  re-drew the whole study?"

For balanced data (every subject rated by every rater) the two are
identical. For sparse data the difference is in the finite-sample variance
of the per-pair distribution. Per-pair resampling is the more standard
choice for per-pair CIs.

The global resampling side product is valid for the mean CI -- that is its
job. Using it for per-pair CIs is a convenience, not a free lunch.


## Total Cost

| Configuration              | Estimator calls                            |
|----------------------------|--------------------------------------------|
| A + B                      | C(n,2) x (1 + n_resamples)                |
| A + B + C independent      | C(n,2) x (1 + n_resamples)                |
| A + B + C correlated       | C(n,2) x (1 + 2 x n_resamples)            |

C(n,2) by number of raters:

| n_raters | C(n,2) |
|----------|--------|
| 2        | 1      |
| 3        | 3      |
| 4        | 6      |
| 5        | 10     |
| 6        | 15     |
| 10       | 45     |
| 20       | 190    |

With n_resamples = 4,000:

| n_raters | A + B          | A + B + C correlated |
|----------|----------------|----------------------|
| 2        | 4,001          | 8,001                |
| 4        | 24,006         | 48,006               |
| 6        | 60,015         | 120,015              |
| 10       | 180,045        | 360,045              |
| 20       | 760,190        | 1,520,190            |

### The cost guard

Before anything runs, the analysis counts these calls and refuses a run above
`max_evaluations` in the config. The default is `GUARD_MAX_EVALUATIONS` in
`contracts`, 2,000,000: at the measured 50 to 260 microseconds per call, about
2 to 9 minutes. Raise it to run a larger analysis on purpose.

With BCa, each interval also runs a jackknife; the guard counts it as one call
per subject in the panel. That is an upper bound -- a pair or subset is
jackknifed only over its own subjects -- so the guard can refuse a little
early, never late. The message says "up to" for that reason.

Subset permutation checks every k before the first one runs: each against
`GUARD_MAX_SUBSETS` (10,000 subsets), and their total against
`max_evaluations`. Without that, `subset_k: "all"` on 20 raters ran k = 2, 3
and 4 in full before k = 5 hit the cap.
