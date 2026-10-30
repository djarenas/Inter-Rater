# Permutation Analysis

Permutation analyses run a coefficient across subsets of the rater panel.
They are standalone -- they are NOT in the Analyzer tree. An Analyzer
computes one coefficient on one data slice. A PermutationAnalysis
enumerates subsets (pairs, leave-one-out, arbitrary size) and calls an
estimator on each, then assembles the results.


## Hierarchy

```
PermutationAnalysis                (parent -- shared validation, estimator resolution)
  ├── PairwisePermutationAnalysis  (every C(n,2) pair)
  └── SubsetPermutationAnalysis    (every C(n,k) subset)
```


## Three Computations

### Computation A -- Observed values (always, no bootstrap)

Run the estimator on each pair using the actual data. No resampling.

- `over_pairs(estimator_input, estimator)` --> pairwise_values
- `subjects_per_subset(estimator_input, 2)` --> subjects_per_pair
- `mean_over_subsets(pairwise_values)` --> mean_over_pairs

Cost: C(n,2) estimator calls.

No config needed beyond the coefficient itself.


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

**Subsets with fewer than two measurable subjects** have nothing to resample.
They get a NaN distribution and an interval whose `method_used` is `None` --
no method ran, which is a different statement from an interval that could not
be computed.


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

Cost: 0 estimator calls (numpy sampling from B's output). Requires B to
have run first.

Controlled by: `bootstrap_mean_config` (method only). The draw count,
interval method, level and seed are the coefficient's own
(`config_of_coefficient`), so the mean and the per-pair intervals are built
the same way. The independent mean needs `ci_method: percentile` there.

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
