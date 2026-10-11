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


### Computation B -- Per-pair bootstrap (if config_of_coefficient requests bootstrap)

For each pair: slice data to that pair's shared subjects, resample within
those subjects, compute kappa each draw.

Produces: per-pair distributions, per-pair CIs.

Answers: "how stable is *this pair's* kappa?"

Resampling: per-pair. Each pair resamples from its own shared subjects.
The effective sample size is fixed across draws (always N_ij subjects).

Cost: C(n,2) x n_resamples estimator calls.

Controlled by: `config_of_coefficient` (n_resamples, alpha, ci_method).


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

Cost: n_bootstraps x C(n,2) estimator calls.

**Independent:** Sample from B's pre-computed per-pair distributions,
average. Treats pairs as if they were independent. Uncertainties add in
quadrature (Pythagorean -- the hypotenuse), so the CI on the mean is
narrower. On real data this is typically too narrow -- pairs sharing a
rater ARE correlated, and ignoring it underestimates the variance of the
mean. The v1 package's analytic SE made this assumption and was ~28% too
narrow.

Cost: 0 estimator calls (numpy sampling from B's output). Requires B to
have run first.

Controlled by: `bootstrap_mean_config` (method, n_bootstraps, alpha).


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
| A only                     | C(n,2)                                     |
| A + B                      | C(n,2) x (1 + n_resamples)                |
| A + B + C independent      | C(n,2) x (1 + n_resamples)                |
| A + B + C correlated       | C(n,2) x (1 + n_resamples + n_bootstraps) |

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

With n_resamples = 4,000 and n_bootstraps = 4,000:

| n_raters | A + B          | A + B + C correlated |
|----------|----------------|----------------------|
| 2        | 4,001          | 8,001                |
| 4        | 24,006         | 48,006               |
| 6        | 60,015         | 120,015              |
| 10       | 180,045        | 360,045              |
| 20       | 760,190        | 1,520,190            |
