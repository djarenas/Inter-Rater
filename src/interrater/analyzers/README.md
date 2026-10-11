# Analyzers

Analyzers are the computational engine of interrater v2. Each analyzer takes
a Ratings object and a Config, validates both, and provides an `estimate()`
method that returns a result object. Analyzers are mutable -- they hold state,
rebuild internals, and can update their own attributes. The data objects they
consume (Ratings, configs) are frozen.


## Hierarchy

The analyzer tree splits at the second level by how each coefficient treats
rater identity:

```
Analyzer (ABC)
  │
  ├── GroupAnalyzer (ABC)
  │   │   Coefficients that discard rater identity.
  │   │   Read a counts matrix (NijMatrix): how many raters put
  │   │   each subject in each category, not which rater said what.
  │   │
  │   ├── WeightableGroupAnalyzer (ABC)
  │   │   │   Adds distance-based weighting for ordinal data.
  │   │   │
  │   │   ├── FleissKappaAnalyzer
  │   │   ├── BrennanPredigerAnalyzer
  │   │   └── KrippendorffAlphaAnalyzer
  │   │
  │   └── (direct GroupAnalyzer children -- no weighting)
  │       ├── GwetAC1Analyzer
  │       └── ScottsAnalyzer
  │
  └── IndividualMarginalAnalyzer (ABC)
      │   Coefficients that keep rater identity.
      │   Read an EstimatorInput: each rater's own category distribution
      │   enters the chance term separately.
      │
      ├── WeightableIndividualAnalyzer (ABC)
      │   │   Adds distance-based weighting.
      │   │
      │   └── CohenKappaAnalyzer
      │
      └── (direct IndividualMarginalAnalyzer children -- no weighting)
          ├── CongerKappaAnalyzer
          └── LightKappaAnalyzer
```

The split is not an implementation choice. Group coefficients treat all raters
as interchangeable -- they pool the raters into one distribution and ask "did
they agree?" Individual marginal coefficients ask "given that rater A uses
these rates and rater B uses those rates, did they agree more than that
difference would predict?" The same data yields different numbers depending
on which question you ask, and the hierarchy makes the question explicit in
the type.


## What each level does

### Analyzer (base)

- Stores Ratings and Config
- Computes BalanceProfile (how subjects were rated -- constant vs. variable
  raters per subject, fill rate, distribution of subject coverage)
- Declares the abstract `estimate()` interface

### GroupAnalyzer

- Validates Config is GroupCoefficientConfig
- Supports two initialization paths: from Ratings, or from NijMatrix directly
- Enforces balance constraints when `config.balanced_data_only=True`

### WeightableGroupAnalyzer / WeightableIndividualAnalyzer

- Resolves the Weights object from `config.weights_type`
- Validates weight matrix dimensions match the number of categories
- Supports "identity" (unweighted), "linear", "quadratic", and "custom"

### IndividualMarginalAnalyzer

- Requires Ratings (never None -- rater identity is essential)
- Builds EstimatorInput from Ratings (category_index + metadata)
- No NijMatrix path: counts discard exactly what these coefficients need

### Leaf analyzers (e.g. FleissKappaAnalyzer)

- Validate their own config type (FleissKappaConfig, etc.)
- Implement `estimate()` with the coefficient-specific bootstrap loop
- Return a coefficient-specific result object (FleissKappaResult, etc.)


## The estimate() flow

Every leaf analyzer's `estimate()` follows the same pattern:

1. **Point estimate.** Compute the coefficient on the original data.
2. **Bootstrap loop.** Resample subjects with replacement, compute the
   coefficient on each resample, collect the distribution.
3. **Confidence interval.** Compute the CI from the bootstrap distribution
   using the method specified in `config.ci_method`.
4. **Return result.** Package everything into a frozen result object.

```python
analyzer = FleissKappaAnalyzer(ratings, config)
result = analyzer.estimate()

result.point               # 0.430
result.confidence_interval # (0.313, 0.525)
result.ci_method           # "percentile"
result.balance_profile     # BalanceProfile(...)
```


## Confidence interval methods

CI computation lives in `ci_methods.py`. Three methods are available:

### Percentile (`ci_method="percentile"`)

The simplest bootstrap CI: take the alpha/2 and 1-alpha/2 quantiles of the
bootstrap distribution directly. No bias correction. Reliable and stable.

### BCa (`ci_method="bca"`)

Bias-corrected and accelerated (Efron 1987). Adjusts the percentile method
for bias and skewness in the bootstrap distribution. More accurate than
percentile, especially for skewed distributions.

The algorithm:

1. Compute bias correction z0 from the ECDF of the bootstrap distribution
   at the point estimate.
2. Compute acceleration *a* from jackknife influence values on the bootstrap
   distribution.
3. Adjust the lower and upper percentiles using z0 and *a*.
4. Return quantiles at the adjusted percentiles.

#### BCa edge cases and fallback

The BCa adjustment can become unstable when the bootstrap distribution has
extreme bias or skewness. Specifically, the adjustment formula:

```
z_adj = z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha))
```

breaks down when the denominator `1 - a * (z0 + z_alpha)` approaches zero.
This produces inf or NaN in the adjusted z-values, which propagates through
`norm.cdf()` into `np.quantile()` and crashes:

```
ValueError: Quantiles must be in the range [0, 1]
```

The root cause is that `np.clip(NaN, 0, 1)` returns NaN -- clip does not
catch NaN values.

A second path to the same crash: if the bootstrap distribution itself
contains NaN values, they propagate through the jackknife into the
acceleration calculation, making *a* = NaN, which makes every downstream
value NaN.

**The fix uses three guard clauses, each falling back to the percentile
method when the BCa correction is unreliable:**

1. **After computing acceleration:** if `acceleration` is not finite (NaN
   from NaN bootstrap values propagating through the jackknife), fall back
   to percentile.

2. **After computing adjusted z-values:** if either `z_lower_adj` or
   `z_upper_adj` is not finite (inf from denominator crossing zero, or NaN
   from other numerical instability), fall back to percentile.

3. **After clipping percentiles:** belt-and-suspenders check. If `p_lower`
   or `p_upper` is NaN despite the earlier guards, fall back to percentile.

The fallback is silent -- the CI is computed, just with the simpler method.
This matches the behavior of R's `boot` package, which falls back to
percentile when the BCa adjustment fails. The CI is still valid; it is
simply less precise about correcting for bias and skewness.

**When does the BCa adjustment actually break in practice?**

Extreme skewness in the bootstrap distribution. This happens when:

- Almost all bootstrap values are on one side of the point estimate (z0 is
  very large positive or negative)
- Combined with non-zero acceleration from asymmetric influence values
- The product `a * (z0 + z)` crosses 1.0, flipping the sign of the
  denominator

On well-behaved data with 1000+ resamples, the BCa adjustment is stable.
The fallback is a safety net for pathological cases.

### Parametric null (`ci_method="parametric_null"`)

Uses the closed-form null variance from Fleiss, Nee & Landis (1979) Eq. 12.
No bootstrap is performed. **Not recommended for new work** -- the null
variance is derived under H0: kappa = 0, and using it to put an interval
around a non-zero kappa produces intervals that are ~55% too narrow.

Provided for reproducing published results that used this formula.


## Subfolders

### `agreement_functions/`

Pure functions that compute each coefficient from its inputs. Called inside
bootstrap loops thousands of times -- they do no validation, no
resampling, no CI computation. They take a counts matrix (or estimator
input) and return a float.

See `agreement_functions/agreement_functions_README.md`.

### `resampling_functions/`

Generators that yield modified arrays: subject bootstrap, subject
jackknife, rater jackknife, rater subset. All return iterators for memory
efficiency -- one resampled array at a time, computed and discarded.

See `resampling_functions/README.md`.

### `permutation_analysis/`

Orchestrators that run a coefficient across subsets of the rater panel.
They are standalone -- NOT in the Analyzer tree. An Analyzer computes one
coefficient on one data slice. An Orchestrator enumerates subsets (pairs,
leave-one-out, arbitrary k) and calls an Analyzer or estimator on each,
then assembles the results.

See `permutation_analysis/README.md`.

### `configs/`

Frozen dataclass configurations for every analyzer. The config hierarchy
mirrors the analyzer hierarchy and encodes mathematical constraints as
type constraints (e.g., a coefficient that doesn't support weighting has
a config that rejects weight specifications at instantiation time).

See `configs/README.md`.


## Other modules

### `ci_methods.py`

Bootstrap CI computation: `percentile_ci()`, `bca_ci()`, and
`scipy_bca_ci()` (reference implementation for validation). See the
confidence interval methods section above.

### `variance_functions.py`

Closed-form variance formulas from the statistical literature. Used for
parametric null CIs and hypothesis tests. Not needed for bootstrap
approaches.

### `convergence_test.py`

Stub for bootstrap convergence diagnostics (not yet implemented). The type
exists so that `ConvergenceTest | None` fields resolve in configs.


## Design principles

**Analyzers are mutable.** They hold state, rebuild internals when needed,
and can update their own attributes. The `add_random_rater()` method
(planned) will modify `self.ratings` and rebuild internal state. This is
deliberate: the analyzer is a computational workspace, not a value object.

**Data objects are frozen.** Ratings and configs are frozen dataclasses. The
analyzer builds new Ratings internally when it needs to change the data --
Ratings itself never mutates.

**Validation flows upward.** Leaf class validates its own constraints,
calls `super().__init__()`, which validates the intermediate class
constraints, which calls its own parent, up to the base. Each level trusts
that its parent has validated what it is responsible for.

**Agreement functions are hot-path code.** They run inside bootstrap loops
(thousands of calls), so they do zero validation. All checking happens
once, before the loop, in the analyzer.

**Resampling is separated from computation.** Resampling functions yield
modified arrays. Agreement functions compute coefficients on arrays. The
analyzer composes them. None knows about the others' internals.
