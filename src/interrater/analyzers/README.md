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
result.confidence_interval # BootstrapCIResult, not a tuple of two floats
result.ci_method           # "percentile"
result.balance_profile     # BalanceProfile(...)

print(result.confidence_interval)
# [0.313, 0.525] (95% CI, percentile, n=1000)
```


## Confidence interval methods

CI computation lives in `ci_methods.py`. Both bootstrap methods return a
`BootstrapCIResult` rather than a pair of floats, because an interval that
reports only its bounds cannot say which method produced them or how many
resamples it had to work with.

### Percentile (`ci_method="percentile"`)

The simplest bootstrap CI: take the alpha/2 and 1-alpha/2 quantiles of the
bootstrap distribution directly. No bias correction. Reliable and stable.

It is also the floor of the module. BCa can decline and hand its work here;
this has nothing simpler to hand off to in turn, so percentile either runs
or nothing does -- it never reports a fallback.

### BCa (`ci_method="bca"`)

Bias-corrected and accelerated (Efron & Tibshirani 1986, eq. 7.15). Two
corrections sit on top of the percentile method:

- **Bias correction `z0`** asks whether the bootstrap distribution is centred
  on the point estimate, and shifts the interval when it is not. It is read
  from the ECDF of the bootstrap distribution at the point estimate.

- **Acceleration `a`** asks whether the coefficient's standard error changes
  as the coefficient itself changes. Agreement near 1.0 varies less between
  samples than agreement near 0.5, simply because there is less room to move;
  the percentile method assumes that away, and `a` corrects for it.

#### Where acceleration comes from

**`a` is the standardised skewness of the leave-one-subject-out estimates,
not of the bootstrap draws.** Efron & Tibshirani (1986) eq. (8.5):

```
a = sum(d_i^3) / (6 * [sum(d_i^2)]^(3/2)),   d_i = mean(theta_(.)) - theta_(i)
```

where `theta_(i)` is the coefficient recomputed with subject *i* removed.

The distinction is the whole point of the quantity. A bootstrap distribution
describes the spread of the coefficient at the one value the data produced.
Acceleration asks how that spread *changes* as the value changes, which no
single distribution can answer. Measuring it means disturbing the data and
watching the estimate move: drop one subject, recompute, repeat.

So `bca_ci` takes `jackknife_estimates` as a **required** argument. Asking
for BCa without supplying them is a mistake in the calling analyzer rather
than a property of the data, and it raises rather than quietly producing
something else. `resampling_functions.subject_jackknife` yields the datasets
those estimates are computed from.

*(Versions before 2.0.1 computed the skewness of the bootstrap draws instead.
That is a different quantity and it is not the acceleration.)*

#### Checked against the paper

Efron & Tibshirani (1986) Example 3, p. 71, works the law school data: 15
schools, correlation between LSAT score and GPA, `theta_hat = .776`.

| | `a` | `z0` | central 90% interval |
|---|---|---|---|
| published | -.0817 | -.0927 | [.43, .92] |
| this package | -.0757 | -.0960 | [0.43, 0.93] |

The paper prints its fifteen empirical influence values `U_i^0`. Feeding
those directly into our implementation of eq. (8.5) returns **-0.0817**, the
published constant, to four decimals. So the formula is right, and the
remaining difference in `a` is a difference in its *input*: Efron computes
`U_i^0` from the infinitesimal jackknife of eq. (10.11), while this package
uses the delete-one jackknife, which is the form the 1993 book uses. Both
are standard and they do not agree exactly.

Two things to watch if you reproduce this. The published interval is
**90%**, so `alpha=0.10`; comparing a 95% interval against it will look like
a failure. And the paper used B = 100,000 draws, "about 100 times more than
actually necessary" in its own words, so a smaller B moves the upper bound
by about .01.

#### The noise gate on acceleration

When every leave-one-out estimate is identical, no subject moves the
coefficient, there is no skewness to correct, and `a` is undefined -- so
`_acceleration` returns `None` and BCa steps down to percentile.

Detecting that case is less obvious than it looks. The natural test is
whether the deviations from the mean are zero, and they are not: the mean of
identical floats is not reliably one of them. Summing twenty copies of 0.6
and dividing by twenty lands an ulp away, leaving deviations around 1e-17 and
a sum of squares around 1e-31 instead of the exact zero you would expect.

That is not a harmless rounding error. **`a` is a standardised skewness**, so
dividing noise by noise of the same order returns a value in the ordinary
range for a real acceleration -- a coefficient that no single subject moves
would get a fabricated correction, and the interval would report
`method_used='bca'` as though the correction meant something. Whether it
happened at all depended on the bit pattern of the value: 0.5 and 0.75 are
exact in binary and came out clean, 0.6 did not.

The gate is therefore on the **range of the values**, not on their deviations
from the mean. Identical values have a range of exactly zero whatever the
value is, because the noise enters through the mean rather than the values.
The tolerance is a few ulps of the values' own magnitude; real leave-one-out
spread is on the order of the coefficient divided by the number of subjects,
around 1e-3 against a floor near 1e-15, so nothing legitimate approaches it.

#### When BCa declines

Two corrections, two ways to be unavailable:

| situation | what happens |
|-----------|--------------|
| acceleration undefined (identical or too few jackknife values) | fall back to percentile |
| `1 - a * (z0 + z)` crosses zero, sending an adjusted quantile to infinity | fall back to percentile |
| no resample usable at all | nothing runs; `method_used` is `None` |

**The fallback is not silent, and it is not partial.** `method_used` records
what actually ran, so `'percentile'` after a `'bca'` request is visible to
the caller. And it is a real fallback: `z0` is discarded along with `a`, so
an interval reporting `percentile` is literally a pair of quantiles rather
than a relabelled half-corrected interval. Table 6 of Efron & Tibshirani
(1986) makes stepping down the right answer -- each method there is correct
under weaker assumptions than the one above it, so the fallback costs
accuracy, not validity.

The third row is a different statement from the first two. `method_used` of
`None` is an absence, not a substitution: no method was applied, so there are
no bounds, and `degraded` is `False` because nothing was swapped for anything.

#### Failed resamples

A resample fails when its coefficient is undefined. For any coefficient of
the form `(p_o - p_e) / (1 - p_e)`, a resample that happens to contain only
one category drives chance agreement to 1 and the expression has nothing to
say; the agreement functions return NaN, and the draw arrives here as a gap.

Both methods compute their bounds from the usable draws and report
`n_failed`. **Neither thresholds on that count.** How much the failures
matter depends on whether failing correlates with the coefficient's value,
which no count can tell you, so the count is reported and the judgement is
left to whoever holds the interval.

#### How many resamples

Efron & Tibshirani (1986) §9, p. 72, recorded as `MIN_RESAMPLES_BCA` and
`MIN_RESAMPLES_PERCENTILE` on `BootstrapCIResult` and re-exported from
`ci_methods`:

| method | rough minimum B |
|--------|-----------------|
| BCa or BC | 1000 |
| percentile | 250 |
| standard errors only | 25 |

**Nothing enforces these, and the interval says when it is under them.** An
interval whose usable draws fall short of its method's minimum reports
`below_recommended_resamples`, and prints it:

```
[0.339, 0.551] (95% CI, bca, n=50, 50 usable, 1000 recommended)
```

The comparison is against `n_usable`, not `n_resamples` -- 1000 draws with
900 failures leaves 100, and 100 is what the quantiles were taken over.

It does not raise and it does not warn. Too few draws makes an interval
noisy, not wrong, so refusing to produce one would be the wrong response;
and a warning would fire once per interval, which for a pairwise run means
one per rater pair. Saying it in the interval's own description puts it
where someone reading the result will see it, exactly once.

### Parametric null (`ci_method="parametric_null"`)

Uses a closed-form null variance -- Fleiss, Nee & Landis (1979) Eq. 12 for
Fleiss' kappa, Fleiss, Cohen & Everitt (1969) for Cohen's. No bootstrap is
performed.

**Not recommended for new work.** The variance is derived under H0:
kappa = 0, so the interval is centred on your estimate and scaled as though
that estimate were zero. Measured against percentile bootstrap intervals on
the shipped datasets:

| dataset / coefficient | bootstrap width | parametric null width | |
|---|---|---|---|
| Fleiss 1971 / Fleiss κ | 0.2166 | 0.0955 | −56% |
| SR1 prelim / Fleiss κ | 0.2265 | 0.0942 | −58% |
| Gwet cont3x3 / Cohen κ | 0.2309 | 0.2979 | +29% |
| Synthetic / Scott π | 0.2599 | 0.2948 | +13% |

**The error is not one-directional.** Fleiss comes out far too narrow,
Cohen and Scott too wide. So there is no correction factor, and no safe
assumption about which way a given number is off.

Provided for reproducing published results that used these formulas, or
comparing against another implementation that does. Available on three
coefficients only -- Fleiss, Cohen and Scott, the ones with a published
closed-form null variance. The other five reject it at config
construction, each with its own explanation.

Choosing it is disclosed twice, in different words:

- **at config construction**, naming the parameter and its alternatives,
  while the caller still knows where the knob is
- **at interval construction**, naming how far off the width is, which is
  only a useful thing to say once a number exists

and the interval says it again whenever it prints:

```
= 0.430 [0.382, 0.478] (95% CI, null variance -- width assumes kappa = 0, reproduction only)
```


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

Bootstrap CI computation: `percentile_ci()` and `bca_ci()`, both returning
a `BootstrapCIResult`. See the confidence interval methods section above.

`scipy_bca_ci()` is gone. It existed as a reference implementation to
validate against, but scipy's BCa computes acceleration from the data it
resamples, which this module never sees -- it receives an array of numbers
and no estimator. Validation against scipy now happens in the tests, where
the data is in hand, rather than through a wrapper that had to guess.

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
