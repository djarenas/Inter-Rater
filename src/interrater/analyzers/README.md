# Analyzers

Each analyzer takes a Ratings object and a config, validates both, and
provides `estimate()`, which returns a result object. An analyzer is built
once: what it works out at construction (the ratings after any injection,
the profiles, the weights) is what `estimate()` and `hypothesis_test()` read.
The data objects it consumes (Ratings, configs) are frozen.


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
  │   └── (direct GroupAnalyzer children, no weighting)
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
      └── (direct IndividualMarginalAnalyzer children, no weighting)
          ├── CongerKappaAnalyzer
          └── LightKappaAnalyzer
```

The splits are meant to be conceptual: Group coefficients treat all raters
as interchangeable; they pool the raters into one distribution; individual marginal coefficients 
use marginal rates. Therefore, some capabilities such as building from nij matrices would be
available to the former but not the latter. 


## What each level does

### Analyzer (base)

- Stores Ratings and config
- Builds the BalanceProfile and, from ratings, the RaterProfile
  (`NotComputed` on the counts path)
- Writes `estimate()` once for every coefficient (see below), with the one
  `_compute_ci` it uses

### GroupAnalyzer

- Validates that the config is a GroupCoefficientConfig
- Supports two initialization paths: from Ratings, or from NijMatrix directly
- Writes `estimate_from_nij_matrix()` once for every group coefficient: the
  counts-path twin of `estimate()`, calling the same `_result`
- Enforces balance constraints when `config.balanced_data_only=True`

### WeightableGroupAnalyzer / WeightableIndividualAnalyzer

- Reads the Weights object from the semantics of the ratings (`Weights.build_from_categories`)
- Validates weight matrix dimensions match the number of categories
- Supports "identity" (unweighted), "linear", "quadratic", and "custom"

### IndividualMarginalAnalyzer

- Requires Ratings (never None: rater identity is essential)
- Builds EstimatorInput from Ratings (category_index + metadata)
- No NijMatrix path: counts discard exactly what these coefficients need

### Leaf analyzers (e.g. FleissKappaAnalyzer)

- Validate their own config type (FleissKappaConfig, etc.)
- Supply only what is specific to their coefficient: the coefficient itself
  (`_estimate_from_array`), the result it is packaged into (`_result`), and,
  for Cohen, Scott and Fleiss, the null-variance interval
  (`_parametric_null_interval`)


## The estimate() flow

`estimate()` is the same for every coefficient, so it is written once, on
`Analyzer`:

1. **Point estimate.** `_estimate_from_array` on the original data.
2. **Bootstrap loop.** Resample subjects with replacement and call
   `_estimate_from_array` on each resample. The bootstrap uses the same
   function as the point, so the two cannot drift apart.
3. **Confidence interval.** Percentile or BCa from the draws; for
   `parametric_null`, `_parametric_null_interval` instead, with no draws.
4. **Return result.** `_result(point, ci, draws)` packages everything into
   the coefficient's frozen result object.

A leaf analyzer writes only the coefficient-specific pieces: steps 2 and 3
need nothing from it beyond `_estimate_from_array`. `_result` is separate
from the rest so that a permutation analysis, which bootstraps each pair or
subset with its own seeds, can package a unit's point, interval and draws
through the same method and get a result built and checked the same way.

The counts path, `estimate_from_nij_matrix()`, is the same four steps on the
counts matrix, written once on `GroupAnalyzer`. `_result` reads its counts and
agreement profile through `_nij_array()` and `_agreement_profile()`, which
take them from whichever source the analyzer was built from, so one `_result`
serves both paths.

```python
analyzer = FleissKappaAnalyzer(ratings, config)
result = analyzer.estimate()

result.point               # 0.430
result.confidence_interval # BootstrapCIResult, not a tuple of two floats
result.ci_method           # "percentile"
result.balance_profile     # BalanceProfile(...)

print(result.confidence_interval)
# [0.316, 0.534] (95% CI, percentile, n=1000)
```


## Confidence interval methods

The intervals are computed in `ci_methods.py` and returned as a
`BootstrapCIResult` (or a `ParametricCIResult`), which records the bounds and
also the method that actually ran, how many draws failed, and whether there
were enough of them.

- **Percentile** (`ci_method="percentile"`): the alpha/2 and 1 - alpha/2
  quantiles of the usable bootstrap draws. No correction. If percentile
  cannot run, there is no interval.
- **BCa** (`ci_method="bca"`): Efron & Tibshirani (1986). Percentile with two
  corrections: `z0` for a bootstrap distribution not centred on the point,
  and the acceleration `a` for a standard error that changes with the
  coefficient. `a` is computed from the leave-one-subject-out estimates, not
  from the draws. When either correction is unavailable, BCa steps down to
  plain percentile and `method_used` says so.
- **Parametric null** (`ci_method="parametric_null"`): the closed-form null
  variance (Fleiss, Nee & Landis 1979 for Fleiss' kappa; Fleiss, Cohen &
  Everitt 1969 for Cohen's). It assumes kappa = 0, so its width is wrong for
  a non-zero estimate, too narrow or too wide depending on the coefficient.
  It is for reproducing published analyses. Only Cohen, Scott and Fleiss
  have one; the other configs refuse it.

Draws whose coefficient is undefined are left out and counted (`n_failed`).
An interval with fewer usable draws than Efron & Tibshirani's rough minimum
(1000 for BCa, 250 for percentile) says so when printed:

```
[0.339, 0.551] (95% CI, bca, n=50, 50 usable, 1000 recommended)
```

Why the code does it this way (where `a` comes from, when BCa falls back,
the zero-spread test): the `ci_methods.py` docstring. The checks against
published answers, and how the interval objects report themselves: the
`intermediates/results/ci_results` README.


## Subfolders

### `agreement_functions/`

Pure functions that compute each coefficient from its inputs. Called inside
bootstrap loops thousands of times, they do no validation, no resampling and
no CI computation. They take a counts matrix (or estimator input) and return
a float.

See `agreement_functions/README.md`.

### `resampling_functions/`

Generators that yield modified arrays: subject bootstrap and jackknife (on
ratings or on a counts matrix) and the three permutation shuffles. All return
iterators, so one resampled array is held in memory at a time.

See `resampling_functions/README.md`.

### `permutation_analysis/`

Orchestrators that run a coefficient across subsets of the rater panel.
They are standalone, not in the Analyzer tree. An Analyzer computes one
coefficient on one data slice. An Orchestrator enumerates subsets (pairs,
leave-one-out, arbitrary k) and calls an Analyzer or estimator on each,
then assembles the results.

Leave-one-rater-out is subset permutation at k = n - 1: `subset_k: 3` on a
4-rater panel gives the coefficient without each rater in turn, with its
interval and its subject count. The `rater_jackknife` config flag, which
would report this on the analyzer's own result, is refused until 2.2.

See `permutation_analysis/README.md`.

### `configs/`

Frozen dataclass configurations for every analyzer. The config hierarchy
mirrors the analyzer hierarchy and refuses, at construction, settings a
coefficient cannot use (e.g. `parametric_null` for a coefficient with no
null variance).

See `configs/README.md`.


## Other modules

### `ci_methods.py`

`percentile_ci()` and `bca_ci()`, both returning a `BootstrapCIResult`. See
the confidence interval methods section above. BCa is checked against the
worked example in Efron & Tibshirani (1986) in
`tests/test_analyzers/test_ci_methods.py` (`TestEfronLawSchool`).

### `variance_functions.py`

Closed-form variance formulas from the statistical literature. Used for
parametric null CIs and hypothesis tests. Not needed for bootstrap
approaches.

### `convergence_test.py`

Stub for bootstrap convergence diagnostics (not yet implemented). The type
exists so that `ConvergenceTest | None` fields resolve in configs; a config
given one is refused when it is built.


## Design principles

**Built once, then read.** An analyzer applies any injection and builds its
profiles and weights in its constructor. After that nothing changes it;
`estimate()` and `hypothesis_test()` only read what it holds.

**Data objects are frozen.** Ratings and configs are frozen dataclasses. An
analyzer that injects raters builds a new Ratings; the one it was given is
never changed.

**Validation flows upward.** A leaf class validates its own constraints and
calls `super().__init__()`, which validates the intermediate class's, up to
the base. Each level trusts that its parent has validated what it is
responsible for.

**Agreement functions are hot-path code.** They run inside bootstrap loops
(thousands of calls), so they do no validation. All checking happens once,
before the loop, in the analyzer.

**Resampling is separated from computation.** Resampling functions yield
modified arrays. Agreement functions compute coefficients on arrays. The
analyzer composes them.

**A step that drops raters measures only what is left.** A pair, a subset,
or the panel without one rater computes its point, and resamples, only on
the subjects still rated by at least two of its raters. With fewer than two
such subjects, or a rater who rated none of them, it has no coefficient and
says why. The rule is written once (`unit_problem` in
`agreement_functions/rater_combinations.py`) and every such step calls it;
see `permutation_analysis/README.md`.
