# Analyzer Configurations

Every analyzer is configured by a frozen dataclass that declares its
decisions up front. One of their utilities is that they can be built from
a YAML file that can be easily archived, modified, and/or shown to collaborators (or your future self).

The decisions are such as: number of bootstrap resamples, confidence level, CI method, random
seed, hypothesis testing, injection, and convergence diagnostics.

Configurations were designed to load from YAML (`build_from_yaml`), but the plain dicts
(`build_from_dict`) option is also availble. `summary()` reports them.

## Principles taking advantage of the OOP

**Refuse nonsensical choices at construction.** A setting a coefficient cannot use (e.g.
`parametric_null` for a coefficient with no null variance) is refused when
the config is built, with a message saying why. Requirements that depend on
the data, such as weights from the semantics, are checked by the analyzer. More obvious mistakes
such as strings in numbers, or negative numbers of bootstraps, are also refused this early
with easily traceable errors. The object is not built if mistakes are found.

**Validation flows upward in the inheritance** Leaf, then intermediate, then base: each class
calls `super().__post_init__()` and wraps it so its own name prefixes any
error.

**Versioned.** Every `summary()` includes `_package_version`.

## Hierarchy

```
AnalyzerConfig (base)
  ├── GroupCoefficientConfig
  │   ├── NoWeightGroupConfig
  │   │   ├── GwetAC1Config
  │   │   └── ScottsConfig
  │   │
  │   └── WeightableGroupConfig
  │       ├── FleissKappaConfig
  │       ├── BrennanPredigerConfig
  │       └── KrippendorffAlphaConfig
  │
  ├── IndividualMarginalConfig
  │   ├── NoWeightIndividualConfig
  │   │   ├── CongerKappaConfig
  │   │   └── LightKappaConfig
  │   │
  │   └── WeightableIndividualConfig
  │       └── CohenKappaConfig
  │
  └── PermutationConfig (holds a nested AnalyzerConfig)
      ├── PairwisePermutationConfig
      └── SubsetPermutationConfig

Supporting types:
  HypothesisTestConfig
  PermutationBootstrapMeanConfig
```

Weights are not configured here. They are declared in the semantics of the
ratings (categories of kind `ordered` with a `weighting`, or kind `weighted`
with a matrix), and the analyzer reads them from there. `weights_type` and
`custom_weights_path` were removed in 2.1.0; a config that still has them is
refused with a message saying so.

## AnalyzerConfig

Fields every analyzer shares:

- `n_resamples` (int, default 1000): bootstrap draws
- `alpha` (float, default 0.05): 0.05 gives a 95% interval
- `random_state` (int | None, default 42): seed; None for unseeded
- `ci_method` (str, default "percentile"): "percentile", "bca", or "parametric_null"
- `convergence_test` (ConvergenceTest | None, default None): not implemented yet
- `hypothesis_test_config` (HypothesisTestConfig | None, default None): test against H₀: coefficient = 0
- `injection` (tuple[InjectorConfig, ...] | None, default None): injection steps applied to the ratings before analysis; a single config or a list is normalised to a tuple

## GroupCoefficientConfig

Coefficients that read a counts matrix. Adds:

- `balanced_data_only` (bool, default False): require the same number of raters on every subject
- `rater_jackknife` (bool, default False): leave-one-rater-out; not supported yet, so `true` is refused

Both must be a real `bool` (`balanced_data_only=1` is refused).

- **NoWeightGroupConfig** (Gwet's AC1, Scott's pi): no weighted form. Adds
  no fields; the analyzers refuse ordered or weighted semantics.
- **WeightableGroupConfig** (Fleiss, Brennan-Prediger, Krippendorff): adds
  no fields; these analyzers use the semantics' weights.

## IndividualMarginalConfig

Coefficients that keep rater identity and read an `EstimatorInput`. Adds
`rater_jackknife` (default False; `true` is refused).

- **NoWeightIndividualConfig** (Conger, Light): no weighted form.
- **WeightableIndividualConfig** (Cohen): uses the semantics' weights.

## Leaf configs: parametric null and balance

Only coefficients with a published null variance accept `parametric_null`,
as `ci_method` or as the hypothesis test's method. The others refuse it when
the config is built.

| Config | Weights | parametric_null | balanced_data_only |
|---|---|---|---|
| FleissKappaConfig | yes | yes (Fleiss, Nee & Landis 1979, Eq. 12) | default **True** |
| ScottsConfig | no | yes (= Fleiss at n = 2) | default False |
| GwetAC1Config | no | no | default False |
| BrennanPredigerConfig | yes | no | default False |
| KrippendorffAlphaConfig | yes | no | default False |
| CohenKappaConfig | yes | yes (Fleiss, Cohen & Everitt 1969) | **True**; `false` refused |
| CongerKappaConfig | no | no | **True**; `false` refused |
| LightKappaConfig | no | no | **True**; `false` refused |

Cohen is always computed on the subjects both raters rated, and Conger and
Light need complete data, so `balanced_data_only: false` does not apply to
them and is refused with a message saying why.

## Permutation configs

`PermutationConfig` wraps a nested coefficient config to compute that
coefficient over rater combinations. Fields:

- `coefficient_to_compute` (str): e.g. "cohen_kappa", "fleiss_kappa"
- `config_of_coefficient` (AnalyzerConfig): the nested config
- `subset_k` (int | list[int] | str | None, default None): subset size(s), for subset permutation
- `bootstrap_mean_config` (PermutationBootstrapMeanConfig | None, default None): bootstrap of the mean
- `balanced_data_only` (bool, default False): `true` refuses an unbalanced panel before anything runs
- `keep_draws` (bool, default True): keep each pair's or subset's draws; `false` discards them after use, with the same values and intervals
- `max_evaluations` (int, default `GUARD_MAX_EVALUATIONS`, 2,000,000): the most coefficient evaluations a run may make; a larger run is refused before it starts

`ci_method: parametric_null` in `config_of_coefficient` is refused:
combinations sharing a rater are correlated, so no closed-form variance
describes them.

- **PairwisePermutationConfig**: `subset_k` must be None. The coefficient
  must be one of `PAIRWISE_PERMUTATION_COEFFICIENTS`, with a matching config.
- **SubsetPermutationConfig**: `subset_k` is required: an int >= 2, a
  non-empty list of distinct ints >= 2, or "all". The coefficient must be
  one of `SUBSET_PERMUTATION_COEFFICIENTS` (Scott's pi is not), with a
  matching config.

Both check that the independent mean has at least 2 draws per combination
and `ci_method: percentile`.

## Supporting types

### HypothesisTestConfig

- `method` (str, default "permutation"): "permutation" or "parametric_null"
- `scheme` (str, default "within_rater_shuffle"): "within_rater_shuffle" or "pool_and_redistribute". Both keep the missing-data structure: answers move only between rated cells. `"column_shuffle"` (renamed) and `"within_subject"` (cannot test coefficient = 0) are refused with a message saying so
- `alternative` (str, default "two-sided"): "two-sided", "greater" or "less". Two-sided by default because below-chance agreement (a scale read backwards, two labels swapped) is a real finding a one-sided test misses
- `n_permutations` (int, default 5000)
- `alpha` (float, default 0.05)
- `random_state` (int | None, default 42)

Numeric fields check for `bool` first, since `True` is an `int` in Python.

### PermutationBootstrapMeanConfig

- `method` (str): "correlated" (resample subjects, recompute every
  combination) or "independent" (one value from each combination's draws)

Its only setting. The draw count, interval method, level and seed are those
in `config_of_coefficient`. Any other key is refused.

## Using configurations

```python
config = FleissKappaConfig.build_from_yaml("fleiss_config.yaml")
config = FleissKappaConfig.build_from_dict({"n_resamples": 2000, "alpha": 0.05})
config = FleissKappaConfig(n_resamples=2000, alpha=0.05)
```

YAML:

```yaml
n_resamples: 2000
alpha: 0.05
random_state: 42
ci_method: percentile
balanced_data_only: true
hypothesis_test:                     # alias for hypothesis_test_config
  method: permutation
  scheme: within_rater_shuffle
  alternative: two-sided
  n_permutations: 5000
```

`build_from_dict` resolves key aliases and turns nested dicts into config
objects. Unknown keys are refused, with the valid keys listed.

`summary()`:

```python
FleissKappaConfig().summary()
# {'n_resamples': 1000, 'alpha': 0.05, 'random_state': 42,
#  'ci_method': 'percentile', '_package_version': '2.1.0',
#  'balanced_data_only': True, 'rater_jackknife': False}
```

## Design notes

- **Frozen dataclasses.** A config cannot change after validation.
  `__post_init__` uses `object.__setattr__` only to normalise `injection` to
  a tuple.
- **`_ERROR_PREFIX`** in every module (e.g. `"[CohenKappaConfig] Error"`)
  starts every error message.
- **TypeError vs ValueError.** Wrong type raises TypeError; a wrong value of
  the right type raises ValueError.
- **Warnings.** The only warning is for `parametric_null`, raised with
  `warnings.warn(..., UserWarning, stacklevel=2)` so it can be caught or
  silenced.

## Krippendorff's interval alpha

Krippendorff's papers name versions of alpha by scale: nominal, interval,
ratio, ordinal. Here they are chosen through the semantics instead. Interval
alpha is alpha with `kind: ordered, weighting: quadratic` (for equally
spaced categories); nominal alpha is alpha with `kind: unordered`. Ratio and
ordinal distances can be given as a `weighted` matrix.

Gwet's AC2, the weighted AC1, is not implemented. When it is, it will
inherit from `WeightableGroupConfig`. See Gwet (2014).
