# Analyzer Configurations

Every analyzer is configured via a frozen dataclass that declares all analytical decisions upfront: bootstrap resamples, confidence level, CI method, random seed, hypothesis testing, injection pipeline, and optional convergence diagnostics. Configurations load from YAML via `build_from_yaml` or from plain dicts via `build_from_dict`, and are fully auditable via `summary()`.

## Philosophy

**Configuration enforces constraints at the boundary.** Each coefficient has mathematical requirements -- balance, weighting support, closed-form null variance. Rather than letting the analyzer discover these at runtime, the config hierarchy encodes them as type constraints that reject invalid combinations at instantiation time with clear error messages.

**Validation flows upward.** Leaf class → intermediate → parent → grandparent. Each class calls `super().__post_init__()` to inherit parent validation, wrapping it in `try/except` to prefix its own class name on errors.

**Version for audit.** Every config's `summary()` includes `_package_version` so every analysis knows which version of interrater created it.

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
  └── PermutationConfig (holds nested AnalyzerConfig)
      ├── PairwisePermutationConfig
      └── SubsetPermutationConfig

Supporting types:
  HypothesisTestConfig
  PermutationBootstrapMeanConfig
```

Weights are not configured here. They are declared in the semantics of the
ratings (categories of kind `ordered` with a `weighting`, or kind `weighted`
with a matrix), and the analyzer reads them from there. `weights_type` and
`custom_weights_path` were removed in 2.1.0; a YAML config that still has them
is refused with a message saying so.

## Base: AnalyzerConfig

All analyzers share these fields:

- `n_resamples` (int, default 1000) -- number of bootstrap draws
- `alpha` (float, default 0.05) -- significance level; 0.05 gives 95% CI via [0.025, 0.975]
- `random_state` (int | None, default 42) -- random seed for reproducibility; None for unseeded
- `ci_method` (str, default "percentile") -- "percentile", "bca", or "parametric_null"
- `convergence_test` (ConvergenceTest | None, default None) -- optional diagnostics verifying n_resamples is large enough for stable CI bounds
- `hypothesis_test_config` (HypothesisTestConfig | None, default None) -- optional hypothesis testing against H₀: coefficient = 0
- `injection` (tuple[InjectorConfig, ...] | None, default None) -- optional pipeline of injection steps applied to ratings before analysis; accepts a single config, list, or tuple (normalized to tuple in `__post_init__`)

Validates all fields in `__post_init__`. Provides `build_from_dict`, `build_from_yaml`, and `summary`.

## GroupCoefficientConfig

Base for coefficients that read counts matrices (discard rater identity). Adds:

- `balanced_data_only` (bool, default False) -- whether to enforce constant raters per subject
- `rater_jackknife` (bool, default False) -- whether to run leave-one-rater-out diagnostics

Both are validated as strict `bool` (rejects `int` to catch `balanced_data_only=1`).

### NoWeightGroupConfig

For coefficients whose chance terms are defined independently of how raters use categories. Weighting would create an inconsistency: weighted observed agreement but unweighted chance term.

Adds no fields. The analyzers for these coefficients refuse ratings whose semantics declares ordered or weighted categories.

Leaf configs: **GwetAC1Config**, **ScottsConfig**.

### WeightableGroupConfig

For coefficients that support distance-based weighting. Adds no fields: the weights come from the semantics of the ratings, and these are the coefficients whose analyzers accept them.

Leaf configs: **FleissKappaConfig**, **BrennanPredigerConfig**, **KrippendorffAlphaConfig**.

## IndividualMarginalConfig

Base for coefficients that keep rater identity (pairwise or per-rater marginals). Adds:

- `rater_jackknife` (bool, default False) -- currently validated to stay False

These coefficients read EstimatorInput (not NijMatrix) and maintain which rater said what.

### NoWeightIndividualConfig

For individual-marginal coefficients with no weighted forms. Same validation as NoWeightGroupConfig.

Leaf configs: **CongerKappaConfig**, **LightKappaConfig**.

### WeightableIndividualConfig

For individual-marginal coefficients that support weighting. Same fields and validation as WeightableGroupConfig, plus enforces `rater_jackknife` stays False.

Leaf configs: **CohenKappaConfig**.

## Leaf Configs: Parametric Null and Balance

Not every coefficient has a closed-form null variance. Only those with a derivable null distribution support `ci_method="parametric_null"` and `hypothesis_test_config.method="parametric_null"`. All others reject these at instantiation.

| Coefficient | Parent | Weights | parametric_null | balanced_data_only default |
|---|---|---|---|---|
| FleissKappaConfig | WeightableGroup | identity, linear, quadratic, custom | yes (Fleiss 1979, Eq. 12) | **True** |
| ScottsConfig | NoWeightGroup | identity only | yes (= Fleiss at n=2) | False |
| GwetAC1Config | NoWeightGroup | identity only | no | False |
| BrennanPredigerConfig | WeightableGroup | identity, linear, quadratic, custom | no | False |
| KrippendorffAlphaConfig | WeightableGroup | identity, linear, quadratic, custom | no | False |
| CohenKappaConfig | WeightableIndividual | identity, linear, quadratic, custom | yes (Fleiss, Cohen & Everitt 1969) | **True** |
| CongerKappaConfig | NoWeightIndividual | identity only | no | **True** |
| LightKappaConfig | NoWeightIndividual | identity only | no | False |

**Auto-correction with warning:** CohenKappaConfig and CongerKappaConfig both require balanced data. If `balanced_data_only=False` is passed, they auto-correct to `True` with a `warnings.warn()` explaining why.

## Permutation Configs

PermutationConfig is a separate hierarchy (not a child of GroupCoefficientConfig or IndividualMarginalConfig). It wraps a nested AnalyzerConfig to compute a coefficient across multiple rater combinations.

### PermutationConfig (base)

- `coefficient_to_compute` (str) -- the coefficient name ("cohen_kappa", "fleiss_kappa", etc.)
- `config_of_coefficient` (AnalyzerConfig) -- the nested config for that coefficient
- `subset_k` (int | list[int] | str | None, default None) -- subset size(s) for subset permutation
- `bootstrap_mean_config` (PermutationBootstrapMeanConfig | None, default None) -- mean-level bootstrap configuration

### PairwisePermutationConfig

For all C(n, 2) pairwise analyses. Validates:
- `subset_k` must be None (pairwise is always k=2)
- `coefficient_to_compute` must be in PAIRWISE_PERMUTATION_COEFFICIENTS
- `config_of_coefficient` must be the correct type for the coefficient
- Cross-validates: "independent" bootstrap method requires per-pair bootstrap enabled (n_resamples >= 2)

### SubsetPermutationConfig

For all C(n, k) subset analyses. Validates:
- `subset_k` is required (not None)
- Accepts int (>= 2), list[int] (non-empty, all >= 2, no duplicates), or "all"
- `coefficient_to_compute` must be in SUBSET_PERMUTATION_COEFFICIENTS
- `config_of_coefficient` must match the coefficient
- Cross-validates bootstrap method against n_resamples

## Supporting Types

### HypothesisTestConfig

Configuration for hypothesis testing against H₀: coefficient = 0.

- `method` (str, default "permutation") -- "permutation" or "parametric_null"
- `scheme` (str, default "column_shuffle") -- "column_shuffle", "within_subject", or "pool_and_redistribute"
- `n_permutations` (int, default 5000) -- number of permutation draws
- `alpha` (float, default 0.05) -- significance level
- `random_state` (int | None, default 42) -- random seed

All numeric fields have bool-first type guards (bool subclasses int in Python, so `True` must be caught before `isinstance(x, int)`).

### PermutationBootstrapMeanConfig

Controls bootstrapping of the mean coefficient across combinations.

- `method` (str) -- "correlated" (resample subjects, recompute all combinations) or "independent" (assemble mean from independently-drawn per-combination values)
- `n_bootstraps` (int) -- number of bootstrap draws for the mean
- `ci_method` (str, default "percentile") -- "percentile" or "bca" (parametric_null rejected)
- `alpha` (float, default 0.05) -- significance level
- `random_state` (int | None, default None) -- random seed

## Using Configurations

### From YAML

```python
config = FleissKappaConfig.build_from_yaml("fleiss_config.yaml")
```

YAML structure:

```yaml
n_resamples: 2000
alpha: 0.05
random_state: 42
ci_method: percentile
balanced_data_only: true
hypothesis_test:                     # alias for hypothesis_test_config
  method: permutation
  scheme: column_shuffle
  n_permutations: 5000
```

Unknown keys are rejected with a clear error listing valid keys (aliases included).

### From dict

```python
config = FleissKappaConfig.build_from_dict({
    "n_resamples": 2000,
    "alpha": 0.05,
})
```

`build_from_dict` normalizes YAML key aliases and hydrates nested objects (sub-dicts become typed config objects).

### Programmatically

```python
config = FleissKappaConfig(
    n_resamples=2000,
    alpha=0.05,
)
```

### Audit trail

Every config's `summary()` returns a dict with all fields plus `_package_version`:

```python
config.summary()
# {
#     "n_resamples": 2000,
#     "alpha": 0.05,
#     "random_state": 42,
#     "ci_method": "percentile",
#     "balanced_data_only": True,
#     "rater_jackknife": False,
#     "_package_version": "2.1.3",
# }
```

## Design Notes

**Frozen dataclasses.** All configs are `@dataclass(frozen=True)`. Immutability guarantees that a config cannot be modified after validation. When auto-correction is needed (e.g., forcing `balanced_data_only=True`), the class uses `object.__setattr__` to bypass the frozen guard in `__post_init__`.

**`_ERROR_PREFIX` at module level.** Every config module defines a module-level `_ERROR_PREFIX` string (e.g., `"[CohenKappaConfig] Error"`) used in all error messages for grep-friendly traceability.

**TypeError vs ValueError separation.** Wrong type → TypeError. Wrong value of the right type → ValueError. Bool-first guards catch `True`/`False` before `isinstance(x, int)` since bool subclasses int.

**`warnings.warn()` not `print()`.** Library code uses `warnings.warn(..., UserWarning, stacklevel=2)` for auto-corrections so users can catch or silence them.

## Future: Extended Weight Types

Several coefficients (Gwet's AC2, Krippendorff's alpha) support ordinal, interval, and ratio weighting schemes beyond what the semantics can declare today (unordered, ordered linear or quadratic, or a weighted matrix):

- Ordinal weights (Krippendorff-style distance-based)
- Interval weights (continuous measurement distances)
- Ratio weights (ratio-scaled measurement distances)

When AC2 is implemented, it will inherit from WeightableGroupConfig. The categories in the semantics, and the Weights class built from them, will be extended to support these additional schemes. See Gwet (2014) for the AC2 formulation.
