# Analysis Results

The output objects that analyzers produce. Each result is a frozen
dataclass holding a point estimate, optional confidence interval, and
enough metadata to interpret, serialize, and visualize the outcome
without going back to the raw data.

## Hierarchy

```
AnalysisResult (ABC)
├── PointEstimateResult
│   ├── FleissKappaResult          ci_method + weights_type, is_exact
│   ├── ScottsPiResult             ci_method, is_exact
│   ├── GwetsAC1Result             ci_method, is_exact
│   ├── BrennanPredigerResult      ci_method + weights_type, is_exact
│   ├── KrippendorffAlphaResult    ci_method + weights_type, is_exact (always True)
│   ├── IndividualMarginalResult   (internal -- validates balance_profile + n_raters)
│   │   ├── CohenKappaResult       ci_method + weights_type
│   │   └── CongerKappaResult      ci_method
│   └── LightKappaResult           n_pairs + ci_method
├── PermutationResult
│   ├── PairwisePermutationResult
│   ├── SubsetPermutationResult
│   └── MultiKSubsetResult

HypothesisTestResult               (standalone, not on the tree)
```

## Classes

| Class                      | What It Holds                                                       | Why It Exists                                                                 |
|----------------------------|---------------------------------------------------------------------|-------------------------------------------------------------------------------|
| `AnalysisResult`           | ABC with `summary() -> dict`                                       | Common contract: every result can describe itself in plain types.             |
| `PointEstimateResult`      | point, CI, bootstrap distribution, balance profile                  | Shared base for all single-coefficient results.                               |
| `IndividualMarginalResult` | (inherits PointEstimateResult)                                      | Internal base for coefficients requiring identified raters (balance_profile and n_raters must not be None). Not exported. |
| `FleissKappaResult`        | ci_method, weights_type, is_exact                                   | Fleiss' kappa for multiple raters. Exact when balanced.                       |
| `ScottsPiResult`           | ci_method, is_exact                                                 | Scott's pi. No weights. Exact when balanced.                                  |
| `GwetsAC1Result`           | ci_method, is_exact                                                 | Gwet's AC1. No weights. Rejects parametric_null (not implemented).            |
| `BrennanPredigerResult`    | ci_method, weights_type, is_exact                                   | Brennan-Prediger coefficient. Exact when balanced.                            |
| `KrippendorffAlphaResult`  | ci_method, weights_type, is_exact (always True)                     | Krippendorff's alpha. Handles unbalanced data natively -- always exact.       |
| `CohenKappaResult`         | ci_method, weights_type                                             | Cohen's kappa for two identified raters. Weighted or unweighted.              |
| `CongerKappaResult`        | ci_method                                                           | Conger's kappa. No weights. Rejects parametric_null (not available).          |
| `LightKappaResult`         | n_pairs, ci_method                                                  | Light's kappa -- mean of C(n,2) pairwise Cohen's kappas. Always unweighted.  |
| `PermutationResult`        | point, CI, bootstrap distribution, balance profile, k_subsets       | Base for permutation-based analyses (leave-k-out).                            |
| `PairwisePermutationResult`| per-pair kappas, spread, pair labels                                | All C(n,2) pairwise kappas with spread diagnostics.                           |
| `SubsetPermutationResult`  | subset results, rater labels                                        | Leave-one-out or leave-k-out at a single k.                                  |
| `MultiKSubsetResult`       | results across multiple k values                                    | Stability curve: how the coefficient changes as raters are removed.           |
| `HypothesisTestResult`     | p_value, alpha, reject_null, method, point, optional z/null dist    | Test against H0: coefficient = 0. Parametric or permutation. Not on the tree.|

## Design decisions

- **Frozen dataclasses.** A result is a fact about a dataset at a moment
  in time. It should not change after construction.

- **`is_exact` varies by coefficient.** Group coefficients (Fleiss, Scott,
  Gwet, Brennan-Prediger) are exact when balanced and approximate otherwise.
  Krippendorff is always exact (it handles missing data natively). Pairwise
  coefficients (Cohen, Conger, Light) do not carry `is_exact` -- they
  require balanced paired data by construction.

- **`ci_method` allowed values vary.** Fleiss, Scott, and Cohen accept
  `"parametric_null"`. Gwet and Conger explicitly reject it (not
  implemented / not available). Brennan-Prediger, Krippendorff, and
  Light allow only `"percentile"` or `"bca"`.

- **`weights_type` is on the result, not just the config.** Weighted
  Cohen's kappa and unweighted Cohen's kappa are different statistics.
  A result must be self-describing. Present on Fleiss, Brennan-Prediger,
  Krippendorff, and Cohen. Absent from Scott, Gwet, Conger (always
  unweighted), and Light (defined as unweighted by Light 1971).

- **`IndividualMarginalResult` is internal.** It validates that
  `balance_profile` is not None and `n_raters` is not None -- the
  invariant shared by Cohen and Conger, which require identified raters.
  Not exported from `__init__.py`.

- **`HypothesisTestResult` is standalone.** A hypothesis test is a
  fundamentally different kind of output from a point estimate with
  confidence interval. It does not inherit from `AnalysisResult`.

- **Error chaining across inheritance.** Each class has its own
  `_ERROR_PREFIX`. Child `__post_init__` wraps `super().__post_init__()`
  in `try/except` so the child prefix appears in the message even when
  the parent validation fails.

## Dependencies

- `data_profiles` (BalanceProfile)
- `numpy` (bootstrap distributions)
