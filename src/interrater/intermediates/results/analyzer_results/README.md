# Analyzer Results

One coefficient, one estimate, one dataset. Every class here descends from
`PointEstimateResult`, which is the shape an analyzer returns when it has
estimated a single coefficient: a point estimate, an optional confidence
interval, the bootstrap distribution behind it, and a balance profile.

```
PointEstimateResult
├── FleissKappaResult          ci_method + weights_type + weight_matrix, is_exact
├── ScottsPiResult             ci_method, is_exact
├── GwetsAC1Result             ci_method, is_exact
├── BrennanPredigerResult      ci_method + weights_type + weight_matrix, is_exact
├── KrippendorffAlphaResult    ci_method + weights_type + weight_matrix, is_exact (always True)
├── IndividualMarginalResult   (internal -- validates balance_profile + n_raters)
│   ├── CohenKappaResult       ci_method + weights_type + weight_matrix
│   └── CongerKappaResult      ci_method
└── LightKappaResult           n_pairs + ci_method
```

| Class                      | What It Holds                                   | Why It Exists                                                                 |
|----------------------------|--------------------------------------------------|-------------------------------------------------------------------------------|
| `PointEstimateResult`      | point, CI, bootstrap distribution, balance profile | Shared base for all single-coefficient results.                              |
| `IndividualMarginalResult` | (inherits `PointEstimateResult`)                  | Internal base for coefficients requiring identified raters. Not exported.     |
| `FleissKappaResult`        | ci_method, weights_type, weight_matrix, is_exact  | Fleiss' kappa for multiple raters. Exact when balanced.                       |
| `ScottsPiResult`           | ci_method, is_exact                               | Scott's pi. No weights. Exact when balanced.                                  |
| `GwetsAC1Result`           | ci_method, is_exact                               | Gwet's AC1. No weights. Rejects `parametric_null` (not implemented).          |
| `BrennanPredigerResult`    | ci_method, weights_type, weight_matrix, is_exact  | Brennan-Prediger coefficient. Exact when balanced.                            |
| `KrippendorffAlphaResult`  | ci_method, weights_type, weight_matrix, is_exact (always True) | Krippendorff's alpha. Handles unbalanced data natively -- always exact.       |
| `CohenKappaResult`         | ci_method, weights_type, weight_matrix            | Cohen's kappa for two identified raters. Weighted or unweighted.              |
| `CongerKappaResult`        | ci_method                                         | Conger's kappa. No weights. Rejects `parametric_null` (not available).        |
| `LightKappaResult`         | n_pairs, ci_method                                | Light's kappa -- mean of C(n,2) pairwise Cohen's kappas. Always unweighted.   |

## Why the fields differ by coefficient

**`is_exact`.** Group coefficients (Fleiss, Scott, Gwet, Brennan-Prediger)
are exact when the data are balanced and approximate otherwise. Krippendorff
is always exact -- it handles missing data natively. The pairwise
coefficients (Cohen, Conger, Light) carry no `is_exact` at all: they require
balanced paired data by construction, so the question does not arise.

**`ci_method` allowed values.** Fleiss, Scott and Cohen accept
`"parametric_null"`, because a closed-form null variance exists for them.
Gwet and Conger reject it explicitly -- not implemented, and not available,
respectively. Brennan-Prediger, Krippendorff and Light allow only
`"percentile"` or `"bca"`.

**`weights_type` and `weight_matrix` are on the result.** Weighted Cohen's
kappa and unweighted Cohen's kappa are different statistics, and a result
must be able to say which one it is. The weights are declared in the
semantics of the ratings, not in the analyzer config, so the result is where
a reader finds them. `weights_type` names the scheme ("identity", "linear",
"quadratic" or "custom"); `weight_matrix` holds the matrix itself, as nested
lists in `summary()`, whenever the result is weighted, and is `None` when it
is not. The matrix is recorded every time, not only for "custom", because a
scheme's name pins down the numbers only together with the category order.
Present on Fleiss, Brennan-Prediger, Krippendorff and Cohen. Absent from
Scott, Gwet and Conger (always unweighted) and from Light (defined as
unweighted by Light 1971).

## `IndividualMarginalResult` is internal

It validates the invariant Cohen and Conger share -- `balance_profile` and
`n_raters` must both be present, because both coefficients require
identified raters rather than counts. It is a base, not an outcome, so it is
deliberately absent from `__all__`.

## Dependencies

- `..analysis_result` (AnalysisResult)
- `...data_profiles` (BalanceProfile)
- `numpy` (bootstrap distributions)
