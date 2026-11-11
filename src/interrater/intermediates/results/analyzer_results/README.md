# Analyzer Results

These classes are organized so that each coefficient gets its own result.
Every class descends from
`PointEstimateResult`, which is the shape an analyzer returns when it has
estimated a single coefficient: a point estimate, an optional confidence
interval, the bootstrap distribution behind it, the observed and chance
terms it was built from (`parts`), its ceiling (`kappa_max`, if appropriate), and three
profiles of the data (balance, rater(what is this), agreement). The parts must rebuild
the point; a result whose parts do not (do not what) is refused.

```
PointEstimateResult
├── FleissKappaResult          ci_method + weights_type + weight_matrix, is_exact
├── ScottsPiResult             ci_method, is_exact
├── GwetsAC1Result             ci_method, is_exact
├── BrennanPredigerResult      ci_method + weights_type + weight_matrix, is_exact
├── KrippendorffAlphaResult    ci_method + weights_type + weight_matrix, is_exact (always True)
├── IndividualMarginalResult   (internal -- requires a RaterProfile and per-rater marginals)
│   ├── CohenKappaResult       ci_method + weights_type + weight_matrix
│   └── CongerKappaResult      ci_method
└── LightKappaResult           n_pairs + ci_method + pair_parts + pair_kappa_max
```

| Class                      | What It Holds                                   | Why It Exists                                                                 |
|----------------------------|--------------------------------------------------|-------------------------------------------------------------------------------|
| `PointEstimateResult`      | point, CI, bootstrap distribution, parts, kappa_max, balance, rater and agreement profiles | Shared base for all single-coefficient results.                              |
| `IndividualMarginalResult` | (inherits `PointEstimateResult`)                  | Internal base for coefficients built from each rater's own marginals. Not exported. |
| `FleissKappaResult`        | ci_method, weights_type, weight_matrix, is_exact  | Fleiss' kappa for multiple raters. Exact when balanced.                       |
| `ScottsPiResult`           | ci_method, is_exact                               | Scott's pi. No weights. Exact when balanced.                                  |
| `GwetsAC1Result`           | ci_method, is_exact                               | Gwet's AC1. No weights. Rejects `parametric_null` (not implemented).          |
| `BrennanPredigerResult`    | ci_method, weights_type, weight_matrix, is_exact  | Brennan-Prediger coefficient. Exact when balanced.                            |
| `KrippendorffAlphaResult`  | ci_method, weights_type, weight_matrix, is_exact (always True) | Krippendorff's alpha. Handles unbalanced data natively -- always exact.       |
| `CohenKappaResult`         | ci_method, weights_type, weight_matrix            | Cohen's kappa for two identified raters. Weighted or unweighted.              |
| `CongerKappaResult`        | ci_method                                         | Conger's kappa. No weights. Rejects `parametric_null` (not available).        |
| `LightKappaResult`         | n_pairs, ci_method, pair_parts, pair_kappa_max    | Light's kappa -- mean of C(n,2) pairwise Cohen's kappas, each pair's Po, Pe and ceiling keyed by rater names. Always unweighted. |

## Why the fields differ by coefficient

The information that is printed, or visualized, has context that 
depends on the coefficient. For example:

**`is_exact`.** A boolean on whether the data is balanced. This is
relevant because for balanced data different literature references
give the same formula; unbalanced data is not as straightforward. Therefore,
in interrater Fleiss kappa uses an approximation (each subject weighted equally) out
of other other options. Interrater philosophy is that this choice should live in the 
result object, and it therefore prints the approximation information. Contrast, 
Gwet's AC1 and Brennan-Prediger print Gwet's (2014) formula for missing ratings. Contrast, Scott is two raters with every subject rated twice, so it is always balanced. Krippendorff is 
originally built to always handle incomplete and unbalanced data so it 
it works in pairs of ratings; missing data does not change the formula.
Other coefficients, Cohen, Conger and Light carry no `is_exact`. As they refuse
incomplete data, the result objects should know not to take is_exact false. 

**`ci_method` allowed values.** The result objects should have the logic
of what ci methods they can report. For example, Fleiss, Scott and Cohen accept
`"parametric_null"`, because a closed-form null variance exists for them.
Gwet and Conger reject it explicitly, not implemented, and not available,
respectively. Brennan-Prediger, Krippendorff and Light allow only
`"percentile"` or `"bca"`. Clean result objects prevent future implementations or extended
features of feeding the wrong parameters to a result cconcept. 

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

## Draws and the interval

`bootstrap_distribution` is present exactly when the interval is a
`BootstrapCIResult`, and must agree with its counts (see `../README.md`).
`summary()` describes the draws with `bootstrap_draws_mean` and
`bootstrap_se`, computed over the usable (finite) draws, with the standard
error dividing by n - 1. Either is `None` when too few draws are usable. The
draw counts are in the interval's own summary, under `ci`.

## `IndividualMarginalResult` is internal

It validates what Cohen and Conger share: `rater_profile` must be a
`RaterProfile`, not `NotComputed`, because both chance terms need to know
which rater said what, so neither can be built from counts; and the
agreement profile must carry one row of marginals per rater, which the
children check their chance term against. It is a base, not an outcome, so
it is deliberately absent from `__all__`.

## Dependencies

- `..analysis_result` (AnalysisResult), `..agreement_parts` (AgreementParts)
- `...data_profiles` (BalanceProfile, RaterProfile, AgreementProfile)
- `...not_computed` (NotComputed)
- `numpy` (bootstrap distributions)
