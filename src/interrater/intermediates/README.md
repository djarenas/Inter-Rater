# Intermediates

The internal data structures that sit between validated input (`Ratings`)
and final output (reports, visualizations). Nothing here is user-facing
API -- callers import from `interrater`, not from `interrater.intermediates`.

## Sub-packages

| Sub-package        | What It Holds                                           | README                        |
|--------------------|---------------------------------------------------------|-------------------------------|
| `estimator_parts/` | Weights, NijMatrix, EstimatorInput -- bare arrays and integers for the hot path | [estimator_parts/README.md](estimator_parts/README.md) |
| `data_profiles/`   | BalanceProfile, RaterSharedMatrix -- structural metadata about the data shape   | [data_profiles/README.md](data_profiles/README.md)     |
| `analysis_results/`| Point estimates, CIs, permutation results, hypothesis tests -- everything an analyzer produces | [analysis_results/README.md](analysis_results/README.md) |

## Modules

| Module             | What It Does                                            |
|--------------------|---------------------------------------------------------|
| `ci_methods.py`    | `percentile_ci`, `bca_ci`, `scipy_bca_ci` -- confidence interval computation from bootstrap distributions. |
| `_profile_utils.py`| Prevalence, observed/expected agreement, category arithmetic. Free functions: numpy in, numpy out. Private to this package. |

## Data flow

```
Ratings (validated input)
    │
    ├──► EstimatorInput ──► coefficient estimator ──► PointEstimateResult
    │         │                                              │
    │         ├──► NijMatrix (group coefficients)            ├──► ci_methods (bootstrap CI)
    │         └──► Weights                                   └──► HypothesisTestResult
    │
    ├──► BalanceProfile ──► attached to every result
    │
    └──► _profile_utils ──► prevalence, agreement arithmetic
```

## Dependencies

- `array_utils` (root-level counting and masking)
- `numpy`, `scipy` (ci_methods)
