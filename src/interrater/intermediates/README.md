# Intermediates

The internal data structures that sit between validated input (`Ratings`)
and final output (reports, visualizations). Nothing here is user-facing
API -- callers import from `interrater`, not from `interrater.intermediates`.

## Sub-packages

| Sub-package        | What It Holds                                           | README                        |
|--------------------|---------------------------------------------------------|-------------------------------|
| `estimator_parts/` | Weights, NijMatrix, EstimatorInput -- bare arrays and integers for the hot path | [estimator_parts/README.md](estimator_parts/README.md) |
| `data_profiles/`   | BalanceProfile, RaterSharedMatrix -- structural metadata about the data shape   | [data_profiles/README.md](data_profiles/README.md)     |
| `results/`         | Point estimates, CIs, permutation results, hypothesis tests -- everything an analysis produces | [results/README.md](results/README.md) |

## Modules

| Module             | What It Does                                            |
|--------------------|---------------------------------------------------------|
| `_profile_utils.py`| Prevalence, observed/expected agreement, category arithmetic. Free functions: numpy in, numpy out. Private to this package. |
| `not_computed.py`  | `NotComputed(reason)`: stands in for a value that does not apply here, with the reason. Keeps three states apart: a number (computed), NaN (computed but undefined for this data), and NotComputed (not applicable). It has no truth value; check `isinstance`. Exported as `interrater.NotComputed`. |

Confidence intervals are computed in `analyzers/ci_methods.py`, not here.
The result objects that carry an interval live in `results/`; the arithmetic
that produces one belongs with the analyzers that call it.

## Data flow

```
Ratings (validated input)
    │
    ├──► EstimatorInput ──► coefficient estimator ──► PointEstimateResult
    │         │                                              │
    │         ├──► NijMatrix (group coefficients)            ├──► analyzers/ci_methods (bootstrap CI)
    │         └──► Weights                                   └──► HypothesisTestResult
    │
    ├──► BalanceProfile ──► attached to every result
    │
    └──► _profile_utils ──► prevalence, agreement arithmetic
```

## Dependencies

- `array_utils` (root-level counting and masking)
- `numpy`
