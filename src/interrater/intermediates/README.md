# Intermediates

The data structures that sit between validated input (`Ratings`) and final
output (reports, visualizations). Every analysis hands some of them back. Some of them
can be useful for reading results with summarized important information, or for feeding
into visualizers. They are rarely built by users. 
(They import what they need from `interrater`, what does this mean).

## Sub-packages

| Sub-package        | What It Holds                                           | README                        |
|--------------------|---------------------------------------------------------|-------------------------------|
| `estimator_parts/` | Weights, NijMatrix, EstimatorInput -- bare arrays and integers for the hot path | [estimator_parts/README.md](estimator_parts/README.md) |
| `data_profiles/`   | BalanceProfile, RaterProfile, RaterSharedMatrix, AgreementProfile, PrevalenceBiasIndex -- what the data looks like, beside the coefficient | [data_profiles/README.md](data_profiles/README.md)     |
| `results/`         | Point estimates, CIs, permutation results, hypothesis tests -- everything an analysis produces | [results/README.md](results/README.md) |

## Modules

| Module             | What It Does                                            |
|--------------------|---------------------------------------------------------|
| `not_computed.py`  | `NotComputed(reason)`: stands in for a value that does not apply here, with the reason. Keeps three states apart: a number (computed), NaN (computed but undefined for this data), and NotComputed (not applicable). It has no truth value; check `isinstance`. Exported as `interrater.NotComputed`. |

Confidence intervals are computed in `analyzers/ci_methods.py`, not here.
The result objects that carry an interval live in `results/`; the arithmetic
that produces one belongs with the analyzers that call it.

## Data flow

```
Ratings (validated input)          Weights (from the semantics)
    │                                  │
    ├──► EstimatorInput ◄──────────────┘
    │         │
    │         ├──► NijMatrix (group coefficients)
    │         ▼
    │    estimator ──► point ─────────────────────┐
    │         │                                   │
    │    resampling ──► analyzers/ci_methods ──► CIResult
    │                                             ▼
    ├──► BalanceProfile, RaterProfile,      PointEstimateResult
    │    AgreementProfile ───────────────────────► (carries all of them)
    │
    └──► analyzer.hypothesis_test() ──► HypothesisTestResult (separate)
```

A result built from counts (`build_from_nij_matrix`) has no rater identity,
so its `rater_profile` is `NotComputed`.

## Dependencies

- `array_utils` (root-level counting and masking)
- `numpy`
