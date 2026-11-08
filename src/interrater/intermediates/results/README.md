# Results

Everything an analysis produces. Each result is a frozen dataclass holding
what is needed to interpret, serialize and plot an outcome without going back
to the raw data.

## Layout

| Path                    | What It Holds                                                      | README                                                       |
|-------------------------|--------------------------------------------------------------------|--------------------------------------------------------------|
| `analysis_result.py`    | `AnalysisResult`, the root ABC both result trees descend from       | below                                                        |
| `analyzer_results/`     | One coefficient, one estimate: the `PointEstimateResult` family     | [analyzer_results/README.md](analyzer_results/README.md)     |
| `permutation_results/`  | One coefficient recomputed over rater subsets: the `PermutationResult` family | [permutation_results/README.md](permutation_results/README.md) |
| `ci_results/`           | The intervals: the `CIResult` family                                | [ci_results/README.md](ci_results/README.md)                 |
| `agreement_parts.py`    | `AgreementParts`, the two terms a coefficient is built from        | below                                                        |
| `hypothesis_test_result.py` | `HypothesisTestResult`, on no tree                              | below                                                        |

```
AnalysisResult (ABC)            analysis_result.py
├── PointEstimateResult         analyzer_results/
└── PermutationResult           permutation_results/

CIResult (ABC)                  ci_results/
├── BootstrapCIResult
└── ParametricCIResult

HypothesisTestResult            hypothesis_test_result.py, on neither tree
```

`AnalysisResult` is at this level because it belongs to neither subfolder
alone. Each subfolder owns its own base (`PointEstimateResult`,
`PermutationResult`, `CIResult`) and documents its own children.

`CIResult` is a separate tree: an interval has no point estimate and no
coefficient name. A `PointEstimateResult` holds one.

**Draws and their interval must agree.** A result that stores bootstrap
draws also holds the interval made from them, and is refused if the two
disagree: draws come exactly with a `BootstrapCIResult`, one per resample it
counted, with as many non-finite draws as it counted failed. The check is in
`_bootstrap_draws.py`, shared by every result that holds draws.

## The classes at this level

| Class                  | What It Holds                                                    | Why It Exists                                                      |
|------------------------|------------------------------------------------------------------|---------------------------------------------------------------------|
| `AnalysisResult`       | ABC with `summary() -> dict`                                      | Every result can describe itself in plain types.                    |
| `AgreementParts`       | observed, chance, kind ('agreement' for Po/Pe, 'disagreement' for Do/De), weighted | Rebuilds the coefficient, so a result can check its point against them. |
| `HypothesisTestResult` | p_value, alpha, reject_null, method, point, alternative; z (parametric) or null distribution and counts (permutation) | Test against H0: coefficient = 0.                                   |

**`HypothesisTestResult` is standalone.** It has no point estimate or
interval of its own, so it does not descend from `AnalysisResult`.

**Its p-value can be recomputed.** A permutation test records `n_usable`,
`n_at_or_above` and `n_at_or_below`, and is refused if its p-value does not
follow from them. Each tail is (count + 1)/(n_usable + 1), counting the
observed data as one arrangement, so p is never 0 (Phipson & Smyth 2010);
two-sided is twice the smaller tail, capped at 1. Shuffles whose coefficient
is NaN are not counted in `n_usable`. A NaN observed coefficient is refused.

## Across the trees

- **Frozen dataclasses.** A result describes one dataset at one moment and
  does not change after construction.
- **Error prefixes.** Each class has its own `_ERROR_PREFIX`. A child's
  `__post_init__` wraps `super().__post_init__()`, so the child's prefix
  appears even when the parent's check fails.
- **Self-describing.** What is needed to read the number travels with it:
  `ci_method`, `weights_type` and `weight_matrix`, `is_exact` (Scott's pi),
  the balance and rater profiles. A caller never has to ask the config.
- **Imports.** Some callers import from the module file rather than from
  these `__init__`s, so moving a file means updating those imports.

## Dependencies

- `data_profiles` (BalanceProfile, RaterProfile)
- `numpy` (bootstrap distributions)
