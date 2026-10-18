# Results

Everything an analysis produces. Each result is a frozen dataclass holding
enough metadata to interpret, serialize, and visualize an outcome without
going back to the raw data.

Four kinds of thing live here, and the folder layout says which is which.

## Layout

| Path                    | What It Holds                                                      | README                                                       |
|-------------------------|--------------------------------------------------------------------|--------------------------------------------------------------|
| `analysis_result.py`    | `AnalysisResult` -- the root ABC both result trees descend from     | below                                                        |
| `analyzer_results/`     | One coefficient, one estimate: the `PointEstimateResult` family     | [analyzer_results/README.md](analyzer_results/README.md)     |
| `permutation_results/`  | One coefficient recomputed over rater subsets: the `PermutationResult` family | [permutation_results/README.md](permutation_results/README.md) |
| `ci_results/`           | The intervals themselves: the `CIResult` family                     | [ci_results/README.md](ci_results/README.md)                 |
| `hypothesis_test_result.py` | `HypothesisTestResult` -- standalone, on no tree                | below                                                        |

```
AnalysisResult (ABC)            results/analysis_result.py
├── PointEstimateResult         analyzer_results/
└── PermutationResult           permutation_results/

CIResult (ABC)                  ci_results/
├── BootstrapCIResult
└── ParametricCIResult

HypothesisTestResult            results/ -- on neither tree
```

`AnalysisResult` sits at this level rather than inside either subfolder
because it belongs to neither alone. Each subfolder owns its own base --
`PointEstimateResult`, `PermutationResult`, `CIResult` -- and documents its
own children.

**`CIResult` is a separate tree on purpose.** An interval is not an analysis
result: it has no point estimate, no coefficient name and no opinion about
the data. A `PointEstimateResult` *holds* one. Putting intervals beside the
results rather than under them keeps that relationship visible, and keeps
`summary()`'s contract meaning one thing.

## The two classes at this level

| Class                  | What It Holds                                                    | Why It Exists                                                      |
|------------------------|------------------------------------------------------------------|---------------------------------------------------------------------|
| `AnalysisResult`       | ABC with `summary() -> dict`                                      | Common contract: every result can describe itself in plain types.   |
| `HypothesisTestResult` | p_value, alpha, reject_null, method, point, optional z/null dist  | Test against H0: coefficient = 0. Parametric or permutation.        |

**`HypothesisTestResult` is standalone.** A hypothesis test is a
fundamentally different kind of output from a point estimate with a
confidence interval -- there is no point estimate of its own for it to
inherit a contract for. It does not inherit from `AnalysisResult`, and it
is kept at this level rather than filed under a tree it is not part of.

## Design decisions spanning the trees

- **Frozen dataclasses.** A result is a fact about a dataset at a moment
  in time. It should not change after construction.

- **Error chaining across inheritance.** Each class has its own
  `_ERROR_PREFIX`. Child `__post_init__` wraps `super().__post_init__()`
  in `try/except` so the child prefix appears in the message even when
  the parent validation fails.

- **A result is self-describing.** Anything needed to interpret the number
  travels with it -- `ci_method`, `weights_type`, `is_exact`, the balance
  profile. A caller holding a result never has to ask the config what it
  asked for.

- **Imports reach into modules, not packages.** Callers elsewhere in the
  package import from the module file rather than from these `__init__`s,
  which is why moving a file here means updating every caller. Deliberate
  for now: it avoids circular imports between analyzers and results.

## Dependencies

- `data_profiles` (BalanceProfile)
- `numpy` (bootstrap distributions)
