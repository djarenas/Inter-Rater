# CI Results

CI intervals with the information that helps them be more interpretable.

Every class here descends from `CIResult`, which
holds the two bounds and the level they were taken at. The subclasses differ
in what can go wrong on the way to them.

```
CIResult (ABC)                      lower, upper, alpha, methods
├── BootstrapCIResult (ABC)         from resampling: draws and counts
│   ├── PercentileBootstrapCIResult bounds = quantiles of the draws
│   ├── BCABootstrapCIResult        point, jackknife, z0, acceleration
│   └── EmptyBootstrapCIResult      every draw failed: no bounds
└── ParametricCIResult              from a variance formula: no counts, no fallback
```

| Class                | What It Holds                                                                 | Why It Exists                                                             |
|----------------------|-------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| `CIResult`           | lower, upper, alpha, method_intended, method_used; `is_defined`, `degraded`, `confidence_level`, `as_tuple()` | What every interval has, whatever produced it.                            |
| `BootstrapCIResult`  | n_resamples, n_failed, draws; `n_usable`, `recommended_resamples`, `below_recommended_resamples`, `compute_draw_moments()` | A resampled interval, with its draws, what failed, and how many draws it ran on. Checks the draws against its counts. |
| `PercentileBootstrapCIResult` | `compute_bounds()` | Percentile ran (asked for, or BCa fell back). Checks its bounds are the quantiles of its draws. |
| `BCABootstrapCIResult` | point, jackknife_estimates, acceleration, bias_correction; `compute_acceleration()`, `compute_bias_correction()`, `compute_bounds()` | BCa ran. Checks its acceleration, z0 and bounds follow from its jackknife, point and draws. |
| `EmptyBootstrapCIResult` | -- | Every draw failed, so no method ran and there are no bounds. |
| `ParametricCIResult` | variance_source, variance_reference                                           | A closed-form interval. Nothing to fall back to.                          |

## Not panel results

`CIResult` does not descend from `PanelResult`. An interval holds no
coefficient name; a `PointEstimateResult` holds the interval. Only a BCa
interval holds a point, the one its z0 was measured against. No interval
checks that the point lies inside it: a strongly bias-corrected BCa interval
can sit entirely to one side of the estimate, correctly. `PointEstimateResult` reports that case as
`point_outside_interval`.

## How an interval prints

`__str__` is abstract, so every subclass describes itself completely,
including anything that qualifies its bounds. Real output examples:

```
[0.481, 0.565] (95% CI, bca, n=1000)
[0.481, 0.565] (95% CI, percentile [bca requested], n=1000)
[0.481, 0.565] (95% CI, bca, n=1200, 12 failed)
[0.339, 0.551] (95% CI, bca, n=50, 50 usable, 1000 recommended)
[undefined] (95% CI, no method ran [bca requested], n=200, all 200 failed)
[0.401, 0.639] (95% CI, null variance -- width assumes kappa = 0, reproduction only)
```

Each stays on one line, because the holder prints it inside its own line.
The third example uses n=1200: a BCa run of exactly 1000 with one failure
has 999 usable draws, below the minimum, and would also print
`999 usable, 1000 recommended`.

## BootstrapCIResult

- **Failed draws.** A draw fails when its coefficient is undefined (a
  resample with a single category drives Pe to 1). Failed draws are left out
  of the bounds and counted in `n_failed`. There is no threshold: whether
  failures matter depends on whether failing correlates with the
  coefficient's value, which a count cannot show.
- **What ran.** `method_intended` and `method_used` are both kept, so a BCa
  request that fell back says `percentile`. The fallback is complete: the
  bias correction is dropped too, and the bounds are plain quantiles.
  `method_used` is `None` when no draw was usable and nothing ran.
- **Too few draws.** Efron & Tibshirani (1986, §9, p. 72) give rough
  minimums:

  | method | rough minimum B |
  |--------|-----------------|
  | BCa or BC | 1000 |
  | percentile | 250 |
  | standard errors only | 25 |

  BCa needs more because it also estimates z0. An interval below the minimum
  for the method that ran prints `N usable, M recommended`. The comparison is
  against `n_usable`, not `n_resamples` (1000 draws with 900 failures is 100
  draws), and against the method that ran, not the one requested.
  `recommended_resamples` is `None` when nothing ran. Nothing is enforced: an
  interval from 50 draws is noisy, not invalid.

## ParametricCIResult

Its half-width comes from the null variance, the spread of the coefficient
when the true value is zero. That is right for a hypothesis test and wrong
for an interval around a non-zero estimate, so the interval is centred on the
estimate but scaled as if the estimate were zero. It is kept for reproducing
published numbers and comparing with other packages.

It warns once at construction (there is one per analysis, and only when
someone chose `parametric_null`) and prints its caveat every time. The config
warns earlier, about the setting itself. The warning is raised after the
bounds are validated, so an interval that cannot exist raises instead of
being explained.

`BootstrapCIResult` and the analyzers do not use warnings for few draws. This is to prevent
warning burnout. For example, a pairwise run builds one per pair, and fifteen
identical warnings teach people to switch warnings off. A note about too few
draws is instead in the result object printed line and `summary()`.

## Checked against known answers

**BCa against Efron & Tibshirani (1986), Example 3, p. 71.** The law school
data: 15 schools, correlation between LSAT score and GPA. Reproduced in
`tests/test_analyzers/test_ci_methods.py` (`TestEfronLawSchool`), with the
data and published values in `tests/fixtures/efron1986_lawschool_*`.

| | point | `a` | `z0` | central 90% interval |
|---|---|---|---|---|
| published (B = 100,000) | .776 | -.0817 | -.0927 | [.43, .92] |
| this package (B = 10,000) | 0.776 | -0.0757 | -0.0773 | [0.44, 0.93] |

- The paper prints its fifteen empirical influence values. Fed directly into
  this package's formula for `a` (their eq. 8.5), they give -0.0817, the
  published value. The remaining difference in `a` comes from its input:
  Efron uses the infinitesimal jackknife (their eq. 10.11), this package the
  delete-one jackknife, as in Efron & Tibshirani (1993).
- The published interval is 90%, so `alpha=0.10`; a 95% interval will look
  wrong against it.
- The paper used B = 100,000, "about 100 times more than actually
  necessary"; at 10,000 the bounds move by about .01 and z0 by about .015.

**The null-variance width against the bootstrap,** on the shipped datasets
(percentile, B = 1000, `random_state=42`):

| dataset / coefficient | bootstrap width | null-variance width | |
|---|---|---|---|
| Fleiss 1971 / Fleiss κ | 0.2186 | 0.0955 | -56% |
| SR1 prelim / Fleiss κ | 0.2298 | 0.0942 | -59% |
| Gwet cont3x3 / Fleiss κ | 0.2347 | 0.2986 | +27% |
| Gwet cont3x3 / Cohen κ | 0.2332 | 0.2979 | +28% |
| Gwet cont3x3 / Cohen κ, quadratic weights | 0.1482 | 0.3916 | +164% |
| Synthetic / Scott π | 0.2627 | 0.2948 | +12% |

The direction depends on the data, not only on the coefficient: Fleiss' κ
is too narrow on Fleiss 1971 and SR1 prelim and too wide on Gwet cont3x3.
So there is no correction factor and no safe direction to assume.

## Where the minimums live

`MIN_RESAMPLES_PERCENTILE` and `MIN_RESAMPLES_BCA` are defined in
`bootstrap_ci_result.py`, since the interval reports against them, and
re-exported by `analyzers/ci_methods.py` so a caller there does not have to
reach into this package.

## Dependencies

- `numpy` (finiteness checks, numeric coercion). No `scipy`; the arithmetic
  is in `analyzers/ci_methods.py`.

## References

- Efron, B. & Tibshirani, R. 1986. "Bootstrap Methods for Standard Errors,
  Confidence Intervals, and Other Measures of Statistical Accuracy."
  *Statistical Science* 1(1): 54-77. §9 p. 72 for the resample minimums;
  eq. 8.5 for the acceleration; Example 3 p. 71 for the law school check.
- Efron, B. & Tibshirani, R. J. 1993. *An Introduction to the Bootstrap.*
  Chapman & Hall.
