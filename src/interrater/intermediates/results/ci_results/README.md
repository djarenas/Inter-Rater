# CI Results

The intervals themselves. Every class here descends from `CIResult`, which
holds the two bounds and the level they were taken at. The subclasses differ
in what can go wrong on the way to them.

```
CIResult (ABC)             lower, upper, alpha
├── BootstrapCIResult      from resampling: counts, methods, diagnostics
└── ParametricCIResult     from a variance formula: no counts, no fallback
```

| Class                | What It Holds                                                                 | Why It Exists                                                             |
|----------------------|-------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| `CIResult`           | lower, upper, alpha; `is_defined`, `confidence_level`, `as_tuple()`           | What every interval has, whatever produced it.                            |
| `BootstrapCIResult`  | n_resamples, n_failed, method_intended, method_used, acceleration, bias_correction; `n_usable`, `degraded`, `recommended_resamples`, `below_recommended_resamples` | A resampled interval, plus an account of what failed, what actually ran, and how much it had to run on. |
| `ParametricCIResult` | variance_source                                                               | A closed-form interval. Nothing to fall back to, so nothing to disclose about a choice never made. |

## These are not analysis results

`CIResult` does not descend from `AnalysisResult`, and that is deliberate. An
interval holds no point estimate, no coefficient name and no opinion about the
data; a `PointEstimateResult` holds *it*. Keeping the two trees separate means
`summary()` on a result and `summary()` on its interval are answering
different questions rather than the same one at two depths.

It also means an interval never checks itself against a point estimate. A
strongly bias-corrected BCa interval can sit entirely to one side of the
estimate -- correct behaviour, not a contradiction -- so whether the two are
consistent is a question for whoever holds both. `PointEstimateResult` asks
it, and warns.

## `__str__` is abstract, not just `summary()`

Every subclass must describe itself completely, including anything that
qualifies how its bounds should be read. That is what lets a holder print any
interval without inspecting fields or branching on type:

```
[0.481, 0.565] (95% CI, bca, n=1000)
[0.481, 0.565] (95% CI, percentile [bca requested], n=1000)
[0.481, 0.565] (95% CI, bca, n=1200, 12 failed)
[0.339, 0.551] (95% CI, bca, n=50, 50 usable, 1000 recommended)
[undefined] (95% CI, no method ran [bca requested], n=200, all 200 failed)
[0.401, 0.639] (95% CI, null variance -- width assumes kappa = 0, reproduction only)
```

Those are copied from real output, not written by hand. The third line is
`n=1200` rather than `n=1000` for a reason worth knowing before reading the
next section: a BCa run of exactly 1000 draws with even one failure has 999
usable, which is below the minimum, so it would carry the
`999 usable, 1000 recommended` clause too.

That clause is two facts set side by side, with no verb between them. An
earlier wording, `B below 1000 recommended`, read as a recommendation to stay
*under* 1000 -- and it was false besides, since B was 1000 and what fell short
was the usable count.

Leaving `__str__` concrete on the base would let a subclass silently inherit a
description that omits its own caveats, which is exactly what
`ParametricCIResult` must never do.

Each line stays on one line, because whoever holds an interval prints it
inside their own single-line description. What survives that compression is
the fact that decides whether a reader should use the number at all.

## What each subclass has to say

**`BootstrapCIResult` reports failures without judging them.** A resample
fails when its coefficient is undefined -- one containing a single category
drives chance agreement to 1, and `(p_o - p_e) / (1 - p_e)` has nothing to
say. Those draws arrive as NaN, are excluded from the bounds, and are counted.
There is no threshold, because how much the failures matter depends on whether
failing correlates with the coefficient's value, which no count can tell you.

It also records `method_intended` against `method_used`, so an interval that
asked for BCa and got a percentile fallback says so. That is a real fallback,
not a relabelling: when `method_used` is `percentile`, the bias correction was
discarded too and the bounds are plain quantiles. `method_used` is `None` when
no draw was usable and nothing ran at all -- an absence rather than a
substitution.

**It also says when it rests on too few draws.** Efron & Tibshirani (1986)
§9 p. 72 give rough minimums: B = 1000 for BCa, and about 250 for a useful
percentile interval, the difference being that the percentile method need not
compute z0. Those are `MIN_RESAMPLES_BCA` and `MIN_RESAMPLES_PERCENTILE`, and
an interval below the minimum for the method that actually ran says so when it
prints.

Three decisions sit inside that sentence, and all three are easy to get
backwards:

- **It compares against `n_usable`, not `n_resamples`.** A run of 1000 draws
  with 900 failures leaves 100, and 100 is the number the quantiles were taken
  over. Checking the requested count would report adequacy the interval does
  not have.

- **It is measured against the method that ran, not the one requested.** A BCa
  request that fell back to percentile is held to the percentile minimum,
  because percentile is what produced the bounds. `recommended_resamples` is
  `None` when nothing ran -- there is no recommendation for a thing that did
  not happen, and the failure count already says everything.

- **Nothing is enforced.** An interval from 900 draws is not invalid, and an
  interval from 50 is noisy rather than wrong. There is no defensible cutoff,
  so the class reports the count and leaves the judgement with whoever holds
  it. Construction validation stays structural: it rejects only intervals no
  correct computation could have produced.

**`ParametricCIResult` carries a standing caveat.** Its half-width comes from
the *null* variance: the spread of the coefficient under H0, on the assumption
the true value is zero. That is the right quantity for a hypothesis test,
which asks how far from zero an estimate is in units of its spread at zero,
and a mismatched one for a confidence interval, which wants the spread at the
value actually observed. So it is centred on the estimate but scaled as though
the estimate were zero.

**How far off is that?** Measured against percentile bootstrap intervals on
this package's own datasets, the width runs from about **56% too narrow**
(Fleiss' kappa on Fleiss 1971 and on SR1 prelim) to **29% too wide** (Cohen's
kappa on Gwet cont3x3). The direction depends on the coefficient, so there is
no correction factor to apply and no safe direction to assume -- which is why
the class discloses rather than adjusts.

It is kept because that is what several published formulas and several other
packages report, which makes it the right object for reproducing a published
number or comparing implementations -- and the wrong one for reporting the
precision of an estimate. A bootstrap interval is what that calls for.

The warning is raised **last** in `__post_init__`, after the parent has
validated the bounds. An interval that cannot exist should not first be
explained; it should raise.

## Why one warns and the other does not

`ParametricCIResult` warns at construction, and prints its caveat every time.
`BootstrapCIResult` stays silent about a low usable count, disclosing it only
in `__str__` and `summary()`. Stated side by side that looks inconsistent, and
it is not. The deciding question is how many of each object one command
produces.

**There is exactly one `ParametricCIResult` per analysis**, and it only exists
because somebody asked for a non-default `ci_method`. A warning there lands
once, on a choice that was made deliberately, and it is the only thing
standing between a mismatched width and a published table.

**One analysis can build many `BootstrapCIResult`s.** A pairwise run makes one
per rater pair -- C(n,2), so fifteen for six raters -- and a multi-k subset run
makes one per subset. If each warned about its resample count, a single
`estimate()` call would emit fifteen identical warnings. What that teaches is
not "raise B"; it teaches `warnings.simplefilter("ignore")`, and then the next
warning worth reading is gone too. A warning's value is spent on repetition,
so it is spent where it lands once.

The same logic puts a second warning in the config rather than merging it into
one. `AnalyzerConfig` defines `_warn_parametric_null()` but **never calls it**;
each config class that genuinely accepts `parametric_null` calls it at the end
of its own `__post_init__`. The two warnings say different things at different
moments: the config one says *where the knob is*, while the choice is still
cheap to reverse, and the interval one says *how wrong the number is*, at the
point where it would be copied into a table. Someone who gets only the second
has to work backwards to find what to change.

## Where the minimums live

`MIN_RESAMPLES_PERCENTILE` and `MIN_RESAMPLES_BCA` are defined in
`bootstrap_ci_result.py` and re-exported by `analyzers/ci_methods.py`, which
is the direction that looks backwards and is not. The interval is what reports
against the minimums, so the interval owns them; the CI functions re-export
them so that a caller working at the function level does not have to reach
into a results sub-package to find a number the functions' own documentation
mentions.

## Dependencies

- `numpy` (finiteness checks, numeric coercion)
- `scipy` is *not* needed here; the arithmetic lives in `analyzers/ci_methods.py`

## References

- Efron, B. & Tibshirani, R. 1986. "Bootstrap Methods for Standard Errors,
  Confidence Intervals, and Other Measures of Statistical Accuracy."
  *Statistical Science* 1(1): 54--77. §9 p. 72 for the resample minimums;
  eq. 8.5 for the acceleration the bootstrap interval reports.
- Efron, B. & Tibshirani, R. J. 1993. *An Introduction to the Bootstrap.*
  Chapman & Hall.
