## What version 2 adds

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

- **Confidence intervals that estimate what they claim to.** Bootstrap
  intervals from resampled subjects, with the null test kept separate and
  separately named. 

- **Intervals that account for shared subjects.** The mean over rater pairs is
  bootstrapped by resampling subjects once and recomputing every pair, because
  the pairs are not independent of each other. Treating them as independent
  narrows the interval.

- **A declaration of the task, apart from the data.** Raters, categories,
  expected subject count and the measurement scale are declared in a semantics
  file before the data is read, and a `Ratings` is validated against that
  declaration. Weights for ordered or custom scales are declared there too,
  not in the analysis, because how far apart two categories are is a fact
  about the task. See [Why objects](why_objects.md).

- **An auditable cleaning layer.** Which strings meant nothing, which were
  repaired and into what, how thin the data was allowed to get: declared in a
  versioned config file and reported afterwards, rather than reconstructed from
  memory when a reviewer asks. If cleaning changed what the data describes
  (a rater dropped, a category removed), the report says so and hands back a
  corrected declaration.

- **Preconditions checked, not assumed.** Exact Fleiss' kappa requires a
  constant number of raters per subject; Conger's and Light's kappas require
  every rater to rate every subject, as their definitions say. Those are
  properties of the data, so they are checked before the coefficient is
  computed (by a data object or by the analyzer) rather than assumed by it. See [Why objects](why_objects.md).

- **Results that cannot contradict themselves.** A bootstrap result is refused
  if its interval does not follow from its draws; a permutation test is
  refused if its p-value does not follow from its counts. Every result records
  the method requested, the method that ran, how many draws were usable, the
  weight matrix if any, and the package version.

- **Support for unbalanced designs.** Where each subject is rated by a subset
  of the rater pool rather than all of them, which is the usual case outside a
  designed reliability study.

- **Pairwise and subset permutation, with visualization.** Compute a
  coefficient for every pair of raters, or every subset of a given size,
  bootstrap each unit and the mean, and plot the agreement structure: which
  pairs agree, which do not, and how much of the spread one rater accounts
  for. Seven coefficients are available pairwise (Cohen's kappa, Scott's pi,
  Conger's kappa, Gwet's AC1, Fleiss' kappa, Brennan-Prediger, Krippendorff's
  alpha), so the same data can be viewed under different chance models side
  by side. See [Why every pair](why_pairwise.md).

- **Synthetic raters, in eight declared kinds.** Add a reviewer who does not
  exist (a conformist, a contrarian, an echo of the panel with tunable noise, a
  guesser drawing from prevalences you declare) and see what your coefficient
  would have been. The original data is never modified and injected columns are
  named as such. See
  [What if one of your raters were someone else?](rater_substitution.md).

- **Hypothesis tests with a declared alternative.** Permutation tests take
  `two-sided`, `greater` or `less`, the same for every coefficient, and report
  p-values as (count + 1)/(B + 1), so a p-value is never zero.

- **Everything declarable in YAML.** Semantics, cleaning, every coefficient
  config, and the pairwise and subset configs with their nested injections
  load from files, so an analysis is recoverable from what it declared.
  Plotting is the exception, below.

The changes between 2.0 and 2.1, several of which alter numbers on the same
data, are listed in the
[changelog](https://github.com/djarenas/Inter-Rater/blob/master/CHANGELOG.md).

## What is coming next

Version 2 is the pipeline: cleaning, validated data objects, eight
coefficients, bootstrap intervals, permutation tests, pairwise and subset
analyses, simulated raters. What follows is mostly about driving that pipeline
with less ceremony and widening what it can compute. This is a statement of
direction, not a schedule.

- **Declarative configuration for plots.** Plotting ships as standalone
  functions: pass a result, get a figure. That is fine once and awkward for a
  figure in a paper, because the choices that shape it (sort order, which
  intervals are drawn, labels) live in the call rather than in a file.
  Plotting will move to the config-plus-plotter pattern already used for
  cleaning and analysis, and used throughout
  [eventus](https://github.com/djarenas/eventus), so that a figure is
  reproducible from a declared file. The standalone functions stay, with the
  config objects built on top.

- **Gwet's AC2 and extended weighting schemes.** AC2 is AC1 with weights, so
  it waits on extending the weight set beyond identity, linear and quadratic
  to ordinal, interval and ratio, which also widens Krippendorff's alpha to
  the variants its literature defines. Further coefficients follow once that
  foundation is in place.

- **Subset permutation, properly visualised.** The current plots show subsets
  sorted by value and as a strip. Neither yet shows what the analysis is for:
  how agreement degrades as raters are removed, and whose absence accounts
  for it. A panel across subset sizes is planned.

- **Orchestrators, for computing every coefficient at once.** Eight
  coefficients today means eight analyzers and eight configs, with the
  bookkeeping left to you. An orchestrator takes one dataset and one
  resampling budget and returns the set, computed on the same resampled
  subjects so the numbers are comparable.
    - *A warning to my future self.* An analyzer that returns eight
      coefficients makes it easy to report the largest one, and the gap
      between Fleiss' kappa and AC1 on skewed data is information only if
      both numbers are shown. So the orchestrator will ask for a primary
      coefficient in the config before anything is computed and record it in
      the result; it will return the set in a fixed order, with no ranking and
      no sort by value; and it will label each coefficient with the chance
      model it assumes, so a reader can see why two of them differ.

- **Testing a difference.** Whether one pair of raters agrees more than
  another within a pairwise result, where the two share subjects and the test
  must be paired; and whether two estimates from two different sets of
  ratings differ, where they do not. In both cases the tempting shortcut,
  reading non-overlapping 95% intervals as a significant difference, is
  closer to a test at p < 0.005 than at p < 0.05 (Schenker & Gentleman,
  2001), so the difference will be bootstrapped directly.

- **The floor of the standard error.** Inject raters with the panel's own
  error rate, one at a time, and bootstrap at each panel size. The interval
  narrows as raters are added and stops narrowing at a floor set by the
  number of subjects, which no panel size can get below. Reading that floor
  off the curve tells a study whether to recruit raters or subjects. This
  needs an injector calibrated to the panel's observed disagreement, which
  the current eight are not.

- **Rater influence by replacement.** For each rater: drop them, inject a
  synthetic replacement (a contrarian, a draw from a declared distribution, a
  resample of the rest of the panel) and re-run. One coefficient per "what if
  this rater were someone else" scenario, which distinguishes a rater who was
  already noise from one who was holding the panel together. The injectors
  ship now; the loop around them does not.

- **Fewer lines, without losing the account of what was decided.** Building
  an analysis today means naming every config explicitly. That is deliberate
  where the choice matters, and more ceremony than it should be where the
  choice is obvious: a first coefficient on clean balanced data does not need
  four objects. A later version will explore defaults for config values, under
  two conditions: anything assumed rather than declared says so in a warning,
  and the result records which defaults it used, so a figure remains
  reproducible from what the run reports. Where no default can be chosen
  safely, it stays an error rather than a guess. This reduces the ceremony
  around the pipeline, not the pipeline itself: there is still no path that
  takes a bare array and hands back a coefficient, because that is the
  interface where the guarantees are lost.

None of this is scheduled for 2.1.x. Version 2.1.0 is large enough already,
and already does several things no other agreement package does; the next
release can be allowed to be small.

Also planned, less immediately: diagnostics for agreement by category and rater
confusion, convergence checking cheap enough to leave on by default,
leave-one-rater-out as a diagnostic (the `rater_jackknife` flag is refused
today rather than silently ignored, as it was in 2.0), pairwise variants of
the individual-marginal coefficients, a clustered bootstrap for subjects that
are not independent draws, and generating ratings from a declared true
coefficient, so that an estimator can be checked against an answer known in
advance, an interval checked for covering it, and a review powered before it
starts.

One thing is deliberately not planned: rater bootstrap will not be offered as
a confidence interval method, because a rater drawn twice agrees perfectly
with themselves and the interval narrows around agreement that was never
observed. It remains available as an injector, where duplicating a rater is
the declared intent.

Suggestions and bug reports are welcome at
[the issue tracker](https://github.com/djarenas/Inter-Rater/issues). If you
are waiting on something above, say so there.
