## What version 2 adds

- **Confidence intervals that estimate what they claim to.** Bootstrap
  intervals from resampled subjects, with the null test kept separate and
  separately named. See above for why.

- **An auditable cleaning layer.** Which strings meant nothing, which were
  repaired and into what, how thin the data was allowed to get: declared in a
  versioned config file and reported afterwards, rather than reconstructed from
  memory when a reviewer asks. If cleaning changed what the data describes
  (a rater dropped, a category removed), the report says so and hands back a
  corrected declaration.

- **Preconditions in the type system.** Exact Fleiss' kappa requires a constant
  number of raters per subject; that is a property of the data, so it is
  checked before the coefficient is computed rather than assumed by it.

- **Synthetic raters, in eight declared kinds.** Add a reviewer who does not
  exist (a conformist, a contrarian, an echo of the panel with tunable noise, a
  guesser drawing from prevalences you declare) and see what your coefficient
  would have been. The original data is never modified and injected columns are
  named as such. This does not appear to be available in other agreement
  packages, and it is the clearest way to ask how much a reported kappa depends
  on the particular people who produced it. See above.

- **Intervals that account for shared subjects.** The mean over rater pairs is
  bootstrapped by resampling subjects once and recomputing every pair, because
  the pairs are not independent of each other. Treating them as independent
  understates the interval by 23 to 47% on real screening data.

- **Support for unbalanced designs.** Where each subject is rated by a subset
  of the rater pool rather than all of them, which is what most rating tasks
  can actually afford.

- **Pairwise permutation with visualization.** Compute any agreement
  coefficient for every pair of raters, bootstrap each pair and the mean, and
  plot the full agreement structure: which pairs agree, which do not, and
  how much of the spread one rater accounts for. The plot supports seven
  coefficients (Cohen's kappa, Scott's pi, Conger's kappa, Gwet's AC1,
  Fleiss' kappa, Brennan-Prediger, Krippendorff's alpha), so the same data
  can be viewed under different chance models side by side.

## What is coming next

Version 2.0 is the pipeline: cleaning, validated data objects, eight
coefficients, bootstrap intervals, permutation tests, pairwise and subset
analyses, simulated raters. What follows is mostly about driving that pipeline
with less ceremony, and widening what it can compute. This is a statement of
direction, not a schedule.

- **Declarative configuration for plots.** Version 2.0 ships plotting as
  standalone functions: pass a result, get a figure. That is fine once and
  awkward for a figure in a paper, because the choices that shape it (sort
  order, which intervals are drawn, labels) live in the call rather than in a
  file. Plotting moves to the config-plus-plotter pattern already used here for
  cleaning and analysis, and used throughout
  [eventus](https://github.com/djarenas/eventus), so that a figure is
  reproducible from a declared file. The standalone functions stay, with the
  config objects built on top.

- **One YAML story for every config.** A prerequisite for the above, and worth
  having on its own. Cleaning, semantics and the top-level analyzer config load
  from YAML today; the leaf configs do not. When they all do, an entire
  analysis (cleaning, coefficient, resampling, plotting) is declarable in
  files and recoverable from them.

- **Gwet's AC2 and extended weighting schemes.** AC2 is AC1 with weights, so it
  waits on extending the weight set beyond identity, linear and quadratic to
  ordinal, interval and ratio, which also widens Krippendorff's alpha to the
  variants its literature defines. Loading a custom weight matrix from file
  arrives with it. Further coefficients follow once that foundation is in
  place.

- **Subset permutation, properly visualised.** Version 2.0 plots subsets sorted
  by value and as a strip. Neither yet shows what the analysis is for: how
  agreement degrades as raters are removed, and whose absence accounts for it.
  A panel across subset sizes is planned, along with a smaller fix that is
  overdue: results carry rater indices, so plots read `(0, 1)` where they
  should read `(Smith, Jones)`.

- **Orchestrators, for computing every coefficient at once.** Eight
  coefficients today means eight analyzers and eight configs, with the
  bookkeeping left to you. An orchestrator takes one dataset and one resampling
  budget and returns the set, computed on the same resampled subjects so the
  numbers are comparable.

  A function that hands back eight coefficients is also a function that invites
  reporting the largest one, and the position of this package is that the
  comparison across chance models is informative only if the set is reported as
  a set: the gap between Fleiss' kappa and AC1 on skewed data is a fact about
  the data, not a menu. So the design is meant to make that the easy path: a
  primary coefficient declared in the config before anything is computed and
  recorded in the result, output in a fixed order with no ranking and no
  sort-by-value convenience, and each coefficient labelled with the chance
  model it assumes so a reader can see why two numbers differ. Alongside it, a
  utility for bootstrapping the difference between two estimates directly,
  because the tempting shortcut, reading non-overlapping 95% intervals as a
  significant difference, is closer to a p < 0.005 test than a p < 0.05 one
  (Schenker & Gentleman, 2001).

- **Rater influence by replacement.** For each rater: drop them, inject a
  synthetic replacement (a contrarian, a draw from a declared distribution, a
  resample of the rest of the panel) and re-run. One coefficient per "what if
  this rater were someone else" scenario, which distinguishes a rater who was
  already noise from one who was holding the panel together. The injectors ship
  in 2.0; the loop around them does not.

- **Fewer lines, without losing the account of what was decided.** Building an
  analysis today means naming every config explicitly. That is deliberate where
  the choice matters, and it is more ceremony than it should be where the choice
  is obvious: a first coefficient on clean balanced data does not need four
  objects. 2.1 will explore defaults for config values, under two conditions:
  anything assumed rather than declared says so in a warning, and the result
  records which defaults it used, so a figure remains reproducible from what the
  run reports. Where no default can be chosen safely, it stays an error rather
  than a guess. Note that this reduces the ceremony around the pipeline, not the
  pipeline itself: there is still no path that takes a bare array and hands back
  a coefficient, because that is the interface where the guarantees are lost.
  More to come.

Also planned, less immediately: diagnostics for agreement by category and rater
confusion, a visualiser for the balance profile, convergence checking cheap
enough to leave on by default, pairwise variants of the individual-marginal
coefficients, a clustered bootstrap for subjects that are not independent
draws, and generating ratings from a declared true coefficient, so that an
estimator can be checked against an answer known in advance, an interval can be
checked for actually covering it, and a review can be powered before it starts
rather than regretted after.

One thing is deliberately **not** planned: rater bootstrap will not be offered
as a confidence interval method, for the reason given above. It remains
available as an *injector*, where duplicating a rater is the declared intent
rather than an accident of the procedure.

Suggestions and bug reports are welcome at
[the issue tracker](https://github.com/djarenas/Inter-Rater/issues). If you are
waiting on something above, saying so is the most useful way to move it up.
