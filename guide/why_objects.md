## Why objects

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

Objects are routine in software engineering and uncommon in statistics
packages, where a coefficient is usually a function that takes an array and
returns a number. `interrater` is built from objects instead. This page lists
what that buys.

### Containers hold data; semantics give it meaning

- NumPy arrays and pandas DataFrames are awesome because they are so
  flexible, and the same flexibility means they cannot carry meaning. A column
  of integers does not know whether it is a rater, a subject, or a score; a 3
  does not know whether it is a category, a count, or a third rater. That is
  what makes them good containers, and it is why they cannot be the thing an
  analysis validates.
- `RatingsSemantics` is the object that gives a container meaning. It says
  which raters there are, which categories they could choose from, how many
  subjects there should be, and whether a subject may appear twice. With it, a
  `Ratings` is a record of a rating task rather than a table of small
  integers.
- The idea has names in software engineering (value objects, typed wrappers)
  and partial forms in data tools (categorical dtypes, named dimensions), but
  it is uncommon in statistics code, where the function takes an array and the
  meaning stays in the analyst's head. Here the object records what the raters
  were told, which is the same fact the statistical background page says the
  chance term depends on.
- Meaning is where the facts the data cannot show are kept. A category that no
  rater used still exists, because the declaration says so, and a coefficient
  whose chance term depends on the number of categories (Bennett's S,
  Brennan-Prediger) sees the scale the raters had rather than the subset they
  used. How far apart two categories are, nominal, ordinal with linear or
  quadratic distances, or a custom matrix, is declared with the categories,
  and the analyzer configs declare no weights at all.
- It is written before the data is seen and loads from YAML, so a reviewer can
  be shown what the raters were asked to do.

### Objects that cannot be wrong

- A `Ratings` is validated once, when it is built. Nothing downstream re-checks
  that a value is in the declared categories, that no rater is idle, or that
  every subject was rated at least twice.
- Construction is also where an error can be specific. The object knows which
  argument was wrong and what the caller asked for, so the message names the
  class, the field, and the alternatives. Compare this to an error inside a
  loop: it can describe the bad array but not where it went bad.
- The same holds for configuration. A YAML file can be wrong in many ways; a
  `CleanerConfig` or a `FleissKappaConfig` either exists and is sound or was
  never built.
- Validation runs before the loop, so the formulas inside ten thousand
  bootstrap draws do no checking at all.

### Inheritance, so each rule is written once

- Checks shared by several configurations live on parent classes. A leaf adds
  only its own and calls `super()`; no function has to remember to call
  another.
- Error messages inherit the same way. Each level prefixes its class name as
  the error passes through it, so the message is the chain from the object the
  user built down to the level whose rule was broken. Two errors on the same
  `FleissKappaConfig`:

  ```
  FleissKappaConfig(n_resamples=-5)
  ValueError: [FleissKappaConfig] Error: [WeightableGroupConfig] Error: [GroupCoefficientConfig] Error: [AnalyzerConfig] Error: n_resamples must be >= 1, got -5

  FleissKappaConfig(balanced_data_only=1)
  TypeError: [FleissKappaConfig] Error: [WeightableGroupConfig] Error: [GroupCoefficientConfig] Error: balanced_data_only must be bool, got int
  ```

  The first chain reaches `AnalyzerConfig`, so the rule is one every analyzer
  shares. The second stops at `GroupCoefficientConfig`, where
  `balanced_data_only` is defined.
- The estimate loop (point, resample, interval, package) is written once on
  `Analyzer`. A coefficient supplies one formula, and the bootstrap calls the
  same function as the point estimate, so the two cannot drift apart.

### Parallel hierarchies, so the pieces fit

- The analyzer tree and the config tree have the same shape, so an analyzer can
  only receive the config built for it.
- The branches of both trees are based on statistical concepts rather than on
  code organization. Coefficients that pool rater marginals read a counts
  matrix; coefficients that use individual marginals read who said what. There
  is no need to write functions checking that the Cohen's kappa analyzer takes
  the right input; inheriting from the right parent class already takes care
  of it.

### Results that remember how they were made

- An agreement statistic is more than a point estimate. It matters how many
  subjects and raters there were, how complete the ratings were, whether some
  categories were picked far more than others, and how often raters agreed on
  each category. The result objects in `interrater` carry these alongside the
  number.
- The result objects also carry the confidence level, the method requested,
  the method that actually ran, and how many bootstrap draws were usable. A
  confidence interval is more than a low and a high.
- Results also check themselves. A bootstrap result stores both the resampled
  values and the interval built from them, and it is refused at construction
  if the interval is not the one those values give. A permutation test stores
  the number of shuffles that matched or beat the observed value, and its
  p-value is recomputed from that count at construction; if the two disagree
  the object is refused. A result that contradicts itself cannot be built.
- Some runs fail: a choice that is incompatible with a statistic, a bootstrap
  draw in which only one category appears, a pair of raters with too few
  shared subjects. Instead of a zero or a silent NaN, the result records what
  did not run and why, as a `NotComputed` object that carries the reason.
- Weighted and unweighted kappa are different statistics, so a result carries
  its weight matrix, and every result carries the package version that produced
  it.

### Frozen data, so operations compose

- A `Ratings` never changes. Anything that would change it returns a new one.
- This allows an injected rater, a subset of raters, and any analyzer to chain
  together without a function written for that combination.

### Where this comes from

Version 1 of this package was 745 lines of functions. It read two text files,
printed a report, drew a figure, and was used in published systematic reviews.
It also only had five error messages, each a `print` followed by `return False`; a
pair of raters that failed its check kept the `0.0` in its cell and entered the
per-rater mean; one chance term, the pooled Σ p̄_j², served every coefficient,
so the two-rater result the code called Cohen's kappa was Scott's pi; and the
standard error of a rater's mean kappa treated pairs sharing that rater as
independent. Every function computed what it was asked to. Nothing recorded
what an array was or what had already been decided about it, so each function
assumed, and the assumptions were not written anywhere a reader could check.
Although it took me years to get back to it, the rebuilt version 2 is, I
believe, worth it.
