# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning follows [semver](https://semver.org/).

## [Unreleased]

### Breaking

Several of these break code or configuration that ran on 2.0.0. The minor
version does not warn you, so this section does.

- `result.confidence_interval` is a `CIResult`, not a tuple. Use `.lower` and
  `.upper`. `summary()` still emits `ci_lower` and `ci_upper`.
- `bca_ci()` takes a fourth argument, `jackknife_estimates`.
- `percentile_ci()` and `bca_ci()` return a `BootstrapCIResult`, not a tuple.
- `scipy_bca_ci()` is removed. See `analyzers/README.md`.
- `interrater.intermediates.analysis_results` is now
  `interrater.intermediates.results`. Top-level imports are unchanged.
- `intermediates/ci_methods.py` is removed; the module is
  `analyzers/ci_methods.py`.
- Subset permutation intervals resample subjects rated by at least two of the
  k raters, not only those all k rated. Numbers change on partially crossed
  data; balanced data is unaffected. See
  `intermediates/results/permutation_results/README.md`.
- Weights are declared in the semantics, not the analyzer config.
  `weights_type` and `custom_weights_path` are removed from every config; a
  config that still has them, including one nested in a pairwise or subset
  config, is refused with a message saying so. For a
  weighted coefficient, declare categories `kind: ordered` with `weighting`,
  or `kind: weighted` with `weights` or `weights_file`.
- Fleiss' kappa, Brennan-Prediger, Krippendorff's alpha and Cohen's kappa
  now use the weights a semantics declares. In 2.0.0 an ordered or weighted
  semantics was silently ignored, so results on such data change.
- Scott's pi, Gwet's AC1, Conger's kappa and Light's kappa refuse ordered
  or weighted semantics, since they have no weighted form. Use a semantics
  with `kind: unordered` for them.
- `build_from_nij_matrix` is always unweighted.
- `CustomWeightsInput` is removed. Nothing used it.
- `Weights.from_file` is removed, with its JSON, CSV and text formats.
  Custom matrices live in the semantics, inline or in a YAML `weights_file`.
- `unused_categories: drop` is removed. The cleaner no longer changes the
  declared categories; use `keep`, or remove the category from the
  semantics file.
- `RatingsSemantics.from_dict` rejects unknown keys, with a did-you-mean
  hint. A misspelt key such as `unique_subject` used to be ignored.
- `scheme="within_subject"` is removed from hypothesis tests. Shuffling the
  raters within a subject keeps the agreement the test is meant to break, so
  it cannot properly test whether agreement is above chance. A config that
  names it is refused, with a message pointing to `within_rater_shuffle`.
- `scheme="column_shuffle"` is renamed `within_rater_shuffle`, and the
  function `column_shuffle()` is renamed `within_rater_shuffle()`. The new
  name says what it does: each rater's answers are shuffled across subjects,
  not the raters themselves. It is still the default. The old name is
  refused with a message giving the new one.
- Hypothesis tests take `alternative`: `"two-sided"` (default), `"greater"`
  or `"less"`, for both methods. In 2.0.0 Gwet's AC1, Brennan-Prediger and
  Krippendorff's alpha were one-sided and the other five two-sided, so the
  same data could give p = 0.02 under one coefficient and 0.99 under
  another. Their p-values change.
- Permutation p-values are (count + 1)/(B + 1) (Phipson & Smyth 2010), so
  they are never 0. Two-sided is twice the smaller tail, capped at 1,
  instead of comparing absolute values, which assumed the null was centred
  on 0.
- `HypothesisTestResult` requires `coefficient`, the name of the
  coefficient tested (e.g. `"fleiss_kappa"`). Its printed form reads
  `H₀: fleiss_kappa = 0` instead of `H₀: point estimate = 0`.
- Category names, or `value_map` keys mapping to different categories, that
  differ only in case or surrounding spaces are refused. With the default
  case-insensitive matching they were silently merged: every answer for one
  was coded as the other.
- An analyzer built from a nij matrix refuses a config with `injection`.
  Injection builds a new rater from individual ratings, which a nij matrix
  does not have; it used to be skipped silently, so the result was reported
  as if the rater had been added.
- `percent_agreement` is no longer listed as a subset-permutation or
  subject-bootstrap coefficient. It had no analyzer or config, so it could
  not run; the `percent_agreement` function remains.
- `parametric_null` (interval and test) is refused on unbalanced data, and
  for Fleiss on weighted data, where its formulas do not apply. Two
  example-README rows now read "not computed".
- Each injection step is seeded with `[random_state, step_position]`.
  Chained steps with the same seed used to produce identical raters;
  injected raters now differ from 2.0, single steps included.

### Added

- `BootstrapCIResult` and `ParametricCIResult`: bounds plus `method_used`,
  `n_usable`, `degraded` and `summary()`.
- `recommended_resamples` — 250 for percentile, 1000 for BCa — and
  `below_recommended_resamples`, compared against the draws that survived.
- A `"ci"` key in `summary()`, holding the interval's `summary()`.
- A `README.md` in each source subpackage.
- `weights_file` in the semantics: a YAML file holding `names` and `matrix`,
  read from the semantics file's folder. The names must match the categories
  exactly, in the same order.
- `weight_matrix` on Fleiss, Brennan-Prediger, Krippendorff and Cohen
  results and in their `summary()`: the matrix used, or `None` when unweighted.
- `Weights.from_categories()`: the weights a semantics declares.
- Weighted SR1 examples in `example_data/SR1_OSF2018/`: ordered (linear)
  and custom weights, file and inline, with reference values.
- Tests run on GitHub on every push: Python 3.10 with the minimum supported
  dependencies, and 3.11–3.13 with the latest.
- `HypothesisTestResult` records `alternative`, `n_usable`, `n_at_or_above`
  and `n_at_or_below`, and refuses a p-value that does not follow from
  them. Permutations whose coefficient is undefined are skipped and
  counted.
- `CleanerConfig.to_dict()`, in the shape `from_dict` reads.
- `plot_bootstrap_distribution(result)`: the bootstrap draws with the point
  estimate and interval. For BCa it also draws the percentile interval for
  comparison and states z0 and a in words.
- `plot_null_distribution(test_result)`: a permutation test's shuffles,
  with the observed value, the shuffles counted toward p coloured, and the
  p-value's arithmetic.

### Changed

- `ci_method="parametric_null"` warns twice, from the config and from the
  interval, and both say the width is scaled as though the estimate were zero.
- The resample-count notice no longer reads "B below 1000 recommended", which
  parsed as a recommendation to stay below 1000.
- Cohen's and Conger's `balanced_data_only` warnings state what was set.
- Example reference values are recomputed at `random_state=42`, the package
  default, instead of 32.
- The cleaning report states the drop rule it applied, marks rare and
  unused categories on their own row, and reads "no rating: 0" when
  nothing is missing.
- The cleaning report's saved record (`summary()`, `to_yaml`) includes the
  cleaner config and the package version. The rendered report shows rated
  counts per rater, category labels, what each drop pass removed, and plain
  plurals.

### Fixed

- A subset with fewer than two measurable subjects returns
  `method_used = None` rather than a tuple of NaNs.
- `hypothesis_test()` with `method="parametric_null"` raised `TypeError`
  for Fleiss' kappa, Scott's pi and Cohen's kappa.
- Custom weights raised `AttributeError` on every group coefficient.
- Unquoted numeric category names in YAML (`names: [1, 0, -1]`) crashed the
  cleaner. Numbers are now read as text; unquoted yes/no are refused with a
  hint to quote them.
- The README Quickstart was not valid Python, and one link was broken on
  PyPI.
- A permutation test on an undefined coefficient (NaN, e.g. every rating in
  one category) reported p = 0 and rejected the null. It is now refused with
  a message saying why.
- BCa's bias correction counted bootstrap draws equal to the point estimate
  as below it, which shifted intervals up. Ties now count half below and
  half above, as in `scipy.stats.bootstrap`. Ties are common on small
  discrete data: Brennan-Prediger tied the point in up to 12% of draws in
  the examples, and 12 example intervals changed, the largest from
  [0.745, 0.925] to [0.715, 0.910].
- Conger's kappa took chance and observed agreement from different
  subjects. Both now use those rated at least twice: pairwise Conger equals
  pairwise Cohen; full-panel values are unchanged.
- Subset permutation checked `subset_k` against the raters before
  injection, so `"all"` skipped the injected raters and a k counting them
  was refused.
- The cleaning report's "no rating" total left out unknowns treated as
  empty while its per-rater lines counted them. One definition now, and
  the report refuses counts that do not add up.
- The cleaning report called categories with fewer than 10 ratings "not
  estimable"; they are estimable, just imprecise.
- The drop report gave each removed rater and subject its count before any
  dropping, not when it was removed, so a casualty of a cascade could appear
  to meet the minimum.

## [2.0.0] - 2026-10-11

Complete rewrite in object-oriented Python. This changelog begins at 2.1.0;
for what 2.0.0 contained, see
[What version 2 adds](README.md#what-version-2-adds).

## [1.5]

Final release of [Inter-Rater v1](https://github.com/djarenas/Inter-Rater),
described in [arXiv:1809.05731](https://arxiv.org/abs/1809.05731) and archived
at [Zenodo](https://doi.org/10.5281/zenodo.1227660).
