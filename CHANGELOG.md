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
- `scipy_bca_ci()` is removed; `bca_ci()` is the one implementation, checked
  against Efron & Tibshirani (1986) in `tests/test_analyzers/test_ci_methods.py`.
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
- `rater_jackknife()` and `rater_subset()` are removed from
  `resampling_functions`. Nothing in the package called them: subset
  permutation builds its combinations in `agreement_functions` and takes
  columns with `EstimatorInput.select_raters()`.
- `Weights.build_from_type()` is removed. It read the `weights_type` key,
  which 2.1.0 also removes; use `Weights.build_from_categories()`, or
  `Weights.linear()` / `Weights.quadratic()` directly.
- `PairedRatings` is removed. No analyzer took it: Cohen and Scott check
  for two raters themselves, and `Ratings` already requires every subject
  rated at least twice, so two raters can have no gaps. Use `Ratings`.
- Fleiss' kappa, Brennan-Prediger, Krippendorff's alpha and Cohen's kappa
  now use the weights a semantics declares. In 2.0.0 an ordered or weighted
  semantics was silently ignored, so results on such data change.
- Scott's pi, Gwet's AC1, Conger's kappa and Light's kappa refuse ordered
  or weighted semantics, since they have no weighted form. Use a semantics
  with `kind: unordered` for them.
- Conger's and Light's kappas require complete data, every rater rating
  every subject, as Conger (1980) and Light (1971) define them. Data with
  missing ratings is refused by their analyzers, and Conger's by pairwise
  permutation, with a message naming the coefficients that handle missing
  ratings. 2.0.0 ran on such data: Conger took chance and observed agreement
  from different subjects, and Light averaged pairs over different subjects.
- Subset permutation applies the coefficient's `balanced_data_only` to each
  subset. Fleiss' kappa defaults to `true`, so on a panel whose subsets are
  unbalanced it is refused unless `config_of_coefficient` sets
  `balanced_data_only: false`; 2.0.0 ignored the setting and ran the
  approximation. Pairs are always balanced on their shared subjects.
- Cohen's, Conger's and Light's configs refuse `balanced_data_only: false`,
  saying why it does not apply. Cohen and Conger used to set it to `true`
  with a warning; Light used to accept it. Cohen is always computed on the
  subjects both raters rated; Conger and Light need complete data. `true` is
  accepted, and is now Light's default, so a config moved from Fleiss'
  kappa still loads.
- `Ratings` refuses subjects with fewer than two ratings, and
  `build_from_nij_matrix` refuses count rows that sum to less than two.
  2.0.0 refused only subjects rated by nobody; a subject rated once joins no
  comparison either.
- The cleaner always drops until nothing more falls below the minimums,
  as `method: iterate` did. `method` is removed from `insufficient_data_drops`
  and `InsufficientDataPolicy`, and a config that still has it is refused
  with a message saying so. The default was `raters_first`; the one-pass
  methods could return raters and subjects below the declared minimums,
  including subjects rated once. What iterating leaves is the largest data
  meeting both minimums, the same whichever axis is swept first.
  `MAX_DROP_PASSES` and `DropReport.method` and `.converged` are removed.
- `rater_jackknife: true` is refused by every coefficient config: leave-one-rater-out
  is not supported in this version yet. 2.0.0's group configs (Fleiss, Scott,
  Gwet, Brennan-Prediger, Krippendorff) accepted it and silently did nothing.
- `build_from_nij_matrix` is always unweighted.
- `CustomWeightsInput` is removed. Nothing used it.
- `Weights.from_file` is removed, with its JSON, CSV and text formats.
  Custom matrices live in the semantics, inline or in a YAML `weights_file`.
- `unused_categories: drop` is removed. The cleaner no longer changes the
  declared categories; use `keep`, or remove the category from the
  semantics file.
- `RatingsSemantics.build_from_dict` rejects unknown keys, with a did-you-mean
  hint. A misspelt key such as `unique_subject` used to be ignored.
- `scheme="within_subject"` is removed from hypothesis tests. Shuffling the
  raters within a subject keeps the agreement the test is meant to break, so
  it cannot test whether agreement is above chance: Fleiss-type coefficients
  never change under it, and for Cohen-type coefficients only the chance
  term moves, so it detects raters using categories at different rates,
  not agreement. A config that names it is refused, with a message pointing
  to `within_rater_shuffle`.
- `scheme="column_shuffle"` is renamed `within_rater_shuffle`, and the
  function `column_shuffle()` is renamed `within_rater_shuffle()`. The new
  name says what it does: each rater's answers are shuffled across subjects,
  not the raters themselves. It is still the default. The old names are
  refused with a message giving the new one: the scheme in a config, and the
  function when called.
- Hypothesis tests take `alternative`: `"two-sided"` (default), `"greater"`
  or `"less"`, for both methods. In 2.0.0 Gwet's AC1, Brennan-Prediger and
  Krippendorff's alpha were one-sided and the other five two-sided, so the
  same data could give p = 0.02 under one coefficient and 0.99 under
  another. Their p-values change.
- Permutation p-values are (count + 1)/(B + 1) (Phipson & Smyth 2010), so
  they are never 0. Two-sided is twice the smaller tail, capped at 1,
  instead of comparing absolute values, which assumed the null was centred
  on 0. A shuffle whose value differs from the observed one only by rounding
  (within 100 machine epsilons of it) counts as a tie, at or beyond it in
  both tails, as in `scipy.stats.permutation_test`.
- `HypothesisTestResult` requires `coefficient`, the name of the
  coefficient tested (e.g. `"fleiss_kappa"`). Its printed form reads
  `H₀: Fleiss' κ = 0` instead of `H₀: point estimate = 0`.
- Category names, or `value_map` keys mapping to different categories, that
  differ only in case or surrounding spaces are refused. With the default
  case-insensitive matching they were silently merged: every answer for one
  was coded as the other.
- An analyzer built from a nij matrix refuses a config with `injection`.
  Injection builds a new rater from individual ratings, which a nij matrix
  does not have; it used to be skipped silently, so the result was reported
  as if the rater had been added.
- `percent_agreement` is no longer an active subset-permutation or
  subject-bootstrap coefficient. It had no analyzer or config, so it could
  not run; it is now a commented "Future" entry in both lists, and the
  `percent_agreement` function remains.
- `parametric_null` (interval and test) is refused on unbalanced data, and
  for Fleiss on weighted data, where its formulas do not apply. Two
  example-README rows now read "not computed".
- Each injection step is seeded with `[random_state, step_position]`.
  Chained steps with the same seed used to produce identical raters;
  injected raters now differ from 2.0, single steps included.
- Injector configs default to `random_state=42`, as every other config
  does; it was `None`, so an analysis with injection was unseeded unless
  each step set a seed. Pass `random_state=None` for fresh randomness.
- `summary()` replaces `bootstrap_n`, `bootstrap_mean` and `bootstrap_std`
  with `bootstrap_draws_mean` and `bootstrap_se`, over the usable draws; the
  counts are in the interval's summary. Results refuse draws that disagree
  with their interval.
- Per-pair and per-subset intervals are required on permutation results.
- The cleaner refuses rater columns read as decimals (pandas does this to a
  numeric column with a blank cell), with a hint to read the file with
  `dtype=str`. They used to fail as unknown values, with advice that would
  have emptied the column.
- `PermutationBootstrapMeanConfig` takes only `method`. The mean's draw count,
  `ci_method`, `alpha` and seed come from `config_of_coefficient`, as the
  per-combination intervals' do. `n_bootstraps`, `ci_method`, `alpha` and
  `random_state` are refused with where each moved, from a config or from
  Python; any other key is refused. Example values are unchanged.
- The independent mean requires `config_of_coefficient.ci_method:
  percentile`. BCa's acceleration needs a jackknife over the units the
  bootstrap resamples, and none matches the independent mean; it is refused
  when the config is built.
- `DropReport` takes two required fields, `min_ratings_per_rater` and
  `min_ratings_per_subject`: the drop rule it applied. Code
  that builds one by hand must pass them; the cleaner does.
- Pairwise and subset results keep every pair's or subset's bootstrap draws.
  2.0.0 kept them only up to 20 pairs or 120 subsets, after holding them
  through the whole run; `MAX_PAIRS_FOR_DISTRIBUTIONS` and
  `MAX_SUBSETS_FOR_DISTRIBUTIONS` are removed. `keep_draws: false` (below)
  leaves them out.
- The independent mean is built while the per-pair or per-subset bootstrap
  runs, from one value picked from each unit's draws in unit order. The
  method is unchanged, but the picks come in a different order, so
  independent intervals differ slightly from 2.0.0 at the same seed.
- Each pair or subset is its coefficient's own result (a pair under Cohen's
  kappa is a `CohenKappaResult`, with its interval, draws, parts and
  profiles), or `NotComputed` with the reason. `pair_results` and
  `subset_results` replace `pairwise_values`, `subjects_per_pair`,
  `bootstrap_cis_per_pair` and `bootstrap_distributions_per_pair` (and the
  `subset_*` and `*_per_subset` fields). The container checks its units
  against its own settings and its mean.
- One rule for which pairs and subsets have a coefficient: at least two
  usable subjects (rated by at least two of its raters), and every rater
  rating at least one of them. A subset with an idle rater used to report
  the smaller subset's number; it is now `NotComputed`, and a bootstrap
  draw or jackknife value that leaves a rater idle counts as failed.
- Scott's pi is refused for subset permutation; it is a two-rater
  coefficient. Use pairs.
- `hypothesis_test()` takes no `method` argument; the method comes from
  `hypothesis_test_config`. The argument was declared but raised
  `TypeError` when passed.
- Constructors are verbs: `from_*` is now `build_from_*` on
  `BalanceProfile`, `AgreementProfile`, `PrevalenceBiasIndex`, `Weights`,
  `RatingsSemantics`, `CleanerConfig`, `FrameLayout`,
  `InsufficientDataPolicy` and `PermutationBootstrapMeanConfig` (14 in
  all). The old names are removed.
- `BalanceProfile` describes subjects only. `n_raters`, `n_raters_exact`,
  `n_per_rater`, `fill_rate` and `n_raters_display` are removed; the raters
  are in the new `RaterProfile`, as `result.rater_profile`. In `summary()`
  those keys move from `"balance"` to a new `"raters"` section.
- Results carry `rater_profile`, a `RaterProfile` built from ratings or
  `NotComputed` on the counts path. `PointEstimateResult` requires a
  `balance_profile`; `None` is refused. Code that builds results by hand
  must pass both.
- Pairwise and subset results carry `weights_type` and `weight_matrix`,
  checked against every pair or subset, and both are in `summary()`. Code
  that builds these results by hand must pass them.
- `plot_pairwise_permutation`, `plot_subset_sorted` and
  `plot_subset_strip` take no `rater_names` argument; names come from the
  result.
- `Ratings` is read-only: setting or deleting an attribute raises. Build a
  new `Ratings` to change the data.
- `DropReport` requires `passes` and refuses a report whose last pass
  dropped something.
- Printed coefficient names come from the result classes and use Greek
  symbols, e.g. `Fleiss' κ`, `Scott's π`, `Gwet's AC1`.

### Added

- `keep_draws` on pairwise and subset permutation configs, default `true`.
  `false` tosses each pair's or subset's draws once its interval is built
  and its share of the mean added, so memory holds one unit's draws at a
  time. Values and intervals are the same either way.
- `max_evaluations` on pairwise and subset permutation configs, default
  `GUARD_MAX_EVALUATIONS` (2,000,000, in `contracts`). The analysis counts
  the coefficient evaluations a run would make and refuses, before anything
  runs, a run above it. Raise it to run a larger analysis on purpose.
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
- `Weights.build_from_categories()`: the weights a semantics declares.
- Weighted SR1 examples in `example_data/SR1_OSF2018/`: ordered (linear)
  and custom weights, file and inline, with reference values.
- Tests run on GitHub on every push: Python 3.10 with the minimum supported
  dependencies, and 3.11–3.13 with the latest.
- `HypothesisTestResult` records `alternative`, `n_usable`, `n_at_or_above`
  and `n_at_or_below`, and refuses a p-value that does not follow from
  them. Permutations whose coefficient is undefined are skipped and
  counted.
- `CleanerConfig.to_dict()`, in the shape `build_from_dict` reads.
- `plot_bootstrap_distribution(result)`: the bootstrap draws with the point
  estimate and interval. For BCa it also draws the percentile interval for
  comparison and states z0 and a in words.
- `plot_null_distribution(test_result)`: a permutation test's shuffles,
  with the observed value, the shuffles counted toward p coloured, and the
  p-value's arithmetic.
- `balanced_data_only` on pairwise and subset permutation configs, default
  `false`: whether the whole panel must be balanced. `true` refuses an
  unbalanced panel before anything runs.
- `RaterProfile`: rater names, how many subjects each rated, and their
  overlap (`rater_overlap`, a `RaterSharedMatrix`), with `n_raters` and
  `fill_rate` derived. `RaterSharedMatrix.build_from_ratings()`.
- `RaterProfilePlotter` (`draw_overlap()`, `draw_counts()`) and
  `BalanceProfilePlotter` (`draw()`), each built from ratings or a result
  (`BalanceProfilePlotter` also from a nij matrix). See
  `visualizers/README.md`.
- `to_frame()` on pairwise and subset results: one row per pair or subset
  (raters, estimate, interval, `point_outside_interval`, method used,
  failed draws, subjects, reason if not computed), readable back from CSV.
- `COEFFICIENT`, `NAME` and `SYMBOL` on every result class, and
  `result_class_for(coefficient)`.
- Tests build their data through the cleaner (`tests/fixtures/through_cleaner.py`),
  and `test_scruffy.py` checks the cleaner against a damaged file's record
  of its damage.
- `example_data/Statistical_background/`: every table on
  `guide/statistical_background.md` as ratings, `generate_tables.py` to
  build them (seed 42 for the simulated ones), and `reproduce_page.py` to
  run interrater on them and check the page's numbers. A test runs both.

### Changed

- `ci_method="parametric_null"` warns twice, from the config and from the
  interval, and both say the width is scaled as though the estimate were zero.
- The resample-count notice no longer reads "B below 1000 recommended", which
  parsed as a recommendation to stay below 1000.
- Example reference values are recomputed at `random_state=42`, the package
  default, instead of 32.
- The cleaning report states the drop rule it applied, marks rare and
  unused categories on their own row, and reads "no rating: 0" when
  nothing is missing.
- The cleaning report's saved record (`summary()`, `to_yaml`) includes the
  cleaner config and the package version. The rendered report shows rated
  counts per rater, category labels, what each drop pass removed, and plain
  plurals.
- Pairwise and subset plots note in the legend that an independent mean's
  interval ignores the correlation between combinations sharing raters, so
  tends to be too narrow.
- A point estimate outside its own interval no longer raises a warning; with
  every pair a result, a run printed one per pair. The result says it itself:
  `point_outside_interval`, the same key in `summary()`, and
  `POINT OUTSIDE INTERVAL` in the printed form. Pairwise and subset tables
  show it on the pair's row, and `to_frame()` has a `point_outside_interval`
  column.
- Pairwise and subset tables, `to_frame()`, summaries and plots name raters
  by the names in the semantics instead of their positions.
- Pairwise and subset printouts state the weighting (`quadratic weights`,
  `unweighted`) for coefficients that have one.
- On unbalanced data, Gwet's AC1 and Brennan-Prediger print `unbalanced
  data: Gwet (2014) formula for missing ratings` instead of `approximation
  (unbalanced data)`: they compute Gwet's formula, as irrCAC does. Fleiss'
  kappa prints `approximation (unbalanced data: each subject weighted
  equally)`.
- `example_data/Fleiss1971/Fleiss1971_psychiatric_data_nij.txt` is
  removed; it duplicated the CSV, which is what the tests load.
- `guide/statistical_background.md`: the more-than-two-raters table shows
  theoretical and observed columns, the incomplete-data section ends with
  Krippendorff's alpha, and a closing section shows how to reproduce the
  page. `guide/statistical_background_numbers.py` is removed;
  `reproduce_page.py` replaces it.

### Fixed

- Analyzer, permutation and injector configs refuse unknown keys in
  `build_from_dict`, naming the class, the nearest match and the valid
  keys; before, only `build_from_yaml` checked, and a dict, a nested
  `config_of_coefficient` or `hypothesis_test`, or an injection step failed
  with Python's bare "unexpected keyword argument".
  `hydrate_injection_step` no longer removes `type` from the caller's dict.
- Semantics declaring `n_subjects` no longer fail when the cleaner drops
  subjects. `n_subjects` is checked against the file before any drop; the
  cleaned semantics now clears it whenever a subject was dropped, as it
  already did when a rater or category was. The report no longer counts
  that as a change of semantics.
- Subset permutation checks every k before running any: `subset_k: "all"`
  on 20 raters used to run k = 2, 3 and 4 in full before k = 5 hit the
  10,000-subset cap. The cap is now `GUARD_MAX_SUBSETS` in `contracts`
  (it was `MAX_SUBSETS` in `rater_combinations`).
- A subset with fewer than two measurable subjects holds `NotComputed`
  with the reason, rather than a tuple of NaNs.
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
- The draws' mean and standard error in `summary()` (now `bootstrap_draws_mean`
  and `bootstrap_se`) were NaN if any draw was undefined. They are computed
  over the usable draws, and the standard error divides by n − 1.
- Pairwise and subset results said "bootstrap was not performed" when draws
  were not stored; they now say why (`draws not kept: keep_draws false`).
- A scale with more than 127 categories is refused when declared. It used
  to crash the cleaner with a numpy `OverflowError`.
- The mean's interval was not reproducible by default (its seed defaulted to
  `None`), and a result could report a different `ci_method` from the one its
  mean interval used.
- A permutation config with `ci_method: parametric_null` was accepted and
  failed once the analysis had started; it is refused when the config is
  built.
- A bootstrap draw of the mean over pairs or subsets that lost one of them
  (undefined in that resample) averaged the rest, a different quantity
  counted as an ordinary draw. It is now counted as failed, as for a single
  coefficient; the jackknife for BCa follows the same rule. Example values
  are unchanged.
- BCa crashed with inverted bounds when `1 - a(z0 + z)` went below zero
  (the point outside nearly all the draws, with a large acceleration). It
  now falls back to percentile, as it did when the quantile was infinite.
- `Ratings` checked category codes after casting to int8, so an
  out-of-range code could wrap into a valid one. Codes are checked first.
- `light_kappa` emitted a numpy "Mean of empty slice" warning when no pair
  was defined; it returns NaN quietly.
- `NijMatrix` refused unsigned integer arrays. Any integer type is
  accepted, and it keeps its own copy.

## [2.0.0] - 2026-10-11

Complete rewrite in object-oriented Python. This changelog begins at 2.1.0;
for what 2.0.0 contained, see
[What version 2 adds](https://djarenas.github.io/packages/interrater/roadmap/).

## [1.5]

Final release of [Inter-Rater v1](https://github.com/djarenas/Inter-Rater),
described in [arXiv:1809.05731](https://arxiv.org/abs/1809.05731) and archived
at [Zenodo](https://doi.org/10.5281/zenodo.1227660).
