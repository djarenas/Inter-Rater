# Tests

1482 tests across eight packages, run with `pytest tests/`.

Every count in this file is a **collected** count -- what `pytest` prints, not
what you get by counting `def test_` in the source. The two differ because
`@pytest.mark.parametrize` runs one function once per input and pytest counts
each run separately; right now 1441 functions expand to 1482 collected tests.
Collected is the number used throughout because it is the number you can check
in one command, so it cannot drift unnoticed:

```
pytest tests/ --co -q | tail -1
```

The test suite is organised bottom-up: lower layers are tested first so that
any failure points at the layer that broke, not at something that depends on it.
`test_example_data/` comes last because it is the only package that tests
everything at once.


## test_semantics/ (63 tests)

Tests for the semantic metadata layer: categories and ratings semantics.
Built bottom-up along the dependency chain:

    Categories (ABC) -> UnorderedCategories / OrderedCategories / WeightedCategories
    RatingsSemantics

| File | What it tests | Tests |
|------|---------------|-------|
| test_categories.py | CategoriesBaseValidation (min names, empty, duplicates, unknown label), UnorderedCategories (construction, n_categories, index_of, label_for, to_dict round-trip, equality), OrderedCategories (quadratic/linear, invalid weighting, weight matrix properties, to_dict, equality), WeightedCategories (construction, wrong shape, not symmetric, diagonal not one, out of range, to_dict, equality), build_categories factory (dispatches unordered/ordered/weighted) | 36 |
| test_ratings_semantics.py | Construction (from names, from count, properties, n_subjects, unique_subjects), type guards (string/bool raters, non-string names, non-categories, bool n_subjects, non-bool unique_subjects), validation (min raters, duplicates, blanks, n_subjects floor), public API (index_of_rater, to_dict/from_dict round-trip, from_dict missing keys), equality | 27 |
| conftest.py | Shared fixtures: UnorderedCategories (K=2), RatingsSemantics instances | -- |


## test_data_objects/ (46 tests)

Tests for the validated data-object chain that sits between raw data and
the agreement functions. Built bottom-up along the dependency chain:

    Ratings -> BalancedRatings -> PairedRatings

Each layer is tested for construction, type guards, validation rejects, public
API, and equality. The two frozen datasets in conftest (SIMPLE and GAPPY) are
small enough that every property can be verified by counting on fingers.

| File | What it tests | Tests |
|------|---------------|-------|
| test_ratings.py | Construction (subject IDs, properties, dtype cast, frozen arrays), validation (dtype, shape, semantics mismatch, value range), subject ID validation (length, duplicates, unique flag), presence checks (idle rater, unrated subject), check methods (balance, crossing, pairing), to_nij_matrix, equality | 33 |
| test_balanced_ratings.py | Balanced validation (accepted, rejected), public API (raters_per_subject, build_from_ratings, preserves subject IDs), error prefix propagation | 8 |
| test_paired_ratings.py | Paired validation (accepted, three raters rejected, gaps rejected), error prefix (own checks, inherited checks) | 5 |
| conftest.py | Shared fixtures: SIMPLE (3x2, fully crossed) and GAPPY (4x3, with empties) datasets, category objects, semantics, Ratings instances | -- |


## test_cleaners/ (92 tests)

Tests for the configuration objects and the cleaning pipeline that turns a
messy DataFrame into a validated Ratings with a CleaningReport receipt.

Built bottom-up: FrameLayout and InsufficientDataPolicy first (they are
delegates of CleanerConfig), then CleanerConfig, then the full
RatingsCleaner pipeline.

| File | What it tests | Tests |
|------|---------------|-------|
| test_cleaners.py | FrameLayout (defaults, validation, from_dict, equality), InsufficientDataPolicy (defaults, validation, bool trap, below-floor, invalid method, from_dict, equality), CleanerConfig (defaults, validation, type guards, empty key, double-claim, from_dict flat/nested, YAML, equality, repr), RatingsCleaner (type guards, happy path, config vs semantics, column matching, encoding, unused categories, drops, report) | 92 |
| conftest.py | Shared fixtures: UnorderedCategories (K=2), RatingsSemantics (ann, bob), CLEAN_DF (4x2 no problems), MESSY_DF (4x2 with null/missing/typo), default CleanerConfig | -- |


## test_intermediates/ (461 tests)

Tests for the intermediates layer: estimator parts, data profiles, utility
functions, and the full analysis-result hierarchy. One test file per
sub-package, not one per class, so each file covers multiple classes in
bottom-up order.

Built bottom-up along the dependency chain:

    NijMatrix -> Weights -> EstimatorInput
    BalanceProfile -> RaterSharedMatrix
    _profile_utils (prevalence, prevalence_per_rater, coverage_per_rater)
    CIResult -> BootstrapCIResult / ParametricCIResult
    PointEstimateResult -> 8 coefficient results
    PermutationResult -> PairwisePermutationResult -> SubsetPermutationResult -> MultiKSubsetResult
    HypothesisTestResult

**The `ci_results` sub-package has no test file of its own in this package.**
`CIResult` and its two subclasses are exercised where they are produced:
`BootstrapCIResult` in `test_analyzers/test_ci_methods.py`, which is the only
thing that builds one, and `ParametricCIResult` in
`test_analyzers/test_analyzer_classes.py`, which owns its two warnings. They
appear here only as inputs, built by the `_ci()` helper at the top of
`test_analysis_results.py`.

| File | What it tests | Tests |
|------|---------------|-------|
| test_estimator_parts.py | NijMatrix (construction, validation, properties, from_ratings, equality, error prefix), Weights (identity/linear/quadratic/custom, validation, equality, error prefix), EstimatorInput (construction, validation, properties, subsetting, equality, error prefix) | 104 |
| test_data_profiles.py | BalanceProfile (from_nij_matrix, n_subjects/n_categories/n_raters, is_balanced, rater_counts, summary, equality, error prefix), RaterSharedMatrix (from_nij_matrix, shape, values, symmetry, self-counts, equality, error prefix) | 59 |
| test_profile_utils.py | prevalence (shape, values, sums-to-1), prevalence_per_rater (shape, values, row-sums), coverage_per_rater (shape, values, range) | 11 |
| test_analysis_results.py | PointEstimateResult (construction, validation, point_outside_interval, __str__, summary), HypothesisTestResult (construction, validation, __str__, summary), PermutationResult (construction, inheritance, validation, summary, error prefix), PairwisePermutationResult (construction, validation, accessor, __str__, summary), SubsetPermutationResult (construction, validation, accessor, __str__, summary), MultiKSubsetResult (construction, validation, container protocol, __str__, summary), FleissKappaResult, ScottsPiResult, GwetsAC1Result, BrennanPredigerResult, KrippendorffAlphaResult (balance required, ci_method sets, is_exact logic), CohenKappaResult (weights_type validation), CongerKappaResult (always unweighted), LightKappaResult (n_pairs = C(n,2)), IndividualMarginalResult, error prefixes for all | 287 |
| conftest.py | Shared fixtures: BALANCED_NIJ (4x3x3, all rows sum to 3), UNBALANCED_NIJ (4x2, variable row sums), SIMPLE_INDEX (2x3, for index-based tests), SIMPLE_SHARED (3x3 symmetric, for shared-matrix tests) | -- |


## test_injectors/ (141 tests)

Tests for the injection pipeline: configs and injectors that add synthetic
raters to a Ratings object for sensitivity analysis.

Built bottom-up: InjectorConfig and its subtypes first, then the injector
classes that use them.

| File | What it tests | Tests |
|------|---------------|-------|
| test_injector_configs.py | InjectorConfig (defaults, validation, build_from_dict, build_from_yaml, summary), DistributionInjectorConfig (prevalences validation, uniform/custom modes), WithinSubjectInjectorConfig (copy_prob validation, noise_mode validation) | 80 |
| test_injectors.py | RaterInjector (input validation, name building, inject output shape/columns/subject IDs), DistributionInjector (config type, prevalences mismatch, uniform/custom generation, reproducibility), WithinSubjectInjector (config type, pure copy, unanimous rows, noise fallback, reproducibility) | 61 |
| conftest.py | Shared fixtures: simple Ratings, category objects, injector configs | -- |


## test_reports/ (77 tests)

Tests for the reporting layer: structured reports that describe what happened
during cleaning and analysis.

| File | What it tests | Tests |
|------|---------------|-------|
| test_reports.py | DropPass (construction, dropped_nothing logic, frozen), DropReport (defaults, properties, summary, render variants), MatchingReport (n_missing, category counts, render), CleaningReport (construction, shape tracking, survived, semantics_changed, render) | 77 |
| conftest.py | Shared fixtures: DropPass instances, MatchingReport instances | -- |


## test_analyzers/ (569 tests)

Tests for the analyzer layer: the classes that wire agreement functions,
CI methods, bootstrap loops, and hypothesis tests into a single `estimate()`
call. Built bottom-up along the inheritance tree:

    Analyzer (ABC)
    ├── GroupAnalyzer (ABC)
    │   ├── WeightableGroupAnalyzer (ABC)
    │   │   ├── FleissKappaAnalyzer
    │   │   ├── BrennanPredigerAnalyzer
    │   │   └── KrippendorffAlphaAnalyzer
    │   ├── GwetAC1Analyzer
    │   └── ScottsAnalyzer
    └── IndividualMarginalAnalyzer (ABC)
        ├── WeightableIndividualAnalyzer (ABC)
        │   └── CohenKappaAnalyzer
        ├── CongerKappaAnalyzer
        └── LightKappaAnalyzer

| File | What it tests | Tests |
|------|---------------|-------|
| test_configs.py | AnalyzerConfig (defaults, custom, CI methods, random_state, frozen, validation, error prefix, summary, from_dict), HypothesisTestConfig (defaults, custom, methods/schemes, frozen, validation, error prefix), GroupCoefficientConfig (defaults, type guards, error prefix), NoWeightGroupConfig / WeightableGroupConfig (validation, error prefix), FleissKappaConfig / ScottsConfig / GwetAC1Config / BrennanPredigerConfig / KrippendorffAlphaConfig (defaults, parametric_null accepted or rejected per coefficient, parametric_null warning, error prefix), IndividualMarginalConfig / NoWeightIndividualConfig / WeightableIndividualConfig (defaults, validation, error prefix), CohenKappaConfig / CongerKappaConfig / LightKappaConfig (defaults, parametric_null, auto-correction, error prefix), PermutationBootstrapMeanConfig (methods, validation, error prefix), PermutationConfig / PairwisePermutationConfig / SubsetPermutationConfig (construction, validation, from_dict, cross-validation, error prefix) | 139 |
| test_agreement_functions.py | Fleiss' kappa (observed agreement, category proportions, chance, kappa value, unweighted/linear/quadratic weights, NijMatrix/EstimatorInput paths), Krippendorff's alpha, Cohen's kappa, Brennan-Prediger, Scott's pi, Gwet's AC1, Conger's kappa (point estimates against hand-computed expected values from Fleiss 1971 and 2-rater datasets) | 108 |
| test_analyzer_classes.py | Error prefixes (all 8 leaf analyzers), class hierarchy (isinstance checks), config type validation (wrong config rejected for each analyzer), constructor paths (from Ratings, from NijMatrix), estimate guards (missing ratings/nij), rater constraints (2-rater enforcement for Scott's/Cohen), hypothesis test guards, smoke estimates (point estimate + CI for all 8), balance profile propagation, CI method in result, ParametricCIResult disclosure (config warning, construction warning, validation-before-warning, bootstrap path stays silent) | 92 |
| test_ci_methods.py | Error prefixes, percentile_ci (normal, edge alpha, sorted, empty, scalar), _acceleration (formula, sign, degenerate input, noise gate), bca_ci (normal, edge cases, jackknife required, fallback to percentile and its disclosure), TestEfronLawSchool (reproduces the published acceleration, z0, and central 90% interval from Efron & Tibshirani 1986 Example 3), cross-method consistency | 51 |
| test_permutation_analysis.py | PermutationAnalysis (construction, validation, weights, _compute_ci, all 7 estimators), PairwisePermutationAnalysis (construction, observed values, per-pair bootstrap, mean bootstrap correlated/independent, result fields, Fleiss 1971 15-pair), SubsetPermutationAnalysis (construction, observed values, measurable pool per subset, per-subset bootstrap, mean bootstrap, result fields, distribution threshold), MultiKSubsetResult (construction, validation, container API, __str__, summary) | 70 |
| test_resampling_functions.py | subject_bootstrap (count, shape, dtype, rows-from-original, reproducibility, seeds), subject_jackknife (count, shape, leave-one-out), subject_bootstrap_from_nij_matrix (count, shape, row sums, reproducibility), rater_jackknife (count, shape, leave-one-out), rater_subset (C(n,k) count, shape, all combos, validation), within_rater_shuffle / within_subject_shuffle / pool_and_redistribute (count, shape, multiset preservation, missing values, reproducibility, changes data) | 82 |
| test_variance_functions.py | fleiss_null_variance (K=2 equal/skewed, K=3 uniform, hand-computed, validation), cohen_null_variance unweighted/weighted (symmetric, asymmetric, identity=None, positive, validation), cross-consistency (Fleiss=Cohen for 2-rater identical marginals) | 27 |
| conftest.py | Shared fixtures: Fleiss 1971 (NijMatrix + expanded Ratings), 2-rater (10x2, K=2), 3-rater (10x3, K=2), expected-value dicts, minimal configs for all 8 coefficients, PairwisePermutationConfig instances | -- |


## test_example_data/ (33 tests)

Integration tests. Every other package tests one layer with data built for the
purpose; this package loads the datasets shipped in `example_data/` and runs
the whole stack on them, from YAML through the cleaner to a point estimate and
an interval. A failure here means the layers do not fit together, even when
each one passes its own tests.

Both construction pathways are exercised:

    YAML -> RatingsSemantics / CleanerConfig -> RatingsCleaner -> Ratings
    csv/txt -> DataFrame -> ndarray -> NijMatrix

**Resample counts here are the published floors, not the fastest numbers that
pass.** 250 for percentile intervals, 1000 for BCa (Efron & Tibshirani 1993,
§9 p. 72). An integration test that reports an interval should report an
interval the package itself would call adequate. The whole package runs in
about five seconds, which is close to half the suite's total runtime and still
a price worth paying. Raising B is the cheap part: BCa's cost is dominated by
the jackknife, one estimator call per subject regardless of B.

Which datasets carry published reference values, and what each test checks
against them, is in `test_example_data/README.md`. The datasets themselves,
with their provenance and shapes, are in `example_data/README.md`.

| File | What it tests | Tests |
|------|---------------|-------|
| test_reference_data.py | Published values: Fleiss 1971 Table 1 (TestFleiss1971), Gwet's 3x3 continuous example (TestGwetCont3x3), a 2-rater synthetic set (TestSynthetic2Rater), and a BCa run end-to-end through an analyzer (TestBcaEndToEnd) | 11 |
| test_research_data.py | The six research datasets -- SR1 prelim, SR1 main, SR2, SR3, SR4, SR5 -- on real partially crossed data: point estimates, intervals, and the balance profile each one produces | 9 |
| test_weighted.py | Weighted coefficients on SR1 prelim, Gwet cont3x3, and the synthetic 2-rater set: linear and quadratic weights through the full pathway | 4 |
| test_pairwise_permutation.py | PairwisePermutationAnalysis on SR1 prelim: per-pair point estimates, per-pair intervals, and the mean bootstrap | 4 |
| test_subset_permutation.py | SubsetPermutationAnalysis on SR1 main at k=3 and k=2, including the regression that pins subset bootstrapping to the measurable pool rather than to complete cases | 5 |
| test_missing_structure.py | The three shuffles on SR1 main's partially crossed design (half the cells empty): every draw keeps the empty cells in place, and each shuffle keeps what it promises -- each reviewer's answer mix, each abstract's answers, or the overall mix | 5 |
| conftest.py | 11 dataset fixtures (3 nij, 8 ratings), 8 single-coefficient configs at 250 resamples, and builders for the subset, pairwise, weighted, unbalanced and BCa configs | -- |
| README.md | Which datasets are validated against which source, and what each test asserts | -- |


## fixtures/ (test-only data, no tests of its own)

Data that is **not** a ratings dataset and therefore does not belong in
`example_data/`, which a user is meant to be able to load and run an analyzer
on. These files exist to pin a specific piece of arithmetic.

| File | What it is | Used by |
|------|------------|---------|
| efron1986_lawschool_data.csv | 15 law schools, LSAT score and GPA. Efron & Tibshirani (1986) Fig. 1, p. 55. No raters, no categories, no coefficient -- the statistic is a Pearson correlation | test_analyzers/test_ci_methods.py (TestEfronLawSchool) |
| efron1986_lawschool_expected.yaml | The paper's published values: point 0.776, z0 -0.0927, acceleration -0.0817, central 90% interval [0.43, 0.92], and the fifteen empirical influence values it computed the acceleration from. Records the two traps in reproducing it -- the interval is central 90% (alpha = 0.10), and the published acceleration comes from the infinitesimal jackknife, where this package's delete-one form gives -0.0757 | test_analyzers/test_ci_methods.py (TestEfronLawSchool) |
| scruffy_ratings_generator.py | A generator for deliberately damaged ratings, with every corruption recorded so a cleaner can be checked against an answer known in advance | -- |
| scruffy_ratings_example.csv | One 40x5 output of that generator, checked in so the file can be read without running the generator | -- |
| scruffy_cleaner_config.yaml | A config permissive enough to clean that file all the way through. The same file raises under the strict defaults, which is the point of the fixture: what is fatal is decided by the config, not by the data | -- |

**Nothing currently imports the three scruffy files.** They are a complete,
documented fixture with no test behind them -- the cleaner is tested in
`test_cleaners/` against the small inline CLEAN_DF and MESSY_DF instead. Stated
here rather than left to be discovered, because a fixture that looks used and
is not will eventually be trusted by somebody.


## Datasets

Small fixtures defined inline or in a conftest. The datasets loaded from disk
by `test_example_data/` are documented in `example_data/README.md`, with
shapes, provenance and reference values, and are not repeated here.

| Name | Shape | Raters | K | Structure | Used in |
|------|-------|--------|---|-----------|---------|
| Fleiss 1971 Table 1 | 30 x 6 | 6 | 5 | NijMatrix | test_analyzers |
| Fleiss 1971 (expanded) | 30 x 6 | 6 | 5 | category_index expanded from NijMatrix | test_analyzers |
| TWO_RATER_INDEX | 10 x 2 | 2 | 2 | category_index, hand-computed expected values | test_analyzers |
| THREE_RATER_INDEX | 10 x 3 | 3 | 2 | category_index, multi-rater individual tests | test_analyzers |
| Inline 2-rater K=3 | 8 x 2 | 2 | 3 | category_index, near/far misses for weighted Cohen | test_analyzers (inline) |
| MISSING_DATA | 5 x 3 | 3 | 2 | category_index with -1 missing entries | test_analyzers (inline) |
| SIMPLE | 3 x 2 | 2 (positional) | 2 | Fully crossed, no empties | test_data_objects |
| GAPPY | 4 x 3 | 3 (ann, bob, cara) | 2 | Empties, unbalanced | test_data_objects |
| CLEAN_DF | 4 x 2 | 2 (ann, bob) | 2 | DataFrame, no problems | test_cleaners |
| MESSY_DF | 4 x 2 | 2 (ann, bob) | 2 | DataFrame, null + missing + typo | test_cleaners |
| BALANCED_NIJ | 4 x 3 | 3 | 3 | NijMatrix, every row sums to 3 | test_intermediates |
| UNBALANCED_NIJ | 4 x 2 | variable | 2 | NijMatrix, variable row sums | test_intermediates |
| SIMPLE_INDEX | 2 x 3 | 3 | 3 | Minimal NijMatrix for index tests | test_intermediates |
| SIMPLE_SHARED | 3 x 3 | 3 | -- | Symmetric matrix for RaterSharedMatrix | test_intermediates |


## Conventions

- **Tolerance**: 1e-12 for deterministic computations. Bootstrap tests check
  structure (shape, keys, type) and determinism (same seed = same output), not
  statistical properties.
- **Expected values**: exact fractions (13/33, not 0.393939) where possible, so
  rounding drift shows up as a failure rather than hiding inside a tolerance.
- **Frozen arrays**: module-level constants in conftest, copied into fixtures so
  no test can mutate what another test sees.
- **Error prefix**: every validation error names the concrete class the caller
  asked for (`[PairedRatings] Error`, not `[Ratings] Error`), and the tests
  verify this.
- **One file per sub-package**: test_intermediates/ groups tests by sub-package
  (estimator_parts, data_profiles, analysis_results) rather than one file per
  class, keeping the file count manageable as the class count grows.
- **A clean run is silent.** The suite emits no warnings. Two things in the
  package warn on purpose -- a config built with `ci_method="parametric_null"`,
  and the `ParametricCIResult` it produces -- and the tests that own those
  warnings assert them with `pytest.warns`. Anywhere else such an object is
  only scaffolding, the warning is suppressed inside a `warnings.catch_warnings()`
  block scoped to the lines that build it, never with a global filter in
  `pyproject.toml`. A global filter would also hide the next warning worth
  seeing.
