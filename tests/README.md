# Tests

1378 test functions across seven packages, run with `pytest tests/`.
`@pytest.mark.parametrize` expands some of these into multiple cases, so
`pytest --co` reports a higher collected count (~1410).

The test suite is organised bottom-up: lower layers are tested first so that
any failure points at the layer that broke, not at something that depends on it.


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


## test_intermediates/ (455 functions, 473 collected)

Tests for the intermediates layer: estimator parts, data profiles, utility
functions, CI methods, and the full analysis-result hierarchy. One test file
per sub-package, not one per class, so each file covers multiple classes in
bottom-up order.

Built bottom-up along the dependency chain:

    NijMatrix -> Weights -> EstimatorInput
    BalanceProfile -> RaterSharedMatrix
    _profile_utils (prevalence, prevalence_per_rater, coverage_per_rater)
    ci_methods (percentile_ci, bca_ci, scipy_bca_ci)
    PointEstimateResult -> 8 coefficient results
    PermutationResult -> PairwisePermutationResult -> SubsetPermutationResult -> MultiKSubsetResult
    HypothesisTestResult

| File | What it tests | Tests |
|------|---------------|-------|
| test_estimator_parts.py | NijMatrix (construction, validation, properties, from_ratings, equality, error prefix), Weights (identity/linear/quadratic/custom, validation, equality, error prefix), EstimatorInput (construction, validation, properties, subsetting, equality, error prefix) | 104 |
| test_data_profiles.py | BalanceProfile (from_nij_matrix, n_subjects/n_categories/n_raters, is_balanced, rater_counts, summary, equality, error prefix), RaterSharedMatrix (from_nij_matrix, shape, values, symmetry, self-counts, equality, error prefix) | 59 |
| test_profile_utils.py | prevalence (shape, values, sums-to-1), prevalence_per_rater (shape, values, row-sums), coverage_per_rater (shape, values, range) | 11 |
| test_ci_methods.py | percentile_ci (normal, edge alpha, sorted, empty, scalar, parametrized alphas), bca_ci (normal, edge cases, fallback-to-percentile), scipy_bca_ci (delegation, unavailable) | 21 (27 collected) |
| test_analysis_results.py | PointEstimateResult (construction, validation, __str__, summary), HypothesisTestResult (construction, validation, __str__, summary), PermutationResult (construction, inheritance, validation, summary, error prefix), PairwisePermutationResult (construction, validation, accessor, __str__, summary), SubsetPermutationResult (construction, validation, accessor, __str__, summary), MultiKSubsetResult (construction, validation, container protocol, __str__, summary), FleissKappaResult, ScottsPiResult, GwetsAC1Result, BrennanPredigerResult, KrippendorffAlphaResult (balance required, ci_method sets, is_exact logic), CohenKappaResult (weights_type validation), CongerKappaResult (always unweighted), LightKappaResult (n_pairs = C(n,2)), IndividualMarginalResult, error prefixes for all | 260 (272 collected) |
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


## test_analyzers/ (504 functions, 518 collected)

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
| test_configs.py | AnalyzerConfig (defaults, custom, CI methods, random_state, frozen, validation, error prefix, summary, from_dict), HypothesisTestConfig (defaults, custom, methods/schemes, frozen, validation, error prefix), GroupCoefficientConfig (defaults, type guards, error prefix), NoWeightGroupConfig / WeightableGroupConfig (validation, error prefix), FleissKappaConfig / ScottsConfig / GwetAC1Config / BrennanPredigerConfig / KrippendorffAlphaConfig (defaults, parametric_null rejection, weights, error prefix), IndividualMarginalConfig / NoWeightIndividualConfig / WeightableIndividualConfig (defaults, validation, error prefix), CohenKappaConfig / CongerKappaConfig / LightKappaConfig (defaults, parametric_null, weights, auto-correction, error prefix), PermutationBootstrapMeanConfig (methods, validation, error prefix), PermutationConfig / PairwisePermutationConfig / SubsetPermutationConfig (construction, validation, from_dict, cross-validation, error prefix), CustomWeightsInput (construction, validation, error prefix) | 120 |
| test_agreement_functions.py | Fleiss' kappa (observed agreement, category proportions, chance, kappa value, unweighted/linear/quadratic weights, NijMatrix/EstimatorInput paths), Krippendorff's alpha, Cohen's kappa, Brennan-Prediger, Scott's pi, Gwet's AC1, Conger's kappa (point estimates against hand-computed expected values from Fleiss 1971 and 2-rater datasets) | 108 |
| test_analyzer_classes.py | Error prefixes (all 8 leaf analyzers), class hierarchy (isinstance checks), config type validation (wrong config rejected for each analyzer), constructor paths (from Ratings, from NijMatrix), estimate guards (missing ratings/nij), rater constraints (2-rater enforcement for Scott's/Cohen), hypothesis test guards, smoke estimates (point estimate + CI for all 8), balance profile propagation, CI method in result | 74 (86 collected) |
| test_ci_methods.py | percentile_ci, bca_ci, scipy_bca_ci (normal cases, edge alphas, fallback, parametrized) | 29 |
| test_permutation_analysis.py | PermutationAnalysis (construction, validation, weights, _compute_ci, all 7 estimators), PairwisePermutationAnalysis (construction, observed values, per-pair bootstrap, mean bootstrap correlated/independent, result fields, Fleiss 1971 15-pair), SubsetPermutationResult / SubsetPermutationAnalysis (construction, observed values, per-subset bootstrap, mean bootstrap, result fields, distribution threshold), MultiKSubsetResult (construction, validation, container API, __str__, summary) | 64 (66 collected) |
| test_resampling_functions.py | subject_bootstrap (count, shape, dtype, rows-from-original, reproducibility, seeds), subject_jackknife (count, shape, leave-one-out), subject_bootstrap_from_nij_matrix (count, shape, row sums, reproducibility), rater_jackknife (count, shape, leave-one-out), rater_subset (C(n,k) count, shape, all combos, validation), column_shuffle / within_subject_shuffle / pool_and_redistribute (count, shape, multiset preservation, missing values, reproducibility, changes data) | 82 |
| test_variance_functions.py | fleiss_null_variance (K=2 equal/skewed, K=3 uniform, hand-computed, validation), cohen_null_variance unweighted/weighted (symmetric, asymmetric, identity=None, positive, validation), cross-consistency (Fleiss=Cohen for 2-rater identical marginals) | 27 |
| conftest.py | Shared fixtures: Fleiss 1971 (NijMatrix + expanded Ratings), 2-rater (10x2, K=2), 3-rater (10x3, K=2), expected-value dicts, minimal configs for all 8 coefficients, PairwisePermutationConfig instances | -- |


## Datasets

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
