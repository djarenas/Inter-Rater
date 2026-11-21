# test_example_data/

Integration tests that run analyzers on the shipped example datasets and assert
point estimates against the validated reference values in `example_data/README.md`.

These sit above the unit test suite: a unit test checks that `FleissKappaAnalyzer`
handles a small array correctly; an integration test checks that the full
pipeline -- load file, clean, build Ratings, configure analyzer,
estimate-- reproduces a published result end-to-end.


## Design decisions

**Point estimates only (deterministic).** Each test asserts the point estimate
to 3 decimal places. This is deterministic, no bootstrap needed, and runs in
milliseconds. CI structure (lower < point < upper) is also asserted, not to
test the interval but as a guard on the resampling: it fails only when the
bootstrap drew from different subjects than the point was computed from (next
item).

**The guard uses percentile at B=250, and the method is not incidental.**
- Under percentile the point estimate always falls inside its own interval, so
asserting it is a real check. A failure means the bootstrap drew from a
different set of subjects than the point estimate was computed from, which is
exactly the bug `test_subset_permutation.py` exists for. 
- Under BCa, strong bias correction can place the interval entirely to
one side of the point, legitimately, which is why `PointEstimateResult` stopped
requiring it. B=250 is Efron & Tibshirani (1986) section 9's minimum for a
useful percentile interval; below it every interval these tests build would
report `below_recommended_resamples` and say so when printed.

**BCa is still covered end-to-end**, once, in `test_reference_data.py` at
B=1000 on the 30-subject Fleiss fixture. See `bca_fleiss_config`.

**One conftest, flat test files.** `conftest.py` provides dataset-loading
fixtures. Each test file covers one logical group (reference data, research
data, permutations).

**Tolerance.** `pytest.approx(expected, abs=5e-4)` for 3-decimal-place
agreement. They are looser than the unit test tolerance (1e-12) because these
are end-to-end checks.

**No slow marks for now.** Everything uses low `n_resamples`. A future
`@pytest.mark.slow` tier could run full `n_resamples=4000` assertions against
the README's CI bounds.


## Test matrix

### Coverage goals

1. All 8 leaf analyzers appear at least once
2. All 5 systematic review datasets represented
3. A dataset published only as counts runs through the ratings route (Fleiss 1971)
4. Both balance conditions covered (balanced and unbalanced)
5. Both pairwise bootstrap methods tested (correlated and independent)
6. Weighted analysis (linear and quadratic) tested on multiple datasets
   for both group and individual-marginal analyzers

### Analyzer coverage map

| Analyzer | Unweighted | Weighted | Total |
|---|---|---|---|
| FleissKappaAnalyzer | Fleiss 1971, SR1 prelim, SR1 main | SR1 prelim (quadratic) | 4 |
| GwetAC1Analyzer | Fleiss 1971, SR2, SR4, Gwet cont3x3 | Gwet cont3x3, SR1 prelim (linear, quadratic; AC2) | 8 |
| BrennanPredigerAnalyzer | Fleiss 1971, SR3, SR5 | Synthetic 2-rater (linear) | 4 |
| KrippendorffAlphaAnalyzer | Gwet cont3x3, SR1 main | SR1 prelim (linear) | 3 |
| CohenKappaAnalyzer | Gwet cont3x3, Synthetic 2-rater | Gwet cont3x3 (quadratic) | 3 |
| ScottsAnalyzer | Synthetic 2-rater | — (not weightable) | 1 |
| CongerKappaAnalyzer | SR1 prelim | — (not weightable) | 1 |
| LightKappaAnalyzer | SR1 prelim | — (not weightable) | 1 |
| PairwisePermutationAnalysis | SR1 prelim (corr. + indep.) | SR1 prelim (semantics weights) | 3 |
| SubsetPermutationAnalysis | SR1 main (k=2 and k=3) | SR1 prelim (semantics weights) | 6 |


## Files

| File | What it tests | Tests |
|------|---------------|-------|
| test_reference_data.py | Datasets with published reference values, plus the one end-to-end BCa test | 11 |
| test_research_data.py | Original systematic review datasets (internal validation against README) | 10 |
| test_subset_permutation.py | Subset permutation on partially crossed data -- the measurable-pool regression | 5 |
| test_weighted.py | Weighted analysis (linear + quadratic) for group and individual-marginal analyzers, Gwet's AC2 against irrCAC, and weights declared in the semantics | 20 |
| test_pairwise_permutation.py | Pairwise permutation analysis with both bootstrap methods | 4 |
| test_missing_structure.py | The three shuffles keep SR1 main's missing-data design | 5 |
| test_example_configs.py | Every YAML shipped in `example_data/` loads with the current code, one test per file, found by globbing | 42 |
| test_statistical_background.py | `example_data/Statistical_background/`: the committed tables reproduce, interrater gives the page's numbers, and a wrong page value fails | 3 |
| test_made_up_ratings.py | `example_data/made_up_ratings.py`: counts come back exactly, names and ids are kept, short rows get empty slots, a subject rated once is refused | 4 |
| test_external_validation.py | The package against five other implementations, from saved values | 20 |
| conftest.py | Dataset-loading fixtures, shared analyzer config builders | -- |

**Total: 124 tests.** Counts here are from `pytest --collect-only`; if you
change a file, re-derive rather than adjust by hand.


## test_reference_data.py (11 tests)

Tests against datasets with externally published reference values.

### Fleiss 1971 (6 made-up rater slots, 5 categories, balanced)

| Test | Analyzer | Expected | Source |
|------|----------|----------|--------|
| test_fleiss_kappa | FleissKappaAnalyzer | 0.430 | Fleiss (1971) original |
| test_gwet_ac1 | GwetAC1Analyzer | 0.448 | V2 validated |
| test_brennan_prediger | BrennanPredigerAnalyzer | 0.444 | V2 validated, Pe = 1/5 |

### Gwet cont3x3 (2 raters, 3 categories, balanced)

| Test | Analyzer | Expected | Source |
|------|----------|----------|--------|
| test_cohen_kappa_ratings | CohenKappaAnalyzer | 0.796 | irrCAC: 0.79641 |
| test_gwet_ac1_ratings | GwetAC1Analyzer | 0.849 | irrCAC: 0.84933 |
| test_krippendorff_alpha_ratings | KrippendorffAlphaAnalyzer | 0.797 | V2 validated |

### Synthetic 2-rater (ratings pathway, 2 raters, 3 categories, balanced)

| Test | Analyzer | Expected | Source |
|------|----------|----------|--------|
| test_cohen_kappa | CohenKappaAnalyzer | 0.641 | V2 validated |
| test_scotts_pi | ScottsAnalyzer | 0.641 | V2 validated, = Fleiss at n=2 |


## test_research_data.py (10 tests)

Tests against the original systematic review datasets. Each SR gets one
coefficient chosen to exercise a different analyzer and avoid redundancy.

### SR1 prelim (ratings, 4 raters, K=3, balanced, ~83% skew)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_fleiss_kappa | FleissKappaAnalyzer | 0.493 | Group coefficient on balanced |
| test_conger_kappa | CongerKappaAnalyzer | 0.494 | Group-level, individual marginals |
| test_light_kappa | LightKappaAnalyzer | 0.496 | Mean of C(4,2)=6 pairwise Cohens; n_pairs=6 |
| test_unquoted_category_names | FleissKappaAnalyzer | 0.493 | `names: [1, 0, -1]` unquoted in YAML gives the same ratings and kappa |

Conger's group-level estimate (0.494) differs from the pairwise-mean Cohen
kappa (0.496) because it uses each rater's full marginal distribution, not
pair-by-pair marginals.  Light's kappa (0.496) is by definition the mean of
all pairwise Cohen kappas, so it matches the pairwise mean exactly.

### SR1 main (ratings, 4 raters, K=3, unbalanced, ~90% skew)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_fleiss_kappa | FleissKappaAnalyzer | 0.619 | Unbalanced approximation; `balanced_data_only=False` |
| test_krippendorff_alpha | KrippendorffAlphaAnalyzer | 0.618 | Exact on unbalanced |

Note: `FleissKappaConfig` defaults to `balanced_data_only=True`, which rejects
unbalanced data.  This test builds the config with `balanced_data_only=False`
to exercise the unbalanced-data approximation path.  Other configs (Krippendorff,
AC1, B-P) default to `False`.

### SR2 (ratings, 4 raters, K=3, nearly balanced 80/81)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_gwet_ac1 | GwetAC1Analyzer | 0.782 | AC1 on near-balanced, 3 categories |

### SR3 (ratings, 3 raters, K=2, balanced, ~90% skew)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_brennan_prediger | BrennanPredigerAnalyzer | 0.886 | B-P on binary, 3 raters |

### SR4 (ratings, 4 raters, K=2, balanced, ~82% skew)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_gwet_ac1 | GwetAC1Analyzer | 0.879 | AC1 on binary, 4 raters |

### SR5 (ratings, 3 raters, K=2, balanced, ~88% skew)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_brennan_prediger | BrennanPredigerAnalyzer | 0.788 | B-P on binary, 3 raters |


## test_pairwise_permutation.py (4 tests)

Pairwise permutation analysis on SR1 prelim (4 raters, balanced).
Tests both bootstrap methods to verify the CI width relationship.

### SR1 prelim pairwise (C(4,2) = 6 pairs)

| Test | Coefficient | Bootstrap | Mean | CI | Notes |
|------|-------------|-----------|------|-----|-------|
| test_cohen_correlated | Cohen κ | correlated | 0.496 | ~[0.359, 0.596] | Correct uncertainty; width 0.238 |
| test_cohen_independent | Cohen κ | independent | 0.496 | ~[0.437, 0.558] | Width 0.121 -- 51% of correlated |
| test_correlated_wider_than_independent | -- | -- | -- | -- | Structural: correlated CI width > independent CI width |
| test_n_pairs | -- | -- | -- | -- | Result has n_pairs = 6 |


## test_subset_permutation.py (5 tests)

Subset permutation on SR1 main, where a subset's **measurable pool** and its
**complete cases** differ by hundreds of subjects.

A subject contributes to a subset's agreement only if at least two of that
subset's raters rated it. Before 2.0.1 the bootstrap drew from the complete
cases -- subjects *all* k raters rated. On balanced data the two sets are
identical, which is why no test caught it:

| k=3 subset | measurable (>=2 of 3) | complete case (all 3) |
|---|---|---|
| (0, 1, 2) | 756 | 0 |
| (0, 1, 3) | 760 | 0 |
| (0, 2, 3) | 757 | 2 |
| (1, 2, 3) | 967 | 1 |

| Test | What it pins |
|------|--------------|
| test_subset_sizes_use_the_measurable_pool | The divergence above. **Not a regression test** -- see below |
| test_every_subset_brackets_its_own_estimate | Regression. Fails under the old rule |
| test_no_subset_failed_to_compute | Regression. Fails under the old rule |
| test_measurable_equals_complete_case (k=2) | "At least two of the k" is one rule, not two |
| test_subsets_vary_in_size (k=2) | Pairs rest on 164 to 426 subjects -- not the same study |

**Only two of the five catch the bug.** Verified by reverting `measurable` to
`present.all(axis=1)`: the two marked Regression fail, with
`method_used=None` and 250 of 250 draws failed on three subsets. The first
one passes either way, because `subjects_per_subset` comes from Computation A,
which hands the whole array to the estimator and was never wrong. It earns its
place by establishing that the two pools diverge here at all -- without that,
the other two would prove nothing.


## test_weighted.py (20 tests)

Weighted analysis across multiple datasets. Tests that distance-based weighting
changes the point estimate (weighted != unweighted), for the group analyzers
with a weighted form (Fleiss, Gwet's AC2, B-P, Krippendorff) and for Cohen.

Weights only matter when K >= 3 -- with binary categories, linear and quadratic
collapse to identity (there is only one distance). All datasets below have K=3.

### SR1 prelim (ratings, 4 raters, K=3, balanced, ordinal -1/0/1)

| Test | Analyzer | Weights | Expected | Notes |
|------|----------|---------|----------|-------|
| test_fleiss_kappa_quadratic | FleissKappaAnalyzer | quadratic | 0.593 | -1/0 closer than -1/1 |
| test_krippendorff_alpha_linear | KrippendorffAlphaAnalyzer | linear | 0.555 | Different weight type on same data |

### Gwet cont3x3 (ratings, 2 raters, K=3, balanced)

| Test | Analyzer | Weights | Expected | Notes |
|------|----------|---------|----------|-------|
| test_cohen_kappa_quadratic | CohenKappaAnalyzer | quadratic | 0.892 | irrCAC: 0.89216 |

### Gwet's AC2 against irrCAC

The formula is the one in Gwet's irrCAC 1.4 (`gwet.ac1.raw`), so the
references are its values, rounded to 5 decimals; checked to 5e-6.

| Test | Data | Weights | irrCAC |
|------|------|---------|--------|
| test_matches_irrcac | Gwet cont3x3 | linear | 0.90002 |
| test_matches_irrcac | Gwet cont3x3 | quadratic | 0.94024 |
| test_matches_irrcac | SR1 prelim | linear | 0.85247 |
| test_matches_irrcac | SR1 prelim | quadratic | 0.86749 |
| test_unweighted_is_ac1 | Gwet cont3x3 | unordered: AC1, no AC2 in the printout | 0.84933 |

### Synthetic 2-rater (ratings, 2 raters, K=3, balanced, ordinal -1/0/1)

| Test | Analyzer | Weights | Expected | Notes |
|------|----------|---------|----------|-------|
| test_brennan_prediger_linear | BrennanPredigerAnalyzer | linear | 0.651 | Second group analyzer with weights |

### Weights from the semantics (SR1 prelim)

| Test | What it pins |
|------|--------------|
| test_file_and_inline_give_identical_results (x4: Fleiss, AC2, B-P, Krippendorff) | `weights_file` and an inline matrix give the same ratings and results |
| test_quadratic_matrix_in_a_file_reproduces_ordered_quadratic (x4) | A weighted matrix equal to quadratic reproduces `ordered, quadratic` exactly |
| test_loads_from_any_working_directory | A relative `weights_file` is read from the semantics file's folder |
| test_pairs_use_the_semantics_weights | Pairwise permutation computes every pair with the semantics' weights |
| test_subsets_use_the_semantics_weights | Subset permutation does the same |

*Weighted values established by discovery run (V2 validated). Unweighted
references for comparison: Fleiss 0.493, Krippendorff 0.493, Cohen 0.796,
B-P 0.670, AC1 0.849. The AC2 references are irrCAC's values, and irrCAC
gives 0.89216 for Gwet cont3x3 weighted Cohen.*

### Structural assertions (all weighted tests)

Each weighted test in the dataset sections above also checks:
- `result.point != unweighted_point` (weights changed the estimate)
- CI structure (lower < point < upper)
- `result.weights_type` matches the requested type


## test_missing_structure.py (5 tests)

The three shuffles behind the permutation tests, on SR1 main: 1617 abstracts,
4 reviewers, half of all cells empty, each abstract read by 2 of the 4 (3 by
three). Who was shown what is the design of the study, so every shuffle must
leave the empty cells exactly where they are: no answer may land on an
abstract its reviewer was never shown, and every reviewer keeps their
workload. Each shuffle test runs 200 draws.

| Test | What it pins |
|------|--------------|
| test_half_the_cells_are_empty | The fixture is the design the rest rely on: 1617 x 4, 3231 empty cells |
| test_each_abstract_read_by_two_or_three | Every abstract was read by 2 or 3 reviewers |
| test_within_rater_shuffle | Empty cells stay put; each reviewer keeps their own answer mix |
| test_within_subject_shuffle | Empty cells stay put; each abstract keeps the answers it got |
| test_pool_and_redistribute | Empty cells stay put; the overall mix of answers is unchanged |

Each shuffle test also checks that answers do move, so a shuffle that returned
its input unchanged fails. Verified by replacing the three shuffles with naive
ones that shuffle empty cells along with the answers, as a generic routine
would: all three shuffle tests fail. The small-fixture versions of these
checks are in `tests/test_analyzers/test_resampling_functions.py`.


## test_external_validation.py (20 tests)

`analyses/external_validation/results.csv` holds values that irrCAC and irr
(R), and statsmodels, scikit-learn and krippendorff (Python) computed on the
shipped datasets. None of them are installed in CI; these tests recompute the
package's side with the analysis's own `measure.py` and check it against the
saved values.

| Test | What it checks |
|---|---|
| test_values_match (x9, one per dataset) | Every estimate and null SE another tool computed: within 5e-6 of irrCAC (it rounds to 5 decimals), within 1e-10 of the others |
| test_refused_rows_still_refused | What only another tool computes (weighted Scott/Conger/Light, Conger and Light on incomplete data) is still refused |
| test_irr_kripp_still_differs | The one known disagreement, irr's `kripp.alpha` on complete data with more than 2 raters, is the only one in `results.csv` |
| test_matches_cleaning (x9) | Each file in `data/` is still what cleaning produces from `example_data/` |

Verified by moving one saved value by 1e-6: test_values_match fails and names
the row. When a test here fails after an intended change, rerun the four
scripts in `analyses/external_validation/` (the order is in `measure.py`).


## Not yet covered

Nothing listed.


## conftest.py

### Data loading

**`_load_ratings(data_path, semantics_path, cleaner_config_path)`** runs the
full RatingsCleaner pipeline. It reads the data file into a DataFrame, builds
`RatingsSemantics` and `CleanerConfig` from their YAML files, and returns
`RatingsCleaner(df, semantics, config).ratings`.

**Header detection** uses the file extension: `.csv` → `header=0`,
`.txt` → `header=None`.  This is a property of the file format, NOT of
`CleanerConfig.frame_layout.use_column_names`.  The `use_column_names`
flag controls rater-label matching inside `RatingsCleaner` (whether
DataFrame column labels must equal declared rater names), not whether the
source file has a header row.  The Gwet CSV proves this: it has headers
(`Subject,Rater1,Rater2`) but `use_column_names: false` (raters are
positional, not named).

### Fixtures

**Dataset fixtures** (scope="module"):

| Fixture | What it provides |
|---------|------------------|
| fleiss1971_counts | The published counts, Fleiss1971_psychiatric_data_nij.csv |
| fleiss1971_ratings | Ratings built in memory from those counts by `example_data/made_up_ratings.py` (made-up rater slots) |
| gwet_ratings | Ratings from gwet_cont3x3_data.csv (via cleaner) |
| synthetic_ratings | Ratings from synthetic_2raters_data.txt (via cleaner) |
| sr1_prelim_ratings | Ratings from SR1_prelim_data.csv (via cleaner) |
| sr1_prelim_with_semantics | Callable: the SR1 prelim data read under another semantics file, for checking that two declarations of the same thing agree |
| sr1_main_ratings | Ratings from SR1_main_data.txt (via cleaner) |
| sr2_ratings | Ratings from SR2_data.txt (via cleaner) |
| sr3_ratings | Ratings from SR3_data.txt (via cleaner) |
| sr4_ratings | Ratings from SR4_data.txt (via cleaner) |
| sr5_ratings | Ratings from SR5_data.txt (via cleaner) |

### Config builders

Two private builders, so that the resample count and CI method are defined
once rather than at every call site. `_INTEGRATION_RESAMPLES = 250` is the
single number they share.

| Builder | What it makes |
|---------|---------------|
| `_integration_config(cls, **overrides)` | percentile at B=250 -- the default for everything |
| `_bca_config(cls, **overrides)` | bca at B=1000, for the one end-to-end BCa test |

**Config fixtures** (scope="session"), all percentile at B=250 unless noted:

| Fixture | Config class | Notes |
|---------|--------------|-------|
| fast_fleiss_config | FleissKappaConfig | |
| fast_gwet_config | GwetAC1Config | |
| fast_bp_config | BrennanPredigerConfig | |
| fast_krippendorff_config | KrippendorffAlphaConfig | |
| fast_cohen_config | CohenKappaConfig | |
| fast_scotts_config | ScottsConfig | |
| fast_conger_config | CongerKappaConfig | |
| fast_light_config | LightKappaConfig | |
| bca_fleiss_config | FleissKappaConfig | **bca at B=1000**; the only BCa path here |
| fleiss_unbalanced_config | FleissKappaConfig | `balanced_data_only=False`, for SR1 main |

**Callable fixtures**, where the varying part is an argument rather than a
separate fixture per combination:

| Fixture | Call | Why a callable |
|---------|------|----------------|
| ordered_ratings | `ordered_ratings("sr1_prelim", "linear")` | each (dataset, weighting) pair would otherwise be its own fixture |
| pairwise_config | `pairwise_config("correlated")` | the mean-bootstrap method is what varies |
| subset_config | `subset_config(3)` | the tests compare behaviour across k |

The `fast_*` names predate the move to percentile. They are still accurate in
spirit -- the suite got roughly eight times faster -- but they describe speed
rather than purpose, and `_integration_config` is the better name for what
they are.


## Datasets used

| Dataset | Source | Format | Raters | K | Balance | Tests |
|---------|--------|--------|--------|---|---------|-------|
| Fleiss 1971 | Fleiss (1971) | counts (csv), turned into Ratings in memory | 6 | 5 | Balanced | 3 |
| Gwet cont3x3 | Gwet (2014) Table 4.1 | ratings (csv) | 2 | 3 | Balanced | 3 + 1 weighted |
| Synthetic 2-rater | Synthetic | ratings (txt) | 2 | 3 | Balanced | 2 + 1 weighted |
| SR1 prelim | OSF 2018 | ratings (csv) | 4 | 3 | Balanced | 4 + 11 weighted + 4 pairwise |
| SR1 main | OSF 2018 | ratings (txt) | 4 | 3 | Partially crossed | 2 + 5 subset |
| SR2 | JSRC 2019 | ratings (txt) | 4 | 3 | Nearly balanced | 1 |
| SR3 | JGIM 2019 | ratings (txt) | 3 | 2 | Balanced | 1 |
| SR4 | -- | ratings (txt) | 4 | 2 | Balanced | 1 |
| SR5 | -- | ratings (txt) | 3 | 2 | Balanced | 1 |
| Statistical background | Synthetic, from the page's tables | ratings (csv) | 2 to 32 | 2 | Balanced and unbalanced | 3 |

**Total: 124 tests** (11 reference + 10 research + 5 subset + 20 weighted + 4
pairwise + 5 missing structure + 42 example configs + 3 statistical background
+ 4 made-up ratings + 20 external validation)
