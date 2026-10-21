# test_example_data/

Integration tests that run analyzers on the shipped example datasets and assert
point estimates against the validated reference values in `example_data/README.md`.

These sit above the unit test suite: a unit test checks that `FleissKappaAnalyzer`
handles a 4x3 NijMatrix correctly; an integration test checks that the full
pipeline -- load file, clean, build Ratings/NijMatrix, configure analyzer,
estimate -- reproduces a published result end-to-end.


## Design decisions

**Point estimates only (deterministic).** Each test asserts the point estimate
to 3 decimal places. This is deterministic -- no bootstrap needed, runs in
milliseconds. CI structure (lower < point < upper) is checked as a sanity
guard.

**The guard uses percentile at B=250, and the method is not incidental.**
Under percentile the point estimate always falls inside its own interval, so
asserting it is a real check -- a failure means the bootstrap drew from a
different set of subjects than the point estimate was computed from, which is
exactly the bug `test_subset_permutation.py` exists for. Under BCa it is not
an invariant at all: strong bias correction can place the interval entirely to
one side of the point, legitimately, which is why `PointEstimateResult` stopped
requiring it. B=250 is Efron & Tibshirani (1986) section 9's minimum for a
useful percentile interval; below it every interval these tests build would
report `below_recommended_resamples` and say so when printed.

These tests previously used BCa at B=50, which was under the published minimum
of 1000 for *that* method and paid for a jackknife besides -- BCa needs one
estimator call per subject regardless of B, so on the larger fixtures the
jackknife cost far more than the bootstrap. Percentile at 250 runs the four
largest fixtures in 0.2s against 1.8s.

**BCa is still covered end-to-end**, once, in `test_reference_data.py` at
B=1000 on the 30-subject Fleiss fixture. See `bca_fleiss_config`.

**One conftest, flat test files.** `conftest.py` provides dataset-loading
fixtures. Each test file covers one logical group (reference data, research
data, permutations).

**Tolerance.** `pytest.approx(expected, abs=5e-4)` for 3-decimal-place
agreement. This is looser than the unit test tolerance (1e-12) because these
are end-to-end checks, not formula-level assertions.

**No slow marks for now.** Everything uses low `n_resamples`. A future
`@pytest.mark.slow` tier could run full `n_resamples=4000` assertions against
the README's CI bounds.


## Test matrix

### Coverage goals

1. All 8 leaf analyzers appear at least once
2. All 5 systematic review datasets represented
3. Both data pathways exercised (ratings and nij)
4. Both balance conditions covered (balanced and unbalanced)
5. Both pairwise bootstrap methods tested (correlated and independent)
6. Weighted analysis (linear and quadratic) tested on multiple datasets
   across both WeightableGroupConfig and WeightableIndividualConfig

### Analyzer coverage map

| Analyzer | Unweighted | Weighted | Total |
|---|---|---|---|
| FleissKappaAnalyzer | Fleiss 1971, SR1 prelim, SR1 main | SR1 prelim (quadratic) | 4 |
| GwetAC1Analyzer | Fleiss 1971, SR2, SR4 | — (not weightable) | 3 |
| BrennanPredigerAnalyzer | Fleiss 1971, SR3, SR5 | Synthetic 2-rater (linear) | 4 |
| KrippendorffAlphaAnalyzer | Gwet cont3x3, SR1 main | SR1 prelim (linear) | 3 |
| CohenKappaAnalyzer | Gwet cont3x3, Synthetic 2-rater | Gwet cont3x3 (quadratic) | 3 |
| ScottsAnalyzer | Synthetic 2-rater | — (not weightable) | 1 |
| CongerKappaAnalyzer | SR1 prelim | — (not weightable) | 1 |
| LightKappaAnalyzer | SR1 prelim | — (not weightable) | 1 |
| PairwisePermutationAnalysis | SR1 prelim (corr. + indep.) | — | 2 |
| SubsetPermutationAnalysis | SR1 main (k=2 and k=3) | — | 5 |


## Files

| File | What it tests | Tests |
|------|---------------|-------|
| test_reference_data.py | Datasets with published reference values, plus the one end-to-end BCa test | 11 |
| test_research_data.py | Original systematic review datasets (internal validation against README) | 9 |
| test_subset_permutation.py | Subset permutation on partially crossed data -- the measurable-pool regression | 5 |
| test_weighted.py | Weighted analysis (linear + quadratic) across both weightable config branches | 4 |
| test_pairwise_permutation.py | Pairwise permutation analysis with both bootstrap methods | 4 |
| conftest.py | Dataset-loading fixtures, shared analyzer config builders | -- |

**Total: 33 tests.** Counts here are from `pytest --collect-only`; if you
change a file, re-derive rather than adjust by hand.


## test_reference_data.py (11 tests)

Tests against datasets with externally published reference values.

### Fleiss 1971 (nij pathway, 6 raters, 5 categories, balanced)

| Test | Analyzer | Expected | Source |
|------|----------|----------|--------|
| test_fleiss_kappa | FleissKappaAnalyzer | 0.430 | Fleiss (1971) original |
| test_gwet_ac1 | GwetAC1Analyzer | 0.448 | V2 validated |
| test_brennan_prediger | BrennanPredigerAnalyzer | 0.444 | V2 validated, Pe = 1/5 |

### Gwet cont3x3 (ratings + nij pathways, 2 raters, 3 categories, balanced)

| Test | Analyzer | Expected | Source |
|------|----------|----------|--------|
| test_cohen_kappa_ratings | CohenKappaAnalyzer | 0.796 | irrCAC: 0.79641 |
| test_gwet_ac1_ratings | GwetAC1Analyzer | 0.849 | irrCAC: 0.84933 |
| test_krippendorff_alpha_nij | KrippendorffAlphaAnalyzer | 0.797 | V2 validated |

### Synthetic 2-rater (ratings pathway, 2 raters, 3 categories, balanced)

| Test | Analyzer | Expected | Source |
|------|----------|----------|--------|
| test_cohen_kappa | CohenKappaAnalyzer | 0.641 | V2 validated |
| test_scotts_pi | ScottsAnalyzer | 0.641 | V2 validated, = Fleiss at n=2 |


## test_research_data.py (9 tests)

Tests against the original systematic review datasets. Each SR gets one
coefficient chosen to exercise a different analyzer and avoid redundancy.

### SR1 prelim (ratings, 4 raters, K=3, balanced, ~83% skew)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_fleiss_kappa | FleissKappaAnalyzer | 0.493 | Group coefficient on balanced |
| test_conger_kappa | CongerKappaAnalyzer | 0.494 | Group-level, individual marginals |
| test_light_kappa | LightKappaAnalyzer | 0.496 | Mean of C(4,2)=6 pairwise Cohens; n_pairs=6 |

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

### SR3 FIS Psych (ratings, 3 raters, K=2, balanced, ~90% skew)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_brennan_prediger | BrennanPredigerAnalyzer | 0.886 | B-P on binary, 3 raters |

### SR4 FIS CVD (ratings, 4 raters, K=2, balanced, ~82% skew)

| Test | Analyzer | Expected | Notes |
|------|----------|----------|-------|
| test_gwet_ac1 | GwetAC1Analyzer | 0.879 | AC1 on binary, 4 raters |

### SR5 Cocaine Cardio (ratings, 3 raters, K=2, balanced, ~88% skew)

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

Subset permutation on SR1 main, the only shipped dataset where a subset's
**measurable pool** and its **complete cases** differ.

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


## test_weighted.py (4 tests)

Weighted analysis across multiple datasets. Tests that distance-based weighting
changes the point estimate (weighted != unweighted) and exercises both
WeightableGroupConfig (Fleiss, B-P, Krippendorff) and WeightableIndividualConfig
(Cohen).

Weights only matter when K >= 3 -- with binary categories, linear and quadratic
collapse to identity (there is only one distance). All datasets below have K=3.

### SR1 prelim (ratings, 4 raters, K=3, balanced, ordinal -1/0/1)

| Test | Analyzer | Weights | Expected | Notes |
|------|----------|---------|----------|-------|
| test_fleiss_kappa_quadratic | FleissKappaAnalyzer | quadratic | 0.593 | -1/0 closer than -1/1; WeightableGroupConfig |
| test_krippendorff_alpha_linear | KrippendorffAlphaAnalyzer | linear | 0.555 | Different weight type on same data |

### Gwet cont3x3 (ratings, 2 raters, K=3, balanced)

| Test | Analyzer | Weights | Expected | Notes |
|------|----------|---------|----------|-------|
| test_cohen_kappa_quadratic | CohenKappaAnalyzer | quadratic | 0.892 | WeightableIndividualConfig; irrCAC cross-check possible |

### Synthetic 2-rater (ratings, 2 raters, K=3, balanced, ordinal -1/0/1)

| Test | Analyzer | Weights | Expected | Notes |
|------|----------|---------|----------|-------|
| test_brennan_prediger_linear | BrennanPredigerAnalyzer | linear | 0.651 | Second WeightableGroup analyzer with weights |

*Weighted values established by discovery run (V2 validated). Unweighted
references for comparison: Fleiss 0.493, Krippendorff 0.493, Cohen 0.796,
B-P 0.670. For Gwet cont3x3 weighted Cohen, irrCAC may provide an external
cross-check.*

### Structural assertions (all weighted tests)

Each weighted test also checks:
- `result.point != unweighted_point` (weights changed the estimate)
- CI structure (lower < point < upper)
- `result.weights_type` matches the requested type


## Not yet covered

**Unbalanced nij.** The `unbalanced_nij` fixture exists in `conftest.py` and
**no test uses it.** This section previously listed two tests against it with
`TBD` expected values, as though they were written and awaiting numbers. They
were never written. Writing them means establishing reference values for
`Unbalanced_nij` first, since there are none published.

Keeping the fixture is reasonable -- the nij pathway on unbalanced data is a
real gap -- but a table of tests that do not exist is worse than an admission
that they do not.


## conftest.py

### Data loading

Two helpers handle file → object conversion:

- **`_load_ratings(data_path, semantics_path, cleaner_config_path)`** --
  Full RatingsCleaner pipeline.  Reads the data file into a DataFrame,
  builds `RatingsSemantics` and `CleanerConfig` from their YAML files,
  and returns `RatingsCleaner(df, semantics, config).ratings`.

- **`_load_nij(data_path, *, index_col, sep)`** -- Reads a counts file
  into `NijMatrix(array)`.

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
| fleiss1971_nij | NijMatrix from Fleiss1971_psychiatric_data_nij.csv |
| gwet_ratings | Ratings from gwet_cont3x3_data.csv (via cleaner) |
| gwet_nij | NijMatrix from gwet_cont3x3_data_nij.csv |
| synthetic_ratings | Ratings from synthetic_2raters_data.txt (via cleaner) |
| sr1_prelim_ratings | Ratings from SR1_prelim_data.csv (via cleaner) |
| sr1_main_ratings | Ratings from SR1_main_data.txt (via cleaner) |
| sr2_ratings | Ratings from SR2_data.txt (via cleaner) |
| sr3_ratings | Ratings from SR_FIS_Psych_data.txt (via cleaner) |
| sr4_ratings | Ratings from SR_FIS_CVD_data.txt (via cleaner) |
| sr5_ratings | Ratings from SR_cocaine_cardio_data.txt (via cleaner) |
| unbalanced_nij | NijMatrix from unbalanced_data_nij.txt |

### Config builders

Three private builders, so that the resample count and CI method are defined
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
| Fleiss 1971 | Fleiss (1971) | nij (csv) | 6 | 5 | Balanced | 3 |
| Gwet cont3x3 | Gwet (2014) Table 4.1 | ratings (csv) + nij (csv) | 2 | 3 | Balanced | 3 + 1 weighted |
| Synthetic 2-rater | Synthetic | ratings (txt) | 2 | 3 | Balanced | 2 + 1 weighted |
| SR1 prelim | OSF 2018 | ratings (csv) | 4 | 3 | Balanced | 3 + 2 weighted + 4 pairwise |
| SR1 main | OSF 2018 | ratings (txt) | 4 | 3 | Partially crossed | 2 + 5 subset |
| SR2 | JSRC 2019 | ratings (txt) | 4 | 3 | Nearly balanced | 1 |
| SR3 FIS Psych | JGIM 2019 | ratings (txt) | 3 | 2 | Balanced | 1 |
| SR4 FIS CVD | -- | ratings (txt) | 4 | 2 | Balanced | 1 |
| SR5 CocCardio | -- | ratings (txt) | 3 | 2 | Balanced | 1 |
| Unbalanced nij | Synthetic | nij (txt) | variable | 3 | Unbalanced | **0 -- fixture unused** |

**Total: 33 tests** (11 reference + 9 research + 5 subset + 4 weighted + 4 pairwise)
