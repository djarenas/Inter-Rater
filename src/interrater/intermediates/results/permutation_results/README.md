# Permutation Results

One coefficient, recomputed over subsets of the raters. Where the analyzer
results answer *what is the coefficient*, these answer *how much does it move
when the rater set changes* -- across all pairs, across leave-k-out subsets at
one k, or across several k at once.

```
PermutationResult               coefficient_name, ci_method, balance_profile
├── PairwisePermutationResult   all C(n,2) pairs
├── SubsetPermutationResult     leave-k-out at a single k
└── MultiKSubsetResult          several k values at once
```

| Class                       | What It Holds                                                                 | Why It Exists                                                        |
|-----------------------------|--------------------------------------------------------------------------------|----------------------------------------------------------------------|
| `PermutationResult`         | coefficient_name, ci_method, balance_profile                                   | Shared base. Deliberately carries no point estimate -- see below.    |
| `PairwisePermutationResult` | pairwise_values, subjects_per_pair, n_pairs, mean_over_pairs, bootstrap_cis_per_pair; bootstrap_distributions_per_pair when n_pairs <= 20; optionally bootstrap_mean_method, bootstrap_ci_of_mean, bootstrap_distribution_of_mean | All C(n,2) pairwise coefficients, with per-pair overlap counts.      |
| `SubsetPermutationResult`   | subset_values, subjects_per_subset, n_subsets, subset_k, mean_over_subsets, bootstrap_cis_per_subset; bootstrap_distributions_per_subset when n_subsets <= 120; the same three optional fields for the mean | Leave-one-out or leave-k-out at a single k.                          |
| `MultiKSubsetResult`        | per_k_results (dict[int, SubsetPermutationResult]), ks                          | Stability curve: how the coefficient moves as more raters are removed. |

## No single point estimate

`PermutationResult` has no `point` field, which is the main thing separating
it from `PointEstimateResult`. A permutation analysis produces a *collection*
of coefficients -- one per pair, one per subset -- and the mean over them is
reported as its own field (`mean_over_pairs`, `mean_over_subsets`) rather than
promoted to "the" estimate.

That distinction is not cosmetic. **The values in these collections are not
independent**: two pairs sharing a rater are positively correlated, so the
variance of their mean is not the variance of an average of independent
things. Treating the mean as an ordinary point estimate invites exactly the
analytic standard error that assumption would license, which is how version 1
produced intervals roughly 28% too narrow on real data.

## Every interval here is a `BootstrapCIResult`

`bootstrap_ci_of_mean` holds a `BootstrapCIResult`, and so does every value in
`bootstrap_cis_per_pair` and `bootstrap_cis_per_subset`. They were `(lower,
upper)` tuples before, and the change is not cosmetic either: a tuple can only
say where the bounds are, and these intervals routinely need to say more than
that.

Read bounds through `.lower` and `.upper`, and check `.is_defined` before
using them:

```python
ci = result.bootstrap_cis_per_pair[("ann", "bob")]
if ci.is_defined:
    print(ci.lower, ci.upper)
```

`.as_tuple()` is there when a tuple is genuinely what is wanted, but reaching
for it as a habit throws away exactly the fields that make the object worth
having -- how many draws failed, whether BCa fell back to percentile, whether
the usable count fell short of the method's minimum. `__str__` carries all of
it, which is why the result classes print their intervals rather than
formatting bounds themselves.

A thin pair or subset produces an interval with `method_used = None`: no
method ran, which is a different statement from an interval that was computed
and came out wide. The accessors and `summary()` preserve that distinction
rather than flattening it to a pair of NaNs.

## The bootstrap fields are optional

Every bootstrap field on the pairwise and subset results is `| None`:
`bootstrap_mean_method`, `bootstrap_ci_of_mean`,
`bootstrap_distribution_of_mean`, and the per-pair or per-subset dictionaries.
A permutation analysis is useful without any bootstrap at all -- the spread of
the values is itself the diagnostic -- so resampling is requested, not assumed,
and absence is recorded as `None` rather than as an empty container.

When a bootstrap is run, it resamples subjects and recomputes every pair or
subset within each draw, which preserves the correlation above without having
to model it. `bootstrap_mean_method` records which scheme produced the
distribution.

## `subjects_per_pair` and `subjects_per_subset`

Each value is accompanied by the number of subjects it was computed from.
Under a partially crossed design two raters may overlap on very few subjects,
and a coefficient from a thin overlap should not be read beside one from a
thick overlap without that context travelling with it.

**What counts as a subject here is a rule, not an obvious fact.** A subject
contributes to a k-rater subset when **at least two of those k raters rated
it** -- not when all k did. Fewer than two ratings is not disagreement, it is
absence: Fleiss' per-subject term `(sum_j n_ij^2 - n_i) / (n_i (n_i - 1))` is
0/0 at `n_i = 1`. At k = 2 the rule reduces to both raters, which is the
complete-case pairing Cohen's kappa requires, so it is one rule at every k
rather than a special case for pairs.

That pool is what the bootstrap draws from, and the reason is estimand
matching: an interval has to be around the thing the point estimate measured.
Resampling the complete cases of a 3-rater subset instead silently switches
estimands, and on partially crossed data the symptom is an interval that does
not contain its own point estimate. On this package's SR1 main dataset the two
pools differ enormously -- 756 measurable subjects against 0 complete cases
for one 3-rater subset -- so the mistake is not a rounding difference.

A subset with fewer than two measurable subjects has nothing to resample at
all. It gets a NaN distribution and an interval with `method_used = None`.

## Getting at the pieces

`MultiKSubsetResult` is a container: `result[3]` returns the
`SubsetPermutationResult` for k = 3, `len(result)` counts the k values, and
iterating yields them in sorted order. The `ks` field lists them.

`get_bootstrap_distributions_per_pair()` and
`get_bootstrap_distributions_per_subset()` return the raw draw arrays, for
callers who want to plot or re-analyse the distributions rather than read the
summarised interval. The per-pair and per-subset bootstrap always runs, so
every pair or subset has its interval; the raw draws are kept only up to 20
pairs or 120 subsets, to bound memory, and above that the accessor says so.

`summary()` on each class returns plain types, and each interval appears as
its own nested `summary()` dict -- the full record, not just the two bounds.
When the mean was bootstrapped, `bootstrap_draws_mean` and `bootstrap_se`
describe its usable draws, as in `PointEstimateResult`. The mean's draws must
agree with the mean's interval, and each stored pair's or subset's draws
with that unit's interval.

## Dependencies

- `..analysis_result` (AnalysisResult)
- `..ci_results` (BootstrapCIResult)
- `...data_profiles` (BalanceProfile)
- `numpy` (bootstrap distributions)
