# Resampling Functions

Bootstrap, jackknife, and permutation resampling for confidence intervals, null hypothesis testing, and diagnostics.

## Why resampling, not formulas

Published variance formulas for kappa exist, but they have problems:

- **They answer the wrong question.** Fleiss 1971 Eq. 16 and the 1979 correction are both derived under the hypothesis of no agreement beyond chance. Using them to put an interval around a *non-zero* kappa is a category error, and the interval comes out about 43% too narrow.

- **They assume more than your data holds.** Analytic variances assume multinomial ratings and independent subjects. If your subjects cluster or your subjects are not exchangeable, the formula is wrong and nothing says so.

- **They need different derivations per coefficient.** Getting it right took three papers and a decade for Cohen's kappa alone.

**Resampling gives one procedure for all coefficients.** Subject bootstrap works the same way for every coefficient. The cost is computation; the benefit is correctness, and one framework that applies to all of them.

## The eight functions

### Subject-Level Resampling

#### `subject_bootstrap(category_index, n_resamples=1000, random_state=42)`

Resample *subjects* with replacement. Yields n_resamples modified arrays, each drawn with replacement from the original subject rows.

**Supports:** All coefficients (group and pairwise).

**Use for:** Confidence intervals around point estimates.

**Why this works:** Subjects are the unit of independence. Resampling them preserves the correlation structure between raters across subjects. Works at any sample size.

#### `subject_jackknife(category_index)`

Leave-one-subject-out. Yields n_subjects modified arrays, each with one subject removed.

**Supports:** All coefficients.

**Use for:** Bias diagnostics, stability analysis. How much does the estimate change when we drop each subject?

#### `subject_bootstrap_from_nij_matrix(nij_matrix, n_resamples=1000, random_state=42)`

Resample *subjects* with replacement from the nij_matrix (the category-counts matrix, shape n_subjects x n_categories). Yields n_resamples resampled nij matrices.

**Supports:** Group coefficients that operate on nij matrices (Fleiss' kappa, Scott's pi, Gwet's AC1).

**Use for:** Confidence intervals when working directly with nij matrices rather than raw ratings. Same bootstrap logic as `subject_bootstrap`, but on the aggregated representation.

### Rater-Level Resampling

#### `rater_jackknife(category_index)`

Leave-one-rater-out. Yields n_raters modified arrays, each with one rater's column removed.

**Supports:** Pairwise coefficients (Cohen's kappa) and any coefficient that works on subsets.

**Use for:** Rater diagnostics. How critical is each rater to the agreement? Which rater matters most?

#### `rater_subset(category_index, k)`

All combinations of k raters from n. Yields C(n, k) modified arrays, each containing only k raters' columns.

**Supports:** All coefficients.

**Use for:** Scaling analysis. How does agreement change with panel size? What is the smallest panel that still gives stable results?

### Permutation Resampling

These generate null distributions for hypothesis testing. Conceptually different from bootstrap: they shuffle without replacement to destroy structure, testing whether the observed agreement exceeds what chance would produce.

#### `column_shuffle(category_index, n_permutations=5000, random_state=42)`

Shuffle each rater's column independently. Preserves each rater's marginal distribution (how often they use each category) but breaks the subject-rater pairing.

**Null hypothesis:** "Raters agree no more than expected given their individual category-usage biases."

**Use for:** The strongest permutation test that respects rater-level marginals.

#### `within_subject_shuffle(category_index, n_permutations=5000, random_state=42)`

Shuffle category assignments within each subject's row. Preserves per-subject category counts but breaks the association between rater identity and category choice.

**Null hypothesis:** "The pattern of which rater chose which category within a subject is random."

**Note:** For group coefficients that discard rater identity (Fleiss, Scott, AC1), this shuffle does not change the nij matrix and therefore does not change the coefficient. More meaningful for pairwise coefficients (Cohen's kappa) where rater identity matters.

#### `pool_and_redistribute(category_index, n_permutations=5000, random_state=42)`

Pool all non-missing ratings and redistribute randomly across the matrix, preserving only global category proportions.

**Null hypothesis:** "There is no structure in the data beyond the overall category proportions."

**Use for:** The strongest null -- breaks both subject-rater pairing and rater-level marginals.

## Generator Design (Iterators)

All eight functions return **generators** (Iterators), not lists. Each time you iterate, the function yields one resampled array.

**Why generators?**

Memory efficiency at scale:

```python
# List approach (stores everything in RAM):
resamples = list(subject_bootstrap(ratings, n_resamples=10000))  # 10,000 arrays in memory

# Generator approach (one array at a time):
for resample in subject_bootstrap(ratings, n_resamples=10000):
    compute_coefficient(resample)  # Process one, discard, get next
```

At 1000+ resamples and large datasets, generators are essential.

**Speed:** Negligible overhead. The coefficient computation (milliseconds) dwarfs the resampling overhead (microseconds).

**Analyzer usage:**

```python
class FleissKappaAnalyzer:
    def estimate(self, n_resamples=1000):
        point = fleiss_kappa(nij)
        
        distribution = []
        for resample in subject_bootstrap(self.ratings.category_index, n_resamples):
            # Convert resample to NijMatrix
            nij_resample = category_counts(resample, self.semantics.n_categories)
            # Compute coefficient
            distribution.append(fleiss_kappa(nij_resample))
        
        ci = percentile_ci(np.array(distribution))
        return FleissKappaResult(point, ci, balance_profile)
```

## Rater Resampling at Small N

**Rater bootstrap is degenerate at small n and not implemented.** At three raters, 77% of draws duplicate a rater and 11% return exactly 1.0. It would work at crowdsourcing scale (100+ raters) but not at typical panel sizes.

**Rater jackknife and subset are diagnostic tools**, not inferential about a population of raters. They show how agreement depends on which raters you have, not how it would generalize to a population.

## Module Structure

```
resampling_functions/
  __init__.py
  subject_resampling.py       # subject_bootstrap, subject_jackknife, subject_bootstrap_from_nij_matrix
  rater_resampling.py         # rater_jackknife, rater_subset
  permutation_resampling.py   # column_shuffle, within_subject_shuffle, pool_and_redistribute
  README.md
```

## Adding New Resampling Functions

New resampling functions should:

1. **Return a generator** (Iterator[np.ndarray])
2. **Yield modified category_index arrays** (same shape or subset shape)
3. **Document what they yield and when to use them**
4. **Add to __all__ in __init__.py**

Example (leave-one-pair-out):

```python
def pairwise_jackknife(category_index: np.ndarray) -> Iterator[np.ndarray]:
    """Leave-one-pair-out across all rater pairs."""
    n_raters = category_index.shape[1]
    
    for i in range(n_raters):
        for j in range(i + 1, n_raters):
            keep = [k for k in range(n_raters) if k not in (i, j)]
            yield category_index[:, keep]
```

Analyzers decide how to use it. No coupling to coefficients or contracts.

## References

- Efron, Bradley. 1979. "Bootstrap Methods: Another Look at the Jackknife." *The Annals of Statistics* 7(1): 1--26.
- Efron, B. & Tibshirani, R. J. 1993. *An Introduction to the Bootstrap.* Chapman & Hall. Chapter 15: Permutation tests.
- Fleiss, Joseph L. 1971. "Measuring nominal scale agreement among many raters." *Psychological Bulletin* 76(5): 378--382.
- Fleiss, Joseph L., Jacob Cohen, and Brian S. Everitt. 1969. "Large sample standard errors of kappa and weighted kappa." *Psychological Bulletin* 72(5): 323--325.
- Fleiss, Joseph L., Jacob B. Nee, and J. Richard Landis. 1979. "The large sample standard error of kappa and weighted kappa." *Journal of the American Statistical Association* 74(365): 30--35.
