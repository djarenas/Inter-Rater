"""
The array contract.

Every workhorse function in this package operates on one array shape, described
here. Nothing else in the package may invent a second convention.

    category_index : np.ndarray, (n_subjects, n_raters), int8, C-contiguous

Each cell is either a position in ``Categories.names`` or ``EMPTY``.

Why subjects are rows
---------------------
Resampling subjects is then ``category_index[idx]`` -- a contiguous row gather,
the fastest operation numpy offers. "Resample subjects" and "resample rows"
become the same instruction, so a bootstrapper never rebuilds an object.

Why int8 and not strings
------------------------
Integer positions let counts be a single ``bincount`` over a flattened index.
At int8, a 1,617 x 4 dataset is 6.5 KB, so thousands of draws stay in cache.
Names are resolved once, by the cleaner, and never enter the hot path.

Why -1 and not NaN
------------------
NaN forces float64, forbids integer indexing, and turns every count into a
masked operation. ``EMPTY`` keeps the array integral, and ``>= 0`` is a cheap
validity mask.

Where K comes from
------------------
K is ``Categories.n_categories`` and nothing else. It is never read off the
data: ``category_index.max() + 1`` is a lower bound on K and never K itself.

That gap is not academic. A bootstrap draw can miss a rare category entirely,
and a function inferring K from that draw would silently compute with a smaller
K than the rating task had. Kappa would survive it, since a category with no mass
contributes a zero term either way -- but Gwet's AC1 and Brennan-Prediger put
the number of categories in the chance term directly, so their values would
move draw to draw for no reason but the inference. The resulting interval would
look like sampling variability and would not be.

Only an analyser unpacks a Ratings. It reads K from the semantics and passes
it down to the workhorse functions and bootstrappers, which take plain arrays.
Keeping that unpacking in one layer is what stops K from ever being read off a
draw.
"""

from __future__ import annotations

import numpy as np

# ----------------------------------------------------------------- Sentinel #

EMPTY: int = -1
"""No rating exists in this cell.

One sentinel, one meaning: this cell takes no part in any calculation, but it
holds its place so pairwise overlaps and rater permutations stay correct.

I carry no reason. A subject never shown to a rater and a subject shown and
left blank are both simply empty here, because at this level the only fact that
matters is that there is no rating. Why a cell is empty is the cleaner's
business, is recorded in its report, and changes what you can say about a rater
-- never what any coefficient computes.

'Unsure', 'abstain', and anything else a rating task chose to offer are
*categories*, not empties. They were answers. They live in ``Categories.names``
and count like any other response.
"""

# ------------------------------------------------------------------- Dtypes #

INDEX_DTYPE: np.dtype = np.dtype(np.int8)
"""Storage type for the category index array. Caps a rating task at 127 categories."""

LABEL_DTYPE: np.dtype = np.dtype(object)
"""Storage type for subject ids.

``object`` holds real Python strings, so an id of any length survives intact.
A fixed-width dtype such as ``U32`` would be faster, but it truncates anything
longer without complaint, and a silently shortened identifier is a data
corruption nobody would notice.
"""

# ----------------------------------------------------------------- Minimums #

MIN_RATERS: int = 2
"""Fewer than two raters is not an agreement problem."""

MIN_SUBJECTS: int = 2
"""Fewer than two subjects leaves nothing to disagree about."""

MIN_RATINGS_PER_RATER: int = 2
"""Floor on how many subjects a rater must have rated.

A rater with one rating joins one comparison and describes nobody. A cleaner
may set its threshold higher as a deliberate choice, never lower.
"""

MIN_RATINGS_PER_SUBJECT: int = 2
"""Floor on how many raters must have rated a subject.

A subject rated once has nobody to disagree with: it joins no pairwise
comparison and contributes nothing to any group coefficient. A cleaner may set
its threshold higher as a deliberate choice, never lower.
"""

MIN_OVERLAP_HARD: int = 2
"""Mathematical floor on co-rated subjects for a pairwise coefficient.

Below two, expected agreement cannot be estimated and kappa is undefined. This
is a validity threshold, not a reliability one: a pair meeting it can still
produce a number nobody should trust.

The *practical* floor is higher and is not a constant. Under heavy skew,
expected agreement is large and kappa has little room left to work in, so more
co-rated subjects are needed for the same precision. That threshold depends on
prevalence and is not yet in this package; see the roadmap entry on overlap
adequacy.
"""

# --------------------------------------------------- Resampling Coefficients #

SUBJECT_BOOTSTRAP_COEFFICIENTS: frozenset[str] = frozenset([
    "fleiss_kappa",
    "scotts_pi",
    "gwet_ac1",
    "brennan_prediger",
    "percent_agreement",
    "krippendorff_alpha",
    "cohen_kappa",
    "conger_kappa",
    # "gwet_ac2",  # Future: weighted form of AC1, also supports subject bootstrap
])
"""Coefficients that support subject bootstrap (resample subjects, never raters).

A coefficient is eligible for subject bootstrap if:
- It reads counts (NijMatrix) or derives them from Ratings
- It is rater-agnostic (does not depend on which rater said what)
- Or it recomputes per-rater marginals from the resampled Ratings

When a new coefficient is added to agreement_functions/, it must be registered
here before any resampling function will accept it. This is deliberate
enforcement: you cannot bootstrap nonsense silently.

Rater resampling is not supported: at small n, draws duplicate raters and they
agree perfectly with themselves, inflating every coefficient. At three raters,
77% of draws contain a duplicate and 11% return exactly 1.0. Subject bootstrap
is the one that works.
"""

# ------------------------------------------------- Permutation Coefficients #

PAIRWISE_PERMUTATION_COEFFICIENTS: frozenset[str] = frozenset([
    "cohen_kappa",
    "scotts_pi",
    "conger_kappa",
    "gwet_ac1",
    "fleiss_kappa",
    "brennan_prediger",
    "krippendorff_alpha",
])
"""Coefficients that support pairwise permutation (all pairs of raters).

A coefficient is eligible for pairwise permutation if it can be meaningfully
computed on pairs of raters. All agreement coefficients work on 2 raters:
- Cohen's kappa (individual marginals, two raters)
- Scott's pi (pooled marginals, two-rater group coefficient)
- Conger's kappa (individual marginals, computed pairwise)
- Gwet's AC1 (chance-adjusted, no-weight group coefficient)
- Fleiss' kappa (pooled marginals, group coefficient)
- Brennan-Prediger (uniform chance, group coefficient)
- Krippendorff's alpha (reliability coefficient, any number of raters)

Pairwise permutation is useful for detecting raters that may be problematic
or for understanding how each pair of raters agrees.
"""

SUBSET_PERMUTATION_COEFFICIENTS: frozenset[str] = frozenset([
    "fleiss_kappa",
    "scotts_pi",
    "gwet_ac1",
    "brennan_prediger",
    "percent_agreement",
    "krippendorff_alpha",
    # "gwet_ac2",  # Future: weighted form of AC1, also supports subset permutation
])
"""Coefficients that support subset permutation (all k-subsets of raters).

A coefficient is eligible for subset permutation if it can be computed on
subsets of the rater pool. This includes all group coefficients (which read
counts matrices and are rater-agnostic):
- Fleiss' kappa, Scott's pi, Gwet AC1, Brennan-Prediger, Percent agreement
- Krippendorff's alpha

Subset permutation is useful for understanding how agreement scales with panel
size or for analyzing specific subgroups of raters.
"""
