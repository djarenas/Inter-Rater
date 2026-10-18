"""
Integration tests for subset permutation analysis on partially crossed data.

These exist for one reason: **which subjects a rater subset can be measured
on**. A subject contributes to a subset's agreement only if at least two of
that subset's raters rated it -- fewer than two ratings is not disagreement,
it is absence, and Fleiss' per-subject term ``(sum_j n_ij^2 - n_i) /
(n_i (n_i - 1))`` is 0/0 at ``n_i = 1``. The bootstrap has to resample that
same pool, or the interval is around a different quantity than the point
estimate.

Before 2.0.1 the pool was the complete cases -- subjects *all* k raters
rated. On balanced data that is the same set, which is why no test caught
it. SR1 main is partially crossed, and there the two differ enormously:

    subset       measurable (>=2 of 3)    complete case (all 3)
    (0, 1, 2)              756                      0
    (0, 1, 3)              760                      0
    (0, 2, 3)              757                      2
    (1, 2, 3)              967                      1

Three of those four subsets had nothing to resample at all, and the fourth
rested on two subjects. This is the only place in the suite where the
distinction is visible: every other shipped dataset is balanced or near
enough that the pools coincide.

**Which tests here actually catch the bug.** Only the two that look at the
bootstrap: ``test_every_subset_brackets_its_own_estimate`` and
``test_no_subset_failed_to_compute``. Reverting the pool to the complete
cases makes both fail, with ``method_used=None`` and 250 of 250 draws
failed on three subsets.

``test_subset_sizes_use_the_measurable_pool`` does *not* fail under the old
rule, and is not a regression test despite reading like one.
``subjects_per_subset`` comes from Computation A, which hands the whole
array to the estimator and was never wrong. It is here to establish the
numbers the other tests rest on, and to record the divergence.

SR1 main: ratings, 1617 subjects, 4 raters, K=3, partially crossed.
"""

import itertools

import pytest

from interrater import SubsetPermutationAnalysis
from interrater.contracts import EMPTY

# Tolerance: 3 decimal places.
TOL = 5e-4

# SR1 main, k=3. Measured from the data; see the module docstring for what
# the complete-case counts would have been.
MEASURABLE_K3 = {
    (0, 1, 2): 756,
    (0, 1, 3): 760,
    (0, 2, 3): 757,
    (1, 2, 3): 967,
}


def _complete_case_counts(ratings, k: int) -> dict:
    """How many subjects every k-subset's raters ALL rated.

    The rule this package used before 2.0.1, computed here so the tests can
    assert against it directly rather than against remembered numbers.
    """
    category_index = ratings.category_index
    n_raters = category_index.shape[1]
    counts = {}
    for subset in itertools.combinations(range(n_raters), k):
        present = category_index[:, list(subset)] != EMPTY
        counts[subset] = int(present.all(axis=1).sum())
    return counts


class TestSR1MainSubsetsK3:
    """k=3 on SR1 main -- where measurable and complete-case diverge."""

    def test_subset_sizes_use_the_measurable_pool(
        self, sr1_main_ratings, subset_config
    ):
        """Record the divergence, and pin Computation A's pool.

        This is **not** a regression test for the bootstrap bug, though it
        reads like one. ``subjects_per_subset`` comes from Computation A,
        which was always correct; reverting the bootstrap's pool leaves
        this passing. The tests below are the ones that fail.

        What it does do is establish that measurable and complete-case
        differ on this dataset, which is what makes the tests below mean
        anything -- on balanced data they would pass either way. The
        complete-case half turns four numbers into evidence: without it,
        756 is a value someone could update to whatever the code happens
        to produce.
        """
        result = SubsetPermutationAnalysis(
            sr1_main_ratings, subset_config(3)
        ).run()
        per_k = result.per_k_results[3]

        assert per_k.subjects_per_subset == MEASURABLE_K3

        complete = _complete_case_counts(sr1_main_ratings, 3)
        assert complete == {
            (0, 1, 2): 0,
            (0, 1, 3): 0,
            (0, 2, 3): 2,
            (1, 2, 3): 1,
        }
        # Every subset must be measured on more subjects than the old rule
        # allowed -- that difference is the bug.
        for subset, n_measurable in per_k.subjects_per_subset.items():
            assert n_measurable > complete[subset]

    def test_every_subset_brackets_its_own_estimate(
        self, sr1_main_ratings, subset_config
    ):
        """Regression test. Fails under the complete-case rule.

        An interval around a narrower pool than the point estimate was
        computed from is an interval around a different quantity, and it
        can land nowhere near its own estimate. Under percentile on a
        matched pool, containment holds.

        Verified by reverting ``measurable`` to ``present.all(axis=1)``:
        three of these four subsets then produce no interval at all, and
        this fails on the first of them.
        """
        result = SubsetPermutationAnalysis(
            sr1_main_ratings, subset_config(3)
        ).run()
        per_k = result.per_k_results[3]

        for subset, ci in per_k.bootstrap_cis_per_subset.items():
            point = per_k.subset_values[subset]
            assert ci.is_defined, f"subset {subset} produced no interval: {ci}"
            assert ci.lower <= point <= ci.upper, (
                f"subset {subset}: point {point:.4f} outside its own "
                f"interval {ci}"
            )

    def test_no_subset_failed_to_compute(self, sr1_main_ratings, subset_config):
        """Regression test. Fails under the complete-case rule.

        All four intervals must be real, and none a fallback or an
        absence. Verified by reverting ``measurable`` to
        ``present.all(axis=1)``, which gives

            BootstrapCIResult(lower=nan, upper=nan, n_resamples=250,
                              n_failed=250, method_used=None)

        on three of the four -- a correct report of nothing having
        happened, and nothing is not the answer here.
        """
        result = SubsetPermutationAnalysis(
            sr1_main_ratings, subset_config(3)
        ).run()
        per_k = result.per_k_results[3]

        assert per_k.n_subsets == 4
        for subset, ci in per_k.bootstrap_cis_per_subset.items():
            assert ci.method_used == "percentile", f"subset {subset}: {ci}"
            assert ci.n_failed == 0, f"subset {subset}: {ci}"
            assert ci.below_recommended_resamples is False


class TestSR1MainSubsetsK2:
    """k=2 on SR1 main -- where the rule has to agree with pairing."""

    def test_measurable_equals_complete_case(
        self, sr1_main_ratings, subset_config
    ):
        """"At least two of the k" is one rule, not two.

        At k=2 it says both raters, which is the complete-case pairing
        Cohen's kappa requires; above k=2 it admits unbalanced subsets.
        The pairwise case is not an exception to the rule, it is the rule
        at k=2 -- and on this dataset, where k=3 diverges wildly, k=2 must
        not diverge at all.
        """
        result = SubsetPermutationAnalysis(
            sr1_main_ratings, subset_config(2)
        ).run()
        per_k = result.per_k_results[2]

        complete = _complete_case_counts(sr1_main_ratings, 2)
        assert per_k.subjects_per_subset == complete

    def test_subsets_vary_in_size(self, sr1_main_ratings, subset_config):
        """Partially crossed data means pairs rest on different studies.

        Spelled out because it is the reason subjects_per_subset exists:
        a pair measured on 164 subjects and one measured on 426 are not
        estimating the same thing, however close their coefficients look.
        """
        result = SubsetPermutationAnalysis(
            sr1_main_ratings, subset_config(2)
        ).run()
        sizes = result.per_k_results[2].subjects_per_subset

        assert len(set(sizes.values())) > 1
        assert min(sizes.values()) < max(sizes.values()) / 2
