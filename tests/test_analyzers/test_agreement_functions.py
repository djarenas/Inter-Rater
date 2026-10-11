"""
Tests for the agreement_functions sub-package.

I cover the four modules -- group_agreement, individual_marginals,
krippendorff, and rater_combinations -- with hand-computed expected
values verified against the published formulas.

Every expected value is a fraction rather than a rounded decimal, so a
failing assertion points to a logic error rather than a tolerance
question. Tolerance is 1e-12: these are deterministic computations
with no resampling.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from interrater.contracts import EMPTY
from interrater.intermediates.estimator_parts import NijMatrix, Weights
from interrater.intermediates.estimator_parts.estimator_input import EstimatorInput

# -- group_agreement -------------------------------------------------------- #
from interrater.analyzers.agreement_functions.group_agreement import (
    agreement_corrected_by_chance,
    brennan_prediger,
    category_proportions,
    chance_brennan_prediger,
    chance_fleiss,
    chance_gwet_ac1,
    fleiss_kappa,
    gwet_ac1,
    observed_agreement,
    percent_agreement,
    scotts_pi,
)

# -- individual_marginals --------------------------------------------------- #
from interrater.analyzers.agreement_functions.individual_marginals import (
    _observed_pairwise,
    chance_cohen,
    chance_conger,
    cohen_kappa,
    conger_kappa,
    light_kappa,
    rater_marginals,
    shared_rows,
)

# -- krippendorff ----------------------------------------------------------- #
from interrater.analyzers.agreement_functions.krippendorff import (
    expected_disagreement,
    krippendorff_alpha,
    observed_disagreement,
    pairable_ratings,
)

# -- rater_combinations ----------------------------------------------------- #
from interrater.analyzers.agreement_functions.rater_combinations import (
    MAX_SUBSETS,
    leave_one_rater_out,
    mean_over_subsets,
    over_pairs,
    over_subsets,
    rater_subsets,
    subjects_per_leave_one_out,
    subjects_per_subset,
)


TOL = 1e-12


# ========================================================================== #
#  Fixtures                                                                    #
# ========================================================================== #

# -- NijMatrix fixtures ---------------------------------------------------- #

# 4 subjects, 2 categories, 3 raters each (balanced).
BALANCED_NIJ = NijMatrix(np.array([
    [3, 0],
    [3, 0],
    [2, 1],
    [0, 3],
], dtype=np.int64))

# Perfect agreement: 3 subjects, 3 categories, 2 raters.
PERFECT_NIJ = NijMatrix(np.array([
    [2, 0, 0],
    [0, 2, 0],
    [0, 0, 2],
], dtype=np.int64))

# Degenerate: every rating in one category.
DEGENERATE_NIJ = NijMatrix(np.array([
    [3, 0],
    [3, 0],
    [3, 0],
], dtype=np.int64))

# 3 subjects, 3 categories, 2 raters (for weighted tests).
THREE_CAT_NIJ = NijMatrix(np.array([
    [2, 0, 0],
    [0, 2, 0],
    [1, 1, 0],
], dtype=np.int64))

# Unbalanced: subject 2 has a single rater (skipped in all computations).
UNBALANCED_NIJ = NijMatrix(np.array([
    [2, 1],
    [1, 1],
    [1, 0],
], dtype=np.int64))


# -- EstimatorInput fixtures ----------------------------------------------- #

# Two raters, 2 categories, 5 subjects, no missingness.
TWO_RATER_INPUT = EstimatorInput(
    category_index=np.array([
        [0, 0],
        [0, 0],
        [1, 1],
        [0, 1],
        [0, 0],
    ], dtype=np.int64),
    n_categories=2,
)

# Three raters, 2 categories, 4 subjects, no missingness.
THREE_RATER_INPUT = EstimatorInput(
    category_index=np.array([
        [0, 0, 0],
        [1, 1, 1],
        [0, 1, 0],
        [0, 0, 1],
    ], dtype=np.int64),
    n_categories=2,
)

# Two raters, 3 categories, 4 subjects (for weighted Cohen).
TWO_RATER_K3_INPUT = EstimatorInput(
    category_index=np.array([
        [0, 0],
        [1, 1],
        [0, 1],
        [2, 2],
    ], dtype=np.int64),
    n_categories=3,
)

# Sparse: two raters, only one shared subject.
SPARSE_INPUT = EstimatorInput(
    category_index=np.array([
        [0, 0],
        [1, EMPTY],
        [EMPTY, 1],
    ], dtype=np.int64),
    n_categories=2,
)


# ========================================================================== #
#  group_agreement                                                             #
# ========================================================================== #


class TestAgreementCorrectedByChance:
    """I test the shared (Po - Pe) / (1 - Pe) formula."""

    def test_formula(self):
        assert abs(agreement_corrected_by_chance(0.8, 0.5) - 0.6) < TOL

    def test_perfect_agreement(self):
        assert abs(agreement_corrected_by_chance(1.0, 0.5) - 1.0) < TOL

    def test_agreement_at_chance(self):
        assert abs(agreement_corrected_by_chance(0.5, 0.5) - 0.0) < TOL

    def test_chance_is_one_returns_nan(self):
        assert math.isnan(agreement_corrected_by_chance(1.0, 1.0))


class TestObservedAgreement:
    """I test Po: proportion of agreeing rater pairs, averaged over subjects."""

    def test_balanced(self):
        # S0: 1, S1: 1, S2: 1/3, S3: 1. Mean = 5/6.
        assert abs(observed_agreement(BALANCED_NIJ) - 5 / 6) < TOL

    def test_perfect(self):
        assert abs(observed_agreement(PERFECT_NIJ) - 1.0) < TOL

    def test_degenerate(self):
        assert abs(observed_agreement(DEGENERATE_NIJ) - 1.0) < TOL

    def test_unbalanced_skips_single_rater(self):
        # S2 has n=1, skipped. S0: 1/3, S1: 0. Mean = 1/6.
        assert abs(observed_agreement(UNBALANCED_NIJ) - 1 / 6) < TOL

    def test_three_rater(self):
        # nij = [[3,0],[0,3],[2,1],[2,1]]. Po = 2/3.
        nij = THREE_RATER_INPUT.nij
        assert abs(observed_agreement(nij) - 2 / 3) < TOL


class TestCategoryProportions:
    """I test pj: per-subject proportions averaged over subjects."""

    def test_balanced(self):
        pj = category_proportions(BALANCED_NIJ)
        expected = np.array([2 / 3, 1 / 3])
        np.testing.assert_allclose(pj, expected, atol=TOL)

    def test_degenerate(self):
        pj = category_proportions(DEGENERATE_NIJ)
        expected = np.array([1.0, 0.0])
        np.testing.assert_allclose(pj, expected, atol=TOL)

    def test_preserves_k(self):
        # THREE_CAT_NIJ has K=3; third category is unused but must be present.
        pj = category_proportions(THREE_CAT_NIJ)
        assert pj.shape == (3,)

    def test_unbalanced_skips_single_rater(self):
        # S0: [2/3, 1/3], S1: [1/2, 1/2]. Mean: [7/12, 5/12].
        pj = category_proportions(UNBALANCED_NIJ)
        expected = np.array([7 / 12, 5 / 12])
        np.testing.assert_allclose(pj, expected, atol=TOL)


class TestGroupChanceTerms:
    """I test chance terms for the group coefficients."""

    # -- Fleiss / Scott ---------------------------------------------------- #

    def test_chance_fleiss_balanced(self):
        # pj = [2/3, 1/3]. Pe = 4/9 + 1/9 = 5/9.
        assert abs(chance_fleiss(BALANCED_NIJ) - 5 / 9) < TOL

    def test_chance_fleiss_degenerate(self):
        # pj = [1, 0]. Pe = 1.
        assert abs(chance_fleiss(DEGENERATE_NIJ) - 1.0) < TOL

    # -- Gwet AC1 ---------------------------------------------------------- #

    def test_chance_gwet_balanced(self):
        # pj = [2/3, 1/3]. (pj*(1-pj)).sum() = 4/9. / (K-1=1) = 4/9.
        assert abs(chance_gwet_ac1(BALANCED_NIJ) - 4 / 9) < TOL

    def test_chance_gwet_degenerate(self):
        # pj = [1, 0]. (1*0 + 0*1) / 1 = 0.
        assert abs(chance_gwet_ac1(DEGENERATE_NIJ) - 0.0) < TOL

    # -- Brennan-Prediger -------------------------------------------------- #

    def test_chance_bp_two_categories(self):
        assert abs(chance_brennan_prediger(BALANCED_NIJ) - 0.5) < TOL

    def test_chance_bp_three_categories(self):
        assert abs(chance_brennan_prediger(THREE_CAT_NIJ) - 1 / 3) < TOL


class TestGroupCoefficients:
    """I test Fleiss, AC1, Brennan-Prediger, and percent agreement."""

    # -- Fleiss' kappa ----------------------------------------------------- #

    def test_fleiss_balanced(self):
        # (5/6 - 5/9) / (1 - 5/9) = 5/8.
        assert abs(fleiss_kappa(BALANCED_NIJ) - 5 / 8) < TOL

    def test_fleiss_degenerate_is_nan(self):
        # Pe = 1 -- denominator = 0.
        assert math.isnan(fleiss_kappa(DEGENERATE_NIJ))

    def test_fleiss_negative(self):
        # UNBALANCED_NIJ: Po=1/6, Pe_fleiss=37/72. fleiss = -5/7.
        assert abs(fleiss_kappa(UNBALANCED_NIJ) - (-5 / 7)) < TOL

    def test_fleiss_three_rater(self):
        nij = THREE_RATER_INPUT.nij
        assert abs(fleiss_kappa(nij) - 11 / 35) < TOL

    # -- Gwet AC1 ---------------------------------------------------------- #

    def test_ac1_balanced(self):
        # (5/6 - 4/9) / (1 - 4/9) = 7/10.
        assert abs(gwet_ac1(BALANCED_NIJ) - 7 / 10) < TOL

    def test_ac1_degenerate_is_one(self):
        # Po=1, Pe=0 -- ac1 = 1.
        assert abs(gwet_ac1(DEGENERATE_NIJ) - 1.0) < TOL

    def test_ac1_three_rater(self):
        nij = THREE_RATER_INPUT.nij
        assert abs(gwet_ac1(nij) - 13 / 37) < TOL

    # -- Brennan-Prediger -------------------------------------------------- #

    def test_bp_balanced(self):
        # (5/6 - 1/2) / (1 - 1/2) = 2/3.
        assert abs(brennan_prediger(BALANCED_NIJ) - 2 / 3) < TOL

    def test_bp_degenerate_is_one(self):
        # Pe = 1/2, Po = 1 -- bp = 1.
        assert abs(brennan_prediger(DEGENERATE_NIJ) - 1.0) < TOL

    def test_bp_three_rater(self):
        nij = THREE_RATER_INPUT.nij
        assert abs(brennan_prediger(nij) - 1 / 3) < TOL

    # -- percent agreement ------------------------------------------------- #

    def test_percent_agreement_is_po(self):
        assert abs(percent_agreement(BALANCED_NIJ) - 5 / 6) < TOL

    # -- perfect agreement ------------------------------------------------- #

    def test_perfect_agreement_all_one(self):
        assert abs(fleiss_kappa(PERFECT_NIJ) - 1.0) < TOL
        assert abs(gwet_ac1(PERFECT_NIJ) - 1.0) < TOL
        assert abs(brennan_prediger(PERFECT_NIJ) - 1.0) < TOL


class TestScottsPi:
    """I test Scott's pi: unweighted only, rejects non-identity weights."""

    def test_scotts_pi_balanced(self):
        # Same chance term as Fleiss (unweighted). fleiss = 5/8.
        assert abs(scotts_pi(BALANCED_NIJ) - 5 / 8) < TOL

    def test_scotts_pi_equals_fleiss_unweighted(self):
        # For unweighted data, Scott's pi = Fleiss' kappa.
        assert abs(
            scotts_pi(BALANCED_NIJ) - fleiss_kappa(BALANCED_NIJ)
        ) < TOL

    def test_scotts_pi_with_identity_weights(self):
        wi = Weights.identity(2)
        assert abs(scotts_pi(BALANCED_NIJ, wi) - 5 / 8) < TOL

    def test_scotts_pi_rejects_non_identity_weights(self):
        wq = Weights.quadratic(3)
        with pytest.raises(ValueError, match="nominal"):
            scotts_pi(THREE_CAT_NIJ, wq)

    def test_scotts_pi_degenerate_is_nan(self):
        assert math.isnan(scotts_pi(DEGENERATE_NIJ))


class TestGroupWeighted:
    """I test weighted variants of group coefficients using quadratic weights."""

    @pytest.fixture()
    def wq(self):
        return Weights.quadratic(3)

    def test_observed_agreement_weighted(self, wq):
        # THREE_CAT_NIJ, quadratic K=3. Po_w = 11/12.
        assert abs(observed_agreement(THREE_CAT_NIJ, wq) - 11 / 12) < TOL

    def test_chance_fleiss_weighted(self, wq):
        # pj = [1/2, 1/2, 0]. Pe_w = pj . W . pj = 7/8.
        assert abs(chance_fleiss(THREE_CAT_NIJ, wq) - 7 / 8) < TOL

    def test_fleiss_kappa_weighted(self, wq):
        # (11/12 - 7/8) / (1 - 7/8) = 1/3.
        assert abs(fleiss_kappa(THREE_CAT_NIJ, wq) - 1 / 3) < TOL

    def test_chance_bp_weighted(self, wq):
        # sum(W) / K^2 = 6/9 = 2/3.
        assert abs(chance_brennan_prediger(THREE_CAT_NIJ, wq) - 2 / 3) < TOL

    def test_bp_weighted(self, wq):
        # (11/12 - 2/3) / (1 - 2/3) = 3/4.
        assert abs(brennan_prediger(THREE_CAT_NIJ, wq) - 3 / 4) < TOL

    def test_percent_agreement_weighted(self, wq):
        assert abs(percent_agreement(THREE_CAT_NIJ, wq) - 11 / 12) < TOL


# ========================================================================== #
#  individual_marginals                                                        #
# ========================================================================== #


class TestRaterMarginals:
    """I test rater_marginals and shared_rows."""

    def test_marginals_two_rater(self):
        m = rater_marginals(TWO_RATER_INPUT)
        # Rater 0: [4/5, 1/5], Rater 1: [3/5, 2/5].
        expected = np.array([
            [4 / 5, 1 / 5],
            [3 / 5, 2 / 5],
        ])
        np.testing.assert_allclose(m, expected, atol=TOL)

    def test_marginals_three_rater(self):
        m = rater_marginals(THREE_RATER_INPUT)
        expected = np.array([
            [3 / 4, 1 / 4],
            [1 / 2, 1 / 2],
            [1 / 2, 1 / 2],
        ])
        np.testing.assert_allclose(m, expected, atol=TOL)

    def test_marginals_restricted_to_rows(self):
        # Restrict to shared rows of raters 0 and 1. All 4 rows shared,
        # so result is the same as unrestricted.
        rows = shared_rows(THREE_RATER_INPUT, 0, 1)
        m = rater_marginals(THREE_RATER_INPUT, rows)
        assert m.shape == (3, 2)
        np.testing.assert_allclose(m[0], [3 / 4, 1 / 4], atol=TOL)

    def test_shared_rows_complete(self):
        # TWO_RATER_INPUT: no missingness -- all rows shared.
        rows = shared_rows(TWO_RATER_INPUT, 0, 1)
        assert len(rows) == 5

    def test_shared_rows_sparse(self):
        # SPARSE_INPUT: only S0 shared.
        rows = shared_rows(SPARSE_INPUT, 0, 1)
        assert len(rows) == 1
        assert rows[0] == 0


class TestIndividualChanceTerms:
    """I test chance_cohen and chance_conger."""

    def test_chance_cohen_two_rater(self):
        m = rater_marginals(TWO_RATER_INPUT)
        pe = chance_cohen(m[0], m[1])
        # (4/5)(3/5) + (1/5)(2/5) = 14/25.
        assert abs(pe - 14 / 25) < TOL

    def test_chance_cohen_identical_marginals(self):
        # When marginals are the same, Cohen's Pe = Scott's Pe.
        p = np.array([0.6, 0.4])
        pe = chance_cohen(p, p)
        assert abs(pe - (0.36 + 0.16)) < TOL

    def test_chance_cohen_weighted(self):
        wq = Weights.quadratic(3)
        m = rater_marginals(TWO_RATER_K3_INPUT)
        pe = chance_cohen(m[0], m[1], wq)
        # m[0] = [1/2, 1/4, 1/4], m[1] = [1/4, 1/2, 1/4]. Pe_w = 11/16.
        assert abs(pe - 11 / 16) < TOL

    def test_chance_conger_three_rater(self):
        m = rater_marginals(THREE_RATER_INPUT)
        pe = chance_conger(m)
        # pooled = 37/72, var correction = 1/72. Pe = 1/2.
        assert abs(pe - 1 / 2) < TOL

    def test_chance_conger_identical_marginals_equals_fleiss(self):
        # When all raters have the same marginals, variance is zero,
        # so Conger's chance term equals Fleiss'.
        m = np.array([[0.6, 0.4], [0.6, 0.4], [0.6, 0.4]])
        pe_conger = chance_conger(m)
        pe_fleiss = 0.6**2 + 0.4**2
        assert abs(pe_conger - pe_fleiss) < TOL

    def test_chance_conger_rejects_non_identity_weights(self):
        wq = Weights.quadratic(3)
        m = rater_marginals(THREE_RATER_INPUT)
        with pytest.raises(NotImplementedError, match="weighted Conger"):
            chance_conger(m, wq)


class TestObservedPairwise:
    """I test _observed_pairwise: Po for exactly two raters."""

    def test_unweighted(self):
        index = TWO_RATER_INPUT.category_index
        assert abs(_observed_pairwise(index, 2) - 4 / 5) < TOL

    def test_weighted(self):
        wq = Weights.quadratic(3)
        index = TWO_RATER_K3_INPUT.category_index
        # W[a,b] per subject: [1, 1, 3/4, 1]. Mean = 15/16.
        assert abs(_observed_pairwise(index, 3, wq) - 15 / 16) < TOL

    def test_perfect(self):
        index = np.array([[0, 0], [1, 1], [2, 2]], dtype=np.int64)
        assert abs(_observed_pairwise(index, 3) - 1.0) < TOL


class TestCohenKappa:
    """I test Cohen's kappa: two raters, individual marginals."""

    def test_cohen_two_rater(self):
        # Po = 4/5, Pe = 14/25. kappa = 6/11.
        assert abs(cohen_kappa(TWO_RATER_INPUT) - 6 / 11) < TOL

    def test_cohen_three_categories_unweighted(self):
        # Po = 3/4, Pe = 5/16. kappa = 7/11.
        assert abs(cohen_kappa(TWO_RATER_K3_INPUT) - 7 / 11) < TOL

    def test_cohen_weighted(self):
        wq = Weights.quadratic(3)
        # Po_w = 15/16, Pe_w = 11/16. kappa_w = 4/5.
        assert abs(cohen_kappa(TWO_RATER_K3_INPUT, wq) - 4 / 5) < TOL

    def test_cohen_sparse_too_few_shared(self):
        # SPARSE_INPUT: 1 shared subject < 2 -- NaN.
        assert math.isnan(cohen_kappa(SPARSE_INPUT))

    def test_cohen_perfect(self):
        perfect = EstimatorInput(
            category_index=np.array([[0, 0], [1, 1], [0, 0]], dtype=np.int64),
            n_categories=2,
        )
        assert abs(cohen_kappa(perfect) - 1.0) < TOL


class TestCongerKappa:
    """I test Conger's kappa: many raters, individual marginals."""

    def test_conger_three_rater(self):
        # Po = 2/3 (group), Pe = 1/2. conger = 1/3.
        assert abs(conger_kappa(THREE_RATER_INPUT) - 1 / 3) < TOL

    def test_conger_approximates_light_on_complete_data(self):
        # Conger is meant to approximate Light's kappa; on this complete,
        # balanced dataset they agree exactly.
        assert abs(
            conger_kappa(THREE_RATER_INPUT) - light_kappa(THREE_RATER_INPUT)
        ) < TOL


class TestLightKappa:
    """I test Light's kappa: mean of pairwise Cohen's kappas."""

    def test_light_three_rater(self):
        # Pairwise Cohens: (0,1)=1/2, (0,2)=1/2, (1,2)=0. Mean = 1/3.
        assert abs(light_kappa(THREE_RATER_INPUT) - 1 / 3) < TOL


# ========================================================================== #
#  krippendorff                                                                #
# ========================================================================== #


class TestKrippendorff:
    """I test Krippendorff's alpha, built from disagreement rather than
    agreement."""

    # -- building blocks --------------------------------------------------- #

    def test_pairable_ratings_balanced(self):
        # 4 subjects x 3 raters, all pairable -- 12.
        assert pairable_ratings(BALANCED_NIJ) == 12

    def test_pairable_ratings_unbalanced(self):
        # S0: 3, S1: 2, S2: 1 (not pairable). Total = 5.
        assert pairable_ratings(UNBALANCED_NIJ) == 5

    def test_observed_disagreement_balanced(self):
        # Do = 1/6.
        assert abs(observed_disagreement(BALANCED_NIJ) - 1 / 6) < TOL

    def test_observed_disagreement_unbalanced(self):
        # Do = 4/5.
        assert abs(observed_disagreement(UNBALANCED_NIJ) - 4 / 5) < TOL

    def test_observed_disagreement_three_rater(self):
        nij = THREE_RATER_INPUT.nij
        assert abs(observed_disagreement(nij) - 1 / 3) < TOL

    def test_expected_disagreement_balanced(self):
        # De = 16/33.
        assert abs(expected_disagreement(BALANCED_NIJ) - 16 / 33) < TOL

    def test_expected_disagreement_unbalanced(self):
        # De = 3/5.
        assert abs(expected_disagreement(UNBALANCED_NIJ) - 3 / 5) < TOL

    def test_expected_disagreement_three_rater(self):
        nij = THREE_RATER_INPUT.nij
        assert abs(expected_disagreement(nij) - 35 / 66) < TOL

    # -- alpha ------------------------------------------------------------- #

    def test_alpha_balanced(self):
        # 1 - (1/6) / (16/33) = 21/32.
        assert abs(krippendorff_alpha(BALANCED_NIJ) - 21 / 32) < TOL

    def test_alpha_perfect(self):
        assert abs(krippendorff_alpha(PERFECT_NIJ) - 1.0) < TOL

    def test_alpha_degenerate_is_nan(self):
        # De = 0 -- NaN.
        assert math.isnan(krippendorff_alpha(DEGENERATE_NIJ))

    def test_alpha_negative(self):
        # UNBALANCED_NIJ: alpha = -1/3.
        assert abs(krippendorff_alpha(UNBALANCED_NIJ) - (-1 / 3)) < TOL

    def test_alpha_three_rater(self):
        nij = THREE_RATER_INPUT.nij
        assert abs(krippendorff_alpha(nij) - 13 / 35) < TOL

    # -- weighted ---------------------------------------------------------- #

    def test_observed_disagreement_weighted(self):
        wq = Weights.quadratic(3)
        # THREE_CAT_NIJ: Do_w = 1/12.
        assert abs(observed_disagreement(THREE_CAT_NIJ, wq) - 1 / 12) < TOL

    def test_expected_disagreement_weighted(self):
        wq = Weights.quadratic(3)
        # THREE_CAT_NIJ: De_w = 3/20.
        assert abs(expected_disagreement(THREE_CAT_NIJ, wq) - 3 / 20) < TOL

    def test_alpha_weighted(self):
        wq = Weights.quadratic(3)
        # THREE_CAT_NIJ: alpha_w = 4/9.
        assert abs(krippendorff_alpha(THREE_CAT_NIJ, wq) - 4 / 9) < TOL


# ========================================================================== #
#  rater_combinations                                                          #
# ========================================================================== #


class TestRaterSubsets:
    """I test enumeration of rater subsets."""

    def test_pairs_of_four(self):
        result = rater_subsets(4, 2)
        assert len(result) == 6
        assert result == [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    def test_leave_one_out_of_four(self):
        result = rater_subsets(4, 3)
        assert len(result) == 4
        assert (0, 1, 2) in result
        assert (1, 2, 3) in result

    def test_size_equals_count(self):
        result = rater_subsets(3, 3)
        assert result == [(0, 1, 2)]

    def test_size_below_minimum_raises(self):
        with pytest.raises(ValueError, match="at least"):
            rater_subsets(4, 1)

    def test_size_exceeds_count_raises(self):
        with pytest.raises(ValueError, match="asked for subsets of"):
            rater_subsets(3, 4)

    def test_cap_exceeded_raises(self):
        # C(20, 10) = 184,756 > MAX_SUBSETS (10,000).
        with pytest.raises(ValueError, match="past the cap"):
            rater_subsets(20, 10)

    def test_subsets_are_sorted(self):
        for subset in rater_subsets(5, 3):
            assert subset == tuple(sorted(subset))


class TestSubjectsPerSubset:
    """I test how many subjects survive each rater subset."""

    def test_complete_data(self):
        # THREE_RATER_INPUT: all 4 subjects rated by all 3 raters. Every pair -- 4.
        result = subjects_per_subset(THREE_RATER_INPUT, 2)
        assert result == {(0, 1): 4, (0, 2): 4, (1, 2): 4}

    def test_sparse_data(self):
        # SPARSE_INPUT: only S0 rated by both -- 1.
        result = subjects_per_subset(SPARSE_INPUT, 2)
        assert result == {(0, 1): 1}


class TestSubjectsPerLeaveOneOut:
    """I test subjects_per_leave_one_out: keyed by dropped rater."""

    def test_complete_data(self):
        result = subjects_per_leave_one_out(THREE_RATER_INPUT)
        assert result == {0: 4, 1: 4, 2: 4}

    def test_sparse_data(self):
        # SPARSE_INPUT: 2 raters. Dropping either leaves 1 rater per subject.
        # No subject has >= 2 ratings from 1 rater -- 0 each.
        result = subjects_per_leave_one_out(SPARSE_INPUT)
        assert result == {0: 0, 1: 0}


# -- Test estimator for rater_combinations --------------------------------- #


def _po_estimator(data: EstimatorInput) -> float:
    """I am a simple observed-agreement estimator for testing."""
    return observed_agreement(data.nij)


class TestOverPairs:
    """I test over_pairs: run an estimator on every rater pair."""

    def test_three_rater_po(self):
        # Raters (0,1): Po=3/4. (0,2): Po=3/4. (1,2): Po=1/2.
        result = over_pairs(THREE_RATER_INPUT, _po_estimator)
        assert abs(result[(0, 1)] - 3 / 4) < TOL
        assert abs(result[(0, 2)] - 3 / 4) < TOL
        assert abs(result[(1, 2)] - 1 / 2) < TOL

    def test_result_keys_match_subsets(self):
        result = over_pairs(THREE_RATER_INPUT, _po_estimator)
        expected_keys = {(0, 1), (0, 2), (1, 2)}
        assert set(result.keys()) == expected_keys


class TestOverSubsets:
    """I test over_subsets: generalisation of over_pairs to any subset size."""

    def test_size_two_equals_over_pairs(self):
        pairs = over_pairs(THREE_RATER_INPUT, _po_estimator)
        subsets = over_subsets(THREE_RATER_INPUT, _po_estimator, 2)
        assert pairs.keys() == subsets.keys()
        for key in pairs:
            assert abs(pairs[key] - subsets[key]) < TOL

    def test_full_panel(self):
        result = over_subsets(THREE_RATER_INPUT, _po_estimator, 3)
        # One subset: all 3 raters. Po = 2/3.
        assert len(result) == 1
        assert abs(result[(0, 1, 2)] - 2 / 3) < TOL


class TestLeaveOneRaterOut:
    """I test leave_one_rater_out: keyed by the dropped rater."""

    def test_three_rater_po(self):
        result = leave_one_rater_out(THREE_RATER_INPUT, _po_estimator)
        # Drop 0 -- raters 1,2 -- Po = 1/2.
        assert abs(result[0] - 1 / 2) < TOL
        # Drop 1 -- raters 0,2 -- Po = 3/4.
        assert abs(result[1] - 3 / 4) < TOL
        # Drop 2 -- raters 0,1 -- Po = 3/4.
        assert abs(result[2] - 3 / 4) < TOL

    def test_result_keys_are_rater_indices(self):
        result = leave_one_rater_out(THREE_RATER_INPUT, _po_estimator)
        assert set(result.keys()) == {0, 1, 2}

    def test_two_raters_raises(self):
        # With 2 raters, leave-one-out needs subsets of size 1,
        # which is below MIN_RATERS -- ValueError.
        with pytest.raises(ValueError, match="at least"):
            leave_one_rater_out(TWO_RATER_INPUT, _po_estimator)


class TestMeanOverSubsets:
    """I test mean_over_subsets: average, skipping NaN."""

    def test_basic(self):
        values = {(0, 1): 0.5, (0, 2): 0.5, (1, 2): 0.0}
        assert abs(mean_over_subsets(values) - 1 / 3) < TOL

    def test_skips_nan(self):
        values = {(0, 1): 0.5, (0, 2): float("nan")}
        assert abs(mean_over_subsets(values) - 0.5) < TOL

    def test_all_nan_returns_nan(self):
        values = {(0, 1): float("nan"), (0, 2): float("nan")}
        assert math.isnan(mean_over_subsets(values))

    def test_leave_one_out_mean(self):
        result = leave_one_rater_out(THREE_RATER_INPUT, _po_estimator)
        # mean(1/2, 3/4, 3/4) = 2/3.
        assert abs(mean_over_subsets(result) - 2 / 3) < TOL


# ========================================================================== #
#  Edge cases                                                                  #
# ========================================================================== #


class TestEdgeCases:
    """I test cross-cutting edge cases: NaN returns, degenerate inputs."""

    def test_all_one_category_nij(self):
        # Single-category data: Po=1, Pe_fleiss=1 -- NaN, Pe_gwet=0 -- 1.
        nij = NijMatrix(np.array([[4, 0], [4, 0]], dtype=np.int64))
        assert math.isnan(fleiss_kappa(nij))
        assert abs(gwet_ac1(nij) - 1.0) < TOL
        assert math.isnan(krippendorff_alpha(nij))

    def test_no_usable_subjects(self):
        # Every subject rated by at most 1 rater.
        nij = NijMatrix(np.array([[1, 0], [0, 1]], dtype=np.int64))
        assert math.isnan(observed_agreement(nij))
        assert math.isnan(krippendorff_alpha(nij))

    def test_identity_weights_equal_unweighted(self):
        wi = Weights.identity(2)
        assert abs(
            observed_agreement(BALANCED_NIJ, wi)
            - observed_agreement(BALANCED_NIJ)
        ) < TOL
        assert abs(
            fleiss_kappa(BALANCED_NIJ, wi)
            - fleiss_kappa(BALANCED_NIJ)
        ) < TOL
        assert abs(
            krippendorff_alpha(BALANCED_NIJ, wi)
            - krippendorff_alpha(BALANCED_NIJ)
        ) < TOL

    def test_gwet_ac1_rejects_weights(self):
        # AC1 takes no weights argument -- TypeError.
        wq = Weights.quadratic(3)
        with pytest.raises(TypeError):
            gwet_ac1(THREE_CAT_NIJ, wq)

    def test_conger_rejects_weights(self):
        # conger_kappa takes no weights -- TypeError at call site.
        wq = Weights.quadratic(3)
        with pytest.raises(TypeError):
            conger_kappa(THREE_RATER_INPUT, wq)

    def test_light_rejects_weights(self):
        # light_kappa takes no weights -- TypeError at call site.
        wq = Weights.quadratic(3)
        with pytest.raises(TypeError):
            light_kappa(THREE_RATER_INPUT, wq)
