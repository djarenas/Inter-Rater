"""
Shared fixtures for test_analyzers/.

I hold data fixtures (Ratings, NijMatrix) and minimal config instances
that the analyzer test files share. The raw numeric arrays and expected
values are inlined here -- they are small, frozen constants duplicated
from test_agreement_functions/conftest.py so that this conftest is
self-contained (no cross-directory imports).

Construction recipe (bottom-up):
    UnorderedCategories -> RatingsSemantics -> Ratings
    np.ndarray          -> NijMatrix
"""

import numpy as np
import pytest

from interrater.contracts import INDEX_DTYPE
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.data_objects.ratings import Ratings
from interrater.intermediates.estimator_parts import NijMatrix

from interrater.analyzers.configs import (
    AnalyzerConfig,
    FleissKappaConfig,
    GwetAC1Config,
    BrennanPredigerConfig,
    KrippendorffAlphaConfig,
    ScottsConfig,
    CohenKappaConfig,
    CongerKappaConfig,
    LightKappaConfig,
)
from interrater.analyzers.configs.hypothesis_test_config import HypothesisTestConfig
from interrater.analyzers.configs.pairwise_permutation_config import PairwisePermutationConfig
from interrater.analyzers.configs.subset_permutation_config import SubsetPermutationConfig
from interrater.analyzers.configs.permutation_bootstrap_mean_config import PermutationBootstrapMeanConfig
from interrater.analyzers.convergence_test import ConvergenceTest


# -- Raw arrays (inlined from test_agreement_functions/conftest.py) --------- #
# Source of truth: Fleiss (1971) Table 1 and a hand-computed 2-rater example.
# See test_agreement_functions/conftest.py for full provenance commentary.

FLEISS_1971_NIJ = np.array([
    [0, 0, 0, 6, 0],
    [0, 3, 0, 0, 3],
    [0, 1, 4, 0, 1],
    [0, 0, 0, 0, 6],
    [0, 3, 0, 3, 0],
    [2, 0, 4, 0, 0],
    [0, 0, 4, 0, 2],
    [2, 0, 3, 1, 0],
    [2, 0, 0, 4, 0],
    [0, 0, 0, 0, 6],
    [1, 0, 0, 5, 0],
    [1, 1, 0, 4, 0],
    [0, 3, 3, 0, 0],
    [1, 0, 0, 5, 0],
    [0, 2, 0, 3, 1],
    [0, 0, 5, 0, 1],
    [3, 0, 0, 1, 2],
    [5, 1, 0, 0, 0],
    [0, 2, 0, 4, 0],
    [1, 0, 2, 0, 3],
    [0, 0, 0, 0, 6],
    [0, 1, 0, 5, 0],
    [0, 2, 0, 1, 3],
    [2, 0, 0, 4, 0],
    [1, 0, 0, 4, 1],
    [0, 5, 0, 1, 0],
    [4, 0, 0, 0, 2],
    [0, 2, 0, 4, 0],
    [1, 0, 5, 0, 0],
    [0, 0, 0, 0, 6],
], dtype=np.int32)

TWO_RATER_INDEX = np.array([
    [0, 0],  # agree on 0
    [0, 0],  # agree on 0
    [0, 0],  # agree on 0
    [0, 1],  # disagree
    [1, 1],  # agree on 1
    [1, 1],  # agree on 1
    [1, 1],  # agree on 1
    [1, 1],  # agree on 1
    [1, 0],  # disagree
    [1, 0],  # disagree
], dtype=INDEX_DTYPE)

TWO_RATER_K = 2


# -- Expected values ------------------------------------------------------- #
# Exact fractions where possible. See test_agreement_functions/conftest.py
# for the arithmetic derivations.

FLEISS_1971_EXPECTED = {
    "fleiss_kappa":        0.430244520060141,
    "scotts_pi":           0.430244520060141,
    "gwet_ac1":            0.447884515844564,
    "brennan_prediger":    4.0 / 9.0,
    "krippendorff_alpha":  0.433409828282029,
}

TWO_RATER_EXPECTED = {
    "cohen_kappa":         2.0 / 5.0,
    "conger_kappa":        2.0 / 5.0,
    "fleiss_kappa":        13.0 / 33.0,
    "scotts_pi":           13.0 / 33.0,
    "gwet_ac1":            41.0 / 101.0,
    "brennan_prediger":    2.0 / 5.0,
    "krippendorff_alpha":  14.0 / 33.0,
}


# -- Helper: expand nij matrix to category_index array --------------------- #

def _nij_to_category_index(nij: np.ndarray) -> np.ndarray:
    """
    Expand an nij matrix (n_subjects x n_categories) into a dense
    category_index array (n_subjects x n_raters_per_subject).

    Each row of nij gives category counts that sum to n (the number of
    raters for that subject). I expand each count into that many copies
    of the category index. For example, [0, 3, 0, 3, 0] becomes
    [1, 1, 1, 3, 3, 3].

    This is deterministic -- categories are laid out in ascending order
    within each subject row.
    """
    n_subjects = nij.shape[0]
    n_raters = int(nij[0].sum())  # assumes balanced
    result = np.empty((n_subjects, n_raters), dtype=INDEX_DTYPE)

    for i in range(n_subjects):
        col = 0
        for cat_idx, count in enumerate(nij[i]):
            for _ in range(count):
                result[i, col] = cat_idx
                col += 1
    return result


# -- Semantics ------------------------------------------------------------- #

TWO_CAT_NAMES = ("cat0", "cat1")
FIVE_CAT_NAMES = ("Depression", "Personality_disorder", "Schizophrenia",
                   "Neurosis", "Other")

TWO_CAT = UnorderedCategories(TWO_CAT_NAMES)
FIVE_CAT = UnorderedCategories(FIVE_CAT_NAMES)

SEM_2RATER_2CAT = RatingsSemantics(raters=2, categories=TWO_CAT)
SEM_6RATER_5CAT = RatingsSemantics(raters=6, categories=FIVE_CAT)
SEM_3RATER_2CAT = RatingsSemantics(raters=3, categories=TWO_CAT)


# -- Expanded Fleiss 1971 array -------------------------------------------- #

FLEISS_1971_INDEX = _nij_to_category_index(FLEISS_1971_NIJ)


# -- Three-rater dataset: 10 subjects, 3 raters, K=2 ---------------------- #
# Hand-picked so Conger and Light have meaningful multi-rater data.
#
#   subject | r0 | r1 | r2
#   --------+----+----+----
#        0  |  0 |  0 |  0    unanimous 0
#        1  |  1 |  1 |  1    unanimous 1
#        2  |  0 |  0 |  1    majority 0
#        3  |  1 |  1 |  0    majority 1
#        4  |  0 |  1 |  0    majority 0
#        5  |  1 |  0 |  1    majority 1
#        6  |  0 |  0 |  0    unanimous 0
#        7  |  1 |  1 |  1    unanimous 1
#        8  |  0 |  1 |  1    majority 1
#        9  |  1 |  0 |  0    majority 0

THREE_RATER_INDEX = np.array([
    [0, 0, 0],
    [1, 1, 1],
    [0, 0, 1],
    [1, 1, 0],
    [0, 1, 0],
    [1, 0, 1],
    [0, 0, 0],
    [1, 1, 1],
    [0, 1, 1],
    [1, 0, 0],
], dtype=INDEX_DTYPE)


# -- Reusable config instances --------------------------------------------- #

# A valid HypothesisTestConfig with defaults, used by AnalyzerConfig tests later.
HYPOTHESIS_TEST_DEFAULTS = HypothesisTestConfig()

# A valid ConvergenceTest stub, used by AnalyzerConfig tests later.
CONVERGENCE_TEST_STUB = ConvergenceTest()

# A small identity weight matrix for CustomWeightsInput tests.
IDENTITY_3x3 = np.eye(3)


# -- Fixtures: Ratings objects --------------------------------------------- #

@pytest.fixture
def two_rater_ratings():
    """Ratings from TWO_RATER_INDEX: 10 subjects, 2 raters, K=2."""
    return Ratings(TWO_RATER_INDEX.copy(), SEM_2RATER_2CAT)


@pytest.fixture
def fleiss_1971_ratings():
    """Ratings expanded from Fleiss (1971) Table 1: 30 subjects, 6 raters, K=5."""
    return Ratings(FLEISS_1971_INDEX.copy(), SEM_6RATER_5CAT)


@pytest.fixture
def three_rater_ratings():
    """Three-rater dataset: 10 subjects, 3 raters, K=2."""
    return Ratings(THREE_RATER_INDEX.copy(), SEM_3RATER_2CAT)


# -- Fixtures: NijMatrix objects ------------------------------------------- #

@pytest.fixture
def fleiss_1971_nij_matrix():
    """Fleiss (1971) Table 1 as NijMatrix: 30 subjects x 5 categories."""
    return NijMatrix(FLEISS_1971_NIJ.copy())


@pytest.fixture
def two_rater_nij_matrix():
    """TWO_RATER_INDEX converted to NijMatrix: 10 subjects x 2 categories."""
    from interrater.intermediates.estimator_parts.estimator_input import EstimatorInput
    ei = EstimatorInput(TWO_RATER_INDEX.copy(), TWO_RATER_K)
    return ei.nij


# -- Fixtures: minimal configs (deterministic, for point-estimate tests) --- #
# These use the smallest viable bootstrap or parametric_null CI so that
# estimate() runs fast and deterministically. We only check .point in
# analyzer tests -- CI bounds belong to test_ci_methods.py.

@pytest.fixture
def fleiss_config_parametric():
    """FleissKappaConfig with ci_method='parametric_null' -- no bootstrap."""
    return FleissKappaConfig(ci_method="parametric_null")


@pytest.fixture
def scotts_config_parametric():
    """ScottsConfig with ci_method='parametric_null' -- no bootstrap."""
    return ScottsConfig(ci_method="parametric_null")


@pytest.fixture
def cohen_config_parametric():
    """CohenKappaConfig with ci_method='parametric_null' -- no bootstrap."""
    return CohenKappaConfig(ci_method="parametric_null")


@pytest.fixture
def gwet_config_fast():
    """GwetAC1Config with tiny bootstrap for speed. Deterministic via seed."""
    return GwetAC1Config(n_resamples=50, random_state=42)


@pytest.fixture
def bp_config_fast():
    """BrennanPredigerConfig with tiny bootstrap for speed."""
    return BrennanPredigerConfig(n_resamples=50, random_state=42)


@pytest.fixture
def krippendorff_config_fast():
    """KrippendorffAlphaConfig with tiny bootstrap for speed."""
    return KrippendorffAlphaConfig(n_resamples=50, random_state=42)


@pytest.fixture
def conger_config_fast():
    """CongerKappaConfig with tiny bootstrap for speed."""
    return CongerKappaConfig(n_resamples=50, random_state=42)


@pytest.fixture
def light_config_fast():
    """LightKappaConfig with tiny bootstrap for speed."""
    return LightKappaConfig(n_resamples=50, random_state=42)


# -- Fixtures: PairwisePermutationConfig instances ------------------------- #
# Tiny bootstrap sizes for speed. Deterministic via seed.

@pytest.fixture
def pairwise_cohen_config_fast():
    """PairwisePermutationConfig: cohen_kappa, tiny per-pair bootstrap, no mean bootstrap."""
    return PairwisePermutationConfig(
        coefficient_to_compute="cohen_kappa",
        config_of_coefficient=CohenKappaConfig(
            n_resamples=50, random_state=42
        ),
    )


@pytest.fixture
def pairwise_cohen_config_correlated_mean():
    """PairwisePermutationConfig: cohen_kappa, per-pair + correlated mean bootstrap."""
    return PairwisePermutationConfig(
        coefficient_to_compute="cohen_kappa",
        config_of_coefficient=CohenKappaConfig(
            n_resamples=50, random_state=42
        ),
        bootstrap_mean_config=PermutationBootstrapMeanConfig(
            method="correlated",
            n_bootstraps=30,
            random_state=99,
        ),
    )


@pytest.fixture
def pairwise_cohen_config_independent_mean():
    """PairwisePermutationConfig: cohen_kappa, per-pair + independent mean bootstrap."""
    return PairwisePermutationConfig(
        coefficient_to_compute="cohen_kappa",
        config_of_coefficient=CohenKappaConfig(
            n_resamples=50, random_state=42
        ),
        bootstrap_mean_config=PermutationBootstrapMeanConfig(
            method="independent",
            n_bootstraps=30,
            random_state=99,
        ),
    )


# -- Fixtures: SubsetPermutationConfig instances ----------------------------- #
# Use fleiss_kappa (group coefficient, naturally supports k>=2).
# subset_k=2 with 3 raters gives C(3,2)=3 subsets -- small, fast, meaningful.

@pytest.fixture
def subset_fleiss_config_fast():
    """SubsetPermutationConfig: fleiss_kappa, k=2, tiny per-subset bootstrap, no mean bootstrap."""
    return SubsetPermutationConfig(
        coefficient_to_compute="fleiss_kappa",
        config_of_coefficient=FleissKappaConfig(
            n_resamples=50, random_state=42
        ),
        subset_k=2,
    )


@pytest.fixture
def subset_fleiss_config_correlated_mean():
    """SubsetPermutationConfig: fleiss_kappa, k=2, per-subset + correlated mean bootstrap."""
    return SubsetPermutationConfig(
        coefficient_to_compute="fleiss_kappa",
        config_of_coefficient=FleissKappaConfig(
            n_resamples=50, random_state=42
        ),
        subset_k=2,
        bootstrap_mean_config=PermutationBootstrapMeanConfig(
            method="correlated",
            n_bootstraps=30,
            random_state=99,
        ),
    )


@pytest.fixture
def subset_fleiss_config_independent_mean():
    """SubsetPermutationConfig: fleiss_kappa, k=2, per-subset + independent mean bootstrap."""
    return SubsetPermutationConfig(
        coefficient_to_compute="fleiss_kappa",
        config_of_coefficient=FleissKappaConfig(
            n_resamples=50, random_state=42
        ),
        subset_k=2,
        bootstrap_mean_config=PermutationBootstrapMeanConfig(
            method="independent",
            n_bootstraps=30,
            random_state=99,
        ),
    )
