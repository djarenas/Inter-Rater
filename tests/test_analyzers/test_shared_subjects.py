"""
A point estimate and its interval describe the same subjects.

The pairwise and subset bootstraps resample a pool: for a pair, the subjects
both raters rated; for a subset of k raters, the subjects at least two of
them rated. A point estimate computed from the raters' full columns must
therefore equal the same coefficient computed on that pool alone. If it does
not, the coefficient is using subjects the interval never sees, and the
interval is built around a different quantity from the number it reports.

Conger's kappa did exactly that before 2.1.0: its observed agreement skipped
subjects rated once, its chance agreement did not. In one 60-subject
example a pair's point sat 0.03 away from the centre of its own bootstrap.

The data is partially crossed on purpose -- 4 raters, about 40% of cells
empty -- so every pair and subset has subjects only one of its raters saw.
The coefficients come from the package's own lists in ``contracts.py``, so
a coefficient added there is tested here without being listed; one this
file cannot build a config for fails rather than being skipped.
"""

from itertools import combinations

import numpy as np
import pytest

import interrater.analyzers.configs as C
from interrater.analyzers.permutation_analysis.pairwise_permutation_analysis import (
    PairwisePermutationAnalysis,
)
from interrater.analyzers.permutation_analysis.subset_permutation_analysis import (
    SubsetPermutationAnalysis,
)
from interrater.contracts import (
    PAIRWISE_PERMUTATION_COEFFICIENTS,
    SUBSET_PERMUTATION_COEFFICIENTS,
)
from interrater.data_objects.ratings import Ratings
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics

CONFIGS = {
    "cohen_kappa": C.CohenKappaConfig,
    "scotts_pi": C.ScottsConfig,
    "conger_kappa": C.CongerKappaConfig,
    "gwet_ac1": C.GwetAC1Config,
    "fleiss_kappa": C.FleissKappaConfig,
    "brennan_prediger": C.BrennanPredigerConfig,
    "krippendorff_alpha": C.KrippendorffAlphaConfig,
}

_CATEGORIES = UnorderedCategories(
    names=("A", "B", "C"), labels={"A": "A", "B": "B", "C": "C"},
)


def _codes() -> np.ndarray:
    """80 subjects, 4 raters, 3 categories; 60% copy a true answer; ~40% of cells empty."""
    rng = np.random.default_rng(42)
    truth = rng.integers(0, 3, 80)
    codes = np.where(
        rng.random((80, 4)) < 0.6, truth[:, None], rng.integers(0, 3, (80, 4)),
    ).astype(np.int8)
    codes[rng.random((80, 4)) < 0.4] = -1
    return codes[(codes >= 0).sum(axis=1) >= 2]  # every subject compared somewhere


def _ratings(codes: np.ndarray) -> Ratings:
    return Ratings(
        category_index=codes,
        semantics=RatingsSemantics(raters=codes.shape[1], categories=_CATEGORIES),
    )


def _config(name: str):
    if name not in CONFIGS:
        pytest.fail(f"{name}: no config known here; add it to CONFIGS")
    # Two resamples: only the point estimate is under test.
    return CONFIGS[name](n_resamples=2, random_state=42, ci_method="percentile")


def _same(a: float, b: float) -> bool:
    return (np.isnan(a) and np.isnan(b)) or abs(a - b) < 1e-12


@pytest.mark.parametrize("name", sorted(PAIRWISE_PERMUTATION_COEFFICIENTS))
def test_pairwise_point_uses_only_shared_subjects(name):
    codes = _codes()

    def run(data):
        config = C.PairwisePermutationConfig(
            coefficient_to_compute=name, config_of_coefficient=_config(name),
        )
        return PairwisePermutationAnalysis(_ratings(data), config).run().pairwise_values

    full = run(codes)
    for a, b in combinations(range(codes.shape[1]), 2):
        pair = codes[:, [a, b]]
        shared = pair[(pair >= 0).all(axis=1)]
        alone = run(shared)[(0, 1)]
        assert _same(full[(a, b)], alone), (
            f"{name}, pair {(a, b)}: {full[(a, b)]:.6f} on full columns, "
            f"{alone:.6f} on the {len(shared)} shared subjects"
        )


@pytest.mark.parametrize("name", sorted(SUBSET_PERMUTATION_COEFFICIENTS))
def test_subset_point_uses_only_subjects_rated_twice(name):
    codes = _codes()

    def run(data):
        config = C.SubsetPermutationConfig(
            coefficient_to_compute=name, config_of_coefficient=_config(name),
            subset_k=3,
        )
        result = SubsetPermutationAnalysis(_ratings(data), config).run()
        return result.per_k_results[3].subset_values

    full = run(codes)
    for trio in combinations(range(codes.shape[1]), 3):
        sub = codes[:, list(trio)]
        pool = sub[(sub >= 0).sum(axis=1) >= 2]
        alone = run(pool)[(0, 1, 2)]
        assert _same(full[trio], alone), (
            f"{name}, raters {trio}: {full[trio]:.6f} on full columns, "
            f"{alone:.6f} on the {len(pool)} subjects rated at least twice"
        )


def test_conger_equals_cohen_on_a_pair():
    """Conger's kappa on two raters is Cohen's: same agreement, same marginals."""
    codes = _codes()

    def run(name):
        config = C.PairwisePermutationConfig(
            coefficient_to_compute=name, config_of_coefficient=_config(name),
        )
        return PairwisePermutationAnalysis(_ratings(codes), config).run().pairwise_values

    conger, cohen = run("conger_kappa"), run("cohen_kappa")
    for pair in cohen:
        assert _same(conger[pair], cohen[pair]), pair
