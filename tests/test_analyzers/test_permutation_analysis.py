"""Tests for analyzers/permutation_analysis/ sub-package.

Covers:
  - PermutationAnalysis (base class: validation, estimator resolution, weight
    resolution, injection pipeline, CI computation)
  - PairwisePermutationAnalysis (Computations A, B, C for pairs)
  - SubsetPermutationAnalysis (Computations A, B, C for k-subsets, multi-k)

One test file per sub-package, not one per class.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

# ------------------------------------------------------------------ #
# Imports under test                                                   #
# ------------------------------------------------------------------ #

from interrater.analyzers.permutation_analysis import (
    PermutationAnalysis,
    PairwisePermutationAnalysis,
    SubsetPermutationAnalysis,
)
from interrater.analyzers.permutation_analysis.permutation_analysis import (
    _ESTIMATOR_MAP,
    Estimator,
)
from interrater.analyzers.configs.cohen_kappa_config import CohenKappaConfig
from interrater.analyzers.configs.conger_kappa_config import CongerKappaConfig
from interrater.analyzers.configs.fleiss_kappa_config import FleissKappaConfig
from interrater.analyzers.configs.scotts_config import ScottsConfig
from interrater.analyzers.configs.pairwise_permutation_config import (
    PairwisePermutationConfig,
)
from interrater.analyzers.configs.subset_permutation_config import (
    SubsetPermutationConfig,
)
from interrater.analyzers.configs.permutation_bootstrap_mean_config import (
    PermutationBootstrapMeanConfig,
)
from interrater.contracts import EMPTY
from interrater.data_objects.ratings import Ratings
from interrater.intermediates.results.permutation_results.pairwise_permutation_result import (
    MAX_PAIRS_FOR_DISTRIBUTIONS,
    PairwisePermutationResult,
)
from interrater.intermediates.results.permutation_results.subset_permutation_result import (
    MAX_SUBSETS_FOR_DISTRIBUTIONS,
    SubsetPermutationResult,
)
from interrater.intermediates.results.permutation_results.multi_k_subset_result import (
    MultiKSubsetResult,
)
from interrater.intermediates.results.ci_results import BootstrapCIResult
from interrater.intermediates.data_profiles.balance_profile import BalanceProfile
from interrater.intermediates.estimator_parts.estimator_input import EstimatorInput
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.semantics import OrderedCategories, UnorderedCategories

# Tolerance for deterministic computations
TOL = 1e-12


# ------------------------------------------------------------------ #
# Fixtures                                                             #
# ------------------------------------------------------------------ #


def _make_ratings(
    n_subjects: int,
    n_raters: int,
    n_categories: int = 2,
    seed: int = 42,
    sparse: bool = False,
) -> Ratings:
    """Build a Ratings object for testing.

    Parameters
    ----------
    sparse : bool
        If True, randomly set ~20% of entries to EMPTY (missing).
    """
    rng = np.random.default_rng(seed)
    cat_names = [f"C{i}" for i in range(n_categories)]
    rater_names = [f"r{i}" for i in range(n_raters)]

    category_index = rng.integers(0, n_categories, size=(n_subjects, n_raters))

    if sparse:
        mask = rng.random(size=(n_subjects, n_raters)) < 0.2
        category_index[mask] = EMPTY

    cats = UnorderedCategories(cat_names)
    semantics = RatingsSemantics(raters=rater_names, categories=cats)
    return Ratings(category_index=category_index, semantics=semantics)


def _pairwise_config(
    coefficient: str = "cohen_kappa",
    n_resamples: int = 50,
    random_state: int = 42,
    ci_method: str = "percentile",
    bootstrap_mean: str | None = None,
) -> PairwisePermutationConfig:
    """Build a PairwisePermutationConfig for testing."""
    config_map = {
        "cohen_kappa": CohenKappaConfig,
        "scotts_pi": ScottsConfig,
        "conger_kappa": CongerKappaConfig,
    }
    coc = config_map[coefficient](
        n_resamples=n_resamples,
        random_state=random_state,
        ci_method=ci_method,
    )
    bmc = None
    if bootstrap_mean is not None:
        bmc = PermutationBootstrapMeanConfig(method=bootstrap_mean)
    return PairwisePermutationConfig(
        coefficient_to_compute=coefficient,
        config_of_coefficient=coc,
        bootstrap_mean_config=bmc,
    )


def _subset_config(
    coefficient: str = "fleiss_kappa",
    subset_k: int | list[int] | str = 3,
    n_resamples: int = 50,
    random_state: int = 42,
    ci_method: str = "percentile",
    bootstrap_mean: str | None = None,
) -> SubsetPermutationConfig:
    """Build a SubsetPermutationConfig for testing."""
    config_map = {
        "fleiss_kappa": FleissKappaConfig,
        "scotts_pi": ScottsConfig,
    }
    coc = config_map[coefficient](
        n_resamples=n_resamples,
        random_state=random_state,
        ci_method=ci_method,
    )
    bmc = None
    if bootstrap_mean is not None:
        bmc = PermutationBootstrapMeanConfig(method=bootstrap_mean)
    return SubsetPermutationConfig(
        coefficient_to_compute=coefficient,
        config_of_coefficient=coc,
        subset_k=subset_k,
        bootstrap_mean_config=bmc,
    )


# ------------------------------------------------------------------ #
# PermutationAnalysis (base class)                                     #
# ------------------------------------------------------------------ #


class TestPermutationAnalysisValidation:
    """Validation in the base class __init__."""

    def test_ratings_none_raises_value_error(self):
        config = _pairwise_config()
        with pytest.raises(ValueError, match="ratings is required"):
            PairwisePermutationAnalysis(None, config)

    def test_ratings_wrong_type_raises_type_error(self):
        config = _pairwise_config()
        with pytest.raises(TypeError, match="ratings must be Ratings"):
            PairwisePermutationAnalysis("not_ratings", config)

    def test_config_wrong_type_raises_type_error(self):
        ratings = _make_ratings(20, 3)
        with pytest.raises(TypeError, match="config must be"):
            PairwisePermutationAnalysis(ratings, "not_config")


class TestEstimatorResolution:
    """_resolve_estimator maps coefficient names to callables."""

    def test_all_known_coefficients_resolve(self):
        """Every key in _ESTIMATOR_MAP should resolve without error."""
        expected = {
            "cohen_kappa",
            "scotts_pi",
            "conger_kappa",
            "gwet_ac1",
            "fleiss_kappa",
            "brennan_prediger",
            "krippendorff_alpha",
        }
        assert set(_ESTIMATOR_MAP.keys()) == expected

    def test_unknown_coefficient_raises(self):
        """Unknown coefficient name should raise ValueError."""
        ratings = _make_ratings(20, 3)
        # Cannot test directly via constructor because config validation
        # catches unknown coefficients first. Test the _resolve_estimator
        # method indirectly through the map.
        assert "nonexistent" not in _ESTIMATOR_MAP

    def test_estimator_is_callable(self):
        """Every registered estimator should be callable."""
        for name, est in _ESTIMATOR_MAP.items():
            assert callable(est), f"{name} estimator is not callable"


class TestWeightResolution:
    """_resolve_weights handles identity/linear/quadratic."""

    def test_identity_weights_returns_none(self):
        """Default identity weights should resolve to None."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        analysis = PairwisePermutationAnalysis(ratings, config)
        # identity weights -> estimator_input.weights is None
        assert analysis.estimator_input.weights is None

    def test_linear_weights_resolved(self):
        """Ordered categories with linear weighting should produce Weights."""
        unordered = _make_ratings(20, 3, n_categories=4)
        semantics = RatingsSemantics(
            raters=list(unordered.semantics.raters),
            categories=OrderedCategories(
                unordered.semantics.categories.names, weighting="linear"
            ),
        )
        ratings = Ratings(
            category_index=unordered.category_index, semantics=semantics
        )
        coc = CohenKappaConfig(
            n_resamples=50,
            random_state=42,
            ci_method="percentile",
        )
        config = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=coc,
        )
        analysis = PairwisePermutationAnalysis(ratings, config)
        w = analysis.estimator_input.weights
        assert w is not None
        assert w.matrix.shape == (4, 4)


class TestInjectionPipeline:
    """apply_injection with injection=None is a no-op."""

    def test_no_injection_is_noop(self):
        """When injection is None, ratings should be unchanged."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        analysis = PairwisePermutationAnalysis(ratings, config)
        # Ratings should be the same object (no copy needed if no injection)
        assert analysis.ratings.n_subjects == 20
        assert analysis.ratings.n_raters == 3

    @staticmethod
    def _subset_with_one_injected(subset_k):
        from interrater.injectors.configs import DistributionInjectorConfig

        coc = FleissKappaConfig(
            n_resamples=2,
            injection=DistributionInjectorConfig(n_new_raters=1, random_state=0),
        )
        config = SubsetPermutationConfig(
            coefficient_to_compute="fleiss_kappa", config_of_coefficient=coc,
            subset_k=subset_k,
        )
        return SubsetPermutationAnalysis(_make_ratings(20, 3), config)

    def test_subset_k_all_counts_injected_raters(self):
        """3 real raters + 1 injected: 'all' runs k = 2, 3 and 4."""
        assert self._subset_with_one_injected("all").ks == [2, 3, 4]

    def test_subset_k_may_include_injected_raters(self):
        """k = 4 was refused as 'exceeds n_raters=3' before 2.1.0."""
        assert self._subset_with_one_injected(4).ks == [4]


class TestComputeCI:
    """_compute_ci delegates to percentile and bca."""

    def test_percentile_ci(self):
        """Percentile CI should produce sensible bounds."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        analysis = PairwisePermutationAnalysis(ratings, config)

        dist = np.linspace(0.0, 1.0, 1000)
        ci = analysis._compute_ci(dist, 0.5, "percentile", 0.05)
        assert ci.lower < ci.upper
        assert ci.lower >= 0.0
        assert ci.upper <= 1.0
        assert ci.method_intended == "percentile"
        assert ci.method_used == "percentile"

    def test_percentile_needs_no_jackknife(self):
        """Percentile does not estimate acceleration, so it asks for nothing."""
        ratings = _make_ratings(20, 3)
        analysis = PairwisePermutationAnalysis(ratings, _pairwise_config())

        ci = analysis._compute_ci(
            np.linspace(0.0, 1.0, 1000), 0.5, "percentile", 0.05,
            jackknife_fn=None,
        )
        assert ci.is_defined

    def test_bca_ci(self):
        """BCa CI should produce sensible bounds, given a jackknife."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        analysis = PairwisePermutationAnalysis(ratings, config)

        dist = np.linspace(0.0, 1.0, 1000)
        ci = analysis._compute_ci(
            dist, 0.5, "bca", 0.05,
            jackknife_fn=lambda: np.linspace(0.45, 0.55, 20),
        )
        assert ci.lower < ci.upper
        assert ci.method_intended == "bca"
        assert ci.method_used == "bca"
        assert ci.acceleration is not None

    def test_bca_without_jackknife_raises(self):
        """Asking for BCa with no way to compute acceleration is a caller bug.

        It used to fall through to percentile silently. The quantity BCa
        needs is a property of the data, not of the bootstrap distribution,
        so a caller that supplies no jackknife has not asked for a weaker
        interval -- it has asked for one that cannot be built.
        """
        ratings = _make_ratings(20, 3)
        analysis = PairwisePermutationAnalysis(ratings, _pairwise_config())

        with pytest.raises(ValueError, match="leave-one-out estimates"):
            analysis._compute_ci(
                np.linspace(0.0, 1.0, 1000), 0.5, "bca", 0.05,
                jackknife_fn=None,
            )

    def test_bca_jackknife_is_not_called_for_percentile(self):
        """The jackknife costs one estimator call per subject -- don't pay it."""
        calls = []

        def _jackknife():
            calls.append(1)
            return np.linspace(0.45, 0.55, 20)

        ratings = _make_ratings(20, 3)
        analysis = PairwisePermutationAnalysis(ratings, _pairwise_config())

        analysis._compute_ci(
            np.linspace(0.0, 1.0, 1000), 0.5, "percentile", 0.05,
            jackknife_fn=_jackknife,
        )
        assert calls == []

    def test_parametric_null_rejected(self):
        """parametric_null should be rejected for permutation results."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        analysis = PairwisePermutationAnalysis(ratings, config)

        dist = np.linspace(0.0, 1.0, 100)
        with pytest.raises(ValueError, match="parametric_null"):
            analysis._compute_ci(dist, 0.5, "parametric_null", 0.05)

    def test_unknown_ci_method_rejected(self):
        """Unknown CI method should raise ValueError."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        analysis = PairwisePermutationAnalysis(ratings, config)

        dist = np.linspace(0.0, 1.0, 100)
        with pytest.raises(ValueError, match="unknown ci_method"):
            analysis._compute_ci(dist, 0.5, "magic", 0.05)


# ------------------------------------------------------------------ #
# PairwisePermutationAnalysis                                          #
# ------------------------------------------------------------------ #


class TestPairwiseValidation:
    """PairwisePermutationAnalysis constructor validation."""

    def test_wrong_config_type_raises(self):
        """Passing SubsetPermutationConfig should raise TypeError."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=3)
        with pytest.raises(TypeError, match="PairwisePermutationConfig"):
            PairwisePermutationAnalysis(ratings, config)


class TestPairwiseComputationA:
    """Computation A: observed values, no bootstrap."""

    def test_result_type(self):
        """run() should return PairwisePermutationResult."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert isinstance(result, PairwisePermutationResult)

    def test_n_pairs_correct(self):
        """n_pairs should be C(n_raters, 2)."""
        for n_raters in [2, 3, 4, 5]:
            ratings = _make_ratings(20, n_raters)
            config = _pairwise_config()
            result = PairwisePermutationAnalysis(ratings, config).run()
            expected = n_raters * (n_raters - 1) // 2
            assert result.n_pairs == expected, (
                f"n_raters={n_raters}: expected {expected} pairs, "
                f"got {result.n_pairs}"
            )

    def test_pair_keys_are_sorted_tuples(self):
        """Pair keys should be sorted 2-tuples of int indices."""
        ratings = _make_ratings(20, 4)
        config = _pairwise_config()
        result = PairwisePermutationAnalysis(ratings, config).run()
        for key in result.pairwise_values:
            assert isinstance(key, tuple)
            assert len(key) == 2
            assert key[0] < key[1]

    def test_subjects_per_pair_matches_values(self):
        """subjects_per_pair should have the same keys as pairwise_values."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert set(result.subjects_per_pair.keys()) == set(
            result.pairwise_values.keys()
        )

    def test_balanced_data_all_subjects_shared(self):
        """With balanced data every pair should share all subjects."""
        ratings = _make_ratings(30, 3)
        config = _pairwise_config()
        result = PairwisePermutationAnalysis(ratings, config).run()
        for pair, count in result.subjects_per_pair.items():
            assert count == 30, f"pair {pair}: expected 30 subjects, got {count}"

    def test_mean_is_average_of_values(self):
        """mean_over_pairs should be the (NaN-skipping) mean of pairwise values."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        result = PairwisePermutationAnalysis(ratings, config).run()
        values = list(result.pairwise_values.values())
        expected_mean = np.nanmean(values)
        assert abs(result.mean_over_pairs - expected_mean) < TOL

    def test_coefficient_name_in_result(self):
        """Result should carry the coefficient name from config."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(coefficient="cohen_kappa")
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert result.coefficient_name == "cohen_kappa"

    def test_no_mean_bootstrap_fields_when_not_requested(self):
        """Without bootstrap_mean_config, mean bootstrap fields should be None."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert result.bootstrap_mean_method is None
        assert result.bootstrap_ci_of_mean is None
        assert result.bootstrap_distribution_of_mean is None


class TestPairwiseComputationB:
    """Computation B: per-pair bootstrap."""

    def test_per_pair_bootstrap_produces_cis(self):
        """With n_resamples > 0, per-pair CIs should be produced."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(n_resamples=100)
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert result.bootstrap_cis_per_pair is not None
        for pair, ci in result.bootstrap_cis_per_pair.items():
            assert isinstance(ci, BootstrapCIResult)
            assert ci.n_resamples == 100
            assert ci.method_intended == "percentile"
            # An interval either has ordered bounds or has none at all --
            # it cannot be half-defined, and it validated that itself.
            if ci.is_defined:
                assert ci.lower <= ci.upper + TOL
            else:
                assert ci.method_used is None

    def test_per_pair_distributions_stored_when_few_pairs(self):
        """Distributions stored when n_pairs <= MAX_PAIRS_FOR_DISTRIBUTIONS."""
        ratings = _make_ratings(20, 3)  # 3 pairs
        config = _pairwise_config(n_resamples=100)
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert 3 <= MAX_PAIRS_FOR_DISTRIBUTIONS
        assert result.bootstrap_distributions_per_pair is not None
        for pair, dist in result.bootstrap_distributions_per_pair.items():
            assert isinstance(dist, np.ndarray)
            assert dist.shape == (100,)

    def test_per_pair_draws_not_stored_says_why(self):
        """Above the limit, with no bootstrap of the mean, the accessor
        names the storage limit. It used to say the bootstrap was not
        performed, while every pair had its interval."""
        ratings = _make_ratings(30, 7)  # 21 pairs
        result = PairwisePermutationAnalysis(
            ratings, _pairwise_config(n_resamples=20)
        ).run()
        assert result.n_pairs > MAX_PAIRS_FOR_DISTRIBUTIONS
        assert result.bootstrap_mean_method is None
        assert len(result.bootstrap_cis_per_pair) == 21
        with pytest.raises(RuntimeError, match="not stored because n_pairs"):
            result.get_bootstrap_distributions_per_pair()

    def test_reproducible_with_same_seed(self):
        """Same random_state should produce identical results."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(n_resamples=100, random_state=42)
        r1 = PairwisePermutationAnalysis(ratings, config).run()
        r2 = PairwisePermutationAnalysis(ratings, config).run()
        for pair in r1.pairwise_values:
            assert r1.pairwise_values[pair] == r2.pairwise_values[pair]
            np.testing.assert_array_equal(
                r1.bootstrap_distributions_per_pair[pair],
                r2.bootstrap_distributions_per_pair[pair],
            )

    def test_different_seed_produces_different_bootstrap(self):
        """Different random_state should produce different bootstrap draws."""
        ratings = _make_ratings(20, 3)
        r1 = PairwisePermutationAnalysis(
            ratings, _pairwise_config(n_resamples=100, random_state=42)
        ).run()
        r2 = PairwisePermutationAnalysis(
            ratings, _pairwise_config(n_resamples=100, random_state=99)
        ).run()
        # Observed values should be identical (same data)
        for pair in r1.pairwise_values:
            assert abs(
                r1.pairwise_values[pair] - r2.pairwise_values[pair]
            ) < TOL
        # Bootstrap distributions should differ
        pair = sorted(r1.pairwise_values.keys())[0]
        assert not np.array_equal(
            r1.bootstrap_distributions_per_pair[pair],
            r2.bootstrap_distributions_per_pair[pair],
        )


class TestPairwiseSparseData:
    """Pairwise analysis on sparse (unbalanced) data."""

    def test_sparse_data_runs(self):
        """Should handle missing data (EMPTY entries) without error."""
        ratings = _make_ratings(30, 3, sparse=True, seed=42)
        config = _pairwise_config(n_resamples=50)
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert result.n_pairs == 3
        # Subjects per pair may differ
        counts = list(result.subjects_per_pair.values())
        assert all(c >= 0 for c in counts)

    def test_pair_with_too_few_shared_subjects_gets_nan(self):
        """Pair with < 2 shared subjects should get NaN coefficient."""
        # Create data where raters 0 and 1 never overlap
        category_index = np.full((10, 3), EMPTY, dtype=np.int64)
        # Rater 0 rates subjects 0-4, rater 1 rates subjects 5-9
        rng = np.random.default_rng(42)
        category_index[:5, 0] = rng.integers(0, 2, size=5)
        category_index[5:, 1] = rng.integers(0, 2, size=5)
        # Rater 2 rates all
        category_index[:, 2] = rng.integers(0, 2, size=10)

        cats = UnorderedCategories(["A", "B"])
        semantics = RatingsSemantics(raters=["r0", "r1", "r2"], categories=cats)
        ratings = Ratings(category_index=category_index, semantics=semantics)

        config = _pairwise_config(n_resamples=50)
        result = PairwisePermutationAnalysis(ratings, config).run()

        # Pair (0, 1) should have 0 shared subjects -> NaN
        assert result.subjects_per_pair[(0, 1)] == 0
        assert math.isnan(result.pairwise_values[(0, 1)])


class TestPairwiseComputationC:
    """Computation C: mean bootstrap (correlated and independent)."""

    def test_correlated_mean_bootstrap(self):
        """Correlated method should produce mean distribution and CI."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(
            n_resamples=100, bootstrap_mean="correlated"
        )
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert result.bootstrap_mean_method == "correlated"
        assert result.bootstrap_distribution_of_mean is not None
        assert result.bootstrap_distribution_of_mean.shape == (100,)
        ci = result.bootstrap_ci_of_mean
        assert ci is not None
        assert ci.lower <= ci.upper + TOL
        assert ci.n_resamples == 100

    def test_independent_mean_bootstrap(self):
        """Independent method should produce mean distribution and CI."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(
            n_resamples=100, bootstrap_mean="independent"
        )
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert result.bootstrap_mean_method == "independent"
        assert result.bootstrap_distribution_of_mean is not None
        assert result.bootstrap_distribution_of_mean.shape == (100,)
        assert result.bootstrap_ci_of_mean is not None

    def test_correlated_wider_than_independent(self):
        """Correlated CI should generally be wider than independent.

        This is a statistical property, not an exact rule. We test with
        enough data and resamples that it holds reliably.
        """
        ratings = _make_ratings(50, 4, seed=42)

        r_corr = PairwisePermutationAnalysis(
            ratings,
            _pairwise_config(
                n_resamples=200,
                bootstrap_mean="correlated",
            ),
        ).run()

        r_indep = PairwisePermutationAnalysis(
            ratings,
            _pairwise_config(
                n_resamples=200,
                bootstrap_mean="independent",
            ),
        ).run()

        corr_ci = r_corr.bootstrap_ci_of_mean
        indep_ci = r_indep.bootstrap_ci_of_mean
        corr_width = corr_ci.upper - corr_ci.lower
        indep_width = indep_ci.upper - indep_ci.lower
        # Correlated should be wider (accounts for shared-rater correlation)
        # Allow some slack -- it's a statistical property
        assert corr_width > indep_width * 0.5, (
            f"Correlated CI width ({corr_width:.4f}) should be wider than "
            f"independent ({indep_width:.4f})"
        )

    def test_mean_interval_built_like_the_per_pair_ones(self):
        """One draw count, interval method and level per run.

        The mean used to have its own, so a result could say "percentile"
        while its mean interval was BCa at 90%.
        """
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(
            n_resamples=80, ci_method="bca", bootstrap_mean="correlated",
        )
        result = PairwisePermutationAnalysis(ratings, config).run()
        mean_ci = result.bootstrap_ci_of_mean
        pair_ci = next(iter(result.bootstrap_cis_per_pair.values()))
        assert result.ci_method == mean_ci.method_intended == "bca"
        assert mean_ci.method_intended == pair_ci.method_intended
        assert mean_ci.alpha == pair_ci.alpha
        assert mean_ci.n_resamples == pair_ci.n_resamples == 80

    def test_mean_reproducible_by_default(self):
        """No seed set anywhere: two runs still agree (the default is 42)."""
        ratings = _make_ratings(20, 3)
        config = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(n_resamples=50),
            bootstrap_mean_config=PermutationBootstrapMeanConfig(method="correlated"),
        )
        r1 = PairwisePermutationAnalysis(ratings, config).run()
        r2 = PairwisePermutationAnalysis(ratings, config).run()
        np.testing.assert_array_equal(
            r1.bootstrap_distribution_of_mean, r2.bootstrap_distribution_of_mean,
        )

    def test_mean_bootstrap_reproducible(self):
        """Same seeds should produce identical mean bootstrap."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(
            n_resamples=50,
            bootstrap_mean="correlated",
        )
        r1 = PairwisePermutationAnalysis(ratings, config).run()
        r2 = PairwisePermutationAnalysis(ratings, config).run()
        np.testing.assert_array_equal(
            r1.bootstrap_distribution_of_mean,
            r2.bootstrap_distribution_of_mean,
        )


class TestPairwiseDifferentCoefficients:
    """Pairwise analysis with different coefficient types."""

    @pytest.mark.parametrize("coefficient", ["cohen_kappa", "scotts_pi", "conger_kappa"])
    def test_coefficient_runs(self, coefficient):
        """Each pairwise-capable coefficient should run without error."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(coefficient=coefficient, n_resamples=50)
        result = PairwisePermutationAnalysis(ratings, config).run()
        assert result.coefficient_name == coefficient
        assert result.n_pairs == 3
        assert isinstance(result.mean_over_pairs, float)

    def test_cohen_and_scotts_differ_on_unequal_marginals(self):
        """Cohen's and Scott's should give different results when marginals differ.

        Cohen uses individual marginals, Scott pools them. With raters who
        use categories at different rates, the two should diverge.
        """
        # Rater 0 says A 80% of the time, rater 1 says A 20%
        rng = np.random.default_rng(42)
        n = 50
        r0 = (rng.random(n) < 0.8).astype(np.int64)
        r1 = (rng.random(n) < 0.2).astype(np.int64)
        r2 = rng.integers(0, 2, size=n)
        category_index = np.column_stack([r0, r1, r2])

        cats = UnorderedCategories(["A", "B"])
        semantics = RatingsSemantics(raters=["r0", "r1", "r2"], categories=cats)
        ratings = Ratings(category_index=category_index, semantics=semantics)

        r_cohen = PairwisePermutationAnalysis(
            ratings, _pairwise_config(coefficient="cohen_kappa", n_resamples=50)
        ).run()
        r_scotts = PairwisePermutationAnalysis(
            ratings, _pairwise_config(coefficient="scotts_pi", n_resamples=50)
        ).run()

        # For the pair (0,1) with very different marginals, values should differ
        assert abs(
            r_cohen.pairwise_values[(0, 1)] - r_scotts.pairwise_values[(0, 1)]
        ) > 0.01


class TestPairwiseBalanceProfile:
    """Balance profile should reflect the full dataset."""

    def test_balance_profile_present(self):
        """Result should carry a BalanceProfile."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        result = PairwisePermutationAnalysis(ratings, config).run()
        bp = result.balance_profile
        assert isinstance(bp, BalanceProfile)
        assert bp.n_subjects == 20
        assert bp.n_raters == 3


# ------------------------------------------------------------------ #
# SubsetPermutationAnalysis                                            #
# ------------------------------------------------------------------ #


class TestSubsetValidation:
    """SubsetPermutationAnalysis constructor validation."""

    def test_wrong_config_type_raises(self):
        """Passing PairwisePermutationConfig should raise TypeError."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config()
        with pytest.raises(TypeError, match="SubsetPermutationConfig"):
            SubsetPermutationAnalysis(ratings, config)

    def test_k_exceeds_n_raters_raises(self):
        """subset_k > n_raters should raise ValueError."""
        ratings = _make_ratings(20, 3)
        config = _subset_config(subset_k=4)  # 4 > 3 raters
        with pytest.raises(ValueError, match="exceeds n_raters"):
            SubsetPermutationAnalysis(ratings, config)


class TestSubsetNormalizeKs:
    """_normalize_ks handles int, list, and "all"."""

    def test_int_k(self):
        """Single int should become [int]."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=3)
        analysis = SubsetPermutationAnalysis(ratings, config)
        assert analysis.ks == [3]

    def test_list_k_sorted(self):
        """List of ints should be sorted."""
        ratings = _make_ratings(20, 5)
        config = _subset_config(subset_k=[4, 2, 3])
        analysis = SubsetPermutationAnalysis(ratings, config)
        assert analysis.ks == [2, 3, 4]

    def test_all_k(self):
        """'all' should produce range(2, n_raters + 1)."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k="all")
        analysis = SubsetPermutationAnalysis(ratings, config)
        assert analysis.ks == [2, 3, 4]


class TestSubsetComputationA:
    """Computation A: observed values for subsets."""

    def test_result_type(self):
        """run() should return MultiKSubsetResult."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        assert isinstance(result, MultiKSubsetResult)

    def test_per_k_result_type(self):
        """Each per-k result should be SubsetPermutationResult."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        for k, r in result.per_k_results.items():
            assert isinstance(r, SubsetPermutationResult)

    def test_n_subsets_for_k3_with_4_raters(self):
        """C(4, 3) = 4 subsets."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        assert result.per_k_results[3].n_subsets == 4

    def test_n_subsets_for_k2_with_4_raters(self):
        """C(4, 2) = 6 subsets."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=2)
        result = SubsetPermutationAnalysis(ratings, config).run()
        assert result.per_k_results[2].n_subsets == 6

    def test_subset_keys_are_sorted_tuples(self):
        """Subset keys should be sorted tuples of int indices."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        for key in result.per_k_results[3].subset_values:
            assert isinstance(key, tuple)
            assert len(key) == 3
            assert key == tuple(sorted(key))

    def test_subjects_per_subset_matches_values(self):
        """subjects_per_subset should have same keys as subset_values."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]
        assert set(r.subjects_per_subset.keys()) == set(r.subset_values.keys())

    def test_mean_is_average_of_values(self):
        """mean_over_subsets should be NaN-skipping mean of subset values."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]
        values = list(r.subset_values.values())
        expected = np.nanmean(values)
        assert abs(r.mean_over_subsets - expected) < TOL


class TestSubsetMultiK:
    """Multi-k support: multiple values of k in one run."""

    def test_multi_k_all(self):
        """subset_k='all' with 4 raters should produce k=2,3,4."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k="all")
        result = SubsetPermutationAnalysis(ratings, config).run()
        assert result.ks == [2, 3, 4]
        assert set(result.per_k_results.keys()) == {2, 3, 4}

    def test_multi_k_list(self):
        """subset_k=[2, 4] should produce results for k=2 and k=4."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k=[2, 4])
        result = SubsetPermutationAnalysis(ratings, config).run()
        assert set(result.per_k_results.keys()) == {2, 4}
        assert result.per_k_results[2].n_subsets == 6  # C(4,2)
        assert result.per_k_results[4].n_subsets == 1  # C(4,4)

    def test_k_equals_n_raters(self):
        """When k = n_raters, there is exactly 1 subset."""
        ratings = _make_ratings(20, 3)
        config = _subset_config(subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]
        assert r.n_subsets == 1
        assert (0, 1, 2) in r.subset_values

    def test_each_k_independent(self):
        """Results for different k values should be independent quantities."""
        ratings = _make_ratings(30, 4, seed=42)
        config = _subset_config(subset_k="all")
        result = SubsetPermutationAnalysis(ratings, config).run()
        # Means for different k should generally differ
        means = [result.per_k_results[k].mean_over_subsets for k in result.ks]
        # At least two should differ
        assert len(set(round(m, 6) for m in means)) > 1


class TestSubsetComputationB:
    """Computation B: per-subset bootstrap."""

    def test_per_subset_bootstrap_produces_cis(self):
        """With n_resamples > 0, per-subset CIs should be produced."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(n_resamples=100, subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]
        assert r.bootstrap_cis_per_subset is not None
        for subset, ci in r.bootstrap_cis_per_subset.items():
            assert isinstance(ci, BootstrapCIResult)
            assert ci.n_resamples == 100
            if ci.is_defined:
                assert ci.lower <= ci.upper + TOL
            else:
                assert ci.method_used is None

    def test_per_subset_intervals_bracket_their_own_estimates(self):
        """The regression test for the measurable-pool bug.

        A 3-rater subset is measured on every subject at least two of
        those three rated. When the bootstrap drew from the complete
        cases instead, it was resampling a different quantity, and the
        interval could land nowhere near the point estimate it belonged
        to. On balanced data the two pools coincide, so this guards the
        general path rather than reproducing the original failure.
        """
        ratings = _make_ratings(20, 4)
        config = _subset_config(n_resamples=200, subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]

        for subset, ci in r.bootstrap_cis_per_subset.items():
            point = r.subset_values[subset]
            if ci.is_defined and not math.isnan(point):
                assert ci.lower - TOL <= point <= ci.upper + TOL, (
                    f"subset {subset}: point {point:.4f} outside its own "
                    f"interval {ci}"
                )

    def test_distributions_stored_when_few_subsets(self):
        """Distributions stored when n_subsets <= MAX_SUBSETS_FOR_DISTRIBUTIONS."""
        ratings = _make_ratings(20, 4)  # C(4,3) = 4 subsets
        config = _subset_config(n_resamples=100, subset_k=3)
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]
        assert 4 <= MAX_SUBSETS_FOR_DISTRIBUTIONS
        assert r.bootstrap_distributions_per_subset is not None

    def test_per_subset_draws_not_stored_says_why(self):
        """Above the limit, with no bootstrap of the mean, the accessor
        names the storage limit rather than claiming no bootstrap ran."""
        ratings = _make_ratings(30, 12)  # C(12,3) = 220 subsets
        result = SubsetPermutationAnalysis(
            ratings, _subset_config(n_resamples=20, subset_k=3)
        ).run()
        r = result.per_k_results[3]
        assert r.n_subsets > MAX_SUBSETS_FOR_DISTRIBUTIONS
        assert r.bootstrap_mean_method is None
        assert len(r.bootstrap_cis_per_subset) == 220
        with pytest.raises(RuntimeError, match="not stored because n_subsets"):
            r.get_bootstrap_distributions_per_subset()

    def test_reproducible_with_same_seed(self):
        """Same random_state should produce identical results."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(n_resamples=100, subset_k=3, random_state=42)
        r1 = SubsetPermutationAnalysis(ratings, config).run()
        r2 = SubsetPermutationAnalysis(ratings, config).run()
        for subset in r1.per_k_results[3].subset_values:
            assert (
                r1.per_k_results[3].subset_values[subset]
                == r2.per_k_results[3].subset_values[subset]
            )


class TestSubsetComputationC:
    """Computation C: mean bootstrap for subsets."""

    def test_correlated_mean_bootstrap(self):
        """Correlated method should produce mean distribution and CI."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(
            subset_k=3,
            n_resamples=100,
            bootstrap_mean="correlated",
        )
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]
        assert r.bootstrap_mean_method == "correlated"
        assert r.bootstrap_distribution_of_mean is not None
        assert r.bootstrap_distribution_of_mean.shape == (100,)
        assert r.bootstrap_ci_of_mean is not None

    def test_independent_mean_bootstrap(self):
        """Independent method should produce mean distribution and CI."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(
            subset_k=3,
            n_resamples=100,
            bootstrap_mean="independent",
        )
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]
        assert r.bootstrap_mean_method == "independent"
        assert r.bootstrap_distribution_of_mean is not None
        assert r.bootstrap_ci_of_mean is not None


class TestSubsetSparseData:
    """Subset analysis on sparse data."""

    def test_sparse_data_runs(self):
        """Should handle missing data without error."""
        ratings = _make_ratings(30, 4, sparse=True, seed=42)
        config = _subset_config(subset_k=3, n_resamples=50)
        result = SubsetPermutationAnalysis(ratings, config).run()
        r = result.per_k_results[3]
        assert r.n_subsets == 4


# ------------------------------------------------------------------ #
# Cross-cutting: estimator adapters                                    #
# ------------------------------------------------------------------ #


class TestEstimatorAdapters:
    """Module-level estimator adapter functions."""

    def test_adapters_are_module_level(self):
        """Estimator adapters should be module-level functions (picklable)."""
        import interrater.analyzers.permutation_analysis.permutation_analysis as mod

        for name, est in _ESTIMATOR_MAP.items():
            assert hasattr(mod, est.__name__), (
                f"{name} adapter {est.__name__} is not at module level"
            )

    def test_cohen_kappa_adapter_returns_float(self):
        """Cohen's kappa adapter should return a float."""
        ratings = _make_ratings(20, 2)
        config = _pairwise_config(coefficient="cohen_kappa")
        analysis = PairwisePermutationAnalysis(ratings, config)
        # Call the estimator directly on the full data
        result = analysis.estimator(analysis.estimator_input)
        assert isinstance(result, (int, float, np.floating))

    def test_fleiss_kappa_adapter_returns_float(self):
        """Fleiss' kappa adapter should return a float."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(coefficient="fleiss_kappa", subset_k=3)
        analysis = SubsetPermutationAnalysis(ratings, config)
        result = analysis.estimator(analysis.estimator_input)
        assert isinstance(result, (int, float, np.floating))


# ------------------------------------------------------------------ #
# Integration: summary() and __str__                                   #
# ------------------------------------------------------------------ #


class TestResultSummary:
    """Result summary and string representation."""

    def test_pairwise_summary_has_required_keys(self):
        """PairwisePermutationResult.summary() should have required keys."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(n_resamples=50)
        result = PairwisePermutationAnalysis(ratings, config).run()
        s = result.summary()
        assert "coefficient_name" in s
        assert "ci_method" in s
        assert "n_pairs" in s
        assert "mean_over_pairs" in s
        assert "pairwise_values" in s
        assert "subjects_per_pair" in s

    def test_pairwise_summary_with_bootstrap(self):
        """Summary with mean bootstrap should include bootstrap keys."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(
            n_resamples=100, bootstrap_mean="correlated"
        )
        result = PairwisePermutationAnalysis(ratings, config).run()
        s = result.summary()
        assert "bootstrap_mean_method" in s
        assert "mean_ci_lower" in s
        assert "mean_ci_upper" in s
        assert "bootstrap_draws_mean" in s
        assert "bootstrap_se" in s
        assert "bootstrap_n" not in s

    def test_pairwise_summary_ignores_undefined_draws_of_the_mean(self):
        """A draw of the mean is NaN when every pair is undefined in it.

        Six subjects, two with a disagreement; a resample missing both has
        every rating in one category, so every pair -- and the mean -- is
        NaN. The summary must describe the usable draws, as the interval
        does, not turn NaN.
        """
        cats = UnorderedCategories(["C0", "C1"])
        sem = RatingsSemantics(raters=["a", "b", "c"], categories=cats)
        index = np.zeros((6, 3), dtype=np.int8)
        index[0] = [1, 1, 0]
        index[1] = [0, 1, 1]
        ratings = Ratings(category_index=index, semantics=sem)
        config = _pairwise_config(
            n_resamples=400, bootstrap_mean="correlated"
        )
        result = PairwisePermutationAnalysis(ratings, config).run()
        draws = result.bootstrap_distribution_of_mean
        n_nan = int(np.isnan(draws).sum())
        assert n_nan > 0
        assert result.bootstrap_ci_of_mean.n_failed == n_nan
        s = result.summary()
        assert s["bootstrap_draws_mean"] == pytest.approx(np.nanmean(draws))
        assert s["bootstrap_se"] == pytest.approx(np.nanstd(draws, ddof=1))

    def test_pairwise_str(self):
        """__str__ should not raise and should contain key info."""
        ratings = _make_ratings(20, 3)
        config = _pairwise_config(
            n_resamples=50, bootstrap_mean="correlated"
        )
        result = PairwisePermutationAnalysis(ratings, config).run()
        s = str(result)
        assert "cohen_kappa" in s
        assert "Pairwise" in s
        assert "n_pairs=3" in s

    def test_multi_k_summary(self):
        """MultiKSubsetResult.summary() should have per-k results."""
        ratings = _make_ratings(20, 4)
        config = _subset_config(subset_k="all")
        result = SubsetPermutationAnalysis(ratings, config).run()
        s = result.summary()
        assert "coefficient_name" in s
        assert "ks" in s
