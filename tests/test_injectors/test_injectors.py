"""
Tests for injector classes and the registry.

Bottom-up: base RaterInjector validation, then each concrete injector,
then the registry, then chaining.

Injectors are tested with real Ratings objects from conftest. Every
generate() test seeds the RNG so results are deterministic.
"""

import numpy as np
import pytest

from interrater.contracts import INDEX_DTYPE, EMPTY
from interrater.data_objects.ratings import Ratings
from interrater.injectors import (
    RaterInjector,
    DistributionInjector,
    WithinSubjectInjector,
    PooledInjector,
    RaterBootstrapInjector,
    ExtremalInjector,
)
from interrater.injectors.configs import (
    InjectorConfig,
    DistributionInjectorConfig,
    WithinSubjectInjectorConfig,
    PooledInjectorConfig,
    RaterBootstrapInjectorConfig,
    ExtremalInjectorConfig,
)
from interrater.injectors.registry import (
    get_registry,
    hydrate_injection_step,
    get_injector_class,
)


# ============================================================================ #
#  RaterInjector (base ABC)                                                    #
# ============================================================================ #


class TestRaterInjectorValidation:
    """Base class input validation."""

    def test_none_ratings_rejected(self, base_config):
        with pytest.raises(ValueError, match="ratings cannot be None"):
            DistributionInjector(None, DistributionInjectorConfig())

    def test_non_ratings_rejected(self, base_config):
        with pytest.raises(TypeError, match="ratings must be Ratings"):
            DistributionInjector("not_ratings", DistributionInjectorConfig())

    def test_none_config_rejected(self, simple_ratings):
        with pytest.raises(ValueError, match="config cannot be None"):
            DistributionInjector(simple_ratings, None)

    def test_non_config_rejected(self, simple_ratings):
        with pytest.raises(TypeError, match="config must be InjectorConfig"):
            DistributionInjector(simple_ratings, "not_config")


class TestRaterInjectorNameBuilding:
    """Injected rater name generation."""

    def test_names_start_at_zero(self, simple_ratings):
        injector = DistributionInjector(
            simple_ratings, DistributionInjectorConfig(
                n_new_raters=2, random_state=0,
            )
        )
        result = injector.inject()
        raters = result.semantics.raters
        assert raters[-2] == "injected_rater0"
        assert raters[-1] == "injected_rater1"

    def test_chained_names_no_collision(self, simple_ratings):
        """Inject twice -- second batch starts where first left off."""
        inj1 = DistributionInjector(
            simple_ratings, DistributionInjectorConfig(
                n_new_raters=2, random_state=0,
            )
        )
        after_first = inj1.inject()

        inj2 = DistributionInjector(
            after_first, DistributionInjectorConfig(
                n_new_raters=1, random_state=1,
            )
        )
        after_second = inj2.inject()

        raters = after_second.semantics.raters
        # Original: alice, bob, carol
        # First injection: injected_rater0, injected_rater1
        # Second injection: injected_rater2
        assert raters[-1] == "injected_rater2"
        assert len(set(raters)) == len(raters)  # no duplicates


class TestRaterInjectorInject:
    """inject() output shape and semantics."""

    def test_output_shape(self, simple_ratings, distribution_config_uniform):
        injector = DistributionInjector(
            simple_ratings, distribution_config_uniform,
        )
        result = injector.inject()
        # Original: (4, 3). Adding 1 rater -> (4, 4)
        assert result.category_index.shape == (4, 4)

    def test_output_is_new_ratings(self, simple_ratings,
                                   distribution_config_uniform):
        injector = DistributionInjector(
            simple_ratings, distribution_config_uniform,
        )
        result = injector.inject()
        assert isinstance(result, Ratings)
        assert result is not simple_ratings

    def test_original_columns_preserved(self, simple_ratings,
                                        distribution_config_uniform):
        injector = DistributionInjector(
            simple_ratings, distribution_config_uniform,
        )
        result = injector.inject()
        np.testing.assert_array_equal(
            result.category_index[:, :3],
            simple_ratings.category_index,
        )

    def test_subject_ids_preserved(self, simple_ratings,
                                   distribution_config_uniform):
        injector = DistributionInjector(
            simple_ratings, distribution_config_uniform,
        )
        result = injector.inject()
        np.testing.assert_array_equal(
            result.subject_ids, simple_ratings.subject_ids,
        )


# ============================================================================ #
#  DistributionInjector                                                        #
# ============================================================================ #


class TestDistributionInjectorValidation:
    """Config type and prevalences length checks."""

    def test_wrong_config_type_rejected(self, simple_ratings):
        with pytest.raises(TypeError, match="DistributionInjectorConfig"):
            DistributionInjector(simple_ratings, InjectorConfig())

    def test_prevalences_length_mismatch(self, simple_ratings):
        cfg = DistributionInjectorConfig(
            prevalences=(0.3, 0.3, 0.4),  # 3, but ratings has 2 categories
            n_new_raters=1,
        )
        with pytest.raises(ValueError, match="prevalences"):
            DistributionInjector(simple_ratings, cfg)


class TestDistributionInjectorGenerate:
    """Generation behavior."""

    def test_uniform_shape_and_dtype(self, simple_ratings,
                                     distribution_config_uniform):
        injector = DistributionInjector(
            simple_ratings, distribution_config_uniform,
        )
        cols = injector.generate()
        assert cols.shape == (4, 1)
        assert cols.dtype == INDEX_DTYPE

    def test_uniform_values_in_range(self, simple_ratings,
                                     distribution_config_uniform):
        injector = DistributionInjector(
            simple_ratings, distribution_config_uniform,
        )
        cols = injector.generate()
        assert cols.min() >= 0
        assert cols.max() < 2  # 2 categories

    def test_custom_prevalences_values_in_range(self, simple_ratings,
                                                distribution_config):
        injector = DistributionInjector(
            simple_ratings, distribution_config,
        )
        cols = injector.generate()
        assert cols.shape == (4, 2)
        assert cols.min() >= 0
        assert cols.max() < 2

    def test_reproducible_with_seed(self, simple_ratings):
        cfg = DistributionInjectorConfig(
            prevalences=(0.5, 0.5), n_new_raters=3, random_state=123,
        )
        cols1 = DistributionInjector(simple_ratings, cfg).generate()
        cols2 = DistributionInjector(simple_ratings, cfg).generate()
        np.testing.assert_array_equal(cols1, cols2)


# ============================================================================ #
#  WithinSubjectInjector                                                       #
# ============================================================================ #


class TestWithinSubjectInjectorValidation:
    """Config type check."""

    def test_wrong_config_type_rejected(self, simple_ratings):
        with pytest.raises(TypeError, match="WithinSubjectInjectorConfig"):
            WithinSubjectInjector(simple_ratings, InjectorConfig())


class TestWithinSubjectInjectorGenerate:
    """Generation behavior: pure copy and noise."""

    def test_pure_copy_output_subset_of_real(self, simple_ratings,
                                             within_subject_config):
        """With noise_rate=0, every generated value must appear in
        the real ratings at that subject."""
        injector = WithinSubjectInjector(
            simple_ratings, within_subject_config,
        )
        cols = injector.generate()
        assert cols.shape == (4, 2)

        real = simple_ratings.category_index
        for i in range(4):
            real_vals = set(real[i, :])
            for j in range(2):
                assert cols[i, j] in real_vals

    def test_unanimous_row_always_copies_same(self, simple_ratings):
        """Subject s2 is unanimous (0, 0, 0) -- copy must produce 0."""
        cfg = WithinSubjectInjectorConfig(
            noise_rate=0.0, n_new_raters=5, random_state=42,
        )
        injector = WithinSubjectInjector(simple_ratings, cfg)
        cols = injector.generate()
        # Row 2 is all zeros
        assert np.all(cols[2, :] == 0)

    def test_noise_can_produce_unused_category(self, three_cat_ratings):
        """With noise_rate=1.0 on 3-cat data where some subjects don't
        use all categories, the injector should pick unused ones."""
        cfg = WithinSubjectInjectorConfig(
            noise_rate=1.0, n_new_raters=50, random_state=42,
        )
        injector = WithinSubjectInjector(three_cat_ratings, cfg)
        cols = injector.generate()

        # Subject s1: all raters used cat 1. Unused: {0, 2}.
        # With noise=1.0 and 50 tries, we should see 0 or 2.
        s1_vals = set(cols[1, :])
        assert len(s1_vals & {0, 2}) > 0

    def test_noise_fallback_when_all_used(self, three_cat_ratings):
        """Subject s2 uses all 3 categories. noise_rate=1.0 should
        fall back to copying (no unused category available)."""
        cfg = WithinSubjectInjectorConfig(
            noise_rate=1.0, n_new_raters=10, random_state=42,
        )
        injector = WithinSubjectInjector(three_cat_ratings, cfg)
        cols = injector.generate()

        # Row 2 uses {0, 1, 2} -- all categories, so fallback copies
        real_vals = set(three_cat_ratings.category_index[2, :])
        for j in range(10):
            assert cols[2, j] in real_vals

    def test_reproducible_with_seed(self, simple_ratings):
        cfg = WithinSubjectInjectorConfig(
            noise_rate=0.3, n_new_raters=5, random_state=99,
        )
        cols1 = WithinSubjectInjector(simple_ratings, cfg).generate()
        cols2 = WithinSubjectInjector(simple_ratings, cfg).generate()
        np.testing.assert_array_equal(cols1, cols2)


# ============================================================================ #
#  PooledInjector                                                              #
# ============================================================================ #


class TestPooledInjectorValidation:
    """Config type and data-dependent checks."""

    def test_wrong_config_type_rejected(self, simple_ratings):
        with pytest.raises(TypeError, match="PooledInjectorConfig"):
            PooledInjector(simple_ratings, InjectorConfig())

    def test_exclude_out_of_bounds_rejected(self, simple_ratings):
        cfg = PooledInjectorConfig(exclude_raters=(10,))
        with pytest.raises(ValueError, match="exclude_raters"):
            PooledInjector(simple_ratings, cfg)

    def test_exclude_all_rejected(self, simple_ratings):
        cfg = PooledInjectorConfig(exclude_raters=(0, 1, 2))
        with pytest.raises(ValueError, match="removes all"):
            PooledInjector(simple_ratings, cfg)


class TestPooledInjectorGenerate:
    """Generation behavior: pool sampling and smoothing."""

    def test_no_exclusion_shape(self, simple_ratings, pooled_config):
        injector = PooledInjector(simple_ratings, pooled_config)
        cols = injector.generate()
        assert cols.shape == (4, 2)
        assert cols.dtype == INDEX_DTYPE

    def test_values_from_pool(self, simple_ratings, pooled_config):
        injector = PooledInjector(simple_ratings, pooled_config)
        cols = injector.generate()
        assert cols.min() >= 0
        assert cols.max() < 2

    def test_exclude_rater(self, simple_ratings, pooled_config_exclude):
        """Excluding rater 2 (carol) -- pool should only contain
        alice and bob's ratings."""
        injector = PooledInjector(simple_ratings, pooled_config_exclude)
        cols = injector.generate()
        assert cols.shape == (4, 1)
        assert cols.min() >= 0
        assert cols.max() < 2

    def test_smoothing_values_in_range(self, simple_ratings,
                                       pooled_config_smooth):
        injector = PooledInjector(simple_ratings, pooled_config_smooth)
        cols = injector.generate()
        assert cols.min() >= 0
        assert cols.max() < 2

    def test_reproducible_with_seed(self, simple_ratings):
        cfg = PooledInjectorConfig(
            smoothing=0.2, n_new_raters=3, random_state=55,
        )
        cols1 = PooledInjector(simple_ratings, cfg).generate()
        cols2 = PooledInjector(simple_ratings, cfg).generate()
        np.testing.assert_array_equal(cols1, cols2)


# ============================================================================ #
#  RaterBootstrapInjector                                                      #
# ============================================================================ #


class TestRaterBootstrapInjectorValidation:
    """Config type and source rater checks."""

    def test_wrong_config_type_rejected(self, simple_ratings):
        with pytest.raises(TypeError, match="RaterBootstrapInjectorConfig"):
            RaterBootstrapInjector(simple_ratings, InjectorConfig())

    def test_source_out_of_bounds_rejected(self, simple_ratings):
        cfg = RaterBootstrapInjectorConfig(source_raters=(10,))
        with pytest.raises(ValueError, match="source_raters"):
            RaterBootstrapInjector(simple_ratings, cfg)


class TestRaterBootstrapInjectorGenerate:
    """Generation behavior."""

    def test_shape_and_dtype(self, simple_ratings, bootstrap_config):
        injector = RaterBootstrapInjector(simple_ratings, bootstrap_config)
        cols = injector.generate()
        assert cols.shape == (4, 2)
        assert cols.dtype == INDEX_DTYPE

    def test_values_from_source(self, simple_ratings):
        """Single source rater (alice, col 0): all generated values
        must come from alice's ratings {0, 1}."""
        cfg = RaterBootstrapInjectorConfig(
            source_raters=(0,), n_new_raters=5, random_state=42,
        )
        injector = RaterBootstrapInjector(simple_ratings, cfg)
        cols = injector.generate()
        alice_vals = set(simple_ratings.category_index[:, 0])
        for val in cols.ravel():
            assert val in alice_vals

    def test_reproducible_with_seed(self, simple_ratings):
        cfg = RaterBootstrapInjectorConfig(
            source_raters=(0, 1), n_new_raters=3, random_state=77,
        )
        cols1 = RaterBootstrapInjector(simple_ratings, cfg).generate()
        cols2 = RaterBootstrapInjector(simple_ratings, cfg).generate()
        np.testing.assert_array_equal(cols1, cols2)


# ============================================================================ #
#  ExtremalInjector                                                            #
# ============================================================================ #


class TestExtremalInjectorValidation:
    """Config type check."""

    def test_wrong_config_type_rejected(self, simple_ratings):
        with pytest.raises(TypeError, match="ExtremalInjectorConfig"):
            ExtremalInjector(simple_ratings, InjectorConfig())


class TestExtremalInjectorGenerate:
    """Generation behavior for all four strategies."""

    def test_majority_picks_most_common(self, simple_ratings,
                                        extremal_config_majority):
        """
        Subject layout (2 categories):
          s0: [0, 0, 1] -> majority = 0
          s1: [1, 1, 0] -> majority = 1
          s2: [0, 0, 0] -> majority = 0
          s3: [1, 0, 1] -> majority = 1
        """
        injector = ExtremalInjector(
            simple_ratings, extremal_config_majority,
        )
        cols = injector.generate()
        assert cols[0, 0] == 0
        assert cols[1, 0] == 1
        assert cols[2, 0] == 0
        assert cols[3, 0] == 1

    def test_minority_picks_least_common(self, simple_ratings,
                                         extremal_config_minority):
        """
          s0: [0, 0, 1] -> minority = 1 (count 1 vs 2)
          s1: [1, 1, 0] -> minority = 0 (count 1 vs 2)
          s2: [0, 0, 0] -> minority = 0 (unanimous, only option)
          s3: [1, 0, 1] -> minority = 0 (count 1 vs 2)
        """
        injector = ExtremalInjector(
            simple_ratings, extremal_config_minority,
        )
        cols = injector.generate()
        assert cols[0, 0] == 1
        assert cols[1, 0] == 0
        assert cols[2, 0] == 0  # unanimous: min == max
        assert cols[3, 0] == 0

    def test_contrarian_picks_unused(self, three_cat_ratings,
                                     extremal_config_contrarian):
        """
        3-cat data:
          s0: [0, 0, 1] -> unused = {2}
          s1: [1, 1, 1] -> unused = {0, 2}
          s2: [0, 1, 2] -> unused = {} -> fallback to minority
          s3: [2, 2, 0] -> unused = {1}
        """
        injector = ExtremalInjector(
            three_cat_ratings, extremal_config_contrarian,
        )
        cols = injector.generate()

        # s0: must be 2 (only unused)
        assert cols[0, 0] == 2
        # s1: must be 0 or 2 (both unused)
        assert cols[1, 0] in {0, 2}
        # s2: all used, fallback to minority.
        # counts: [1, 1, 1] -- all tied, so any is valid
        assert cols[2, 0] in {0, 1, 2}
        # s3: must be 1 (only unused)
        assert cols[3, 0] == 1

    def test_contrarian_fallback_to_minority(self, simple_ratings,
                                             extremal_config_contrarian):
        """With only 2 categories and all subjects using both (except s2),
        contrarian falls back to minority at covered subjects."""
        injector = ExtremalInjector(
            simple_ratings, extremal_config_contrarian,
        )
        cols = injector.generate()

        # s0: [0, 0, 1] -> both used, fallback minority -> 1
        assert cols[0, 0] == 1
        # s2: [0, 0, 0] -> only 0 used, unused = {1}
        assert cols[2, 0] == 1

    def test_median_picks_middle_prevalence(self, three_cat_ratings,
                                            extremal_config_median):
        """
        3-cat data:
          s0: [0, 0, 1] -> counts [2, 1, 0], used sorted [1, 2],
              median_idx=0 -> count=1 -> cat 1
          s1: [1, 1, 1] -> counts [0, 3, 0], only 1 used -> cat 1
          s2: [0, 1, 2] -> counts [1, 1, 1], all tied at 1,
              median_idx=1 -> count=1 -> all tied -> any of {0,1,2}
          s3: [2, 2, 0] -> counts [1, 0, 2], used sorted [1, 2],
              median_idx=0 -> count=1 -> cat 0
        """
        injector = ExtremalInjector(
            three_cat_ratings, extremal_config_median,
        )
        cols = injector.generate()
        assert cols[0, 0] == 1
        assert cols[1, 0] == 1
        assert cols[2, 0] in {0, 1, 2}  # all tied
        assert cols[3, 0] == 0

    def test_shape_and_dtype(self, simple_ratings,
                             extremal_config_majority):
        injector = ExtremalInjector(
            simple_ratings, extremal_config_majority,
        )
        cols = injector.generate()
        assert cols.shape == (4, 1)
        assert cols.dtype == INDEX_DTYPE

    def test_reproducible_with_seed(self, three_cat_ratings):
        cfg = ExtremalInjectorConfig(
            strategy="contrarian", n_new_raters=3, random_state=42,
        )
        cols1 = ExtremalInjector(three_cat_ratings, cfg).generate()
        cols2 = ExtremalInjector(three_cat_ratings, cfg).generate()
        np.testing.assert_array_equal(cols1, cols2)


# ============================================================================ #
#  Registry                                                                    #
# ============================================================================ #


class TestRegistry:
    """Registry mapping and lookup."""

    def test_all_types_registered(self):
        registry = get_registry()
        expected = {
            "distribution", "within_subject", "pooled",
            "rater_bootstrap", "extremal",
        }
        assert set(registry.keys()) == expected

    def test_registry_values_are_pairs(self):
        registry = get_registry()
        for key, (cfg_cls, inj_cls) in registry.items():
            assert issubclass(cfg_cls, InjectorConfig)
            assert issubclass(inj_cls, RaterInjector)

    def test_registry_cached(self):
        r1 = get_registry()
        r2 = get_registry()
        assert r1 is r2


class TestHydrateInjectionStep:
    """hydrate_injection_step() conversion."""

    def test_valid_step(self):
        step = {"type": "distribution", "n_new_raters": 3}
        cfg = hydrate_injection_step(step)
        assert isinstance(cfg, DistributionInjectorConfig)
        assert cfg.n_new_raters == 3

    def test_type_key_popped(self):
        step = {"type": "extremal", "strategy": "minority"}
        hydrate_injection_step(step)
        assert "type" not in step

    def test_missing_type_rejected(self):
        with pytest.raises(ValueError, match="type"):
            hydrate_injection_step({"n_new_raters": 1})

    def test_unknown_type_rejected(self):
        with pytest.raises(ValueError, match="unknown injection type"):
            hydrate_injection_step({"type": "bogus"})

    def test_non_dict_rejected(self):
        with pytest.raises(TypeError, match="must be a dict"):
            hydrate_injection_step("not_a_dict")


class TestGetInjectorClass:
    """get_injector_class() reverse lookup."""

    def test_distribution(self):
        cfg = DistributionInjectorConfig()
        cls = get_injector_class(cfg)
        assert cls is DistributionInjector

    def test_within_subject(self):
        cfg = WithinSubjectInjectorConfig()
        cls = get_injector_class(cfg)
        assert cls is WithinSubjectInjector

    def test_pooled(self):
        cfg = PooledInjectorConfig()
        cls = get_injector_class(cfg)
        assert cls is PooledInjector

    def test_rater_bootstrap(self):
        cfg = RaterBootstrapInjectorConfig()
        cls = get_injector_class(cfg)
        assert cls is RaterBootstrapInjector

    def test_extremal(self):
        cfg = ExtremalInjectorConfig()
        cls = get_injector_class(cfg)
        assert cls is ExtremalInjector

    def test_non_config_rejected(self):
        with pytest.raises(TypeError, match="InjectorConfig"):
            get_injector_class("not_a_config")

    def test_base_config_not_registered(self):
        """InjectorConfig itself has no injector."""
        with pytest.raises(ValueError, match="no injector registered"):
            get_injector_class(InjectorConfig())


# ============================================================================ #
#  Chaining                                                                    #
# ============================================================================ #


class TestChaining:
    """Inject twice to verify chaining semantics."""

    def test_double_injection_shape(self, simple_ratings):
        """Two rounds of injection produce correct total columns."""
        inj1 = DistributionInjector(
            simple_ratings, DistributionInjectorConfig(
                n_new_raters=2, random_state=0,
            )
        )
        after_first = inj1.inject()
        assert after_first.category_index.shape == (4, 5)

        inj2 = WithinSubjectInjector(
            after_first, WithinSubjectInjectorConfig(
                n_new_raters=1, random_state=1,
            )
        )
        after_second = inj2.inject()
        assert after_second.category_index.shape == (4, 6)

    def test_double_injection_names(self, simple_ratings):
        """Chained injection names are sequential."""
        inj1 = PooledInjector(
            simple_ratings, PooledInjectorConfig(
                n_new_raters=1, random_state=0,
            )
        )
        after_first = inj1.inject()

        inj2 = PooledInjector(
            after_first, PooledInjectorConfig(
                n_new_raters=1, random_state=1,
            )
        )
        after_second = inj2.inject()

        raters = after_second.semantics.raters
        assert "injected_rater0" in raters
        assert "injected_rater1" in raters
        assert len(set(raters)) == len(raters)

    def test_original_data_untouched(self, simple_ratings):
        """Original Ratings object is not modified by injection."""
        original_copy = simple_ratings.category_index.copy()
        inj = DistributionInjector(
            simple_ratings, DistributionInjectorConfig(
                n_new_raters=3, random_state=42,
            )
        )
        _ = inj.inject()
        np.testing.assert_array_equal(
            simple_ratings.category_index, original_copy,
        )


# ============================================================================ #
#  Injection through an analyzer config                                        #
# ============================================================================ #


class TestInjectionNeedsRatings:
    """
    Injection builds a new rater from the real raters' answers. A nij matrix
    holds only counts per subject, so there is nothing to build from: the
    analyzer refuses at construction rather than silently skipping the step
    and reporting an un-augmented coefficient as an augmented one.
    """

    @staticmethod
    def _config(config_cls):
        return config_cls(injection=DistributionInjectorConfig(n_new_raters=1, random_state=0))

    @pytest.mark.parametrize("name", ["FleissKappa", "GwetAC1"])
    def test_nij_path_refused(self, name):
        import interrater.analyzers as A
        import interrater.analyzers.configs as C
        from interrater.intermediates.estimator_parts import NijMatrix

        nij = NijMatrix(matrix=np.array([[3, 0], [2, 1], [1, 2], [0, 3]]))
        with pytest.raises(ValueError, match="injection needs individual ratings"):
            getattr(A, f"{name}Analyzer").build_from_nij_matrix(
                nij, self._config(getattr(C, f"{name}Config")),
            )

    def test_ratings_path_injects(self, simple_ratings):
        from interrater.analyzers import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        analyzer = FleissKappaAnalyzer(simple_ratings, self._config(FleissKappaConfig))
        assert analyzer.ratings.n_raters == simple_ratings.n_raters + 1


# ============================================================================ #
#  Seeds in a chain                                                            #
# ============================================================================ #


class TestChainedSeeds:
    """
    Each step is seeded with [random_state, step_position]. Before 2.1.0 a
    step was seeded with random_state alone, so two chained steps given the
    same seed produced identical raters.
    """

    @staticmethod
    def _injected(simple_ratings, steps):
        from interrater.analyzers import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        analyzer = FleissKappaAnalyzer(
            simple_ratings, FleissKappaConfig(injection=steps, n_resamples=2),
        )
        return analyzer.ratings.category_index[:, simple_ratings.n_raters:]

    def test_same_seed_in_two_steps_gives_different_raters(self, simple_ratings):
        step = DistributionInjectorConfig(n_new_raters=5, random_state=42)
        injected = self._injected(simple_ratings, [step, step])
        assert not np.array_equal(injected[:, :5], injected[:, 5:])

    def test_pipeline_reruns_identically(self, simple_ratings):
        steps = [
            DistributionInjectorConfig(n_new_raters=2, random_state=42),
            PooledInjectorConfig(n_new_raters=2, random_state=42),
        ]
        np.testing.assert_array_equal(
            self._injected(simple_ratings, steps), self._injected(simple_ratings, steps),
        )

    def test_seed_is_the_pair(self, simple_ratings):
        """The rule itself: the generator is default_rng([random_state, step_position])."""
        config = DistributionInjectorConfig(n_new_raters=1, random_state=42)
        injector = DistributionInjector(simple_ratings, config, step_position=3)
        expected = np.random.default_rng([42, 3]).random(5)
        np.testing.assert_array_equal(injector._rng().random(5), expected)

    @pytest.mark.parametrize("bad, error", [(-1, ValueError), (True, TypeError), (1.0, TypeError)])
    def test_step_position_validated(self, simple_ratings, bad, error):
        config = DistributionInjectorConfig(n_new_raters=1, random_state=42)
        with pytest.raises(error, match="step_position"):
            DistributionInjector(simple_ratings, config, step_position=bad)
