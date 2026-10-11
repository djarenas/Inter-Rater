"""
Tests for injector configuration classes.

Bottom-up: base InjectorConfig first, then each leaf config.
Each config is a frozen dataclass that validates its own fields
at construction time, with no data dependency.
"""

import pytest
import tempfile
import os

from interrater.injectors.configs import (
    InjectorConfig,
    DistributionInjectorConfig,
    WithinSubjectInjectorConfig,
    PooledInjectorConfig,
    RaterBootstrapInjectorConfig,
    ExtremalInjectorConfig,
)


# ============================================================================ #
#  InjectorConfig (base)                                                       #
# ============================================================================ #


class TestInjectorConfigDefaults:
    """Default construction and field values."""

    def test_default_n_new_raters(self):
        cfg = InjectorConfig()
        assert cfg.n_new_raters == 1

    def test_default_random_state(self):
        cfg = InjectorConfig()
        assert cfg.random_state is None

    def test_explicit_values(self, base_config):
        assert base_config.n_new_raters == 2
        assert base_config.random_state == 42

    def test_frozen(self, base_config):
        with pytest.raises(AttributeError):
            base_config.n_new_raters = 5


class TestInjectorConfigValidation:
    """Type and value validation for base fields."""

    # -- n_new_raters -------------------------------------------------------- #

    def test_n_new_raters_bool_rejected(self):
        with pytest.raises(TypeError, match="n_new_raters"):
            InjectorConfig(n_new_raters=True)

    def test_n_new_raters_float_rejected(self):
        with pytest.raises(TypeError, match="n_new_raters"):
            InjectorConfig(n_new_raters=1.0)

    def test_n_new_raters_string_rejected(self):
        with pytest.raises(TypeError, match="n_new_raters"):
            InjectorConfig(n_new_raters="2")

    def test_n_new_raters_zero_rejected(self):
        with pytest.raises(ValueError, match="n_new_raters"):
            InjectorConfig(n_new_raters=0)

    def test_n_new_raters_negative_rejected(self):
        with pytest.raises(ValueError, match="n_new_raters"):
            InjectorConfig(n_new_raters=-1)

    # -- random_state -------------------------------------------------------- #

    def test_random_state_none_accepted(self):
        cfg = InjectorConfig(random_state=None)
        assert cfg.random_state is None

    def test_random_state_zero_accepted(self):
        cfg = InjectorConfig(random_state=0)
        assert cfg.random_state == 0

    def test_random_state_bool_rejected(self):
        with pytest.raises(TypeError, match="random_state"):
            InjectorConfig(random_state=True)

    def test_random_state_float_rejected(self):
        with pytest.raises(TypeError, match="random_state"):
            InjectorConfig(random_state=1.0)

    def test_random_state_negative_rejected(self):
        with pytest.raises(ValueError, match="random_state"):
            InjectorConfig(random_state=-1)


class TestInjectorConfigBuildFromDict:
    """build_from_dict validation and construction."""

    def test_valid_dict(self):
        cfg = InjectorConfig.build_from_dict(
            {"n_new_raters": 3, "random_state": 99}
        )
        assert cfg.n_new_raters == 3
        assert cfg.random_state == 99

    def test_empty_dict_uses_defaults(self):
        cfg = InjectorConfig.build_from_dict({})
        assert cfg.n_new_raters == 1
        assert cfg.random_state is None

    def test_none_rejected(self):
        with pytest.raises(ValueError, match="data_dict"):
            InjectorConfig.build_from_dict(None)

    def test_non_dict_rejected(self):
        with pytest.raises(TypeError, match="data_dict"):
            InjectorConfig.build_from_dict([1, 2])


class TestInjectorConfigBuildFromYaml:
    """build_from_yaml validation and construction."""

    def test_valid_yaml(self, tmp_path):
        yaml_file = tmp_path / "config.yaml"
        yaml_file.write_text("n_new_raters: 5\nrandom_state: 7\n")
        cfg = InjectorConfig.build_from_yaml(str(yaml_file))
        assert cfg.n_new_raters == 5
        assert cfg.random_state == 7

    def test_unknown_key_rejected(self, tmp_path):
        yaml_file = tmp_path / "config.yaml"
        yaml_file.write_text("n_new_raters: 1\nbogus: 42\n")
        with pytest.raises(TypeError, match="Unknown YAML"):
            InjectorConfig.build_from_yaml(str(yaml_file))

    def test_missing_file_rejected(self):
        with pytest.raises(FileNotFoundError):
            InjectorConfig.build_from_yaml("/nonexistent/config.yaml")


class TestInjectorConfigSummary:
    """summary() returns a plain dict with version."""

    def test_contains_base_fields(self, base_config):
        s = base_config.summary()
        assert s["n_new_raters"] == 2
        assert s["random_state"] == 42

    def test_contains_package_version(self, base_config):
        s = base_config.summary()
        assert "_package_version" in s


# ============================================================================ #
#  DistributionInjectorConfig                                                  #
# ============================================================================ #


class TestDistributionInjectorConfigDefaults:
    """Default construction."""

    def test_default_prevalences_is_none(self):
        cfg = DistributionInjectorConfig()
        assert cfg.prevalences is None

    def test_explicit_prevalences(self, distribution_config):
        assert distribution_config.prevalences == (0.6, 0.4)

    def test_inherits_base_defaults(self):
        cfg = DistributionInjectorConfig()
        assert cfg.n_new_raters == 1
        assert cfg.random_state is None


class TestDistributionInjectorConfigValidation:
    """Prevalences validation."""

    def test_list_rejected(self):
        with pytest.raises(TypeError, match="prevalences"):
            DistributionInjectorConfig(prevalences=[0.5, 0.5])

    def test_empty_tuple_rejected(self):
        with pytest.raises(ValueError, match="prevalences"):
            DistributionInjectorConfig(prevalences=())

    def test_negative_element_rejected(self):
        with pytest.raises(ValueError, match="prevalences"):
            DistributionInjectorConfig(prevalences=(0.5, -0.5))

    def test_bool_element_rejected(self):
        with pytest.raises(TypeError, match="prevalences"):
            DistributionInjectorConfig(prevalences=(True, False))

    def test_string_element_rejected(self):
        with pytest.raises(TypeError, match="prevalences"):
            DistributionInjectorConfig(prevalences=("0.5", "0.5"))

    def test_sum_not_one_rejected(self):
        with pytest.raises(ValueError, match="sum to 1.0"):
            DistributionInjectorConfig(prevalences=(0.5, 0.3))

    def test_sum_within_tolerance_accepted(self):
        # 0.1 + 0.2 + 0.7 in float: close enough to 1.0
        cfg = DistributionInjectorConfig(prevalences=(0.1, 0.2, 0.7))
        assert cfg.prevalences == (0.1, 0.2, 0.7)

    def test_inherits_base_validation(self):
        with pytest.raises(TypeError, match="n_new_raters"):
            DistributionInjectorConfig(n_new_raters=True)


class TestDistributionInjectorConfigBuildFromDict:
    """build_from_dict hydrates list -> tuple."""

    def test_list_to_tuple(self):
        cfg = DistributionInjectorConfig.build_from_dict(
            {"prevalences": [0.5, 0.5]}
        )
        assert isinstance(cfg.prevalences, tuple)
        assert cfg.prevalences == (0.5, 0.5)

    def test_none_stays_none(self):
        cfg = DistributionInjectorConfig.build_from_dict(
            {"prevalences": None}
        )
        assert cfg.prevalences is None


class TestDistributionInjectorConfigSummary:
    """summary() includes prevalences."""

    def test_contains_prevalences(self, distribution_config):
        s = distribution_config.summary()
        assert s["prevalences"] == (0.6, 0.4)
        assert "n_new_raters" in s

    def test_none_prevalences_in_summary(self, distribution_config_uniform):
        s = distribution_config_uniform.summary()
        assert s["prevalences"] is None


# ============================================================================ #
#  WithinSubjectInjectorConfig                                                 #
# ============================================================================ #


class TestWithinSubjectInjectorConfigDefaults:
    """Default construction."""

    def test_default_noise_rate(self):
        cfg = WithinSubjectInjectorConfig()
        assert cfg.noise_rate == 0.0


class TestWithinSubjectInjectorConfigValidation:
    """noise_rate validation."""

    def test_bool_rejected(self):
        with pytest.raises(TypeError, match="noise_rate"):
            WithinSubjectInjectorConfig(noise_rate=True)

    def test_string_rejected(self):
        with pytest.raises(TypeError, match="noise_rate"):
            WithinSubjectInjectorConfig(noise_rate="0.5")

    def test_negative_rejected(self):
        with pytest.raises(ValueError, match="noise_rate"):
            WithinSubjectInjectorConfig(noise_rate=-0.1)

    def test_above_one_rejected(self):
        with pytest.raises(ValueError, match="noise_rate"):
            WithinSubjectInjectorConfig(noise_rate=1.1)

    def test_boundary_zero_accepted(self):
        cfg = WithinSubjectInjectorConfig(noise_rate=0.0)
        assert cfg.noise_rate == 0.0

    def test_boundary_one_accepted(self):
        cfg = WithinSubjectInjectorConfig(noise_rate=1.0)
        assert cfg.noise_rate == 1.0

    def test_int_accepted(self):
        cfg = WithinSubjectInjectorConfig(noise_rate=0)
        assert cfg.noise_rate == 0

    def test_inherits_base_validation(self):
        with pytest.raises(ValueError, match="n_new_raters"):
            WithinSubjectInjectorConfig(n_new_raters=0)


class TestWithinSubjectInjectorConfigSummary:
    """summary() includes noise_rate."""

    def test_contains_noise_rate(self, within_subject_config):
        s = within_subject_config.summary()
        assert s["noise_rate"] == 0.0
        assert "n_new_raters" in s


# ============================================================================ #
#  PooledInjectorConfig                                                        #
# ============================================================================ #


class TestPooledInjectorConfigDefaults:
    """Default construction."""

    def test_default_smoothing(self):
        cfg = PooledInjectorConfig()
        assert cfg.smoothing == 0.0

    def test_default_exclude_raters(self):
        cfg = PooledInjectorConfig()
        assert cfg.exclude_raters == ()


class TestPooledInjectorConfigValidation:
    """Smoothing and exclude_raters validation."""

    # -- smoothing ----------------------------------------------------------- #

    def test_smoothing_bool_rejected(self):
        with pytest.raises(TypeError, match="smoothing"):
            PooledInjectorConfig(smoothing=True)

    def test_smoothing_string_rejected(self):
        with pytest.raises(TypeError, match="smoothing"):
            PooledInjectorConfig(smoothing="0.5")

    def test_smoothing_negative_rejected(self):
        with pytest.raises(ValueError, match="smoothing"):
            PooledInjectorConfig(smoothing=-0.1)

    def test_smoothing_above_one_rejected(self):
        with pytest.raises(ValueError, match="smoothing"):
            PooledInjectorConfig(smoothing=1.1)

    def test_smoothing_boundary_values(self):
        PooledInjectorConfig(smoothing=0.0)
        PooledInjectorConfig(smoothing=1.0)

    # -- exclude_raters ------------------------------------------------------ #

    def test_exclude_raters_list_rejected(self):
        with pytest.raises(TypeError, match="exclude_raters"):
            PooledInjectorConfig(exclude_raters=[0, 1])

    def test_exclude_raters_bool_element_rejected(self):
        with pytest.raises(TypeError, match="exclude_raters"):
            PooledInjectorConfig(exclude_raters=(True,))

    def test_exclude_raters_negative_rejected(self):
        with pytest.raises(ValueError, match="exclude_raters"):
            PooledInjectorConfig(exclude_raters=(-1,))

    def test_exclude_raters_duplicates_rejected(self):
        with pytest.raises(ValueError, match="duplicates"):
            PooledInjectorConfig(exclude_raters=(0, 0))


class TestPooledInjectorConfigBuildFromDict:
    """build_from_dict hydrates list -> tuple."""

    def test_list_to_tuple(self):
        cfg = PooledInjectorConfig.build_from_dict(
            {"exclude_raters": [0, 1]}
        )
        assert isinstance(cfg.exclude_raters, tuple)
        assert cfg.exclude_raters == (0, 1)


class TestPooledInjectorConfigSummary:
    """summary() includes smoothing and exclude_raters."""

    def test_contains_pooled_fields(self, pooled_config):
        s = pooled_config.summary()
        assert s["smoothing"] == 0.0
        assert s["exclude_raters"] == ()
        assert "n_new_raters" in s


# ============================================================================ #
#  RaterBootstrapInjectorConfig                                                #
# ============================================================================ #


class TestRaterBootstrapInjectorConfigDefaults:
    """Default construction."""

    def test_default_source_raters(self):
        cfg = RaterBootstrapInjectorConfig()
        assert cfg.source_raters == (0,)


class TestRaterBootstrapInjectorConfigValidation:
    """source_raters validation."""

    def test_list_rejected(self):
        with pytest.raises(TypeError, match="source_raters"):
            RaterBootstrapInjectorConfig(source_raters=[0, 1])

    def test_empty_rejected(self):
        with pytest.raises(ValueError, match="source_raters"):
            RaterBootstrapInjectorConfig(source_raters=())

    def test_bool_element_rejected(self):
        with pytest.raises(TypeError, match="source_raters"):
            RaterBootstrapInjectorConfig(source_raters=(True,))

    def test_negative_rejected(self):
        with pytest.raises(ValueError, match="source_raters"):
            RaterBootstrapInjectorConfig(source_raters=(-1,))

    def test_duplicates_rejected(self):
        with pytest.raises(ValueError, match="duplicates"):
            RaterBootstrapInjectorConfig(source_raters=(0, 0))


class TestRaterBootstrapInjectorConfigBuildFromDict:
    """build_from_dict hydrates list -> tuple."""

    def test_list_to_tuple(self):
        cfg = RaterBootstrapInjectorConfig.build_from_dict(
            {"source_raters": [0, 1, 2]}
        )
        assert isinstance(cfg.source_raters, tuple)
        assert cfg.source_raters == (0, 1, 2)


class TestRaterBootstrapInjectorConfigSummary:
    """summary() includes source_raters."""

    def test_contains_source_raters(self, bootstrap_config):
        s = bootstrap_config.summary()
        assert s["source_raters"] == (0, 1)
        assert "n_new_raters" in s


# ============================================================================ #
#  ExtremalInjectorConfig                                                      #
# ============================================================================ #


class TestExtremalInjectorConfigDefaults:
    """Default construction."""

    def test_default_strategy(self):
        cfg = ExtremalInjectorConfig()
        assert cfg.strategy == "majority"


class TestExtremalInjectorConfigValidation:
    """Strategy validation."""

    def test_non_string_rejected(self):
        with pytest.raises(TypeError, match="strategy"):
            ExtremalInjectorConfig(strategy=42)

    def test_invalid_strategy_rejected(self):
        with pytest.raises(ValueError, match="strategy"):
            ExtremalInjectorConfig(strategy="random")

    def test_all_valid_strategies_accepted(self):
        for strategy in ("majority", "minority", "median", "contrarian"):
            cfg = ExtremalInjectorConfig(strategy=strategy)
            assert cfg.strategy == strategy


class TestExtremalInjectorConfigSummary:
    """summary() includes strategy."""

    def test_contains_strategy(self, extremal_config_majority):
        s = extremal_config_majority.summary()
        assert s["strategy"] == "majority"
        assert "n_new_raters" in s


# ============================================================================ #
#  Repr                                                                        #
# ============================================================================ #


class TestRepr:
    """All config types have a sensible repr."""

    def test_base_config(self, base_config):
        r = repr(base_config)
        assert "InjectorConfig" in r
        assert "n_new_raters=2" in r

    def test_distribution_config(self, distribution_config):
        r = repr(distribution_config)
        assert "DistributionInjectorConfig" in r
        assert "prevalences" in r

    def test_within_subject_config(self, within_subject_config):
        r = repr(within_subject_config)
        assert "WithinSubjectInjectorConfig" in r
        assert "noise_rate" in r

    def test_pooled_config(self, pooled_config):
        r = repr(pooled_config)
        assert "PooledInjectorConfig" in r
        assert "smoothing" in r

    def test_bootstrap_config(self, bootstrap_config):
        r = repr(bootstrap_config)
        assert "RaterBootstrapInjectorConfig" in r
        assert "source_raters" in r

    def test_extremal_config(self, extremal_config_majority):
        r = repr(extremal_config_majority)
        assert "ExtremalInjectorConfig" in r
        assert "majority" in r
