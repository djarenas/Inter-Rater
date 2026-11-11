"""Tests for analyzers/configs/ sub-package.

One test file for the entire sub-package (DandC convention).
Tests are organized by class, following the inheritance hierarchy.
"""

from __future__ import annotations

import dataclasses
import warnings

import pytest

from interrater.analyzers import configs
from interrater.analyzers.configs import (
    AnalyzerConfig,
    BrennanPredigerConfig,
    CohenKappaConfig,
    CongerKappaConfig,
    FleissKappaConfig,
    GroupCoefficientConfig,
    GwetAC1Config,
    HypothesisTestConfig,
    IndividualMarginalConfig,
    KrippendorffAlphaConfig,
    LightKappaConfig,
    NoWeightGroupConfig,
    NoWeightIndividualConfig,
    PairwisePermutationConfig,
    PermutationBootstrapMeanConfig,
    PermutationConfig,
    ScottsConfig,
    SubsetPermutationConfig,
    WeightableGroupConfig,
    WeightableIndividualConfig,
)


# ============================================================== #
#  AnalyzerConfig (base)                                         #
# ============================================================== #


class TestAnalyzerConfig:
    """Tests for the base AnalyzerConfig class."""

    def test_default_construction(self):
        cfg = AnalyzerConfig()
        assert cfg.n_resamples == 1000
        assert cfg.alpha == 0.05
        assert cfg.random_state == 42
        assert cfg.ci_method == "percentile"
        assert cfg.convergence_test is None
        assert cfg.hypothesis_test_config is None
        assert cfg.injection is None

    def test_frozen(self):
        cfg = AnalyzerConfig()
        with pytest.raises(AttributeError):
            cfg.n_resamples = 999

    def test_n_resamples_zero_rejected(self):
        with pytest.raises(ValueError, match="n_resamples"):
            AnalyzerConfig(n_resamples=0)

    def test_n_resamples_negative_rejected(self):
        with pytest.raises(ValueError, match="n_resamples"):
            AnalyzerConfig(n_resamples=-1)

    def test_alpha_zero_rejected(self):
        with pytest.raises(ValueError, match="alpha"):
            AnalyzerConfig(alpha=0.0)

    def test_alpha_one_rejected(self):
        with pytest.raises(ValueError, match="alpha"):
            AnalyzerConfig(alpha=1.0)

    def test_alpha_negative_rejected(self):
        with pytest.raises(ValueError, match="alpha"):
            AnalyzerConfig(alpha=-0.1)

    def test_ci_method_invalid_rejected(self):
        with pytest.raises(ValueError, match="ci_method"):
            AnalyzerConfig(ci_method="invalid")

    def test_ci_method_percentile_accepted(self):
        cfg = AnalyzerConfig(ci_method="percentile")
        assert cfg.ci_method == "percentile"

    def test_ci_method_bca_accepted(self):
        cfg = AnalyzerConfig(ci_method="bca")
        assert cfg.ci_method == "bca"

    def test_ci_method_parametric_null_accepted(self):
        """The base accepts it -- and says nothing, deliberately.

        AnalyzerConfig does not know which coefficient it is configuring,
        so it cannot know whether parametric_null is even permitted. Five
        of the eight leaf configs reject it, and they do so *after* calling
        super().__post_init__(). A warning here would therefore fire for
        configs about to raise, telling someone how to read an interval
        they are not allowed to have. The three that accept the method warn
        from their own __post_init__ instead.
        """
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cfg = AnalyzerConfig(ci_method="parametric_null")
        assert cfg.ci_method == "parametric_null"
        assert len(caught) == 0

    def test_injection_single_normalized_to_tuple(self):
        """A single InjectorConfig should be normalized to a 1-tuple."""
        from interrater.injectors.configs import InjectorConfig

        single = InjectorConfig()
        cfg = AnalyzerConfig(injection=single)
        assert isinstance(cfg.injection, tuple)
        assert len(cfg.injection) == 1

    def test_injection_list_normalized_to_tuple(self):
        """A list of InjectorConfigs should be normalized to a tuple."""
        from interrater.injectors.configs import InjectorConfig

        lst = [InjectorConfig(), InjectorConfig()]
        cfg = AnalyzerConfig(injection=lst)
        assert isinstance(cfg.injection, tuple)
        assert len(cfg.injection) == 2

    def test_injection_empty_list_rejected(self):
        with pytest.raises(ValueError, match="injection.*cannot be empty"):
            AnalyzerConfig(injection=[])

    def test_build_from_dict_none_rejected(self):
        with pytest.raises(ValueError, match="data_dict"):
            AnalyzerConfig.build_from_dict(None)

    def test_build_from_dict_non_dict_rejected(self):
        with pytest.raises(TypeError, match="data_dict"):
            AnalyzerConfig.build_from_dict("not a dict")

    def test_build_from_dict_empty(self):
        """Empty dict should use all defaults."""
        cfg = AnalyzerConfig.build_from_dict({})
        assert cfg.n_resamples == 1000
        assert cfg.alpha == 0.05

    def test_build_from_dict_overrides(self):
        cfg = AnalyzerConfig.build_from_dict(
            {"n_resamples": 5000, "alpha": 0.01}
        )
        assert cfg.n_resamples == 5000
        assert cfg.alpha == 0.01

    @pytest.mark.parametrize("key", ["weights_type", "custom_weights_path"])
    def test_build_from_dict_removed_weights_key_refused(self, key):
        """Weights moved to the semantics in 2.1.0; the old keys say so."""
        with pytest.raises(TypeError, match="removed in 2.1.0"):
            FleissKappaConfig.build_from_dict({key: None})

    def test_summary_has_version(self):
        cfg = AnalyzerConfig()
        s = cfg.summary()
        assert "_package_version" in s

    def test_summary_keys(self):
        cfg = AnalyzerConfig()
        s = cfg.summary()
        assert s["n_resamples"] == 1000
        assert s["alpha"] == 0.05
        assert s["random_state"] == 42
        assert s["ci_method"] == "percentile"


# ============================================================== #
#  HypothesisTestConfig                                          #
# ============================================================== #


class TestHypothesisTestConfig:
    """Tests for HypothesisTestConfig (Finding 1b: bool-first guards)."""

    def test_default_construction(self):
        cfg = HypothesisTestConfig()
        assert cfg.method == "permutation"
        assert cfg.scheme == "within_rater_shuffle"
        assert cfg.n_permutations == 5000
        assert cfg.alpha == 0.05
        assert cfg.random_state == 42

    def test_method_invalid_rejected(self):
        with pytest.raises(ValueError, match="method"):
            HypothesisTestConfig(method="invalid")

    def test_method_parametric_null_accepted(self):
        cfg = HypothesisTestConfig(method="parametric_null")
        assert cfg.method == "parametric_null"

    def test_scheme_invalid_rejected(self):
        with pytest.raises(ValueError, match="scheme"):
            HypothesisTestConfig(scheme="invalid")

    def test_alternative_default_two_sided_and_all_accepted(self):
        """Two-sided by default: below-chance agreement is a real finding."""
        assert HypothesisTestConfig().alternative == "two-sided"
        for alt in ("two-sided", "greater", "less"):
            assert HypothesisTestConfig(alternative=alt).alternative == alt

    def test_alternative_invalid_rejected(self):
        with pytest.raises(ValueError, match="alternative"):
            HypothesisTestConfig(alternative="two_sided")

    def test_scheme_within_subject_refused(self):
        with pytest.raises(ValueError, match="removed in 2.1.0"):
            HypothesisTestConfig(scheme="within_subject")

    def test_scheme_column_shuffle_refused(self):
        """The old name is refused, and the message names the new one."""
        with pytest.raises(ValueError, match="renamed to 'within_rater_shuffle' in 2.1.0"):
            HypothesisTestConfig(scheme="column_shuffle")

    # Finding 1b: bool-first guards on n_permutations
    def test_n_permutations_bool_rejected(self):
        """bool subclasses int; True must be caught before isinstance(int)."""
        with pytest.raises(TypeError, match="n_permutations"):
            HypothesisTestConfig(n_permutations=True)

    def test_n_permutations_false_rejected(self):
        with pytest.raises(TypeError, match="n_permutations"):
            HypothesisTestConfig(n_permutations=False)

    def test_n_permutations_string_rejected(self):
        with pytest.raises(TypeError, match="n_permutations"):
            HypothesisTestConfig(n_permutations="1000")

    def test_n_permutations_zero_rejected(self):
        with pytest.raises(ValueError, match="n_permutations"):
            HypothesisTestConfig(n_permutations=0)

    def test_n_permutations_negative_rejected(self):
        with pytest.raises(ValueError, match="n_permutations"):
            HypothesisTestConfig(n_permutations=-1)

    # Finding 1b: bool-first guards on alpha
    def test_alpha_bool_rejected(self):
        """bool subclasses int/float; True must be caught."""
        with pytest.raises(TypeError, match="alpha"):
            HypothesisTestConfig(alpha=True)

    def test_alpha_string_rejected(self):
        with pytest.raises(TypeError, match="alpha"):
            HypothesisTestConfig(alpha="0.05")

    def test_alpha_zero_rejected(self):
        with pytest.raises(ValueError, match="alpha"):
            HypothesisTestConfig(alpha=0.0)

    def test_alpha_one_rejected(self):
        with pytest.raises(ValueError, match="alpha"):
            HypothesisTestConfig(alpha=1.0)

    def test_random_state_string_rejected(self):
        with pytest.raises(TypeError, match="random_state"):
            HypothesisTestConfig(random_state="42")

    def test_summary_keys(self):
        cfg = HypothesisTestConfig()
        s = cfg.summary()
        assert set(s.keys()) == {
            "method", "scheme", "alternative", "n_permutations", "alpha",
            "random_state",
        }


# ============================================================== #
#  Group hierarchy: GroupCoefficientConfig → children             #
# ============================================================== #


class TestGroupCoefficientConfig:
    """Tests for GroupCoefficientConfig."""

    def test_default_construction(self):
        cfg = GroupCoefficientConfig()
        assert cfg.balanced_data_only is False
        assert cfg.rater_jackknife is False

    def test_balanced_data_only_non_bool_rejected(self):
        with pytest.raises(TypeError, match="balanced_data_only"):
            GroupCoefficientConfig(balanced_data_only=1)

    def test_rater_jackknife_non_bool_rejected(self):
        with pytest.raises(TypeError, match="rater_jackknife"):
            GroupCoefficientConfig(rater_jackknife=1)

    def test_rater_jackknife_true_rejected(self):
        """Refused, not silently ignored as in 2.0.0."""
        with pytest.raises(ValueError, match="rater_jackknife.*not supported in this version"):
            GroupCoefficientConfig(rater_jackknife=True)

    def test_summary_includes_group_fields(self):
        cfg = GroupCoefficientConfig()
        s = cfg.summary()
        assert "balanced_data_only" in s
        assert "rater_jackknife" in s


class TestNoWeightGroupConfig:
    """Tests for NoWeightGroupConfig (rejects weights)."""

    def test_default_construction(self):
        cfg = NoWeightGroupConfig()
        assert isinstance(cfg, NoWeightGroupConfig)

class TestWeightableGroupConfig:
    """Tests for WeightableGroupConfig (accepts weights)."""

    def test_default_construction(self):
        cfg = WeightableGroupConfig()
        assert isinstance(cfg, WeightableGroupConfig)

# ============================================================== #
#  Leaf group configs                                            #
# ============================================================== #


class TestFleissKappaConfig:
    """Tests for FleissKappaConfig."""

    def test_default_balanced_data_only_true(self):
        cfg = FleissKappaConfig()
        assert cfg.balanced_data_only is True

    def test_allows_parametric_null_ci(self):
        """Fleiss has a closed-form null variance (Eq. 12).

        Allowed, and warned about: the null variance is the right one for
        a hypothesis test and the wrong one for an interval around a
        non-zero estimate. The message has to name the alternatives, or it
        tells someone they have a problem without telling them the fix.
        """
        with pytest.warns(UserWarning, match="parametric_null") as caught:
            cfg = FleissKappaConfig(ci_method="parametric_null")
        assert cfg.ci_method == "parametric_null"

        message = str(caught[0].message)
        assert "percentile" in message
        assert "bca" in message

    def test_allows_parametric_null_hypothesis(self):
        ht = HypothesisTestConfig(method="parametric_null")
        cfg = FleissKappaConfig(hypothesis_test_config=ht)
        assert cfg.hypothesis_test_config.method == "parametric_null"


class TestScottsConfig:
    """Tests for ScottsConfig (Finding 1c: full dataclass pattern)."""

    def test_construction(self):
        """ScottsConfig should construct properly (not just pass)."""
        cfg = ScottsConfig()
        assert isinstance(cfg, NoWeightGroupConfig)

    def test_allows_parametric_null_ci(self):
        """Scott's pi = Fleiss' kappa at 2 raters, so parametric_null OK.

        Allowed, and warned about: the null variance is the right one for
        a hypothesis test and the wrong one for an interval around a
        non-zero estimate. The message has to name the alternatives, or it
        tells someone they have a problem without telling them the fix.
        """
        with pytest.warns(UserWarning, match="parametric_null") as caught:
            cfg = ScottsConfig(ci_method="parametric_null")
        assert cfg.ci_method == "parametric_null"

        message = str(caught[0].message)
        assert "percentile" in message
        assert "bca" in message

    def test_summary(self):
        cfg = ScottsConfig()
        s = cfg.summary()
        assert "_package_version" in s


class TestGwetAC1Config:
    """Tests for GwetAC1Config (no parametric null)."""

    def test_default_construction(self):
        cfg = GwetAC1Config()
        assert cfg.balanced_data_only is False

    def test_rejects_parametric_null_ci(self):
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            GwetAC1Config(ci_method="parametric_null")

    def test_rejects_parametric_null_hypothesis(self):
        ht = HypothesisTestConfig(method="parametric_null")
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            GwetAC1Config(hypothesis_test_config=ht)

    def test_allows_bca(self):
        cfg = GwetAC1Config(ci_method="bca")
        assert cfg.ci_method == "bca"


class TestBrennanPredigerConfig:
    """Tests for BrennanPredigerConfig (weightable, no parametric null)."""

    def test_default_construction(self):
        cfg = BrennanPredigerConfig()
        assert isinstance(cfg, BrennanPredigerConfig)

    def test_rejects_parametric_null_ci(self):
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            BrennanPredigerConfig(ci_method="parametric_null")

    def test_rejects_parametric_null_hypothesis(self):
        ht = HypothesisTestConfig(method="parametric_null")
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            BrennanPredigerConfig(hypothesis_test_config=ht)

class TestKrippendorffAlphaConfig:
    """Tests for KrippendorffAlphaConfig (weightable, no parametric null)."""

    def test_default_construction(self):
        cfg = KrippendorffAlphaConfig()
        assert isinstance(cfg, KrippendorffAlphaConfig)

    def test_rejects_parametric_null_ci(self):
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            KrippendorffAlphaConfig(ci_method="parametric_null")

    def test_rejects_parametric_null_hypothesis(self):
        ht = HypothesisTestConfig(method="parametric_null")
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            KrippendorffAlphaConfig(hypothesis_test_config=ht)

# ============================================================== #
#  Individual hierarchy: IndividualMarginalConfig → children      #
# ============================================================== #


class TestIndividualMarginalConfig:
    """Tests for IndividualMarginalConfig."""

    def test_default_construction(self):
        cfg = IndividualMarginalConfig()
        assert cfg.rater_jackknife is False

    def test_rater_jackknife_true_rejected(self):
        """Individual marginals don't support jackknife yet."""
        with pytest.raises(ValueError, match="rater_jackknife.*not supported in this version"):
            IndividualMarginalConfig(rater_jackknife=True)

    def test_summary_includes_rater_jackknife(self):
        cfg = IndividualMarginalConfig()
        s = cfg.summary()
        assert "rater_jackknife" in s
        assert s["rater_jackknife"] is False


class TestNoWeightIndividualConfig:
    """Tests for NoWeightIndividualConfig (rejects weights)."""

    def test_default_construction(self):
        cfg = NoWeightIndividualConfig()
        assert isinstance(cfg, NoWeightIndividualConfig)

class TestWeightableIndividualConfig:
    """Tests for WeightableIndividualConfig (accepts weights)."""

    def test_default_construction(self):
        cfg = WeightableIndividualConfig()
        assert isinstance(cfg, WeightableIndividualConfig)

# ============================================================== #
#  Leaf individual configs                                       #
# ============================================================== #


class TestCohenKappaConfig:
    """Tests for CohenKappaConfig (Finding 1a: warnings.warn)."""

    def test_default_construction(self):
        cfg = CohenKappaConfig()
        assert cfg.balanced_data_only is True

    def test_balanced_data_false_refused_with_reason(self):
        """Not rewritten to True: refused, saying why it does not apply."""
        with pytest.raises(ValueError, match="not available for Cohen's kappa.*subjects both raters rated"):
            CohenKappaConfig(balanced_data_only=False)

    def test_balanced_data_true_accepted(self):
        """The key a Fleiss config carries loads here too."""
        assert CohenKappaConfig.build_from_dict({"balanced_data_only": True}).balanced_data_only

    def test_balanced_data_not_bool(self):
        with pytest.raises(TypeError, match="balanced_data_only must be bool"):
            CohenKappaConfig(balanced_data_only="yes")

    def test_allows_parametric_null_ci(self):
        """Cohen has asymptotic null variance (Fleiss, Cohen & Everitt 1969).

        Allowed, and warned about: the null variance is the right one for
        a hypothesis test and the wrong one for an interval around a
        non-zero estimate. The message has to name the alternatives, or it
        tells someone they have a problem without telling them the fix.
        """
        with pytest.warns(UserWarning, match="parametric_null") as caught:
            cfg = CohenKappaConfig(ci_method="parametric_null")
        assert cfg.ci_method == "parametric_null"

        message = str(caught[0].message)
        assert "percentile" in message
        assert "bca" in message

class TestCongerKappaConfig:
    """Tests for CongerKappaConfig (Finding 1a + no parametric null)."""

    def test_default_construction(self):
        cfg = CongerKappaConfig()
        assert cfg.balanced_data_only is True

    def test_balanced_data_false_refused_with_reason(self):
        """Not rewritten to True: refused, pointing to what handles missing ratings."""
        with pytest.raises(ValueError, match="not available for Conger's kappa.*pairwise permutation"):
            CongerKappaConfig(balanced_data_only=False)

    def test_balanced_data_true_accepted(self):
        assert CongerKappaConfig.build_from_dict({"balanced_data_only": True}).balanced_data_only

    def test_balanced_data_not_bool(self):
        with pytest.raises(TypeError, match="balanced_data_only must be bool"):
            CongerKappaConfig(balanced_data_only="yes")

    def test_rejects_parametric_null_ci(self):
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            CongerKappaConfig(ci_method="parametric_null")

    def test_rejects_parametric_null_hypothesis(self):
        ht = HypothesisTestConfig(method="parametric_null")
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            CongerKappaConfig(hypothesis_test_config=ht)

class TestLightKappaConfig:
    """Tests for LightKappaConfig (no parametric null, complete data only)."""

    def test_default_balanced_data_true(self):
        assert LightKappaConfig().balanced_data_only is True

    def test_balanced_data_false_refused_with_reason(self):
        """Light (1971) is for complete data; the message says what to use instead."""
        with pytest.raises(ValueError, match="not available for Light's kappa.*Light 1971.*pairwise permutation"):
            LightKappaConfig(balanced_data_only=False)

    def test_balanced_data_not_bool(self):
        with pytest.raises(TypeError, match="balanced_data_only must be bool"):
            LightKappaConfig(balanced_data_only="no")

    def test_rejects_parametric_null_ci(self):
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            LightKappaConfig(ci_method="parametric_null")

    def test_rejects_parametric_null_hypothesis(self):
        ht = HypothesisTestConfig(method="parametric_null")
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            LightKappaConfig(hypothesis_test_config=ht)

# ============================================================== #
#  Permutation configs                                           #
# ============================================================== #


class TestPermutationBootstrapMeanConfig:
    """Tests for PermutationBootstrapMeanConfig: the method, and nothing else."""

    def test_construction(self):
        cfg = PermutationBootstrapMeanConfig(method="correlated")
        assert cfg.method == "correlated"

    def test_method_invalid_rejected(self):
        with pytest.raises(ValueError, match="method"):
            PermutationBootstrapMeanConfig(method="invalid")

    def test_summary_keys(self):
        cfg = PermutationBootstrapMeanConfig(method="independent")
        assert cfg.summary() == {"method": "independent"}

    def test_removed_key_refused_from_yaml(self):
        """A removed key says where it moved."""
        with pytest.raises(
            TypeError, match=r"\['n_bootstraps'\] were removed in 2\.1\.0.*n_bootstraps -> n_resamples"
        ):
            PermutationBootstrapMeanConfig.build_from_dict(
                {"method": "correlated", "n_bootstraps": 4000}
            )

    @pytest.mark.parametrize("key, value", [
        ("n_bootstraps", 4000), ("ci_method", "bca"), ("alpha", 0.05), ("random_state", 42),
    ])
    def test_removed_key_refused_from_python(self, key, value):
        """The same message when called from Python, not a generic TypeError."""
        with pytest.raises(TypeError, match="removed in 2.1.0.*config_of_coefficient"):
            PermutationBootstrapMeanConfig(method="correlated", **{key: value})

    def test_other_unknown_key_refused(self):
        with pytest.raises(ValueError, match="only key is 'method'"):
            PermutationBootstrapMeanConfig.build_from_dict({"method": "correlated", "n_draws": 10})

    def test_removed_fields_not_in_summary_or_repr(self):
        cfg = PermutationBootstrapMeanConfig(method="correlated")
        assert cfg.summary() == {"method": "correlated"}
        assert repr(cfg) == "PermutationBootstrapMeanConfig(method='correlated')"

    def test_from_yaml_shape(self):
        cfg = PermutationBootstrapMeanConfig.build_from_dict({"method": "correlated"})
        assert cfg == PermutationBootstrapMeanConfig(method="correlated")


class TestPermutationConfig:
    """Tests for PermutationConfig (base)."""

    def test_construction(self):
        cfg = PermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(),
        )
        assert cfg.coefficient_to_compute == "cohen_kappa"
        assert cfg.subset_k is None
        assert cfg.bootstrap_mean_config is None
        assert cfg.keep_draws is True
        from interrater.contracts import GUARD_MAX_EVALUATIONS
        assert cfg.max_evaluations == GUARD_MAX_EVALUATIONS

    @pytest.mark.parametrize("value, error", [
        (True, TypeError), (2.5, TypeError), ("1000", TypeError), (0, ValueError),
    ])
    def test_max_evaluations_refused(self, value, error):
        with pytest.raises(error, match="max_evaluations must be"):
            PermutationConfig(
                coefficient_to_compute="cohen_kappa",
                config_of_coefficient=CohenKappaConfig(),
                max_evaluations=value,
            )

    def test_keep_draws_non_bool_rejected(self):
        with pytest.raises(TypeError, match="keep_draws must be bool"):
            PermutationConfig(
                coefficient_to_compute="cohen_kappa",
                config_of_coefficient=CohenKappaConfig(),
                keep_draws=0,
            )

    def test_invalid_coefficient_rejected(self):
        with pytest.raises(ValueError, match="coefficient_to_compute"):
            PermutationConfig(
                coefficient_to_compute="not_a_real_coefficient",
                config_of_coefficient=CohenKappaConfig(),
            )

    def test_non_analyzer_config_rejected(self):
        with pytest.raises(TypeError, match="config_of_coefficient"):
            PermutationConfig(
                coefficient_to_compute="cohen_kappa",
                config_of_coefficient="not_a_config",
            )

    def test_summary_keys(self):
        cfg = PermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(),
        )
        s = cfg.summary()
        assert "coefficient_to_compute" in s
        assert "config_of_coefficient" in s
        assert "subset_k" in s
        assert "bootstrap_mean_config" in s
        assert s["keep_draws"] is True
        assert "max_evaluations" in s


class TestPairwisePermutationConfig:
    """Tests for PairwisePermutationConfig."""

    def test_construction(self):
        cfg = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(),
        )
        assert cfg.coefficient_to_compute == "cohen_kappa"

    def test_subset_k_rejected(self):
        """Pairwise is always k=2; subset_k must be None."""
        with pytest.raises(ValueError, match="subset_k must be None"):
            PairwisePermutationConfig(
                coefficient_to_compute="cohen_kappa",
                config_of_coefficient=CohenKappaConfig(),
                subset_k=3,
            )

    def test_config_type_mismatch_rejected(self):
        """cohen_kappa requires CohenKappaConfig, not FleissKappaConfig."""
        with pytest.raises(TypeError, match="CohenKappaConfig"):
            PairwisePermutationConfig(
                coefficient_to_compute="cohen_kappa",
                config_of_coefficient=FleissKappaConfig(),
            )

    def test_independent_bootstrap_requires_resamples(self):
        """Independent method needs at least 2 draws per pair."""
        bmc = PermutationBootstrapMeanConfig(method="independent")
        with pytest.raises(ValueError, match="independent.*n_resamples"):
            PairwisePermutationConfig(
                coefficient_to_compute="cohen_kappa",
                config_of_coefficient=CohenKappaConfig(n_resamples=1),
                bootstrap_mean_config=bmc,
            )

    @pytest.mark.parametrize("config_cls, coefficient, extra", [
        (PairwisePermutationConfig, "cohen_kappa", {}),
        (SubsetPermutationConfig, "fleiss_kappa", {"subset_k": 2}),
    ])
    def test_independent_with_bca_refused(self, config_cls, coefficient, extra):
        """No jackknife matches the independent mean, so BCa has no
        acceleration to use. Refused when the config is built."""
        coc_cls = CohenKappaConfig if coefficient == "cohen_kappa" else FleissKappaConfig
        with pytest.raises(ValueError, match="no jackknife matches it"):
            config_cls(
                coefficient_to_compute=coefficient,
                config_of_coefficient=coc_cls(ci_method="bca"),
                bootstrap_mean_config=PermutationBootstrapMeanConfig(method="independent"),
                **extra,
            )

    @pytest.mark.parametrize("config_cls, coefficient, extra", [
        (PairwisePermutationConfig, "cohen_kappa", {}),
        (SubsetPermutationConfig, "fleiss_kappa", {"subset_k": 2}),
    ])
    def test_parametric_null_refused(self, config_cls, coefficient, extra):
        """No closed-form variance describes correlated combinations.
        Refused when the config is built, not once the analysis runs."""
        coc_cls = CohenKappaConfig if coefficient == "cohen_kappa" else FleissKappaConfig
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # the config's own notice
            coc = coc_cls(ci_method="parametric_null")
        with pytest.raises(ValueError, match="'parametric_null' is not available"):
            config_cls(
                coefficient_to_compute=coefficient,
                config_of_coefficient=coc,
                **extra,
            )

    def test_independent_with_percentile_accepted(self):
        cfg = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(ci_method="percentile"),
            bootstrap_mean_config=PermutationBootstrapMeanConfig(method="independent"),
        )
        assert cfg.bootstrap_mean_config.method == "independent"

    def test_correlated_with_bca_accepted(self):
        cfg = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(ci_method="bca"),
            bootstrap_mean_config=PermutationBootstrapMeanConfig(method="correlated"),
        )
        assert cfg.config_of_coefficient.ci_method == "bca"

    def test_extra_key_in_mean_block_refused(self):
        """A removed key in the mean block says where it moved, read from a dict."""
        with pytest.raises(TypeError, match="removed in 2.1.0.*n_bootstraps -> n_resamples"):
            PairwisePermutationConfig.build_from_dict({
                "coefficient_to_compute": "cohen_kappa",
                "config_of_coefficient": {"n_resamples": 100},
                "bootstrap_mean_config": {
                    "method": "correlated", "n_bootstraps": 4000,
                },
            })

    def test_nested_removed_weights_key_refused(self):
        """An old key inside config_of_coefficient gets the 2.1.0 message too."""
        with pytest.raises(TypeError, match="removed in 2.1.0"):
            PairwisePermutationConfig.build_from_dict({
                "coefficient_to_compute": "cohen_kappa",
                "config_of_coefficient": {"weights_type": "identity"},
            })


class TestSubsetPermutationConfig:
    """Tests for SubsetPermutationConfig."""

    def test_construction_int(self):
        cfg = SubsetPermutationConfig(
            coefficient_to_compute="fleiss_kappa",
            config_of_coefficient=FleissKappaConfig(),
            subset_k=3,
        )
        assert cfg.subset_k == 3

    def test_construction_list(self):
        cfg = SubsetPermutationConfig(
            coefficient_to_compute="fleiss_kappa",
            config_of_coefficient=FleissKappaConfig(),
            subset_k=[2, 3, 4],
        )
        assert cfg.subset_k == [2, 3, 4]

    def test_construction_all(self):
        cfg = SubsetPermutationConfig(
            coefficient_to_compute="fleiss_kappa",
            config_of_coefficient=FleissKappaConfig(),
            subset_k="all",
        )
        assert cfg.subset_k == "all"

    def test_subset_k_none_rejected(self):
        with pytest.raises(ValueError, match="subset_k is required"):
            SubsetPermutationConfig(
                coefficient_to_compute="fleiss_kappa",
                config_of_coefficient=FleissKappaConfig(),
            )

    def test_subset_k_below_two_rejected(self):
        with pytest.raises(ValueError, match="subset_k must be >= 2"):
            SubsetPermutationConfig(
                coefficient_to_compute="fleiss_kappa",
                config_of_coefficient=FleissKappaConfig(),
                subset_k=1,
            )

    def test_subset_k_list_empty_rejected(self):
        with pytest.raises(ValueError, match="must not be empty"):
            SubsetPermutationConfig(
                coefficient_to_compute="fleiss_kappa",
                config_of_coefficient=FleissKappaConfig(),
                subset_k=[],
            )

    def test_subset_k_list_duplicates_rejected(self):
        with pytest.raises(ValueError, match="duplicates"):
            SubsetPermutationConfig(
                coefficient_to_compute="fleiss_kappa",
                config_of_coefficient=FleissKappaConfig(),
                subset_k=[3, 3],
            )

    def test_subset_k_invalid_string_rejected(self):
        with pytest.raises(ValueError, match="must be 'all'"):
            SubsetPermutationConfig(
                coefficient_to_compute="fleiss_kappa",
                config_of_coefficient=FleissKappaConfig(),
                subset_k="some",
            )

    def test_config_type_mismatch_rejected(self):
        with pytest.raises(TypeError, match="FleissKappaConfig"):
            SubsetPermutationConfig(
                coefficient_to_compute="fleiss_kappa",
                config_of_coefficient=CohenKappaConfig(),
                subset_k=3,
            )

    def test_nested_removed_weights_key_refused(self):
        """An old key inside config_of_coefficient gets the 2.1.0 message too."""
        with pytest.raises(TypeError, match="removed in 2.1.0"):
            SubsetPermutationConfig.build_from_dict({
                "coefficient_to_compute": "fleiss_kappa",
                "config_of_coefficient": {"custom_weights_path": None},
                "subset_k": 3,
            })

    def test_scotts_pi_refused(self):
        """Scott's pi is two-rater only: on a subset of k >= 3 it was
        Fleiss' kappa printed under Scott's name."""
        with pytest.raises(ValueError, match="Scott's pi is defined for exactly two raters"):
            SubsetPermutationConfig(
                coefficient_to_compute="scotts_pi",
                config_of_coefficient=ScottsConfig(),
                subset_k=2,
            )

    def test_scotts_pi_refused_from_yaml(self):
        """The same message from a YAML dict, before its keys are read."""
        with pytest.raises(ValueError, match="Use coefficient_to_compute: fleiss_kappa for subsets"):
            SubsetPermutationConfig.build_from_dict({
                "coefficient_to_compute": "scotts_pi",
                "config_of_coefficient": {},
                "subset_k": 3,
            })


# ============================================================== #
#  parametric_null: who warns, who raises, who stays silent      #
# ============================================================== #


# The three coefficients with a published closed-form null variance.
_PARAMETRIC_OK = [FleissKappaConfig, CohenKappaConfig, ScottsConfig]

# The five without one. Each rejects parametric_null with its own message.
_PARAMETRIC_REJECTED = [
    BrennanPredigerConfig,
    GwetAC1Config,
    KrippendorffAlphaConfig,
    CongerKappaConfig,
    LightKappaConfig,
]


class TestParametricNullWarningBoundary:
    """
    Where the parametric_null warning fires, and where it must not.

    The warning lives on AnalyzerConfig as a helper that AnalyzerConfig
    never calls. That is unusual enough to look like an oversight, and
    "fixing" it by wiring it into the base __post_init__ would break the
    five configs below: they call super() first and raise second, so they
    would be warned about a method they are about to refuse.

    These tests exist to make that refactor fail loudly instead of quietly
    producing a warning-then-error for five of eight coefficients.
    """

    @pytest.mark.parametrize("config_cls", _PARAMETRIC_OK)
    def test_accepting_configs_warn_exactly_once(self, config_cls):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            config_cls(ci_method="parametric_null")

        parametric = [
            w for w in caught if "parametric_null" in str(w.message)
        ]
        assert len(parametric) == 1
        assert issubclass(parametric[0].category, UserWarning)

    @pytest.mark.parametrize("config_cls", _PARAMETRIC_OK)
    @pytest.mark.parametrize("ci_method", ["percentile", "bca"])
    def test_accepting_configs_silent_on_other_methods(
        self, config_cls, ci_method
    ):
        """The warning is about the method, not about the config."""
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            config_cls(ci_method=ci_method)

        assert not [w for w in caught if "parametric_null" in str(w.message)]

    @pytest.mark.parametrize("config_cls", _PARAMETRIC_REJECTED)
    def test_rejecting_configs_raise_without_warning_first(self, config_cls):
        """Raise, and say nothing before raising.

        A warning here would be actively misleading -- it would explain
        how to interpret an interval that this coefficient cannot produce.
        """
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with pytest.raises(ValueError, match="parametric_null"):
                config_cls(ci_method="parametric_null")

        assert len(caught) == 0

    @pytest.mark.parametrize("config_cls", _PARAMETRIC_REJECTED)
    def test_rejecting_configs_keep_their_own_message(self, config_cls):
        """Each refusal explains itself in its own terms.

        These messages are coefficient-specific on purpose -- Krippendorff's
        says alpha has no published closed-form null variance. Collapsing
        them into one generic message would lose that, which is the usual
        cost of centralising validation.
        """
        with pytest.raises(ValueError) as exc_info:
            config_cls(ci_method="parametric_null")

        message = str(exc_info.value)
        assert config_cls.__name__.replace("Config", "") in message or (
            "not available" in message
        )


# -- Weights live in the semantics, not in any config ------------------------ #


_ALL_CONFIGS = sorted(
    (getattr(configs, name) for name in dir(configs) if name.endswith("Config")),
    key=lambda cls: cls.__name__,
)


class TestConfigsCarryNoWeights:
    """
    Since 2.1.0 the weights come from the semantics. No config declares
    them, every coefficient config builds without them, and the old way of
    passing a weights file is gone.
    """

    @pytest.mark.parametrize("config_cls", _ALL_CONFIGS, ids=lambda cls: cls.__name__)
    def test_no_weight_field(self, config_cls):
        fields = [f.name for f in dataclasses.fields(config_cls)]
        assert [name for name in fields if "weight" in name] == []

    @pytest.mark.parametrize("config_cls", [
        FleissKappaConfig, ScottsConfig, GwetAC1Config, BrennanPredigerConfig,
        KrippendorffAlphaConfig, CohenKappaConfig, CongerKappaConfig, LightKappaConfig,
    ], ids=lambda cls: cls.__name__)
    def test_coefficient_config_builds_without_weights(self, config_cls):
        summary = config_cls().summary()
        assert [key for key in summary if "weight" in key] == []

    def test_custom_weights_input_is_gone(self):
        import interrater

        assert not hasattr(interrater, "CustomWeightsInput")
        assert not hasattr(configs, "CustomWeightsInput")

    def test_weights_from_file_is_gone(self):
        from interrater.intermediates.estimator_parts import Weights

        assert not hasattr(Weights, "from_file")
        assert not hasattr(Weights, "build_from_file")



class TestUnknownKeysRefused:
    """
    build_from_dict refuses a misspelt key itself, naming the config, the
    nearest valid key and every valid key. Before, only build_from_yaml
    checked, and a dict -- or a config nested in a permutation config --
    reached the constructor and failed with Python's bare
    "unexpected keyword argument".
    """

    def test_coefficient_config(self):
        with pytest.raises(TypeError) as err:
            FleissKappaConfig.build_from_dict({"n_resample": 2000})
        message = str(err.value)
        assert message.startswith("[FleissKappaConfig] Error: unknown key(s) ['n_resample']")
        assert "did you mean 'n_resamples'?" in message
        assert "'balanced_data_only'" in message

    def test_nested_config_of_coefficient(self):
        with pytest.raises(TypeError, match=r"\[CohenKappaConfig\] Error: unknown key\(s\) \['ci_metod'\]"):
            PairwisePermutationConfig.build_from_dict({
                "coefficient_to_compute": "cohen_kappa",
                "config_of_coefficient": {"ci_metod": "bca"},
            })

    def test_permutation_config(self):
        with pytest.raises(TypeError, match=r"\[SubsetPermutationConfig\] Error: unknown key\(s\) \['subset_kk'\]"):
            SubsetPermutationConfig.build_from_dict({
                "coefficient_to_compute": "fleiss_kappa",
                "config_of_coefficient": {},
                "subset_kk": 3,
            })

    def test_hypothesis_test(self):
        with pytest.raises(TypeError, match=r"\[HypothesisTestConfig\] Error: unknown key\(s\) \['n_permutation'\]"):
            CohenKappaConfig.build_from_dict({"hypothesis_test": {"n_permutation": 100}})

    def test_aliases_still_accepted(self):
        cfg = CohenKappaConfig.build_from_dict({"hypothesis_test": {"n_permutations": 100}})
        assert cfg.hypothesis_test_config.n_permutations == 100
        perm = PairwisePermutationConfig.build_from_dict({
            "coefficient_to_compute": "cohen_kappa",
            "config_of_coefficient": {},
            "bootstrap_mean": {"method": "correlated"},
        })
        assert perm.bootstrap_mean_config.method == "correlated"

    def test_yaml_gives_the_same_message(self, tmp_path):
        path = tmp_path / "fleiss.yaml"
        path.write_text("n_resample: 2000\n")
        with pytest.raises(TypeError, match=r"\[FleissKappaConfig\] Error: unknown key\(s\) \['n_resample'\]"):
            FleissKappaConfig.build_from_yaml(path)
