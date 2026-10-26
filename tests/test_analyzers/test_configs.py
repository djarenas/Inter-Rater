"""Tests for analyzers/configs/ sub-package.

One test file for the entire sub-package (DandC convention).
Tests are organized by class, following the inheritance hierarchy.
"""

from __future__ import annotations

import warnings

import pytest

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
        with pytest.raises(ValueError, match="rater_jackknife.*not currently"):
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

    def test_balanced_data_forced_true_with_warning(self):
        """Setting balanced_data_only=False should warn and force True."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cfg = CohenKappaConfig(balanced_data_only=False)
            assert cfg.balanced_data_only is True
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)
            assert "balanced data" in str(w[0].message).lower()

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

    def test_balanced_data_forced_true_with_warning(self):
        """Setting balanced_data_only=False should warn and force True."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cfg = CongerKappaConfig(balanced_data_only=False)
            assert cfg.balanced_data_only is True
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)
            assert "balanced data" in str(w[0].message).lower()

    def test_rejects_parametric_null_ci(self):
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            CongerKappaConfig(ci_method="parametric_null")

    def test_rejects_parametric_null_hypothesis(self):
        ht = HypothesisTestConfig(method="parametric_null")
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            CongerKappaConfig(hypothesis_test_config=ht)

class TestLightKappaConfig:
    """Tests for LightKappaConfig (no parametric null, unbalanced OK)."""

    def test_default_balanced_data_false(self):
        cfg = LightKappaConfig()
        assert cfg.balanced_data_only is False

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
    """Tests for PermutationBootstrapMeanConfig."""

    def test_construction(self):
        cfg = PermutationBootstrapMeanConfig(
            method="correlated", n_bootstraps=500
        )
        assert cfg.method == "correlated"
        assert cfg.n_bootstraps == 500
        assert cfg.ci_method == "percentile"
        assert cfg.alpha == 0.05
        assert cfg.random_state is None

    def test_method_invalid_rejected(self):
        with pytest.raises(ValueError, match="method"):
            PermutationBootstrapMeanConfig(
                method="invalid", n_bootstraps=500
            )

    def test_n_bootstraps_zero_rejected(self):
        with pytest.raises(ValueError, match="n_bootstraps"):
            PermutationBootstrapMeanConfig(
                method="correlated", n_bootstraps=0
            )

    def test_ci_method_parametric_null_rejected(self):
        with pytest.raises(ValueError, match="parametric_null.*not available"):
            PermutationBootstrapMeanConfig(
                method="correlated",
                n_bootstraps=500,
                ci_method="parametric_null",
            )

    def test_independent_with_bca_refused(self):
        """No jackknife matches the independent mean, so BCa has no
        acceleration to use. Refused when the config is built."""
        with pytest.raises(ValueError, match="no jackknife matches it"):
            PermutationBootstrapMeanConfig(
                method="independent", n_bootstraps=100, ci_method="bca",
            )

    def test_independent_with_percentile_accepted(self):
        cfg = PermutationBootstrapMeanConfig(
            method="independent", n_bootstraps=100, ci_method="percentile",
        )
        assert cfg.ci_method == "percentile"

    def test_correlated_with_bca_accepted(self):
        cfg = PermutationBootstrapMeanConfig(
            method="correlated", n_bootstraps=100, ci_method="bca",
        )
        assert cfg.ci_method == "bca"

    def test_summary_keys(self):
        cfg = PermutationBootstrapMeanConfig(
            method="independent", n_bootstraps=1000
        )
        s = cfg.summary()
        assert set(s.keys()) == {
            "method", "n_bootstraps", "ci_method", "alpha", "random_state",
        }


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
        """Independent method needs per-pair bootstrap (n_resamples >= 2)."""
        bmc = PermutationBootstrapMeanConfig(
            method="independent", n_bootstraps=100
        )
        with pytest.raises(ValueError, match="independent.*n_resamples"):
            PairwisePermutationConfig(
                coefficient_to_compute="cohen_kappa",
                config_of_coefficient=CohenKappaConfig(n_resamples=1),
                bootstrap_mean_config=bmc,
            )

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
