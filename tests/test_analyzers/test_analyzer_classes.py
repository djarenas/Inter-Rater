"""Tests for analyzer classes -- hierarchy, construction, guards, smoke estimates.

Deterministic tests only where possible. Smoke estimates use n_resamples=10
for speed. No numerical correctness assertions on point estimates -- those
belong in the agreement-function tests. This file exercises the analyzer
layer: config validation, constructor paths, estimate guards, rater
constraints, hypothesis-test guards, and result propagation.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

# ------------------------------------------------------------------ #
# Fixtures                                                             #
# ------------------------------------------------------------------ #

from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.data_objects.ratings import Ratings
from interrater.intermediates.estimator_parts import NijMatrix
from interrater.intermediates.results.ci_results import (
    BootstrapCIResult,
    ParametricCIResult,
)


def _assert_bootstrap_interval(
    result, n_resamples: int, method: str = "percentile"
) -> None:
    """
    Assert the interval is a real one, not merely present.

    These smoke tests used to check ``isinstance(ci, tuple)``, which a pair
    of NaNs passed as happily as a working interval. An interval now states
    what it is, so the smoke test states what it expects: the method that
    was asked for, the number of draws it saw, and that bounds exist at all.

    Containment of the point estimate is deliberately not asserted. It is
    not an invariant -- strong bias correction can legitimately break it --
    and with n_resamples=10 it would be flaky besides.
    """
    ci = result.confidence_interval
    assert isinstance(ci, BootstrapCIResult)
    assert ci.method_intended == method
    assert ci.n_resamples == n_resamples
    assert 0 <= ci.n_failed <= n_resamples
    assert ci.is_defined
    assert ci.lower <= ci.upper


def _assert_parametric_interval(result) -> None:
    """The parametric path returns a different kind of interval entirely."""
    ci = result.confidence_interval
    assert isinstance(ci, ParametricCIResult)
    assert not isinstance(ci, BootstrapCIResult)
    assert ci.variance_source == "null"
    assert ci.is_defined
    assert ci.lower <= ci.upper


def _make_ratings(n_raters: int = 3) -> Ratings:
    """Build a minimal balanced Ratings fixture."""
    cats = UnorderedCategories(
        names=("A", "B"), labels={"A": "A", "B": "B"}
    )
    sem = RatingsSemantics(raters=n_raters, categories=cats)
    if n_raters == 2:
        ci = np.array(
            [[0, 1], [1, 0], [0, 0], [1, 1]],
            dtype=np.int8,
        )
    elif n_raters == 3:
        ci = np.array(
            [[0, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 0]],
            dtype=np.int8,
        )
    else:
        raise ValueError(f"Fixture does not support n_raters={n_raters}")
    return Ratings(category_index=ci, semantics=sem)


def _make_nij(raters_per_subject: int = 3) -> NijMatrix:
    """Build a minimal NijMatrix fixture with given raters per subject."""
    if raters_per_subject == 2:
        arr = np.array([[2, 0], [0, 2], [1, 1], [1, 1]])
    elif raters_per_subject == 3:
        arr = np.array([[2, 1], [1, 2], [1, 2], [2, 1]])
    else:
        raise ValueError(
            f"Fixture does not support raters_per_subject={raters_per_subject}"
        )
    return NijMatrix(matrix=arr)


# ------------------------------------------------------------------ #
# Module-level _ERROR_PREFIX                                           #
# ------------------------------------------------------------------ #


_ANALYZER_MODULES = [
    ("interrater.analyzers.analyzer", "[Analyzer] Error"),
    ("interrater.analyzers.group_analyzer", "[GroupAnalyzer] Error"),
    (
        "interrater.analyzers.weightable_group_analyzer",
        "[WeightableGroupAnalyzer] Error",
    ),
    (
        "interrater.analyzers.individual_marginal_analyzer",
        "[IndividualMarginalAnalyzer] Error",
    ),
    (
        "interrater.analyzers.weightable_individual_analyzer",
        "[WeightableIndividualAnalyzer] Error",
    ),
    (
        "interrater.analyzers.fleiss_kappa_analyzer",
        "[FleissKappaAnalyzer] Error",
    ),
    ("interrater.analyzers.gwet_ac1_analyzer", "[GwetAC1Analyzer] Error"),
    (
        "interrater.analyzers.brennan_prediger_analyzer",
        "[BrennanPredigerAnalyzer] Error",
    ),
    (
        "interrater.analyzers.krippendorff_alpha_analyzer",
        "[KrippendorffAlphaAnalyzer] Error",
    ),
    ("interrater.analyzers.scotts_analyzer", "[ScottsAnalyzer] Error"),
    (
        "interrater.analyzers.cohen_kappa_analyzer",
        "[CohenKappaAnalyzer] Error",
    ),
    (
        "interrater.analyzers.conger_kappa_analyzer",
        "[CongerKappaAnalyzer] Error",
    ),
    (
        "interrater.analyzers.light_kappa_analyzer",
        "[LightKappaAnalyzer] Error",
    ),
]


class TestErrorPrefixes:
    """Every analyzer module exposes _ERROR_PREFIX at module level."""

    @pytest.mark.parametrize(
        "module_path, expected_prefix",
        _ANALYZER_MODULES,
        ids=[m.rsplit(".", 1)[-1] for m, _ in _ANALYZER_MODULES],
    )
    def test_error_prefix_at_module_level(self, module_path, expected_prefix):
        import importlib

        mod = importlib.import_module(module_path)
        assert hasattr(mod, "_ERROR_PREFIX")
        assert mod._ERROR_PREFIX == expected_prefix


# ------------------------------------------------------------------ #
# Class hierarchy                                                      #
# ------------------------------------------------------------------ #


class TestClassHierarchy:
    """Verify inheritance chains match the design."""

    def test_fleiss_extends_weightable_group(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.weightable_group_analyzer import (
            WeightableGroupAnalyzer,
        )

        assert issubclass(FleissKappaAnalyzer, WeightableGroupAnalyzer)

    def test_brennan_prediger_extends_weightable_group(self):
        from interrater.analyzers.brennan_prediger_analyzer import (
            BrennanPredigerAnalyzer,
        )
        from interrater.analyzers.weightable_group_analyzer import (
            WeightableGroupAnalyzer,
        )

        assert issubclass(BrennanPredigerAnalyzer, WeightableGroupAnalyzer)

    def test_krippendorff_extends_weightable_group(self):
        from interrater.analyzers.krippendorff_alpha_analyzer import (
            KrippendorffAlphaAnalyzer,
        )
        from interrater.analyzers.weightable_group_analyzer import (
            WeightableGroupAnalyzer,
        )

        assert issubclass(KrippendorffAlphaAnalyzer, WeightableGroupAnalyzer)

    def test_gwet_extends_group_not_weightable(self):
        """AC1 is unweighted -- extends GroupAnalyzer directly."""
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.group_analyzer import GroupAnalyzer
        from interrater.analyzers.weightable_group_analyzer import (
            WeightableGroupAnalyzer,
        )

        assert issubclass(GwetAC1Analyzer, GroupAnalyzer)
        assert not issubclass(GwetAC1Analyzer, WeightableGroupAnalyzer)

    def test_scott_extends_group_not_weightable(self):
        """Scott's pi is unweighted -- extends GroupAnalyzer directly."""
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.group_analyzer import GroupAnalyzer
        from interrater.analyzers.weightable_group_analyzer import (
            WeightableGroupAnalyzer,
        )

        assert issubclass(ScottsAnalyzer, GroupAnalyzer)
        assert not issubclass(ScottsAnalyzer, WeightableGroupAnalyzer)

    def test_cohen_extends_weightable_individual(self):
        from interrater.analyzers.cohen_kappa_analyzer import CohenKappaAnalyzer
        from interrater.analyzers.weightable_individual_analyzer import (
            WeightableIndividualAnalyzer,
        )

        assert issubclass(CohenKappaAnalyzer, WeightableIndividualAnalyzer)

    def test_conger_extends_individual_not_weightable(self):
        from interrater.analyzers.conger_kappa_analyzer import CongerKappaAnalyzer
        from interrater.analyzers.individual_marginal_analyzer import (
            IndividualMarginalAnalyzer,
        )
        from interrater.analyzers.weightable_individual_analyzer import (
            WeightableIndividualAnalyzer,
        )

        assert issubclass(CongerKappaAnalyzer, IndividualMarginalAnalyzer)
        assert not issubclass(CongerKappaAnalyzer, WeightableIndividualAnalyzer)

    def test_light_extends_individual_not_weightable(self):
        from interrater.analyzers.light_kappa_analyzer import LightKappaAnalyzer
        from interrater.analyzers.individual_marginal_analyzer import (
            IndividualMarginalAnalyzer,
        )
        from interrater.analyzers.weightable_individual_analyzer import (
            WeightableIndividualAnalyzer,
        )

        assert issubclass(LightKappaAnalyzer, IndividualMarginalAnalyzer)
        assert not issubclass(LightKappaAnalyzer, WeightableIndividualAnalyzer)


# ------------------------------------------------------------------ #
# Config type validation                                               #
# ------------------------------------------------------------------ #


class TestConfigTypeValidation:
    """Each leaf analyzer rejects a config from a different analyzer."""

    def test_fleiss_rejects_gwet_config(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import GwetAC1Config

        with pytest.raises(TypeError, match=r"\[FleissKappaAnalyzer\] Error"):
            FleissKappaAnalyzer(_make_ratings(3), GwetAC1Config())

    def test_gwet_rejects_fleiss_config(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import FleissKappaConfig

        with pytest.raises(TypeError, match=r"\[GwetAC1Analyzer\] Error"):
            GwetAC1Analyzer(_make_ratings(3), FleissKappaConfig())

    def test_brennan_prediger_rejects_fleiss_config(self):
        from interrater.analyzers.brennan_prediger_analyzer import (
            BrennanPredigerAnalyzer,
        )
        from interrater.analyzers.configs import FleissKappaConfig

        with pytest.raises(
            TypeError, match=r"\[BrennanPredigerAnalyzer\] Error"
        ):
            BrennanPredigerAnalyzer(_make_ratings(3), FleissKappaConfig())

    def test_krippendorff_rejects_fleiss_config(self):
        from interrater.analyzers.krippendorff_alpha_analyzer import (
            KrippendorffAlphaAnalyzer,
        )
        from interrater.analyzers.configs import FleissKappaConfig

        with pytest.raises(
            TypeError, match=r"\[KrippendorffAlphaAnalyzer\] Error"
        ):
            KrippendorffAlphaAnalyzer(_make_ratings(3), FleissKappaConfig())

    def test_scott_rejects_gwet_config(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import GwetAC1Config

        with pytest.raises(TypeError, match=r"\[ScottsAnalyzer\] Error"):
            ScottsAnalyzer(_make_ratings(2), GwetAC1Config())

    def test_cohen_rejects_conger_config(self):
        from interrater.analyzers.cohen_kappa_analyzer import CohenKappaAnalyzer
        from interrater.analyzers.configs import CongerKappaConfig

        with pytest.raises(TypeError, match=r"\[CohenKappaAnalyzer\] Error"):
            CohenKappaAnalyzer(_make_ratings(2), CongerKappaConfig())

    def test_conger_rejects_cohen_config(self):
        from interrater.analyzers.conger_kappa_analyzer import CongerKappaAnalyzer
        from interrater.analyzers.configs import CohenKappaConfig

        with pytest.raises(TypeError, match=r"\[CongerKappaAnalyzer\] Error"):
            CongerKappaAnalyzer(_make_ratings(3), CohenKappaConfig())

    def test_light_rejects_cohen_config(self):
        from interrater.analyzers.light_kappa_analyzer import LightKappaAnalyzer
        from interrater.analyzers.configs import CohenKappaConfig

        with pytest.raises(TypeError, match=r"\[LightKappaAnalyzer\] Error"):
            LightKappaAnalyzer(_make_ratings(3), CohenKappaConfig())


# ------------------------------------------------------------------ #
# Constructor paths -- __init__ and build_from_nij_matrix              #
# ------------------------------------------------------------------ #


class TestConstructorPaths:
    """Group analyzers support both __init__ and build_from_nij_matrix."""

    def test_fleiss_from_ratings(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        a = FleissKappaAnalyzer(_make_ratings(3), FleissKappaConfig())
        assert a.ratings is not None
        assert a.nij_matrix is None

    def test_fleiss_from_nij_matrix(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        a = FleissKappaAnalyzer.build_from_nij_matrix(
            _make_nij(3), FleissKappaConfig()
        )
        assert a.ratings is None
        assert a.nij_matrix is not None

    def test_gwet_from_ratings(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import GwetAC1Config

        a = GwetAC1Analyzer(_make_ratings(3), GwetAC1Config())
        assert a.ratings is not None

    def test_gwet_from_nij_matrix(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import GwetAC1Config

        a = GwetAC1Analyzer.build_from_nij_matrix(_make_nij(3), GwetAC1Config())
        assert a.ratings is None
        assert a.nij_matrix is not None

    def test_brennan_prediger_from_ratings(self):
        from interrater.analyzers.brennan_prediger_analyzer import (
            BrennanPredigerAnalyzer,
        )
        from interrater.analyzers.configs import BrennanPredigerConfig

        a = BrennanPredigerAnalyzer(_make_ratings(3), BrennanPredigerConfig())
        assert a.ratings is not None

    def test_brennan_prediger_from_nij_matrix(self):
        from interrater.analyzers.brennan_prediger_analyzer import (
            BrennanPredigerAnalyzer,
        )
        from interrater.analyzers.configs import BrennanPredigerConfig

        a = BrennanPredigerAnalyzer.build_from_nij_matrix(
            _make_nij(3), BrennanPredigerConfig()
        )
        assert a.ratings is None
        assert a.nij_matrix is not None

    def test_krippendorff_from_ratings(self):
        from interrater.analyzers.krippendorff_alpha_analyzer import (
            KrippendorffAlphaAnalyzer,
        )
        from interrater.analyzers.configs import KrippendorffAlphaConfig

        a = KrippendorffAlphaAnalyzer(_make_ratings(3), KrippendorffAlphaConfig())
        assert a.ratings is not None

    def test_krippendorff_from_nij_matrix(self):
        from interrater.analyzers.krippendorff_alpha_analyzer import (
            KrippendorffAlphaAnalyzer,
        )
        from interrater.analyzers.configs import KrippendorffAlphaConfig

        a = KrippendorffAlphaAnalyzer.build_from_nij_matrix(
            _make_nij(3), KrippendorffAlphaConfig()
        )
        assert a.ratings is None
        assert a.nij_matrix is not None

    def test_scott_from_ratings(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        a = ScottsAnalyzer(_make_ratings(2), ScottsConfig())
        assert a.ratings is not None

    def test_scott_from_nij_matrix(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        a = ScottsAnalyzer.build_from_nij_matrix(_make_nij(2), ScottsConfig())
        assert a.ratings is None
        assert a.nij_matrix is not None

    def test_individual_analyzers_have_no_nij_path(self):
        """Cohen, Conger, Light -- no build_from_nij_matrix classmethod."""
        from interrater.analyzers.cohen_kappa_analyzer import CohenKappaAnalyzer
        from interrater.analyzers.conger_kappa_analyzer import CongerKappaAnalyzer
        from interrater.analyzers.light_kappa_analyzer import LightKappaAnalyzer

        # These classes should NOT have build_from_nij_matrix as their own method.
        # They do inherit it from Analyzer via GroupAnalyzer, but individual
        # analyzers do not extend GroupAnalyzer, so it should not exist.
        assert not hasattr(CohenKappaAnalyzer, "build_from_nij_matrix")
        assert not hasattr(CongerKappaAnalyzer, "build_from_nij_matrix")
        assert not hasattr(LightKappaAnalyzer, "build_from_nij_matrix")


# ------------------------------------------------------------------ #
# Estimate guards                                                      #
# ------------------------------------------------------------------ #


class TestEstimateGuards:
    """Guard on estimate() when ratings is None, and vice versa."""

    def test_fleiss_estimate_raises_when_built_from_nij(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        a = FleissKappaAnalyzer.build_from_nij_matrix(
            _make_nij(3), FleissKappaConfig()
        )
        with pytest.raises(ValueError, match="ratings"):
            a.estimate()

    def test_fleiss_estimate_from_nij_raises_when_built_from_ratings(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        a = FleissKappaAnalyzer(_make_ratings(3), FleissKappaConfig())
        with pytest.raises(ValueError, match="nij_matrix"):
            a.estimate_from_nij_matrix()

    def test_gwet_estimate_raises_when_built_from_nij(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import GwetAC1Config

        a = GwetAC1Analyzer.build_from_nij_matrix(_make_nij(3), GwetAC1Config())
        with pytest.raises(ValueError, match="ratings"):
            a.estimate()

    def test_gwet_estimate_from_nij_raises_when_built_from_ratings(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import GwetAC1Config

        a = GwetAC1Analyzer(_make_ratings(3), GwetAC1Config())
        with pytest.raises(ValueError, match="nij_matrix"):
            a.estimate_from_nij_matrix()

    def test_brennan_prediger_estimate_raises_when_built_from_nij(self):
        from interrater.analyzers.brennan_prediger_analyzer import (
            BrennanPredigerAnalyzer,
        )
        from interrater.analyzers.configs import BrennanPredigerConfig

        a = BrennanPredigerAnalyzer.build_from_nij_matrix(
            _make_nij(3), BrennanPredigerConfig()
        )
        with pytest.raises(ValueError, match="ratings"):
            a.estimate()

    def test_krippendorff_estimate_raises_when_built_from_nij(self):
        from interrater.analyzers.krippendorff_alpha_analyzer import (
            KrippendorffAlphaAnalyzer,
        )
        from interrater.analyzers.configs import KrippendorffAlphaConfig

        a = KrippendorffAlphaAnalyzer.build_from_nij_matrix(
            _make_nij(3), KrippendorffAlphaConfig()
        )
        with pytest.raises(ValueError, match="ratings"):
            a.estimate()

    def test_scott_estimate_raises_when_built_from_nij(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        a = ScottsAnalyzer.build_from_nij_matrix(_make_nij(2), ScottsConfig())
        with pytest.raises(ValueError, match="ratings"):
            a.estimate()


# ------------------------------------------------------------------ #
# Rater constraints                                                    #
# ------------------------------------------------------------------ #


class TestRaterConstraints:
    """Cohen and Scott require exactly 2 raters."""

    def test_cohen_rejects_three_raters(self):
        from interrater.analyzers.cohen_kappa_analyzer import CohenKappaAnalyzer
        from interrater.analyzers.configs import CohenKappaConfig

        with pytest.raises(ValueError):
            CohenKappaAnalyzer(_make_ratings(3), CohenKappaConfig())

    def test_cohen_accepts_two_raters(self):
        from interrater.analyzers.cohen_kappa_analyzer import CohenKappaAnalyzer
        from interrater.analyzers.configs import CohenKappaConfig

        a = CohenKappaAnalyzer(_make_ratings(2), CohenKappaConfig())
        assert a.ratings is not None

    def test_scott_rejects_three_raters_from_ratings(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        with pytest.raises(ValueError):
            ScottsAnalyzer(_make_ratings(3), ScottsConfig())

    def test_scott_accepts_two_raters_from_ratings(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        a = ScottsAnalyzer(_make_ratings(2), ScottsConfig())
        assert a.ratings is not None

    def test_scott_rejects_three_raters_from_nij(self):
        """Row sums = 3 means 3 raters -- Scott rejects."""
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        with pytest.raises(ValueError):
            ScottsAnalyzer.build_from_nij_matrix(_make_nij(3), ScottsConfig())

    def test_scott_accepts_two_raters_from_nij(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        a = ScottsAnalyzer.build_from_nij_matrix(_make_nij(2), ScottsConfig())
        assert a.nij_matrix is not None


# ------------------------------------------------------------------ #
# Hypothesis test guards                                               #
# ------------------------------------------------------------------ #


class TestHypothesisTestGuards:
    """Hypothesis test guards: no config, wrong method, no ratings."""

    def test_fleiss_no_ht_config_raises(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        a = FleissKappaAnalyzer(_make_ratings(3), FleissKappaConfig())
        with pytest.raises(ValueError, match="hypothesis_test_config"):
            a.hypothesis_test()

    def test_gwet_no_ht_config_raises(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import GwetAC1Config

        a = GwetAC1Analyzer(_make_ratings(3), GwetAC1Config())
        with pytest.raises(ValueError, match="hypothesis_test_config"):
            a.hypothesis_test()

    def test_gwet_rejects_parametric_null_at_config_level(self):
        """AC1 rejects parametric_null in config __post_init__."""
        from interrater.analyzers.configs import GwetAC1Config
        from interrater.analyzers.configs.hypothesis_test_config import (
            HypothesisTestConfig,
        )

        with pytest.raises(ValueError, match="parametric_null"):
            GwetAC1Config(
                hypothesis_test_config=HypothesisTestConfig(
                    method="parametric_null"
                )
            )

    def test_light_rejects_parametric_null_at_config_level(self):
        """Light rejects parametric_null in config __post_init__."""
        from interrater.analyzers.configs import LightKappaConfig
        from interrater.analyzers.configs.hypothesis_test_config import (
            HypothesisTestConfig,
        )

        with pytest.raises(ValueError, match="parametric_null"):
            LightKappaConfig(
                hypothesis_test_config=HypothesisTestConfig(
                    method="parametric_null"
                )
            )

    def test_bp_rejects_parametric_null_at_config_level(self):
        """BP rejects parametric_null in config __post_init__."""
        from interrater.analyzers.configs import BrennanPredigerConfig
        from interrater.analyzers.configs.hypothesis_test_config import (
            HypothesisTestConfig,
        )

        with pytest.raises(ValueError, match="parametric_null"):
            BrennanPredigerConfig(
                hypothesis_test_config=HypothesisTestConfig(
                    method="parametric_null"
                )
            )

    def test_krippendorff_rejects_parametric_null_at_config_level(self):
        """Krippendorff rejects parametric_null in config __post_init__."""
        from interrater.analyzers.configs import KrippendorffAlphaConfig
        from interrater.analyzers.configs.hypothesis_test_config import (
            HypothesisTestConfig,
        )

        with pytest.raises(ValueError, match="parametric_null"):
            KrippendorffAlphaConfig(
                hypothesis_test_config=HypothesisTestConfig(
                    method="parametric_null"
                )
            )

    def test_conger_rejects_parametric_null_at_config_level(self):
        """Conger rejects parametric_null in config __post_init__."""
        from interrater.analyzers.configs import CongerKappaConfig
        from interrater.analyzers.configs.hypothesis_test_config import (
            HypothesisTestConfig,
        )

        with pytest.raises(ValueError, match="parametric_null"):
            CongerKappaConfig(
                hypothesis_test_config=HypothesisTestConfig(
                    method="parametric_null"
                )
            )

    def test_permutation_from_nij_raises(self):
        """Permutation requires raw ratings, not NijMatrix."""
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig
        from interrater.analyzers.configs.hypothesis_test_config import (
            HypothesisTestConfig,
        )

        config = FleissKappaConfig(
            hypothesis_test_config=HypothesisTestConfig(method="permutation")
        )
        a = FleissKappaAnalyzer.build_from_nij_matrix(_make_nij(3), config)
        with pytest.raises(ValueError, match="Ratings"):
            a.hypothesis_test()


# ------------------------------------------------------------------ #
# Hypothesis test on an undefined coefficient                          #
# ------------------------------------------------------------------ #


def _all_same_ratings(n_raters: int) -> Ratings:
    """Every rater gives every subject the same answer."""
    cats = UnorderedCategories(names=("A", "B"), labels={"A": "A", "B": "B"})
    sem = RatingsSemantics(raters=n_raters, categories=cats)
    ci = np.zeros((6, n_raters), dtype=np.int8)
    return Ratings(category_index=ci, semantics=sem)


class TestHypothesisTestUndefinedCoefficient:
    """
    When every rating is the same, six coefficients are undefined (NaN).
    Their permutation test used to report p = 0 and reject; it must refuse.
    Gwet's AC1 and Brennan-Prediger are defined there (1.0) and still run.
    """

    @pytest.mark.parametrize("name, n_raters", [
        ("FleissKappa", 3), ("Scotts", 2), ("KrippendorffAlpha", 3),
        ("CohenKappa", 2), ("CongerKappa", 3), ("LightKappa", 3),
    ])
    # Light averages its pair values; all NaN here, so numpy says so.
    @pytest.mark.filterwarnings("ignore:Mean of empty slice:RuntimeWarning")
    def test_undefined_refused(self, name, n_raters):
        import interrater.analyzers as A
        import interrater.analyzers.configs as C
        from interrater.analyzers.configs import HypothesisTestConfig

        ht = HypothesisTestConfig(method="permutation", n_permutations=20)
        config = getattr(C, f"{name}Config")(hypothesis_test_config=ht)
        a = getattr(A, f"{name}Analyzer")(_all_same_ratings(n_raters), config)
        with pytest.raises(ValueError, match="undefined"):
            a.hypothesis_test()

    @pytest.mark.parametrize("name", ["GwetAC1", "BrennanPrediger"])
    def test_defined_still_runs(self, name):
        import interrater.analyzers as A
        import interrater.analyzers.configs as C
        from interrater.analyzers.configs import HypothesisTestConfig

        ht = HypothesisTestConfig(method="permutation", n_permutations=20)
        config = getattr(C, f"{name}Config")(hypothesis_test_config=ht)
        result = getattr(A, f"{name}Analyzer")(_all_same_ratings(3), config).hypothesis_test()
        assert result.point == 1.0
        assert not result.reject_null


# ------------------------------------------------------------------ #
# parametric_null only where its formula applies                       #
# ------------------------------------------------------------------ #


def _three_cat_ratings(n_raters: int, ordered: bool = False, unbalanced: bool = False) -> Ratings:
    """40 subjects, 3 categories; raters copy a true answer 70% of the time."""
    from interrater.semantics.categories import OrderedCategories

    rng = np.random.default_rng(42)
    truth = rng.integers(0, 3, 40)
    codes = np.where(
        rng.random((40, n_raters)) < 0.7, truth[:, None], rng.integers(0, 3, (40, n_raters)),
    ).astype(np.int8)
    if unbalanced:
        codes[:5, -1] = -1  # five subjects lose their last rater
    names, labels = ("A", "B", "C"), {"A": "A", "B": "B", "C": "C"}
    cats = (OrderedCategories(names=names, labels=labels, weighting="linear") if ordered
            else UnorderedCategories(names=names, labels=labels))
    return Ratings(category_index=codes, semantics=RatingsSemantics(raters=n_raters, categories=cats))


@pytest.mark.filterwarnings("ignore::UserWarning")  # parametric_null warns; not the subject
class TestParametricNullNeedsItsAssumptions:
    """
    The null variances behind parametric_null assume every subject has the
    same number of raters, and Fleiss' (Eq. 12) has no weighted form. Both
    the interval and the hypothesis test are refused at construction when
    the data break that; Cohen's weighted null variance is published, so
    weighted Cohen still runs.
    """

    @staticmethod
    def _configs(config_cls, **kw):
        from interrater.analyzers.configs import HypothesisTestConfig

        return [
            config_cls(ci_method="parametric_null", **kw),
            config_cls(hypothesis_test_config=HypothesisTestConfig(method="parametric_null"), **kw),
        ]

    @pytest.mark.parametrize("which", [0, 1], ids=["interval", "test"])
    def test_fleiss_unbalanced_refused(self, which):
        from interrater.analyzers import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        config = self._configs(FleissKappaConfig, balanced_data_only=False)[which]
        with pytest.raises(ValueError, match="same number of raters; these data have 2 to 3"):
            FleissKappaAnalyzer(_three_cat_ratings(3, unbalanced=True), config)

    def test_fleiss_unbalanced_nij_refused(self):
        from interrater.analyzers import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        nij = NijMatrix(matrix=np.array([[3, 0], [2, 1], [1, 1], [0, 2]]))
        config = self._configs(FleissKappaConfig, balanced_data_only=False)[0]
        with pytest.raises(ValueError, match="same number of raters"):
            FleissKappaAnalyzer.build_from_nij_matrix(nij, config)

    @pytest.mark.parametrize("which", [0, 1], ids=["interval", "test"])
    def test_fleiss_weighted_refused(self, which):
        from interrater.analyzers import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        with pytest.raises(ValueError, match="unweighted null variance"):
            FleissKappaAnalyzer(_three_cat_ratings(3, ordered=True), self._configs(FleissKappaConfig)[which])

    def test_scott_unbalanced_refused_with_the_right_reason(self):
        """It used to fail by dividing by zero, blaming 'all ratings in one category'."""
        from interrater.analyzers import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        config = self._configs(ScottsConfig, balanced_data_only=False)[1]
        with pytest.raises(ValueError, match="same number of raters; these data have 1 to 2"):
            ScottsAnalyzer(_three_cat_ratings(2, unbalanced=True), config)

    @pytest.mark.parametrize("which", [0, 1], ids=["interval", "test"])
    def test_balanced_unweighted_fleiss_still_runs(self, which):
        from interrater.analyzers import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        FleissKappaAnalyzer(_three_cat_ratings(3), self._configs(FleissKappaConfig)[which])

    @pytest.mark.parametrize("which", [0, 1], ids=["interval", "test"])
    def test_weighted_cohen_still_runs(self, which):
        from interrater.analyzers import CohenKappaAnalyzer
        from interrater.analyzers.configs import CohenKappaConfig

        analyzer = CohenKappaAnalyzer(_three_cat_ratings(2, ordered=True), self._configs(CohenKappaConfig)[which])
        result = analyzer.estimate() if which == 0 else analyzer.hypothesis_test()
        assert np.isfinite(result.point)

    def test_other_methods_on_unbalanced_data_unaffected(self):
        from interrater.analyzers import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        FleissKappaAnalyzer(
            _three_cat_ratings(3, unbalanced=True),
            FleissKappaConfig(balanced_data_only=False, n_resamples=20),
        ).estimate()


# ------------------------------------------------------------------ #
# p-values: the shared helpers and what every coefficient reports      #
# ------------------------------------------------------------------ #


def _fleiss_analyzer():
    from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
    from interrater.analyzers.configs import FleissKappaConfig

    return FleissKappaAnalyzer(_make_ratings(3), FleissKappaConfig())


class TestPValueHelpers:
    """The two helpers every analyzer now calls."""

    def test_nan_draws_skipped_and_counted(self):
        """3 of 5 draws usable, none at or above 0.35: greater = 1/4."""
        null = np.array([0.1, 0.2, np.nan, 0.3, np.nan])
        d = _fleiss_analyzer()._permutation_p_value(null, 0.35, "two-sided")
        assert (d["n_usable"], d["n_at_or_above"], d["n_at_or_below"]) == (3, 0, 3)
        assert d["p_value"] == pytest.approx(2 * 1 / 4)

    def test_all_nan_draws_refused(self):
        null = np.array([np.nan, np.nan])
        with pytest.raises(ValueError, match="every permutation"):
            _fleiss_analyzer()._permutation_p_value(null, 0.35, "two-sided")

    @pytest.mark.parametrize("alternative, expected", [
        ("greater", 1 / 100), ("two-sided", 2 / 100), ("less", 1.0),
    ])
    def test_observed_beyond_every_draw(self, alternative, expected):
        """99 draws, all below the observed value."""
        null = np.linspace(-0.1, 0.1, 99)
        d = _fleiss_analyzer()._permutation_p_value(null, 0.5, alternative)
        assert d["p_value"] == pytest.approx(expected)

    @pytest.mark.parametrize("alternative, expected", [
        ("two-sided", 0.0455), ("greater", 0.0228), ("less", 0.9772),
    ])
    def test_parametric_z_of_two(self, alternative, expected):
        from interrater.analyzers.analyzer import Analyzer

        assert Analyzer._parametric_p_value(2.0, alternative) == pytest.approx(expected, abs=1e-4)


class TestPValueAgreesAcrossCoefficients:
    """
    With each rater's answer mix held fixed and no missing cells, Fleiss,
    Gwet, Brennan-Prediger, Krippendorff and Conger are all increasing
    rescalings of raw agreement, so the same shuffles rank them the same
    way and their p-values must match. Before 2.1.0 they did not: on this
    data Fleiss gave p = 0.02, Conger 0.055, and Gwet, Brennan-Prediger and
    Krippendorff (then one-sided) about 0.99. All five now give 0.0398.
    """

    @staticmethod
    def _p(name: str) -> float:
        import interrater.analyzers as A
        import interrater.analyzers.configs as C
        from interrater.analyzers.configs import HypothesisTestConfig

        cats = UnorderedCategories(
            names=("A", "B", "C"), labels={"A": "A", "B": "B", "C": "C"}
        )
        rng = np.random.default_rng(42)
        ratings = Ratings(
            category_index=rng.integers(0, 3, (30, 3)).astype(np.int8),
            semantics=RatingsSemantics(raters=3, categories=cats),
        )
        ht = HypothesisTestConfig(n_permutations=200)
        config = getattr(C, f"{name}Config")(hypothesis_test_config=ht)
        result = getattr(A, f"{name}Analyzer")(ratings, config).hypothesis_test()
        assert result.alternative == "two-sided"
        return result.p_value

    @pytest.mark.parametrize("name, key", [
        ("FleissKappa", "fleiss_kappa"), ("GwetAC1", "gwet_ac1"),
        ("BrennanPrediger", "brennan_prediger"),
        ("KrippendorffAlpha", "krippendorff_alpha"),
        ("CongerKappa", "conger_kappa"), ("LightKappa", "light_kappa"),
    ])
    def test_result_names_its_coefficient(self, name, key):
        import interrater.analyzers as A
        import interrater.analyzers.configs as C
        from interrater.analyzers.configs import HypothesisTestConfig

        ht = HypothesisTestConfig(n_permutations=20)
        config = getattr(C, f"{name}Config")(hypothesis_test_config=ht)
        result = getattr(A, f"{name}Analyzer")(_make_ratings(3), config).hypothesis_test()
        assert result.coefficient == key

    @pytest.mark.parametrize("name, key", [
        ("Scotts", "scotts_pi"), ("CohenKappa", "cohen_kappa"),
    ])
    @pytest.mark.parametrize("method", ["permutation", "parametric_null"])
    def test_two_rater_result_names_its_coefficient(self, name, key, method):
        import interrater.analyzers as A
        import interrater.analyzers.configs as C
        from interrater.analyzers.configs import HypothesisTestConfig

        ht = HypothesisTestConfig(method=method, n_permutations=20)
        config = getattr(C, f"{name}Config")(hypothesis_test_config=ht)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # parametric_null warns; not the subject
            result = getattr(A, f"{name}Analyzer")(_make_ratings(2), config).hypothesis_test()
        assert result.coefficient == key

    def test_same_two_sided_p(self):
        names = ["FleissKappa", "GwetAC1", "BrennanPrediger",
                 "KrippendorffAlpha", "CongerKappa"]
        p = {name: self._p(name) for name in names}
        assert len(set(p.values())) == 1, p


# ------------------------------------------------------------------ #
# Smoke estimates -- all 8 leaf analyzers                              #
# ------------------------------------------------------------------ #


class TestSmokeEstimates:
    """Run estimate() with n_resamples=10 and verify result structure."""

    def test_fleiss_estimate_from_ratings(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        config = FleissKappaConfig(n_resamples=10, random_state=42)
        a = FleissKappaAnalyzer(_make_ratings(3), config)
        result = a.estimate()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None
        assert len(result.bootstrap_distribution) == 10

    def test_fleiss_estimate_from_nij_matrix(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        config = FleissKappaConfig(n_resamples=10, random_state=42)
        a = FleissKappaAnalyzer.build_from_nij_matrix(_make_nij(3), config)
        result = a.estimate_from_nij_matrix()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None
        assert len(result.bootstrap_distribution) == 10

    def test_gwet_estimate_from_ratings(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import GwetAC1Config

        config = GwetAC1Config(n_resamples=10, random_state=42)
        a = GwetAC1Analyzer(_make_ratings(3), config)
        result = a.estimate()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None
        assert len(result.bootstrap_distribution) == 10

    def test_gwet_estimate_from_nij_matrix(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import GwetAC1Config

        config = GwetAC1Config(n_resamples=10, random_state=42)
        a = GwetAC1Analyzer.build_from_nij_matrix(_make_nij(3), config)
        result = a.estimate_from_nij_matrix()

        assert isinstance(result.point, float)
        assert result.bootstrap_distribution is not None

    def test_brennan_prediger_estimate_from_ratings(self):
        from interrater.analyzers.brennan_prediger_analyzer import (
            BrennanPredigerAnalyzer,
        )
        from interrater.analyzers.configs import BrennanPredigerConfig

        config = BrennanPredigerConfig(n_resamples=10, random_state=42)
        a = BrennanPredigerAnalyzer(_make_ratings(3), config)
        result = a.estimate()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None

    def test_brennan_prediger_estimate_from_nij_matrix(self):
        from interrater.analyzers.brennan_prediger_analyzer import (
            BrennanPredigerAnalyzer,
        )
        from interrater.analyzers.configs import BrennanPredigerConfig

        config = BrennanPredigerConfig(n_resamples=10, random_state=42)
        a = BrennanPredigerAnalyzer.build_from_nij_matrix(
            _make_nij(3), config
        )
        result = a.estimate_from_nij_matrix()

        assert isinstance(result.point, float)
        assert result.bootstrap_distribution is not None

    def test_krippendorff_estimate_from_ratings(self):
        from interrater.analyzers.krippendorff_alpha_analyzer import (
            KrippendorffAlphaAnalyzer,
        )
        from interrater.analyzers.configs import KrippendorffAlphaConfig

        config = KrippendorffAlphaConfig(n_resamples=10, random_state=42)
        a = KrippendorffAlphaAnalyzer(_make_ratings(3), config)
        result = a.estimate()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None

    def test_krippendorff_estimate_from_nij_matrix(self):
        from interrater.analyzers.krippendorff_alpha_analyzer import (
            KrippendorffAlphaAnalyzer,
        )
        from interrater.analyzers.configs import KrippendorffAlphaConfig

        config = KrippendorffAlphaConfig(n_resamples=10, random_state=42)
        a = KrippendorffAlphaAnalyzer.build_from_nij_matrix(
            _make_nij(3), config
        )
        result = a.estimate_from_nij_matrix()

        assert isinstance(result.point, float)
        assert result.bootstrap_distribution is not None

    def test_scott_estimate_from_ratings(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        config = ScottsConfig(n_resamples=10, random_state=42)
        a = ScottsAnalyzer(_make_ratings(2), config)
        result = a.estimate()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None

    def test_scott_estimate_from_nij_matrix(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        config = ScottsConfig(n_resamples=10, random_state=42)
        a = ScottsAnalyzer.build_from_nij_matrix(_make_nij(2), config)
        result = a.estimate_from_nij_matrix()

        assert isinstance(result.point, float)
        assert result.bootstrap_distribution is not None

    def test_cohen_estimate(self):
        from interrater.analyzers.cohen_kappa_analyzer import CohenKappaAnalyzer
        from interrater.analyzers.configs import CohenKappaConfig

        config = CohenKappaConfig(n_resamples=10, random_state=42)
        a = CohenKappaAnalyzer(_make_ratings(2), config)
        result = a.estimate()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None
        assert len(result.bootstrap_distribution) == 10

    def test_conger_estimate(self):
        from interrater.analyzers.conger_kappa_analyzer import CongerKappaAnalyzer
        from interrater.analyzers.configs import CongerKappaConfig

        config = CongerKappaConfig(n_resamples=10, random_state=42)
        a = CongerKappaAnalyzer(_make_ratings(3), config)
        result = a.estimate()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None

    def test_light_estimate(self):
        from interrater.analyzers.light_kappa_analyzer import LightKappaAnalyzer
        from interrater.analyzers.configs import LightKappaConfig

        config = LightKappaConfig(n_resamples=10, random_state=42)
        a = LightKappaAnalyzer(_make_ratings(3), config)
        result = a.estimate()

        assert isinstance(result.point, float)
        _assert_bootstrap_interval(result, n_resamples=10)
        assert result.bootstrap_distribution is not None
        assert hasattr(result, "n_pairs")

    def test_cohen_parametric_null_ci(self):
        """Cohen supports ci_method='parametric_null' -- no bootstrap."""
        from interrater.analyzers.cohen_kappa_analyzer import CohenKappaAnalyzer
        from interrater.analyzers.configs import CohenKappaConfig

        # Two warnings now travel with this path -- the config names the
        # parameter, the interval names how far off its width is. Both are
        # asserted properly in TestParametricNullDisclosure below.
        with pytest.warns(UserWarning):
            config = CohenKappaConfig(ci_method="parametric_null", random_state=42)
            a = CohenKappaAnalyzer(_make_ratings(2), config)
            result = a.estimate()

        assert isinstance(result.point, float)
        _assert_parametric_interval(result)
        assert result.bootstrap_distribution is None

    def test_fleiss_parametric_null_ci(self):
        """Fleiss supports ci_method='parametric_null' -- no bootstrap."""
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        # Two warnings now travel with this path -- the config names the
        # parameter, the interval names how far off its width is. Both are
        # asserted properly in TestParametricNullDisclosure below.
        with pytest.warns(UserWarning):
            config = FleissKappaConfig(ci_method="parametric_null", random_state=42)
            a = FleissKappaAnalyzer(_make_ratings(3), config)
            result = a.estimate()

        assert isinstance(result.point, float)
        _assert_parametric_interval(result)
        assert result.bootstrap_distribution is None

    def test_scott_parametric_null_ci(self):
        """Scott supports ci_method='parametric_null' -- no bootstrap."""
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        # Two warnings now travel with this path -- the config names the
        # parameter, the interval names how far off its width is. Both are
        # asserted properly in TestParametricNullDisclosure below.
        with pytest.warns(UserWarning):
            config = ScottsConfig(ci_method="parametric_null", random_state=42)
            a = ScottsAnalyzer(_make_ratings(2), config)
            result = a.estimate()

        assert isinstance(result.point, float)
        _assert_parametric_interval(result)
        assert result.bootstrap_distribution is None


# ------------------------------------------------------------------ #
# Balance profile propagation                                         #
# ------------------------------------------------------------------ #


class TestBalanceProfilePropagation:
    """Result carries the same balance_profile as the analyzer."""

    def test_fleiss_balance_profile_in_result(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        config = FleissKappaConfig(n_resamples=10, random_state=42)
        a = FleissKappaAnalyzer(_make_ratings(3), config)
        result = a.estimate()

        assert result.balance_profile is a.balance_profile

    def test_cohen_balance_profile_in_result(self):
        from interrater.analyzers.cohen_kappa_analyzer import CohenKappaAnalyzer
        from interrater.analyzers.configs import CohenKappaConfig

        config = CohenKappaConfig(n_resamples=10, random_state=42)
        a = CohenKappaAnalyzer(_make_ratings(2), config)
        result = a.estimate()

        assert result.balance_profile is a.balance_profile

    def test_gwet_balance_profile_in_result(self):
        from interrater.analyzers.gwet_ac1_analyzer import GwetAC1Analyzer
        from interrater.analyzers.configs import GwetAC1Config

        config = GwetAC1Config(n_resamples=10, random_state=42)
        a = GwetAC1Analyzer(_make_ratings(3), config)
        result = a.estimate()

        assert result.balance_profile is a.balance_profile

    def test_conger_balance_profile_in_result(self):
        from interrater.analyzers.conger_kappa_analyzer import CongerKappaAnalyzer
        from interrater.analyzers.configs import CongerKappaConfig

        config = CongerKappaConfig(n_resamples=10, random_state=42)
        a = CongerKappaAnalyzer(_make_ratings(3), config)
        result = a.estimate()

        assert result.balance_profile is a.balance_profile

    def test_scott_balance_profile_in_result(self):
        from interrater.analyzers.scotts_analyzer import ScottsAnalyzer
        from interrater.analyzers.configs import ScottsConfig

        config = ScottsConfig(n_resamples=10, random_state=42)
        a = ScottsAnalyzer(_make_ratings(2), config)
        result = a.estimate()

        assert result.balance_profile is a.balance_profile

    def test_nij_matrix_balance_profile_in_result(self):
        """Balance profile propagated even via NijMatrix path."""
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        config = FleissKappaConfig(n_resamples=10, random_state=42)
        a = FleissKappaAnalyzer.build_from_nij_matrix(_make_nij(3), config)
        result = a.estimate_from_nij_matrix()

        assert result.balance_profile is a.balance_profile


# ------------------------------------------------------------------ #
# Bootstrap summary with failed draws                                  #
# ------------------------------------------------------------------ #


class TestBootstrapSummaryWithFailedDraws:
    """A draw whose coefficient is undefined must not poison the summary.

    Six subjects, one of them with a disagreement: a resample that misses
    it has every rating in one category, so its coefficient is NaN. That
    happens in (5/6)^6 = 33% of draws.
    """

    def _result(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs.fleiss_kappa_config import FleissKappaConfig

        cats = UnorderedCategories(names=("A", "B"), labels={"A": "A", "B": "B"})
        sem = RatingsSemantics(raters=3, categories=cats)
        index = np.zeros((6, 3), dtype=np.int8)
        index[0] = [1, 1, 0]
        ratings = Ratings(category_index=index, semantics=sem)
        config = FleissKappaConfig(n_resamples=1000, random_state=42)
        return FleissKappaAnalyzer(ratings, config).estimate()

    def test_summary_uses_the_usable_draws(self):
        result = self._result()
        draws = result.bootstrap_distribution
        assert result.confidence_interval.n_failed == 329
        d = result.summary()
        assert d["bootstrap_draws_mean"] == pytest.approx(np.nanmean(draws))
        assert d["bootstrap_se"] == pytest.approx(np.nanstd(draws, ddof=1))
        assert round(d["bootstrap_draws_mean"], 3) == 0.394
        assert round(d["bootstrap_se"], 3) == 0.062


# ------------------------------------------------------------------ #
# CI method stored in result                                           #
# ------------------------------------------------------------------ #


class TestCiMethodInResult:
    """Result carries the ci_method from config."""

    def test_percentile_method_in_result(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        config = FleissKappaConfig(
            n_resamples=10, random_state=42, ci_method="percentile"
        )
        a = FleissKappaAnalyzer(_make_ratings(3), config)
        result = a.estimate()
        assert result.ci_method == "percentile"

    def test_bca_method_in_result(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        config = FleissKappaConfig(
            n_resamples=10, random_state=42, ci_method="bca"
        )
        a = FleissKappaAnalyzer(_make_ratings(3), config)
        result = a.estimate()
        assert result.ci_method == "bca"

    def test_parametric_null_method_in_result(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        with pytest.warns(UserWarning):
            config = FleissKappaConfig(
                ci_method="parametric_null", random_state=42
            )
            a = FleissKappaAnalyzer(_make_ratings(3), config)
            result = a.estimate()
        assert result.ci_method == "parametric_null"


# ------------------------------------------------------------------ #
# parametric_null: the two disclosure channels                        #
# ------------------------------------------------------------------ #


class TestParametricNullDisclosure:
    """
    Choosing parametric_null is disclosed twice, saying different things.

    The method is kept because it reproduces published analyses, and it is
    wrong for reporting uncertainty: the null variance is the spread of the
    coefficient at zero, not at the value observed. Keeping a thing that is
    easy to misuse means the disclosure is load-bearing, so it is tested
    rather than assumed.

    The two warnings are deliberately different, and must stay different:

    - the **config** warning names the parameter and its alternatives,
      firing at the moment the caller still knows where the knob is
    - the **interval** warning names how far off the width is, which only
      becomes a useful thing to say once a number exists

    If someone later consolidates them into one message, the test below
    that asserts they differ is what should stop them.
    """

    def _estimate_with_parametric_ci(self):
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        config = FleissKappaConfig(ci_method="parametric_null", random_state=42)
        return FleissKappaAnalyzer(_make_ratings(3), config).estimate()

    def test_both_channels_fire(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self._estimate_with_parametric_ci()

        sources = [str(w.message).split(":")[0] for w in caught]
        assert "[AnalyzerConfig] Warning" in sources
        assert "[ParametricCIResult] Warning" in sources

    def test_the_two_messages_say_different_things(self):
        """One points at the knob, the other at the magnitude.

        Two identical warnings for one call is how people learn to filter
        warnings, which would cost us both.
        """
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            self._estimate_with_parametric_ci()

        by_source = {str(w.message).split(":")[0]: str(w.message) for w in caught}
        config_msg = by_source["[AnalyzerConfig] Warning"]
        interval_msg = by_source["[ParametricCIResult] Warning"]

        assert config_msg != interval_msg
        # The config names the alternatives; the interval does not.
        assert "percentile" in config_msg and "bca" in config_msg
        # The interval names the magnitude; the config does not.
        assert "%" in interval_msg
        assert "%" not in config_msg

    def test_interval_says_so_when_printed(self):
        """Channel three: visible to anyone who prints the result.

        A warning reaches whoever runs the analysis. The string reaches
        whoever reads the output afterwards, which is often not the same
        person -- and the number is what gets copied into a table.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = self._estimate_with_parametric_ci()

        rendered = str(result)
        assert "null variance" in rendered
        assert "kappa = 0" in rendered
        assert "reproduction only" in rendered

    def test_direct_construction_warns_too(self):
        """The reproduction recipe is still disclosed.

        A user rebuilding a published interval from variance_functions
        never touches a config, so the config warning never fires for
        them. The interval still tells them what they are holding.
        """
        from interrater.intermediates.results import ParametricCIResult

        with pytest.warns(UserWarning, match="reproduction of the published"):
            ci = ParametricCIResult(
                lower=0.383, upper=0.478, alpha=0.05, variance_source="null",
            )
        assert ci.variance_source == "null"

    def test_invalid_interval_raises_before_warning(self):
        """Validation first, disclosure second.

        An interval that cannot exist should not first be explained.
        """
        from interrater.intermediates.results import ParametricCIResult

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with pytest.raises(ValueError):
                ParametricCIResult(
                    lower=0.9, upper=0.1, alpha=0.05, variance_source="null",
                )

        assert len(caught) == 0

    def test_bootstrap_path_is_silent(self):
        """The warnings belong to the method, not to running an analysis."""
        from interrater.analyzers.fleiss_kappa_analyzer import FleissKappaAnalyzer
        from interrater.analyzers.configs import FleissKappaConfig

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            config = FleissKappaConfig(
                n_resamples=250, random_state=42, ci_method="percentile",
            )
            FleissKappaAnalyzer(_make_ratings(3), config).estimate()

        assert not [w for w in caught if "parametric" in str(w.message).lower()]
