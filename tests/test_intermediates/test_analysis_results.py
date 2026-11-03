"""
Tests for the results sub-package.

One test file for one sub-package. Replaces the four previous test files
(test_point_estimate_result.py, test_hypothesis_test_result.py,
test_permutation_results.py, test_coefficient_results.py) and adds
coverage for SubsetPermutationResult and MultiKSubsetResult.

Tests are grouped by class in hierarchy order:

  PointEstimateResult
  HypothesisTestResult (standalone -- not an AnalysisResult)
  PermutationResult
  PairwisePermutationResult (with accessor method)
  SubsetPermutationResult (with accessor method)
  MultiKSubsetResult (with container protocol)
  Eight coefficient result classes

Each class section: Construction → Validation → __str__ → summary().
Hand-verifiable values throughout. No randomness.
"""

import dataclasses
import itertools
import math
import warnings

import numpy as np
import pytest

from interrater.intermediates.results import (
    AnalysisResult,
    PointEstimateResult,
    HypothesisTestResult,
    FleissKappaResult,
    ScottsPiResult,
    GwetsAC1Result,
    BrennanPredigerResult,
    KrippendorffAlphaResult,
    CohenKappaResult,
    CongerKappaResult,
)
from interrater.intermediates.results.hypothesis_test_result import (
    permutation_p_value_from_counts,
)
from interrater.intermediates.results.analyzer_results.individual_marginal_result import (
    IndividualMarginalResult,
)
from interrater.intermediates.results.analyzer_results.light_kappa_result import (
    LightKappaResult,
)
from interrater.intermediates.results.permutation_results.permutation_result import (
    PermutationResult,
)
from interrater.intermediates.results.permutation_results.pairwise_permutation_result import (
    PairwisePermutationResult,
)
from interrater.intermediates.results.permutation_results.subset_permutation_result import (
    SubsetPermutationResult,
)
from interrater.intermediates.results.permutation_results.multi_k_subset_result import (
    MultiKSubsetResult,
)
from interrater.intermediates.results.ci_results import (
    BootstrapCIResult,
    ParametricCIResult,
)
from interrater.intermediates.estimator_parts.nij_matrix import NijMatrix
from interrater.intermediates.data_profiles.balance_profile import BalanceProfile
from interrater.intermediates.data_profiles.agreement_profile import AgreementProfile
from interrater.intermediates.not_computed import NotComputed
from interrater.intermediates.results import AgreementParts
from interrater.data_objects.ratings import Ratings
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics

from .conftest import BALANCED_NIJ, UNBALANCED_NIJ


# -- Helpers ----------------------------------------------------------------- #


_PROFILE = AgreementProfile.from_nij_matrix(NijMatrix(BALANCED_NIJ.copy()))

# The three fields every point estimate carries that most tests here are not
# about. The parts and kappa_max say they were left out; the profile is a
# small real one. Tests about these fields build their own.
_NOT_UNDER_TEST = dict(
    parts=NotComputed("not under test"),
    kappa_max=NotComputed("not under test"),
    agreement_profile=_PROFILE,
)


def _fields(point: float, pe: float, weighted: bool = False) -> dict:
    """
    The three fields at ``point`` for a coefficient with chance term ``pe``.

    Po is whatever makes ``point`` with ``pe``; kappa_max is NotComputed,
    as for every coefficient whose chance term puts no limit below 1.
    """
    return dict(
        parts=AgreementParts(
            observed=point * (1 - pe) + pe, chance=pe,
            kind="agreement", weighted=weighted,
        ),
        kappa_max=NotComputed("no limit below 1"),
        agreement_profile=_PROFILE,
    )


def _pooled(point: float, weighted: bool = False) -> dict:
    """Fleiss' kappa or Scott's pi: Pe is the sum of squared shares in _PROFILE."""
    pe = sum(p * p for p in _PROFILE.prevalence.values())
    return _fields(point, pe, weighted)


def _gwet(point: float) -> dict:
    """Gwet's AC1: Pe is Sum p_j (1 - p_j) / (K - 1) over the shares in _PROFILE."""
    shares = list(_PROFILE.prevalence.values())
    return _fields(point, sum(p * (1 - p) for p in shares) / (len(shares) - 1))


def _uniform(point: float, weighted: bool = False) -> dict:
    """Brennan-Prediger: Pe is 1/K for the K categories in _PROFILE."""
    return _fields(point, 1 / len(_PROFILE.category_names), weighted)


def _rater_profile(index) -> AgreementProfile:
    """An agreement profile from ratings, so it carries each rater's marginals."""
    index = np.array(index, dtype=np.int8)
    semantics = RatingsSemantics(
        raters=index.shape[1], categories=UnorderedCategories(("a", "b", "c")),
    )
    return AgreementProfile.from_ratings(Ratings(category_index=index, semantics=semantics))


# Two raters with identical marginals (1/4, 1/2, 1/4), so Cohen's ceiling is 1
_PAIR = _rater_profile([[0, 0], [1, 1], [2, 1], [1, 2]])
# Three raters, the BALANCED_NIJ subjects with rater identity
_PANEL = _rater_profile([[0, 0, 1], [1, 1, 1], [2, 1, 2], [0, 0, 0]])

# As _NOT_UNDER_TEST, for the individual-marginal results, which need each
# rater's marginals: three raters, as _bp_with_n_raters(3) says
_PANEL_NOT_UNDER_TEST = dict(_NOT_UNDER_TEST, agreement_profile=_PANEL)


def _marginals(profile: AgreementProfile) -> np.ndarray:
    return np.array([list(r.values()) for r in profile.rater_marginals.values()])


def _cohen(point: float, weighted: bool = False) -> dict:
    """Cohen's kappa for _PAIR: Pe = Sum p_a p_b, and kappa_max = 1."""
    p_a, p_b = _marginals(_PAIR)
    fields = _fields(point, float((p_a * p_b).sum()), weighted)
    return dict(fields, kappa_max=1.0, agreement_profile=_PAIR)


def _conger(point: float) -> dict:
    """Conger's kappa for _PANEL: Conger's Pe on its marginals, no kappa_max."""
    m = _marginals(_PANEL)
    pe = (m.mean(axis=0) ** 2).sum() - m.var(axis=0, ddof=1).sum() / 3
    return dict(_fields(point, float(pe)), agreement_profile=_PANEL)


def _light(point: float) -> dict:
    """Light's kappa for _PANEL: every pair's kappa is ``point``, ceiling 1."""
    pairs = list(itertools.combinations(_PANEL.rater_marginals, 2))
    pair = _fields(point, 0.5)["parts"]
    return dict(
        parts=NotComputed("a mean of pairs"),
        kappa_max=NotComputed("an upper bound only"),
        agreement_profile=_PANEL,
        pair_parts={p: pair for p in pairs},
        pair_kappa_max={p: 1.0 for p in pairs},
    )


def _alpha(point: float, weighted: bool = False) -> dict:
    """Krippendorff's alpha: De = 0.5, and Do whatever makes ``point``."""
    return dict(
        parts=AgreementParts(
            observed=(1 - point) * 0.5, chance=0.5,
            kind="disagreement", weighted=weighted,
        ),
        kappa_max=NotComputed("no limit below 1"),
        agreement_profile=_PROFILE,
    )


def _ci(lower: float = 0.45, upper: float = 0.55, **overrides) -> BootstrapCIResult:
    """
    A plain percentile interval, for tests that need one but are not about it.

    Intervals validate themselves, so the results that hold them no longer
    re-check bounds, ordering or types. Tests that used to pass a tuple of two
    floats pass one of these instead; tests that are *about* the interval build
    their own with the fields they care about.
    """
    base = dict(
        lower=lower,
        upper=upper,
        alpha=0.05,
        n_resamples=1000,
        n_failed=0,
        method_intended="percentile",
        method_used="percentile",
    )
    base.update(overrides)
    return BootstrapCIResult(**base)


def _undefined_ci(**overrides) -> BootstrapCIResult:
    """An interval that was attempted and produced nothing: every draw failed."""
    base = dict(
        lower=float("nan"),
        upper=float("nan"),
        alpha=0.05,
        n_resamples=200,
        n_failed=200,
        method_intended="bca",
        method_used=None,
    )
    base.update(overrides)
    return BootstrapCIResult(**base)


def _draws_for(ci: BootstrapCIResult) -> np.ndarray:
    """Draws that agree with ``ci``: n_resamples of them, n_failed NaN."""
    draws = np.full(ci.n_resamples, np.nan)
    if ci.n_usable:
        low, high = (ci.lower, ci.upper) if ci.is_defined else (0.0, 0.0)
        draws[ci.n_failed:] = np.linspace(low, high, ci.n_usable)
    return draws


def _pe(point, ci=None, balance_profile=None) -> PointEstimateResult:
    """A PointEstimateResult whose draws, if any, agree with its interval."""
    draws = _draws_for(ci) if isinstance(ci, BootstrapCIResult) else None
    return PointEstimateResult(point, ci, draws, balance_profile, **_NOT_UNDER_TEST)


def _bp() -> BalanceProfile:
    """Balanced profile: 4 subjects, 3 categories, 3 raters each."""
    return BalanceProfile.from_nij_matrix(NijMatrix(BALANCED_NIJ.copy()))


def _unbalanced_bp() -> BalanceProfile:
    """Unbalanced profile: 4 subjects, 2 categories, variable raters."""
    return BalanceProfile.from_nij_matrix(NijMatrix(UNBALANCED_NIJ.copy()))


def _complete_bp(n_raters: int = 3) -> BalanceProfile:
    """As _bp_with_n_raters, with every rater rating every subject (Conger's kappa)."""
    return dataclasses.replace(_bp_with_n_raters(n_raters), fill_rate=1.0)


def _bp_with_n_raters(n_raters: int = 3) -> BalanceProfile:
    """Balanced profile with n_raters set (needed for Light's kappa).

    BalanceProfile.from_nij_matrix sets n_raters=None because rater identity
    is lost. Light's kappa requires rater identity, so its BalanceProfile
    always comes from Ratings (which sets n_raters). This helper simulates
    that by constructing a BalanceProfile directly with n_raters set.
    """
    base = _bp()
    return BalanceProfile(
        is_balanced=base.is_balanced,
        n_raters=n_raters,
        n_raters_exact=True,
        n_per_subject=base.n_per_subject,
        n_per_rater=None,
        subject_counts=base.subject_counts,
        fill_rate=None,
        raters_per_subject=base.raters_per_subject,
    )


# -- Per-pair and per-subset bootstrap, as every analysis produces it ------- #
# The per-unit bootstrap always runs (n_resamples >= 1), so every result has
# an interval per unit, and its draws whenever the units fit under the
# storage limit.


def _cis_for(values) -> dict:
    """One interval per unit, keyed like ``values``."""
    if not isinstance(values, dict):
        return {}
    return {key: _ci(0.40, 0.60) for key in values}


def _dists_for(values) -> dict:
    """Draws per unit that agree with ``_cis_for(values)``."""
    return _dists_matching(_cis_for(values))


def _dists_matching(cis) -> dict:
    """Draws per unit that agree with each unit's interval in ``cis``.

    Anything that is not an interval gets a placeholder; the result refuses
    the interval first, so the placeholder is never compared.
    """
    if not isinstance(cis, dict):
        return {}
    return {
        key: _draws_for(ci) if isinstance(ci, BootstrapCIResult)
        else np.array([0.5])
        for key, ci in cis.items()
    }


def _many_units(n_raters: int, k: int) -> dict:
    """Every k-subset of n_raters, each valued 0.5."""
    return {s: 0.5 for s in itertools.combinations(range(n_raters), k)}


# -- Pairwise test data: 3 raters → C(3,2) = 3 pairs ----------------------- #

PAIRS = {(0, 1): 0.6, (0, 2): 0.4, (1, 2): 0.5}
SUBJECTS_PER_PAIR = {(0, 1): 100, (0, 2): 80, (1, 2): 90}
PAIRWISE_MEAN = 0.5  # (0.6 + 0.4 + 0.5) / 3


def _pairwise_kwargs(**overrides):
    """Defaults for PairwisePermutationResult without a bootstrap of the mean."""
    base = dict(
        coefficient_name="cohen_kappa",
        ci_method="percentile",
        balance_profile=_bp(),
        pairwise_values=dict(PAIRS),
        subjects_per_pair=dict(SUBJECTS_PER_PAIR),
        n_pairs=3,
        mean_over_pairs=PAIRWISE_MEAN,
        bootstrap_mean_method=None,
        bootstrap_ci_of_mean=None,
        bootstrap_distribution_of_mean=None,
    )
    base.update(overrides)
    base.setdefault("bootstrap_cis_per_pair", _cis_for(base["pairwise_values"]))
    base.setdefault(
        "bootstrap_distributions_per_pair",
        _dists_matching(base["bootstrap_cis_per_pair"]),
    )
    return base


def _pairwise_with_bootstrap(**overrides):
    """Defaults for PairwisePermutationResult with correlated bootstrap."""
    # Overrides go in before the per-unit draws are built, so draws
    # follow any per-unit intervals a test supplies.
    return _pairwise_kwargs(**{
        "bootstrap_mean_method": "correlated",
        "bootstrap_ci_of_mean": _ci(0.45, 0.55, n_resamples=3),
        "bootstrap_distribution_of_mean": np.array([0.48, 0.50, 0.52]),
        **overrides,
    })


# -- Subset test data: k=2 (same structure as pairs) and k=3 --------------- #

SUBSETS_K2 = {(0, 1): 0.6, (0, 2): 0.4, (1, 2): 0.5}
SUBJECTS_K2 = {(0, 1): 100, (0, 2): 80, (1, 2): 90}
SUBSET_K2_MEAN = 0.5

SUBSETS_K3 = {(0, 1, 2): 0.55}
SUBJECTS_K3 = {(0, 1, 2): 70}
SUBSET_K3_MEAN = 0.55


def _subset_kwargs(**overrides):
    """Defaults for SubsetPermutationResult (k=2) without a bootstrap of the mean."""
    base = dict(
        coefficient_name="fleiss_kappa",
        ci_method="percentile",
        balance_profile=_bp(),
        subset_values=dict(SUBSETS_K2),
        subjects_per_subset=dict(SUBJECTS_K2),
        n_subsets=3,
        subset_k=2,
        mean_over_subsets=SUBSET_K2_MEAN,
        bootstrap_mean_method=None,
        bootstrap_ci_of_mean=None,
        bootstrap_distribution_of_mean=None,
    )
    base.update(overrides)
    base.setdefault("bootstrap_cis_per_subset", _cis_for(base["subset_values"]))
    base.setdefault(
        "bootstrap_distributions_per_subset",
        _dists_matching(base["bootstrap_cis_per_subset"]),
    )
    return base


def _subset_with_bootstrap(**overrides):
    """Defaults for SubsetPermutationResult (k=2) with correlated bootstrap."""
    # Overrides go in before the per-unit draws are built, so draws
    # follow any per-unit intervals a test supplies.
    return _subset_kwargs(**{
        "bootstrap_mean_method": "correlated",
        "bootstrap_ci_of_mean": _ci(0.45, 0.55, n_resamples=3),
        "bootstrap_distribution_of_mean": np.array([0.48, 0.50, 0.52]),
        **overrides,
    })


def _multi_k_kwargs(**overrides):
    """Defaults for MultiKSubsetResult with ks=[2, 3].

    Shares a single balance_profile object across all per-k results
    and the parent (identity check in validation).
    """
    bp = _bp()
    result_k2 = SubsetPermutationResult(
        coefficient_name="fleiss_kappa",
        ci_method="percentile",
        balance_profile=bp,
        subset_values=dict(SUBSETS_K2),
        subjects_per_subset=dict(SUBJECTS_K2),
        n_subsets=3,
        subset_k=2,
        mean_over_subsets=SUBSET_K2_MEAN,
        bootstrap_mean_method=None,
        bootstrap_ci_of_mean=None,
        bootstrap_distribution_of_mean=None,
        bootstrap_cis_per_subset=_cis_for(SUBSETS_K2),
        bootstrap_distributions_per_subset=_dists_for(SUBSETS_K2),
    )
    result_k3 = SubsetPermutationResult(
        coefficient_name="fleiss_kappa",
        ci_method="percentile",
        balance_profile=bp,
        subset_values=dict(SUBSETS_K3),
        subjects_per_subset=dict(SUBJECTS_K3),
        n_subsets=1,
        subset_k=3,
        mean_over_subsets=SUBSET_K3_MEAN,
        bootstrap_mean_method=None,
        bootstrap_ci_of_mean=None,
        bootstrap_distribution_of_mean=None,
        bootstrap_cis_per_subset=_cis_for(SUBSETS_K3),
        bootstrap_distributions_per_subset=_dists_for(SUBSETS_K3),
    )
    base = dict(
        coefficient_name="fleiss_kappa",
        ci_method="percentile",
        balance_profile=bp,
        per_k_results={2: result_k2, 3: result_k3},
        ks=[2, 3],
    )
    base.update(overrides)
    return base


# ============================================================================ #
#                          PointEstimateResult                                 #
# ============================================================================ #


class TestPointEstimateResultConstruction:
    """Valid construction with various combinations of optional fields."""

    def test_point_only(self):
        """Minimal: just a point estimate, everything else None."""
        r = PointEstimateResult(
            point=0.75,
            confidence_interval=None,
            bootstrap_distribution=None,
            balance_profile=None,
            **_NOT_UNDER_TEST,
        )
        assert r.point == 0.75
        assert r.confidence_interval is None
        assert r.bootstrap_distribution is None
        assert r.balance_profile is None

    def test_with_ci(self):
        """Point estimate with a bootstrap interval and its draws."""
        ci = _ci(0.3, 0.7)
        r = PointEstimateResult(
            point=0.5,
            confidence_interval=ci,
            bootstrap_distribution=_draws_for(ci),
            balance_profile=None,
            **_NOT_UNDER_TEST,
        )
        assert r.confidence_interval is ci
        assert r.confidence_interval.as_tuple() == (0.3, 0.7)

    def test_with_bootstrap(self):
        """Point estimate with bootstrap distribution and its interval."""
        boot = np.array([0.4, 0.5, 0.6])
        r = PointEstimateResult(
            point=0.5,
            confidence_interval=_ci(0.4, 0.6, n_resamples=3),
            bootstrap_distribution=boot,
            balance_profile=None,
            **_NOT_UNDER_TEST,
        )
        np.testing.assert_array_equal(r.bootstrap_distribution, boot)

    def test_with_balance_profile(self):
        """Point estimate with balance profile."""
        bp = _bp()
        r = PointEstimateResult(
            point=0.5,
            confidence_interval=None,
            bootstrap_distribution=None,
            balance_profile=bp,
            **_NOT_UNDER_TEST,
        )
        assert r.balance_profile is bp

    def test_all_fields(self):
        """All optional fields populated."""
        boot = np.array([0.3, 0.4, 0.5, 0.6, 0.7])
        bp = _bp()
        r = PointEstimateResult(
            point=0.5,
            confidence_interval=_ci(0.3, 0.7, n_resamples=5),
            bootstrap_distribution=boot,
            balance_profile=bp,
            **_NOT_UNDER_TEST,
        )
        assert r.point == 0.5
        assert r.confidence_interval.as_tuple() == (0.3, 0.7)
        assert len(r.bootstrap_distribution) == 5
        assert r.balance_profile is bp

    def test_with_parametric_ci(self):
        """Any CIResult is accepted, not only a bootstrap one."""
        # The null-variance warning is the interval's own business, and it is
        # pinned in test_analyzer_classes.py. Here the interval is scaffolding,
        # so the filter is suppressed for these three lines only --
        # catch_warnings restores the previous filters on exit.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            ci = ParametricCIResult(
                lower=0.401, upper=0.639, alpha=0.05, variance_source="null",
            )
        r = PointEstimateResult(0.52, ci, None, None, **_NOT_UNDER_TEST)
        assert r.confidence_interval is ci

    def test_with_undefined_ci(self):
        """An interval that produced no bounds is still an interval."""
        r = _pe(0.5, _undefined_ci())
        assert r.confidence_interval.is_defined is False
        assert r.point_outside_interval is False

    def test_is_analysis_result(self):
        """PointEstimateResult is a subclass of AnalysisResult."""
        r = PointEstimateResult(0.5, None, None, None, **_NOT_UNDER_TEST)
        assert isinstance(r, AnalysisResult)

    def test_is_frozen(self):
        """Frozen dataclass: cannot reassign attributes."""
        r = PointEstimateResult(0.5, None, None, None, **_NOT_UNDER_TEST)
        with pytest.raises(AttributeError):
            r.point = 0.9

    def test_int_point_accepted(self):
        """Integer point values are accepted (int is numeric)."""
        r = PointEstimateResult(1, None, None, None, **_NOT_UNDER_TEST)
        assert r.point == 1


class TestPointEstimateResultValidationRejects:
    """Each constraint in __post_init__ tested with one failing input."""

    def test_point_not_numeric(self):
        """Point must be int or float, not str."""
        with pytest.raises(TypeError, match="point must be float"):
            PointEstimateResult("0.5", None, None, None, **_NOT_UNDER_TEST)

    def test_ci_tuple_rejected(self):
        """A bare pair of floats is no longer an interval."""
        with pytest.raises(TypeError, match="confidence_interval must be"):
            PointEstimateResult(0.5, (0.3, 0.7), None, None, **_NOT_UNDER_TEST)

    def test_ci_list_rejected(self):
        with pytest.raises(TypeError, match="confidence_interval must be"):
            PointEstimateResult(0.5, [0.3, 0.7], None, None, **_NOT_UNDER_TEST)

    # Bounds ordering, bounds types and interval length are no longer checked
    # here. An interval validates those when it is built, so a result that
    # re-checked them would be asserting something already guaranteed, in a
    # second place that could drift from the first. See the CIResult tests.

    def test_bootstrap_not_ndarray(self):
        """Bootstrap distribution must be np.ndarray."""
        with pytest.raises(TypeError, match="bootstrap_distribution must be"):
            PointEstimateResult(0.5, None, [0.4, 0.5, 0.6], None, **_NOT_UNDER_TEST)

    def test_bootstrap_wrong_ndim(self):
        """Bootstrap distribution must be 1-dimensional."""
        with pytest.raises(ValueError, match="must be 1-dimensional"):
            PointEstimateResult(0.5, None, np.array([[0.4, 0.5]]), None, **_NOT_UNDER_TEST)

    def test_bootstrap_empty(self):
        """Bootstrap distribution must not be empty."""
        with pytest.raises(ValueError, match="must have at least 1 value"):
            PointEstimateResult(0.5, None, np.array([]), None, **_NOT_UNDER_TEST)

    def test_balance_profile_wrong_type(self):
        """Balance profile must be BalanceProfile, not dict."""
        with pytest.raises(TypeError, match="balance_profile must be"):
            PointEstimateResult(0.5, None, None, {"n_raters": 3}, **_NOT_UNDER_TEST)

    def test_error_prefix(self):
        """Error messages start with [PointEstimateResult] Error."""
        with pytest.raises(TypeError, match=r"\[PointEstimateResult\] Error"):
            PointEstimateResult("bad", None, None, None, **_NOT_UNDER_TEST)


class TestPointEstimateResultPointOutsideInterval:
    """
    A point estimate outside its own interval is allowed, and is never silent.

    It used to raise. Strong bias correction puts a BCa interval entirely to
    one side of the point legitimately, so refusing to exist was the wrong
    response -- but it is also the signature of a degenerate bootstrap, and of
    a bootstrap that resampled a different pool than the point estimate was
    computed from. So it is permitted and disclosed three ways.
    """

    def test_point_below_interval_is_allowed(self):
        with pytest.warns(UserWarning, match="falls outside its own interval"):
            r = _pe(0.2, _ci(0.3, 0.7))
        assert r.point == 0.2
        assert r.point_outside_interval is True

    def test_point_above_interval_is_allowed(self):
        with pytest.warns(UserWarning, match="falls outside its own interval"):
            r = _pe(0.8, _ci(0.3, 0.7))
        assert r.point_outside_interval is True

    def test_point_inside_does_not_warn(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            r = _pe(0.5, _ci(0.3, 0.7))
        assert r.point_outside_interval is False

    def test_boundary_counts_as_inside(self):
        """A point exactly on a bound is in, not out."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            assert (
                _pe(0.3, _ci(0.3, 0.7))
                .point_outside_interval
                is False
            )
            assert (
                _pe(0.7, _ci(0.3, 0.7))
                .point_outside_interval
                is False
            )

    def test_undefined_interval_is_not_outside(self):
        """No bounds means no comparison to make -- not a violation."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            r = _pe(0.5, _undefined_ci())
        assert r.point_outside_interval is False

    def test_no_interval_is_not_outside(self):
        assert (
            PointEstimateResult(0.5, None, None, None, **_NOT_UNDER_TEST).point_outside_interval
            is False
        )

    def test_marker_appears_in_str(self):
        """Channel two: visible to anyone who prints the result."""
        with pytest.warns(UserWarning):
            r = _pe(0.2, _ci(0.3, 0.7))
        assert "POINT OUTSIDE INTERVAL" in str(r)

    def test_flag_appears_in_summary(self):
        """Channel three: visible to anyone who serializes it."""
        with pytest.warns(UserWarning):
            r = _pe(0.2, _ci(0.3, 0.7))
        assert r.summary()["point_outside_interval"] is True

    def test_warning_names_the_two_explanations(self):
        """The warning has to be actionable, not just alarming."""
        with pytest.warns(UserWarning) as caught:
            _pe(0.2, _ci(0.3, 0.7))
        message = str(caught[0].message)
        assert "bias correction" in message
        assert "degenerate" in message


class TestPointEstimateResultStr:
    """The compact display string under different field combinations."""

    def test_point_only(self):
        """Just the point value, no CI or bootstrap."""
        s = str(PointEstimateResult(0.523, None, None, None, **_NOT_UNDER_TEST))
        assert s == "= 0.523"

    def test_with_ci(self):
        """Point + CI. The interval renders itself."""
        s = str(_pe(0.523, _ci(0.481, 0.565)))
        assert s == "= 0.523 [0.481, 0.565] (95% CI, percentile, n=1000)"

    def test_ci_qualifiers_travel_with_the_interval(self):
        """A fallback and a failure count appear without this class's help."""
        fallen = _ci(
            0.481, 0.565,
            method_intended="bca", method_used="percentile", n_failed=12,
        )
        s = str(_pe(0.523, fallen))
        assert "percentile [bca requested]" in s
        assert "12 failed" in s

    def test_with_ci_and_bootstrap(self):
        """The bootstrap count is dropped once an interval says n= itself."""
        boot = np.ones(1000)
        s = str(PointEstimateResult(1.0, _ci(1.0, 1.0), boot, None, **_NOT_UNDER_TEST))
        assert "n=1000" in s
        assert "n_bootstrap" not in s

    def test_undefined_ci_says_so(self):
        s = str(_pe(0.5, _undefined_ci()))
        assert "[undefined]" in s
        assert "no method ran" in s


class TestPointEstimateResultSummary:
    """The plain-types dict under different field combinations."""

    def test_point_only(self):
        """Minimal summary: the point and the three required fields."""
        d = PointEstimateResult(0.75, None, None, None, **_NOT_UNDER_TEST).summary()
        assert set(d) == {"point", "parts", "kappa_max", "agreement_profile"}
        assert d["point"] == 0.75

    def test_with_ci(self):
        """CI adds the interval's own summary plus hoisted bounds."""
        d = _pe(0.5, _ci(0.3, 0.7)).summary()
        assert d["point"] == 0.5
        assert d["ci_lower"] == 0.3
        assert d["ci_upper"] == 0.7
        assert d["point_outside_interval"] is False
        assert d["ci"]["method_used"] == "percentile"
        assert d["ci"]["n_failed"] == 0

    def test_with_bootstrap(self):
        """Draws add their mean and standard error; counts live in 'ci'."""
        boot = np.array([1.0, 2.0, 3.0])
        d = PointEstimateResult(
            2.0, _ci(1.0, 3.0, n_resamples=3), boot, None,
            **_NOT_UNDER_TEST,
        ).summary()
        assert abs(d["bootstrap_draws_mean"] - 2.0) < 1e-12
        assert abs(d["bootstrap_se"] - 1.0) < 1e-12
        assert d["ci"]["n_resamples"] == 3
        assert "bootstrap_n" not in d

    def test_with_balance_profile(self):
        """Balance profile adds nested 'balance' dict."""
        bp = _bp()
        d = PointEstimateResult(0.5, None, None, bp, **_NOT_UNDER_TEST).summary()
        assert "balance" in d
        assert isinstance(d["balance"], dict)

    def test_all_fields(self):
        """Full summary has all keys."""
        boot = np.array([0.4, 0.5, 0.6])
        bp = _bp()
        d = PointEstimateResult(
            0.5, _ci(0.4, 0.6, n_resamples=3), boot, bp,
            **_NOT_UNDER_TEST,
        ).summary()
        expected_keys = {
            "point", "ci", "ci_lower", "ci_upper", "point_outside_interval",
            "bootstrap_draws_mean", "bootstrap_se",
            "balance", "parts", "kappa_max", "agreement_profile",
        }
        assert set(d.keys()) == expected_keys


class TestPointEstimateResultParts:
    """What the coefficient is made of, checked against the point estimate."""

    def _result(self, point, **overrides):
        kwargs = dict(_NOT_UNDER_TEST, **overrides)
        return PointEstimateResult(point, None, None, None, **kwargs)

    def test_parts_that_rebuild_the_point_accepted(self):
        parts = AgreementParts(observed=0.8, chance=0.5, kind="agreement", weighted=False)
        assert self._result(0.6, parts=parts).parts is parts

    def test_parts_that_do_not_rebuild_the_point_refused(self):
        parts = AgreementParts(observed=0.8, chance=0.5, kind="agreement", weighted=False)
        with pytest.raises(ValueError, match="make 0.600000, but the point estimate is 0.610000"):
            self._result(0.61, parts=parts)

    def test_disagreement_parts(self):
        """Krippendorff: 1 - Do / De."""
        parts = AgreementParts(observed=0.2, chance=0.5, kind="disagreement", weighted=False)
        assert self._result(0.6, parts=parts).parts is parts

    def test_nan_point_needs_parts_that_make_nan(self):
        undefined = AgreementParts(observed=1.0, chance=1.0, kind="agreement", weighted=False)
        assert self._result(float("nan"), parts=undefined).parts is undefined
        with pytest.raises(ValueError, match="parts"):
            self._result(0.0, parts=undefined)

    def test_parts_none_refused(self):
        """None cannot say why; NotComputed must."""
        with pytest.raises(TypeError, match="parts must be AgreementParts or NotComputed"):
            self._result(0.6, parts=None)

    def test_kappa_max_bounds_the_point(self):
        assert self._result(0.6, kappa_max=0.6).kappa_max == 0.6
        with pytest.raises(ValueError, match="exceeds kappa_max"):
            self._result(0.61, kappa_max=0.6)

    def test_kappa_max_above_one_refused(self):
        with pytest.raises(ValueError, match=r"kappa_max must be in \[0, 1\]"):
            self._result(0.6, kappa_max=1.2)

    def test_negative_kappa_max_refused(self):
        """Independent ratings are always possible, so the maximum is at least 0."""
        with pytest.raises(ValueError, match=r"kappa_max must be in \[0, 1\]"):
            self._result(-0.2, kappa_max=-0.1)

    def test_kappa_max_none_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be a number or NotComputed"):
            self._result(0.6, kappa_max=None)

    def test_agreement_profile_required(self):
        with pytest.raises(TypeError, match="agreement_profile must be AgreementProfile"):
            self._result(0.6, agreement_profile=None)

    def test_summary_plain_types(self):
        parts = AgreementParts(observed=0.8, chance=0.5, kind="agreement", weighted=False)
        d = self._result(0.6, parts=parts, kappa_max=0.75).summary()
        assert d["parts"] == parts.summary()
        assert d["kappa_max"] == 0.75
        assert d["agreement_profile"]["n_subjects"] == 4

    def test_summary_says_why_not(self):
        d = self._result(0.6).summary()
        assert d["parts"] == {"not_computed": "not under test"}
        assert d["kappa_max"] == {"not_computed": "not under test"}


class TestPointEstimateResultDrawsMatchInterval:
    """Draws and interval are two views of one bootstrap; I refuse a mismatch."""

    def test_draws_without_interval_refused(self):
        with pytest.raises(ValueError, match="not a BootstrapCIResult"):
            PointEstimateResult(0.5, None, np.array([0.4, 0.5, 0.6]), None, **_NOT_UNDER_TEST)

    def test_draws_with_parametric_interval_refused(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            ci = ParametricCIResult(
                lower=0.401, upper=0.639, alpha=0.05, variance_source="null",
            )
        with pytest.raises(ValueError, match="ParametricCIResult"):
            PointEstimateResult(0.52, ci, np.array([0.5, 0.5]), None, **_NOT_UNDER_TEST)

    def test_bootstrap_interval_without_draws_refused(self):
        with pytest.raises(ValueError, match="comes with the draws"):
            PointEstimateResult(0.5, _ci(0.3, 0.7), None, None, **_NOT_UNDER_TEST)

    def test_draw_count_must_match(self):
        with pytest.raises(ValueError, match="has 3 draws but .* resampled 1000"):
            PointEstimateResult(0.5, _ci(0.3, 0.7), np.array([0.4, 0.5, 0.6]), None, **_NOT_UNDER_TEST)

    def test_failed_count_must_match(self):
        draws = np.array([np.nan, 0.4, 0.5, 0.6])
        with pytest.raises(ValueError, match="1 unusable .* counted 0 failed"):
            PointEstimateResult(0.5, _ci(0.3, 0.7, n_resamples=4), draws, None, **_NOT_UNDER_TEST)

    def test_infinite_draw_counts_as_failed(self):
        """Same rule as the interval: usable means finite."""
        draws = np.array([np.inf, 0.4, 0.5, 0.6])
        ci = _ci(0.4, 0.6, n_resamples=4, n_failed=1)
        assert PointEstimateResult(0.5, ci, draws, None, **_NOT_UNDER_TEST).point == 0.5

    def test_every_draw_failed_is_consistent(self):
        ci = _undefined_ci()
        r = PointEstimateResult(0.5, ci, np.full(200, np.nan), None, **_NOT_UNDER_TEST)
        assert r.confidence_interval.is_defined is False

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[PointEstimateResult\] Error"):
            PointEstimateResult(0.5, _ci(0.3, 0.7), None, None, **_NOT_UNDER_TEST)


class TestPointEstimateResultDrawsNotKept:
    """
    A bootstrap whose draws were not kept (keep_draws: false in a permutation
    analysis) holds NotComputed in their place, saying why. None still means
    no bootstrap at all.
    """

    NOT_KEPT = NotComputed("draws not kept: keep_draws false")

    def test_bootstrap_interval_with_draws_not_kept(self):
        r = PointEstimateResult(0.5, _ci(0.3, 0.7), self.NOT_KEPT, None, **_NOT_UNDER_TEST)
        assert r.bootstrap_distribution is self.NOT_KEPT

    def test_none_still_refused_with_a_bootstrap_interval(self):
        with pytest.raises(ValueError, match="or with NotComputed saying why"):
            PointEstimateResult(0.5, _ci(0.3, 0.7), None, None, **_NOT_UNDER_TEST)

    def test_not_kept_needs_a_bootstrap_interval(self):
        with pytest.raises(ValueError, match="only a bootstrap has draws to not keep"):
            PointEstimateResult(0.5, None, self.NOT_KEPT, None, **_NOT_UNDER_TEST)

    def test_summary_says_why(self):
        d = PointEstimateResult(0.5, _ci(0.3, 0.7), self.NOT_KEPT, None, **_NOT_UNDER_TEST).summary()
        assert d["bootstrap_draws_mean"] == {"not_computed": "draws not kept: keep_draws false"}
        assert d["bootstrap_se"] == {"not_computed": "draws not kept: keep_draws false"}

    def test_plot_says_why(self):
        from interrater.visualizers import plot_bootstrap_distribution

        r = PointEstimateResult(0.5, _ci(0.3, 0.7), self.NOT_KEPT, None, **_NOT_UNDER_TEST)
        with pytest.raises(ValueError, match="draws were not kept"):
            plot_bootstrap_distribution(r)


class TestPointEstimateResultDrawMoments:
    """bootstrap_draws_mean and bootstrap_se use the usable draws only."""

    def test_nan_draws_do_not_poison_the_summary(self):
        """One NaN draw used to make both numbers NaN."""
        draws = np.array([np.nan, 1.0, 2.0, 3.0])
        ci = _ci(1.0, 3.0, n_resamples=4, n_failed=1)
        d = PointEstimateResult(2.0, ci, draws, None, **_NOT_UNDER_TEST).summary()
        assert d["bootstrap_draws_mean"] == 2.0
        assert d["bootstrap_se"] == 1.0

    def test_se_divides_by_n_minus_1(self):
        """Efron & Tibshirani (1993) eq. 6.6."""
        draws = np.array([0.0, 1.0])
        d = PointEstimateResult(
            0.5, _ci(0.0, 1.0, n_resamples=2), draws, None,
            **_NOT_UNDER_TEST,
        ).summary()
        assert abs(d["bootstrap_se"] - np.sqrt(0.5)) < 1e-12

    def test_every_draw_failed_gives_none(self):
        d = PointEstimateResult(
            0.5, _undefined_ci(), np.full(200, np.nan), None,
            **_NOT_UNDER_TEST,
        ).summary()
        assert d["bootstrap_draws_mean"] is None
        assert d["bootstrap_se"] is None

    def test_one_usable_draw_has_a_mean_but_no_se(self):
        draws = np.array([np.nan, 0.5])
        ci = _ci(0.5, 0.5, n_resamples=2, n_failed=1)
        d = PointEstimateResult(0.5, ci, draws, None, **_NOT_UNDER_TEST).summary()
        assert d["bootstrap_draws_mean"] == 0.5
        assert d["bootstrap_se"] is None


# ============================================================================ #
#                         HypothesisTestResult                                 #
# ============================================================================ #


class TestHypothesisTestResultConstruction:
    """Valid construction with different method flavours."""

    def test_parametric_null(self):
        """Parametric null: z_statistic and null_variance populated."""
        r = HypothesisTestResult(
            p_value=0.001,
            alpha=0.05,
            reject_null=True,
            method="parametric_null",
            point=0.65,
            z_statistic=3.42,
            null_variance=0.036,
            coefficient="fleiss_kappa",
        )
        assert r.p_value == 0.001
        assert r.reject_null is True
        assert r.z_statistic == 3.42
        assert r.null_variance == 0.036
        assert r.null_distribution is None
        assert r.n_permutations is None

    def test_permutation(self):
        """Permutation: null_distribution and n_permutations populated."""
        null_dist = np.array([0.01, -0.02, 0.03, 0.05, -0.01])
        r = HypothesisTestResult(
            p_value=0.04,
            alpha=0.05,
            reject_null=True,
            method="permutation",
            point=0.5,
            null_distribution=null_dist,
            n_permutations=5000,
            coefficient="fleiss_kappa",
        )
        assert r.method == "permutation"
        assert r.n_permutations == 5000
        np.testing.assert_array_equal(r.null_distribution, null_dist)
        assert r.z_statistic is None

    def test_permutation_minimal(self):
        """Permutation with only required optional fields."""
        r = HypothesisTestResult(
            p_value=0.5,
            alpha=0.05,
            reject_null=False,
            method="permutation",
            point=0.1,
            n_permutations=1000,
            coefficient="fleiss_kappa",
        )
        assert r.null_distribution is None
        assert r.z_statistic is None
        assert r.null_variance is None

    def test_is_frozen(self):
        """Frozen dataclass: cannot reassign attributes."""
        r = HypothesisTestResult(
            0.001, 0.05, True, "parametric_null", 0.65,
            z_statistic=3.42, null_variance=0.036,
            coefficient="fleiss_kappa",
        )
        with pytest.raises(AttributeError):
            r.p_value = 0.01

    def test_boundary_p_equals_alpha_reject(self):
        """At p == alpha, either decision is valid -- reject accepted."""
        r = HypothesisTestResult(
            p_value=0.05, alpha=0.05, reject_null=True,
            method="permutation", point=0.3, n_permutations=100,
            coefficient="fleiss_kappa",
        )
        assert r.reject_null is True

    def test_boundary_p_equals_alpha_no_reject(self):
        """At p == alpha, either decision is valid -- fail to reject accepted."""
        r = HypothesisTestResult(
            p_value=0.05, alpha=0.05, reject_null=False,
            method="permutation", point=0.3, n_permutations=100,
            coefficient="fleiss_kappa",
        )
        assert r.reject_null is False


class TestHypothesisTestResultValidation:
    """__post_init__ validation with TypeError/ValueError separation."""

    # --- p_value ------------------------------------------------------------ #

    def test_p_value_not_numeric(self):
        with pytest.raises(TypeError, match="p_value must be float"):
            HypothesisTestResult(
                p_value="0.5", alpha=0.05, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    def test_p_value_negative(self):
        with pytest.raises(ValueError, match=r"p_value must be in \[0, 1\]"):
            HypothesisTestResult(
                p_value=-0.1, alpha=0.05, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    def test_p_value_above_one(self):
        with pytest.raises(ValueError, match=r"p_value must be in \[0, 1\]"):
            HypothesisTestResult(
                p_value=1.5, alpha=0.05, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    # --- alpha -------------------------------------------------------------- #

    def test_alpha_not_numeric(self):
        with pytest.raises(TypeError, match="alpha must be float"):
            HypothesisTestResult(
                p_value=0.5, alpha="0.05", reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    def test_alpha_zero(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.0, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    def test_alpha_one(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            HypothesisTestResult(
                p_value=0.5, alpha=1.0, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    # --- reject_null -------------------------------------------------------- #

    def test_reject_null_not_bool(self):
        with pytest.raises(TypeError, match="reject_null must be bool"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=0,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    def test_reject_null_false_but_p_below_alpha(self):
        """p < alpha but reject_null=False is inconsistent."""
        with pytest.raises(ValueError, match="reject_null is False but"):
            HypothesisTestResult(
                p_value=0.01, alpha=0.05, reject_null=False,
                method="permutation", point=0.5, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    def test_reject_null_true_but_p_above_alpha(self):
        """p > alpha but reject_null=True is inconsistent."""
        with pytest.raises(ValueError, match="reject_null is True but"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=True,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )

    # --- method ------------------------------------------------------------- #

    def test_method_not_str(self):
        with pytest.raises(TypeError, match="method must be str"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=False,
                method=42, point=0.1,
                coefficient="fleiss_kappa",
            )

    def test_method_invalid(self):
        with pytest.raises(ValueError, match="method must be"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=False,
                method="custom", point=0.1,
                coefficient="fleiss_kappa",
            )

    # --- point -------------------------------------------------------------- #

    def test_point_not_numeric(self):
        with pytest.raises(TypeError, match="point must be float"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=False,
                method="permutation", point="0.1", n_permutations=100,
                coefficient="fleiss_kappa",
            )

    def test_point_nan_refused(self):
        """What a permutation test used to return for an undefined coefficient."""
        with pytest.raises(ValueError, match="undefined"):
            HypothesisTestResult(
                p_value=0.0, alpha=0.05, reject_null=True,
                method="permutation", point=float("nan"), n_permutations=100,
                coefficient="fleiss_kappa",
            )

    # --- Cross-validation: method-specific required fields ------------------ #

    def test_parametric_null_requires_z_statistic(self):
        with pytest.raises(ValueError, match="z_statistic must not be None"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                null_variance=0.036,
                coefficient="fleiss_kappa",
            )

    def test_parametric_null_requires_null_variance(self):
        with pytest.raises(ValueError, match="null_variance must not be None"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                z_statistic=3.42,
                coefficient="fleiss_kappa",
            )

    def test_permutation_requires_n_permutations(self):
        with pytest.raises(ValueError, match="n_permutations must not be None"):
            HypothesisTestResult(
                p_value=0.04, alpha=0.05, reject_null=True,
                method="permutation", point=0.5,
                coefficient="fleiss_kappa",
            )

    # --- Optional field type/value checks ----------------------------------- #

    def test_z_statistic_not_numeric(self):
        with pytest.raises(TypeError, match="z_statistic must be float"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                z_statistic="3.42", null_variance=0.036,
                coefficient="fleiss_kappa",
            )

    def test_null_variance_not_numeric(self):
        with pytest.raises(TypeError, match="null_variance must be float"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                z_statistic=3.42, null_variance="0.036",
                coefficient="fleiss_kappa",
            )

    def test_null_variance_negative(self):
        with pytest.raises(ValueError, match="null_variance must be >= 0"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                z_statistic=3.42, null_variance=-0.01,
                coefficient="fleiss_kappa",
            )

    def test_null_distribution_not_ndarray(self):
        with pytest.raises(TypeError, match="null_distribution must be"):
            HypothesisTestResult(
                p_value=0.04, alpha=0.05, reject_null=True,
                method="permutation", point=0.5,
                null_distribution=[0.1, 0.2], n_permutations=100,
                coefficient="fleiss_kappa",
            )

    def test_n_permutations_not_int(self):
        with pytest.raises(TypeError, match="n_permutations must be int"):
            HypothesisTestResult(
                p_value=0.04, alpha=0.05, reject_null=True,
                method="permutation", point=0.5,
                n_permutations=100.0,
                coefficient="fleiss_kappa",
            )

    def test_n_permutations_zero(self):
        with pytest.raises(ValueError, match="n_permutations must be >= 1"):
            HypothesisTestResult(
                p_value=0.04, alpha=0.05, reject_null=True,
                method="permutation", point=0.5,
                n_permutations=0,
                coefficient="fleiss_kappa",
            )

    # --- Error prefix ------------------------------------------------------- #

    def test_error_prefix(self):
        with pytest.raises(TypeError, match=r"\[HypothesisTestResult\] Error"):
            HypothesisTestResult(
                p_value="bad", alpha=0.05, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
                coefficient="fleiss_kappa",
            )


class TestHypothesisTestResultStr:
    """Two __str__ branches: parametric and permutation."""

    def test_parametric_format(self):
        """Parametric: shows z statistic."""
        r = HypothesisTestResult(
            p_value=0.0006,
            alpha=0.05,
            reject_null=True,
            method="parametric_null",
            point=0.65,
            z_statistic=3.42,
            null_variance=0.036,
            coefficient="fleiss_kappa",
        )
        s = str(r)
        assert "z = 3.42" in s
        assert "reject" in s
        assert "α = 0.05" in s
        assert "H₀: fleiss_kappa = 0" in s

    def test_permutation_format(self):
        """Permutation: shows n_permutations."""
        r = HypothesisTestResult(
            p_value=0.003,
            alpha=0.05,
            reject_null=True,
            method="permutation",
            point=0.5,
            n_permutations=5000,
            coefficient="fleiss_kappa",
        )
        s = str(r)
        assert "5000 permutations" in s
        assert "reject" in s
        assert "H₀: fleiss_kappa = 0" in s

    def test_very_small_p_value(self):
        """p < 0.000001 uses 'p < 0.000001' format."""
        r = HypothesisTestResult(
            1e-10, 0.05, True, "parametric_null", 0.9,
            z_statistic=8.0, null_variance=0.012,
            coefficient="fleiss_kappa",
        )
        s = str(r)
        assert "p < 0.000001" in s

    def test_normal_p_value_format(self):
        """Normal p-values show 6 decimal places."""
        r = HypothesisTestResult(
            0.042, 0.05, True, "permutation", 0.5, n_permutations=1000,
            coefficient="fleiss_kappa",
        )
        s = str(r)
        assert "p = 0.042000" in s


class TestHypothesisTestResultSummary:
    """summary() returns method-specific keys."""

    def test_required_keys_always_present(self):
        """Permutation with only required extras has base keys + n_permutations."""
        d = HypothesisTestResult(
            0.5, 0.05, False, "permutation", 0.1, n_permutations=100,
            coefficient="fleiss_kappa",
        ).summary()
        assert set(d.keys()) == {
            "p_value", "alpha", "reject_null", "method", "point",
            "coefficient", "n_permutations",
        }

    def test_parametric_keys(self):
        """z_statistic and null_variance included when present."""
        d = HypothesisTestResult(
            0.001, 0.05, True, "parametric_null", 0.65,
            z_statistic=3.42, null_variance=0.036,
            coefficient="fleiss_kappa",
        ).summary()
        assert d["z_statistic"] == 3.42
        assert d["null_variance"] == 0.036
        assert "n_permutations" not in d

    def test_permutation_keys(self):
        """n_permutations included when present; null_distribution is not."""
        d = HypothesisTestResult(
            0.04, 0.05, True, "permutation", 0.5,
            null_distribution=np.array([0.0, 0.1]),
            n_permutations=5000,
            coefficient="fleiss_kappa",
        ).summary()
        assert d["n_permutations"] == 5000
        assert "null_distribution" not in d
        assert "z_statistic" not in d


class TestHypothesisTestResultCounts:
    """
    The p-value follows from the counts the result records.

    Hand-worked counts: 9 usable permutations, 2 at or above the observed
    value, 7 at or below. With the observed data counted as one more
    arrangement (Phipson & Smyth 2010): greater = 3/10, less = 8/10,
    two-sided = 2 * 3/10.
    """

    @pytest.mark.parametrize("alternative, expected", [
        ("greater", 0.3), ("less", 0.8), ("two-sided", 0.6),
    ])
    def test_formula(self, alternative, expected):
        assert permutation_p_value_from_counts(2, 7, 9, alternative) == pytest.approx(expected)

    def test_two_sided_capped_at_one(self):
        """Both tails at 0.6; twice the smaller would be 1.2."""
        assert permutation_p_value_from_counts(5, 5, 9, "two-sided") == 1.0

    def test_never_zero(self):
        """No permutation as extreme as the observed: 1/(B+1), not 0."""
        assert permutation_p_value_from_counts(0, 9, 9, "greater") == pytest.approx(0.1)

    def _result(self, **overrides):
        kwargs = dict(
            p_value=0.6, alpha=0.05, reject_null=False, method="permutation",
            point=0.2, coefficient="fleiss_kappa", n_permutations=10,
            alternative="two-sided",
            n_usable=9, n_at_or_above=2, n_at_or_below=7,
        )
        kwargs.update(overrides)
        return HypothesisTestResult(**kwargs)

    def test_consistent_counts_accepted(self):
        assert self._result().n_usable == 9

    def test_p_value_not_following_counts_refused(self):
        with pytest.raises(ValueError, match="does not follow from the counts"):
            self._result(p_value=0.3)

    def test_counts_without_alternative_refused(self):
        with pytest.raises(ValueError, match="alternative must be too"):
            self._result(alternative=None)

    def test_counts_recorded_together(self):
        with pytest.raises(ValueError, match="recorded together"):
            self._result(n_at_or_below=None)

    @pytest.mark.parametrize("overrides, match", [
        (dict(n_usable=11), "cannot exceed n_permutations"),
        (dict(n_usable=0), "n_usable must be >= 1"),
        (dict(n_at_or_above=10), r"n_at_or_above must be in \[0"),
        (dict(n_at_or_below=-1), r"n_at_or_below must be in \[0"),
    ])
    def test_counts_out_of_range_refused(self, overrides, match):
        with pytest.raises(ValueError, match=match):
            self._result(**overrides)

    def test_alternative_invalid_refused(self):
        with pytest.raises(ValueError, match="alternative must be"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=False,
                method="permutation", point=0.1, n_permutations=10,
                alternative="sideways",
                coefficient="fleiss_kappa",
            )

    def test_str_names_alternative_and_usable_draws(self):
        s = str(self._result())
        assert "two-sided p = 0.600000" in s
        assert "9 of 10 permutations usable" in s

    def test_coefficient_required(self):
        """Every coefficient returns this class, so it must say which it is."""
        kwargs = dict(p_value=0.5, alpha=0.05, reject_null=False,
                      method="permutation", point=0.1, n_permutations=10)
        with pytest.raises(TypeError, match="coefficient"):
            HypothesisTestResult(**kwargs)

    def test_coefficient_unknown_refused(self):
        with pytest.raises(ValueError, match="coefficient must be one of"):
            self._result(coefficient="fleiss")

    def test_str_names_coefficient(self):
        assert str(self._result()).startswith("H₀: fleiss_kappa = 0")

    def test_summary_has_counts(self):
        d = self._result().summary()
        assert d["alternative"] == "two-sided"
        assert (d["n_usable"], d["n_at_or_above"], d["n_at_or_below"]) == (9, 2, 7)


# ============================================================================ #
#                           PermutationResult                                  #
# ============================================================================ #


class TestPermutationResult:
    """Base class: coefficient_name, ci_method, balance_profile."""

    def test_valid_construction(self):
        r = PermutationResult("cohen_kappa", "percentile", _bp())
        assert r.coefficient_name == "cohen_kappa"
        assert r.ci_method == "percentile"

    def test_is_analysis_result(self):
        r = PermutationResult("cohen_kappa", "bca", _bp())
        assert isinstance(r, AnalysisResult)

    def test_coefficient_name_not_str(self):
        """Non-str coefficient_name raises TypeError (not ValueError)."""
        with pytest.raises(TypeError, match="coefficient_name must be str"):
            PermutationResult(42, "percentile", _bp())

    def test_coefficient_name_empty(self):
        with pytest.raises(ValueError, match="coefficient_name must be non-empty"):
            PermutationResult("", "percentile", _bp())

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            PermutationResult("cohen_kappa", 42, _bp())

    def test_ci_method_parametric_null(self):
        """parametric_null is not available for permutation results."""
        with pytest.raises(ValueError, match="ci_method must be"):
            PermutationResult("cohen_kappa", "parametric_null", _bp())

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            PermutationResult("cohen_kappa", "invalid", _bp())

    def test_balance_profile_none(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            PermutationResult("cohen_kappa", "percentile", None)

    def test_balance_profile_wrong_type(self):
        with pytest.raises(TypeError, match="balance_profile must be BalanceProfile"):
            PermutationResult("cohen_kappa", "percentile", {"n_raters": 3})

    def test_summary_keys(self):
        d = PermutationResult("fleiss_kappa", "bca", _bp()).summary()
        assert d["coefficient_name"] == "fleiss_kappa"
        assert d["ci_method"] == "bca"
        assert isinstance(d["balance"], dict)

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[PermutationResult\] Error"):
            PermutationResult("", "percentile", _bp())


# ============================================================================ #
#                     PairwisePermutationResult                                #
# ============================================================================ #


class TestPairwisePermutationResultConstruction:
    """Valid construction with and without bootstrap."""

    def test_without_mean_bootstrap(self):
        r = PairwisePermutationResult(**_pairwise_kwargs())
        assert r.n_pairs == 3
        assert r.mean_over_pairs == PAIRWISE_MEAN
        assert r.bootstrap_mean_method is None
        # The per-pair bootstrap always runs
        assert len(r.bootstrap_cis_per_pair) == r.n_pairs

    def test_with_bootstrap(self):
        r = PairwisePermutationResult(**_pairwise_with_bootstrap())
        assert r.bootstrap_mean_method == "correlated"
        assert r.bootstrap_ci_of_mean.as_tuple() == (0.45, 0.55)
        assert len(r.bootstrap_distribution_of_mean) == 3

    def test_independent_method(self):
        r = PairwisePermutationResult(
            **_pairwise_with_bootstrap(bootstrap_mean_method="independent")
        )
        assert r.bootstrap_mean_method == "independent"

    def test_is_permutation_result(self):
        r = PairwisePermutationResult(**_pairwise_kwargs())
        assert isinstance(r, PermutationResult)
        assert isinstance(r, AnalysisResult)


class TestPairwisePermutationResultValidation:
    """Each validation constraint tested with one failing input."""

    def test_pairwise_values_not_dict(self):
        with pytest.raises(TypeError, match="pairwise_values must be dict"):
            PairwisePermutationResult(
                **_pairwise_kwargs(pairwise_values=[(0, 1, 0.6)])
            )

    def test_pairwise_values_bad_keys(self):
        with pytest.raises(ValueError, match="pairwise_values keys must be 2-tuples"):
            PairwisePermutationResult(**_pairwise_kwargs(
                pairwise_values={0: 0.6, 1: 0.5, 2: 0.4},
            ))

    def test_n_pairs_not_int(self):
        """Non-int n_pairs raises TypeError."""
        with pytest.raises(TypeError, match="n_pairs must be int"):
            PairwisePermutationResult(**_pairwise_kwargs(n_pairs=3.0))

    def test_n_pairs_zero(self):
        """n_pairs must be >= 1."""
        with pytest.raises(ValueError, match="n_pairs must be >= 1"):
            PairwisePermutationResult(**_pairwise_kwargs(n_pairs=0))

    def test_n_pairs_mismatch(self):
        with pytest.raises(
            ValueError, match="pairwise_values has .* entries but n_pairs"
        ):
            PairwisePermutationResult(**_pairwise_kwargs(n_pairs=5))

    def test_subjects_per_pair_keys_mismatch(self):
        with pytest.raises(
            ValueError, match="subjects_per_pair must have the same keys"
        ):
            PairwisePermutationResult(**_pairwise_kwargs(
                subjects_per_pair={(0, 1): 100, (0, 2): 80, (2, 3): 90},
            ))

    def test_mean_not_numeric(self):
        with pytest.raises(TypeError, match="mean_over_pairs must be float"):
            PairwisePermutationResult(**_pairwise_kwargs(mean_over_pairs="0.5"))

    def test_bootstrap_mean_method_invalid(self):
        with pytest.raises(ValueError, match="bootstrap_mean_method must be"):
            PairwisePermutationResult(**_pairwise_with_bootstrap(
                bootstrap_mean_method="invalid",
            ))

    def test_bootstrap_fields_without_method(self):
        """Bootstrap data present but bootstrap_mean_method is None."""
        with pytest.raises(ValueError, match="bootstrap_mean_method must be set"):
            PairwisePermutationResult(**_pairwise_kwargs(
                bootstrap_ci_of_mean=(0.4, 0.6),
            ))

    def test_method_without_bootstrap_fields(self):
        """bootstrap_mean_method set but no bootstrap data."""
        with pytest.raises(
            ValueError, match="no mean bootstrap fields are populated"
        ):
            PairwisePermutationResult(**_pairwise_kwargs(
                bootstrap_mean_method="correlated",
            ))

    def test_bootstrap_ci_of_mean_must_be_an_interval(self):
        """Bounds ordering is the interval's business; its type is ours."""
        with pytest.raises(TypeError, match="must be BootstrapCIResult"):
            PairwisePermutationResult(**_pairwise_with_bootstrap(
                bootstrap_ci_of_mean=(0.7, 0.3),
            ))

    def test_per_pair_ci_must_be_an_interval(self):
        with pytest.raises(TypeError, match="bootstrap CI for pair"):
            PairwisePermutationResult(**_pairwise_with_bootstrap(
                bootstrap_cis_per_pair={
                    (0, 1): (0.50, 0.70),
                    (0, 2): _ci(0.30, 0.50),
                    (1, 2): _ci(0.40, 0.60),
                },
            ))

    def test_pairs_may_fail_independently(self):
        """The reason each pair carries its own interval object.

        One pair with few shared subjects can lose draws, or lose all of
        them, while its neighbours return ordinary intervals. A single
        shared summary could not say that.
        """
        r = PairwisePermutationResult(**_pairwise_with_bootstrap(
            bootstrap_cis_per_pair={
                (0, 1): _ci(0.50, 0.70),
                (0, 2): _ci(0.30, 0.50, n_failed=137),
                (1, 2): _undefined_ci(),
            },
        ))
        cis = r.bootstrap_cis_per_pair
        assert cis[(0, 1)].n_failed == 0
        assert cis[(0, 2)].n_failed == 137
        assert cis[(1, 2)].method_used is None
        assert cis[(1, 2)].is_defined is False

    def test_bootstrap_dist_wrong_ndim(self):
        with pytest.raises(ValueError, match="must be 1-dimensional"):
            PairwisePermutationResult(**_pairwise_with_bootstrap(
                bootstrap_distribution_of_mean=np.array([[0.5, 0.6]]),
            ))

    def test_bootstrap_dist_empty(self):
        with pytest.raises(ValueError, match="must have at least 1 value"):
            PairwisePermutationResult(**_pairwise_with_bootstrap(
                bootstrap_distribution_of_mean=np.array([]),
            ))

    def test_per_pair_cis_required(self):
        """The per-pair bootstrap always runs, so None is refused."""
        with pytest.raises(TypeError, match="always runs"):
            PairwisePermutationResult(**_pairwise_kwargs(bootstrap_cis_per_pair=None))

    def test_draws_may_be_left_out(self):
        """None is what keep_draws: false leaves; the intervals stay."""
        r = PairwisePermutationResult(**_pairwise_kwargs(bootstrap_distributions_per_pair=None))
        assert len(r.bootstrap_cis_per_pair) == 3

    def test_draws_kept_past_the_old_limit(self):
        """2.0.0 refused draws above 20 pairs; there is no limit now."""
        values = _many_units(7, 2)
        r = PairwisePermutationResult(**_pairwise_kwargs(
            pairwise_values=values,
            subjects_per_pair={key: 10 for key in values},
            n_pairs=len(values),
        ))
        assert len(r.bootstrap_distributions_per_pair) == 21

    def test_draws_for_every_pair_or_none(self):
        dists = _dists_for(PAIRS)
        del dists[(1, 2)]
        with pytest.raises(ValueError, match="same keys as pairwise_values"):
            PairwisePermutationResult(**_pairwise_kwargs(bootstrap_distributions_per_pair=dists))

    def test_bootstrap_cis_per_pair_keys_mismatch(self):
        with pytest.raises(
            ValueError,
            match="bootstrap_cis_per_pair must have the same keys",
        ):
            PairwisePermutationResult(**_pairwise_with_bootstrap(
                bootstrap_cis_per_pair={
                    (0, 1): (0.5, 0.7),
                    (0, 2): (0.3, 0.5),
                    # Missing (1, 2)
                },
            ))

    def test_error_prefix(self):
        with pytest.raises(
            TypeError, match=r"\[PairwisePermutationResult\] Error"
        ):
            PairwisePermutationResult(
                **_pairwise_kwargs(pairwise_values="not a dict")
            )

    def test_parent_validation_propagated(self):
        """Parent error messages are wrapped with child prefix."""
        with pytest.raises(
            TypeError, match=r"\[PairwisePermutationResult\] Error"
        ):
            PairwisePermutationResult(**_pairwise_kwargs(
                balance_profile={"not": "a profile"},
            ))


class TestPairwisePermutationResultAccessor:
    """get_bootstrap_distributions_per_pair() accessor method."""

    def test_returns_distributions_when_stored(self):
        """When distributions are stored, accessor returns them."""
        dists = _dists_for(PAIRS)
        r = PairwisePermutationResult(**_pairwise_with_bootstrap(
            bootstrap_distributions_per_pair=dists,
        ))
        result = r.get_bootstrap_distributions_per_pair()
        assert result is dists

    def test_raises_when_draws_not_kept(self):
        """The one reason draws can be missing: keep_draws was false.

        It must not claim the bootstrap was not performed -- it always is.
        """
        r = PairwisePermutationResult(**_pairwise_kwargs(bootstrap_distributions_per_pair=None))
        with pytest.raises(RuntimeError) as exc:
            r.get_bootstrap_distributions_per_pair()
        message = str(exc.value)
        assert "keep_draws false" in message
        assert "bootstrap_cis_per_pair" in message
        assert "not performed" not in message


class TestPairwisePermutationResultStr:
    """__str__ output."""

    def test_no_bootstrap(self):
        s = str(PairwisePermutationResult(**_pairwise_kwargs()))
        assert s.startswith("Pairwise Cohen's κ, 3 pairs: mean = 0.500 (no interval for the mean)")
        assert "each pair on the subjects both raters rated" in s

    def test_with_bootstrap(self):
        s = str(PairwisePermutationResult(**_pairwise_with_bootstrap()))
        assert "95% CI" in s
        assert "correlated" in s


class TestPairwisePermutationResultSummary:
    """summary() returns observed + bootstrap keys."""

    def test_no_bootstrap(self):
        d = PairwisePermutationResult(**_pairwise_kwargs()).summary()
        assert d["n_pairs"] == 3
        assert d["mean_over_pairs"] == PAIRWISE_MEAN
        assert d["coefficient_name"] == "cohen_kappa"
        # Tuple keys stringified for JSON
        assert "(0, 1)" in d["pairwise_values"]

    def test_with_bootstrap(self):
        d = PairwisePermutationResult(**_pairwise_with_bootstrap()).summary()
        assert d["bootstrap_mean_method"] == "correlated"
        assert d["mean_ci_lower"] == 0.45
        assert d["mean_ci_upper"] == 0.55
        # Draws of the mean are 0.48, 0.50, 0.52
        assert d["bootstrap_draws_mean"] == pytest.approx(0.50)
        assert d["bootstrap_se"] == pytest.approx(0.02)
        assert "bootstrap_n" not in d
        assert d["pair_distributions_stored"] is True
        assert len(d["bootstrap_cis_per_pair"]) == 3

    def test_pair_distributions_stored_flag(self):
        """Draws not kept: the flag is False."""
        d = PairwisePermutationResult(**_pairwise_kwargs(
            bootstrap_distributions_per_pair=None,
        )).summary()
        assert d["pair_distributions_stored"] is False


# ============================================================================ #
#                      SubsetPermutationResult                                 #
# ============================================================================ #


class TestSubsetPermutationResultConstruction:
    """Valid construction with and without bootstrap."""

    def test_without_mean_bootstrap(self):
        r = SubsetPermutationResult(**_subset_kwargs())
        assert r.subset_k == 2
        assert r.n_subsets == 3
        assert r.mean_over_subsets == SUBSET_K2_MEAN
        assert r.bootstrap_mean_method is None
        # The per-subset bootstrap always runs
        assert len(r.bootstrap_cis_per_subset) == r.n_subsets

    def test_with_bootstrap(self):
        r = SubsetPermutationResult(**_subset_with_bootstrap())
        assert r.bootstrap_mean_method == "correlated"
        assert r.bootstrap_ci_of_mean.as_tuple() == (0.45, 0.55)
        assert len(r.bootstrap_distribution_of_mean) == 3

    def test_independent_method(self):
        r = SubsetPermutationResult(
            **_subset_with_bootstrap(bootstrap_mean_method="independent")
        )
        assert r.bootstrap_mean_method == "independent"

    def test_is_permutation_result(self):
        r = SubsetPermutationResult(**_subset_kwargs())
        assert isinstance(r, PermutationResult)
        assert isinstance(r, AnalysisResult)

    def test_k3_single_subset(self):
        """k=3 with 3 raters gives 1 subset."""
        r = SubsetPermutationResult(**_subset_kwargs(
            subset_values=dict(SUBSETS_K3),
            subjects_per_subset=dict(SUBJECTS_K3),
            n_subsets=1,
            subset_k=3,
            mean_over_subsets=SUBSET_K3_MEAN,
        ))
        assert r.subset_k == 3
        assert r.n_subsets == 1
        assert r.mean_over_subsets == SUBSET_K3_MEAN


class TestSubsetPermutationResultValidation:
    """Each validation constraint tested with one failing input."""

    def test_subset_k_not_int(self):
        with pytest.raises(TypeError, match="subset_k must be int"):
            SubsetPermutationResult(**_subset_kwargs(subset_k=2.0))

    def test_subset_k_too_small(self):
        with pytest.raises(ValueError, match="subset_k must be >= 2"):
            SubsetPermutationResult(**_subset_kwargs(subset_k=1))

    def test_subset_values_not_dict(self):
        with pytest.raises(TypeError, match="subset_values must be dict"):
            SubsetPermutationResult(**_subset_kwargs(
                subset_values=[(0, 1, 0.6)],
            ))

    def test_subset_values_bad_keys(self):
        """Keys must be k-tuples matching subset_k."""
        with pytest.raises(ValueError, match="subset_values keys must be"):
            SubsetPermutationResult(**_subset_kwargs(
                subset_values={0: 0.6, 1: 0.5, 2: 0.4},
            ))

    def test_n_subsets_not_int(self):
        with pytest.raises(TypeError, match="n_subsets must be int"):
            SubsetPermutationResult(**_subset_kwargs(n_subsets=3.0))

    def test_n_subsets_zero(self):
        with pytest.raises(ValueError, match="n_subsets must be >= 1"):
            SubsetPermutationResult(**_subset_kwargs(n_subsets=0))

    def test_n_subsets_mismatch(self):
        with pytest.raises(
            ValueError, match="subset_values has .* entries but n_subsets"
        ):
            SubsetPermutationResult(**_subset_kwargs(n_subsets=5))

    def test_subjects_per_subset_keys_mismatch(self):
        with pytest.raises(
            ValueError, match="subjects_per_subset must have the same keys"
        ):
            SubsetPermutationResult(**_subset_kwargs(
                subjects_per_subset={(0, 1): 100, (0, 2): 80, (2, 3): 90},
            ))

    def test_mean_not_numeric(self):
        with pytest.raises(TypeError, match="mean_over_subsets must be float"):
            SubsetPermutationResult(**_subset_kwargs(mean_over_subsets="0.5"))

    def test_bootstrap_mean_method_invalid(self):
        with pytest.raises(ValueError, match="bootstrap_mean_method must be"):
            SubsetPermutationResult(**_subset_with_bootstrap(
                bootstrap_mean_method="invalid",
            ))

    def test_bootstrap_fields_without_method(self):
        with pytest.raises(ValueError, match="bootstrap_mean_method must be set"):
            SubsetPermutationResult(**_subset_kwargs(
                bootstrap_ci_of_mean=(0.4, 0.6),
            ))

    def test_method_without_bootstrap_fields(self):
        with pytest.raises(
            ValueError, match="no mean bootstrap fields are populated"
        ):
            SubsetPermutationResult(**_subset_kwargs(
                bootstrap_mean_method="correlated",
            ))

    def test_bootstrap_ci_of_mean_must_be_an_interval(self):
        """Bounds ordering is the interval's business; its type is ours."""
        with pytest.raises(TypeError, match="must be BootstrapCIResult"):
            SubsetPermutationResult(**_subset_with_bootstrap(
                bootstrap_ci_of_mean=(0.7, 0.3),
            ))

    def test_per_subset_ci_must_be_an_interval(self):
        with pytest.raises(TypeError, match="bootstrap CI for subset"):
            SubsetPermutationResult(**_subset_with_bootstrap(
                bootstrap_cis_per_subset={
                    (0, 1): (0.50, 0.70),
                    (0, 2): _ci(0.30, 0.50),
                    (1, 2): _ci(0.40, 0.60),
                },
            ))

    def test_subsets_may_fail_independently(self):
        """Routine here rather than exceptional.

        A subset with fewer than two measurable subjects has nothing to
        resample and gets an interval that ran no method, while its
        siblings on the same data come back normally.
        """
        r = SubsetPermutationResult(**_subset_with_bootstrap(
            bootstrap_cis_per_subset={
                (0, 1): _ci(0.50, 0.70),
                (0, 2): _ci(0.30, 0.50, n_failed=137),
                (1, 2): _undefined_ci(),
            },
        ))
        cis = r.bootstrap_cis_per_subset
        assert cis[(0, 1)].n_failed == 0
        assert cis[(0, 2)].n_failed == 137
        assert cis[(1, 2)].method_used is None
        assert cis[(1, 2)].is_defined is False

    def test_bootstrap_dist_wrong_ndim(self):
        with pytest.raises(ValueError, match="must be 1-dimensional"):
            SubsetPermutationResult(**_subset_with_bootstrap(
                bootstrap_distribution_of_mean=np.array([[0.5, 0.6]]),
            ))

    def test_bootstrap_dist_empty(self):
        with pytest.raises(ValueError, match="must have at least 1 value"):
            SubsetPermutationResult(**_subset_with_bootstrap(
                bootstrap_distribution_of_mean=np.array([]),
            ))

    def test_per_subset_cis_required(self):
        """The per-subset bootstrap always runs, so None is refused."""
        with pytest.raises(TypeError, match="always runs"):
            SubsetPermutationResult(**_subset_kwargs(bootstrap_cis_per_subset=None))

    def test_draws_may_be_left_out(self):
        """None is what keep_draws: false leaves; the intervals stay."""
        r = SubsetPermutationResult(**_subset_kwargs(bootstrap_distributions_per_subset=None))
        assert len(r.bootstrap_cis_per_subset) == 3

    def test_draws_kept_past_the_old_limit(self):
        """2.0.0 refused draws above 120 subsets; there is no limit now."""
        values = _many_units(17, 2)
        r = SubsetPermutationResult(**_subset_kwargs(
            subset_values=values,
            subjects_per_subset={key: 10 for key in values},
            n_subsets=len(values),
        ))
        assert len(r.bootstrap_distributions_per_subset) == 136

    def test_draws_for_every_subset_or_none(self):
        dists = _dists_for(SUBSETS_K2)
        del dists[(1, 2)]
        with pytest.raises(ValueError, match="same keys as subset_values"):
            SubsetPermutationResult(**_subset_kwargs(bootstrap_distributions_per_subset=dists))

    def test_bootstrap_cis_per_subset_keys_mismatch(self):
        with pytest.raises(
            ValueError,
            match="bootstrap_cis_per_subset must have the same keys",
        ):
            SubsetPermutationResult(**_subset_with_bootstrap(
                bootstrap_cis_per_subset={
                    (0, 1): (0.5, 0.7),
                    (0, 2): (0.3, 0.5),
                    # Missing (1, 2)
                },
            ))

    def test_error_prefix(self):
        with pytest.raises(
            TypeError, match=r"\[SubsetPermutationResult\] Error"
        ):
            SubsetPermutationResult(
                **_subset_kwargs(subset_values="not a dict")
            )

    def test_parent_validation_propagated(self):
        """Parent error messages are wrapped with child prefix."""
        with pytest.raises(
            TypeError, match=r"\[SubsetPermutationResult\] Error"
        ):
            SubsetPermutationResult(**_subset_kwargs(
                balance_profile={"not": "a profile"},
            ))


class TestSubsetPermutationResultAccessor:
    """get_bootstrap_distributions_per_subset() accessor method."""

    def test_returns_distributions_when_stored(self):
        """When distributions are stored, accessor returns them."""
        dists = _dists_for(SUBSETS_K2)
        r = SubsetPermutationResult(**_subset_with_bootstrap(
            bootstrap_distributions_per_subset=dists,
        ))
        result = r.get_bootstrap_distributions_per_subset()
        assert result is dists

    def test_raises_when_draws_not_kept(self):
        """The one reason draws can be missing: keep_draws was false.

        It must not claim the bootstrap was not performed -- it always is.
        """
        r = SubsetPermutationResult(**_subset_kwargs(bootstrap_distributions_per_subset=None))
        with pytest.raises(RuntimeError) as exc:
            r.get_bootstrap_distributions_per_subset()
        message = str(exc.value)
        assert "keep_draws false" in message
        assert "bootstrap_cis_per_subset" in message
        assert "not performed" not in message


class TestSubsetPermutationResultStr:
    """__str__ output."""

    def test_no_bootstrap(self):
        s = str(SubsetPermutationResult(**_subset_kwargs()))
        assert s.startswith("Subsets of 2 raters, Fleiss' κ, 3 subsets: mean = 0.500")
        assert "each subset on the subjects at least two of its raters rated" in s

    def test_with_bootstrap(self):
        s = str(SubsetPermutationResult(**_subset_with_bootstrap()))
        assert "95% CI" in s
        assert "correlated" in s


class TestSubsetPermutationResultSummary:
    """summary() returns observed + bootstrap keys."""

    def test_no_bootstrap(self):
        d = SubsetPermutationResult(**_subset_kwargs()).summary()
        assert d["subset_k"] == 2
        assert d["n_subsets"] == 3
        assert d["mean_over_subsets"] == SUBSET_K2_MEAN
        assert d["coefficient_name"] == "fleiss_kappa"
        # Tuple keys stringified for JSON
        assert "(0, 1)" in d["subset_values"]

    def test_with_bootstrap(self):
        d = SubsetPermutationResult(**_subset_with_bootstrap()).summary()
        assert d["bootstrap_mean_method"] == "correlated"
        assert d["mean_ci_lower"] == 0.45
        assert d["mean_ci_upper"] == 0.55
        # Draws of the mean are 0.48, 0.50, 0.52
        assert d["bootstrap_draws_mean"] == pytest.approx(0.50)
        assert d["bootstrap_se"] == pytest.approx(0.02)
        assert "bootstrap_n" not in d
        assert d["subset_distributions_stored"] is True
        assert len(d["bootstrap_cis_per_subset"]) == 3

    def test_subset_distributions_stored_flag(self):
        """Draws not kept: the flag is False."""
        d = SubsetPermutationResult(**_subset_kwargs(
            bootstrap_distributions_per_subset=None,
        )).summary()
        assert d["subset_distributions_stored"] is False


# -- Draws agree with intervals: pairwise and subset alike ------------------ #

_PERMUTATION_KINDS = [
    pytest.param(
        (PairwisePermutationResult, _pairwise_kwargs, "pair", PAIRS),
        id="pairwise",
    ),
    pytest.param(
        (SubsetPermutationResult, _subset_kwargs, "subset", SUBSETS_K2),
        id="subset",
    ),
]


@pytest.mark.parametrize("kind", _PERMUTATION_KINDS)
class TestPermutationDrawsMatchIntervals:
    """The mean's draws agree with the mean's interval; each stored unit's
    draws agree with that unit's interval."""

    def test_mean_draws_without_interval_refused(self, kind):
        cls, kwargs, _, _ = kind
        with pytest.raises(ValueError, match="not a BootstrapCIResult"):
            cls(**kwargs(
                bootstrap_mean_method="correlated",
                bootstrap_distribution_of_mean=np.array([0.5, 0.5]),
            ))

    def test_mean_interval_without_draws_refused(self, kind):
        cls, kwargs, _, _ = kind
        with pytest.raises(ValueError, match="comes with the draws"):
            cls(**kwargs(
                bootstrap_mean_method="correlated",
                bootstrap_ci_of_mean=_ci(0.45, 0.55),
            ))

    def test_mean_draw_count_must_match(self, kind):
        cls, kwargs, _, _ = kind
        with pytest.raises(ValueError, match="has 2 draws but .* resampled 1000"):
            cls(**kwargs(
                bootstrap_mean_method="correlated",
                bootstrap_ci_of_mean=_ci(0.45, 0.55),
                bootstrap_distribution_of_mean=np.array([0.5, 0.5]),
            ))

    def test_mean_failed_count_must_match(self, kind):
        cls, kwargs, _, _ = kind
        with pytest.raises(ValueError, match="1 unusable .* counted 0 failed"):
            cls(**kwargs(
                bootstrap_mean_method="correlated",
                bootstrap_ci_of_mean=_ci(0.45, 0.55, n_resamples=3),
                bootstrap_distribution_of_mean=np.array([np.nan, 0.5, 0.5]),
            ))

    def test_every_mean_draw_failed_is_consistent(self, kind):
        cls, kwargs, _, _ = kind
        r = cls(**kwargs(
            bootstrap_mean_method="correlated",
            bootstrap_ci_of_mean=_undefined_ci(),
            bootstrap_distribution_of_mean=np.full(200, np.nan),
        ))
        d = r.summary()
        assert d["bootstrap_draws_mean"] is None
        assert d["bootstrap_se"] is None

    def test_nan_mean_draw_does_not_poison_the_summary(self, kind):
        """One NaN draw of the mean used to make both numbers NaN."""
        cls, kwargs, _, _ = kind
        d = cls(**kwargs(
            bootstrap_mean_method="correlated",
            bootstrap_ci_of_mean=_ci(0.45, 0.55, n_resamples=4, n_failed=1),
            bootstrap_distribution_of_mean=np.array([np.nan, 0.4, 0.5, 0.6]),
        )).summary()
        assert d["bootstrap_draws_mean"] == pytest.approx(0.5)
        assert d["bootstrap_se"] == pytest.approx(0.1)

    def test_unit_draw_count_must_match(self, kind):
        cls, kwargs, unit, values = kind
        draws = _dists_for(values)
        draws[(0, 2)] = np.array([0.4, 0.5, 0.6])
        with pytest.raises(
            ValueError,
            match=rf"{unit} \(0, 2\) has 3 draws but its interval resampled 1000",
        ):
            cls(**kwargs(**{f"bootstrap_distributions_per_{unit}": draws}))

    def test_unit_failed_count_must_match(self, kind):
        cls, kwargs, unit, values = kind
        draws = _dists_for(values)
        draws[(0, 2)] = draws[(0, 2)].copy()
        draws[(0, 2)][0] = np.nan
        with pytest.raises(
            ValueError, match=rf"{unit} \(0, 2\) has 1 unusable"
        ):
            cls(**kwargs(**{f"bootstrap_distributions_per_{unit}": draws}))

    def test_one_unit_failing_entirely_is_consistent(self, kind):
        """A unit with nothing to resample sits beside healthy ones."""
        cls, kwargs, unit, values = kind
        cis = _cis_for(values)
        cis[(1, 2)] = _undefined_ci(n_resamples=1000, n_failed=1000)
        r = cls(**kwargs(**{f"bootstrap_cis_per_{unit}": cis}))
        stored = getattr(r, f"bootstrap_distributions_per_{unit}")
        assert np.isnan(stored[(1, 2)]).all()
        assert np.isfinite(stored[(0, 1)]).all()


# ============================================================================ #
#                         MultiKSubsetResult                                   #
# ============================================================================ #


class TestMultiKSubsetResultConstruction:
    """Valid construction with single and multiple k values."""

    def test_two_ks(self):
        r = MultiKSubsetResult(**_multi_k_kwargs())
        assert r.ks == [2, 3]
        assert len(r.per_k_results) == 2
        assert r.per_k_results[2].subset_k == 2
        assert r.per_k_results[3].subset_k == 3

    def test_single_k(self):
        """Single k produces a dict with one key."""
        bp = _bp()
        result_k2 = SubsetPermutationResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            subset_values=dict(SUBSETS_K2),
            subjects_per_subset=dict(SUBJECTS_K2),
            n_subsets=3,
            subset_k=2,
            mean_over_subsets=SUBSET_K2_MEAN,
            bootstrap_mean_method=None,
            bootstrap_ci_of_mean=None,
            bootstrap_distribution_of_mean=None,
            bootstrap_cis_per_subset=_cis_for(SUBSETS_K2),
            bootstrap_distributions_per_subset=_dists_for(SUBSETS_K2),
        )
        r = MultiKSubsetResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            per_k_results={2: result_k2},
            ks=[2],
        )
        assert r.ks == [2]
        assert len(r.per_k_results) == 1

    def test_is_permutation_result(self):
        r = MultiKSubsetResult(**_multi_k_kwargs())
        assert isinstance(r, PermutationResult)
        assert isinstance(r, AnalysisResult)


class TestMultiKSubsetResultValidation:
    """Each validation constraint tested with one failing input."""

    def test_ks_not_list(self):
        with pytest.raises(TypeError, match="ks must be list"):
            MultiKSubsetResult(**_multi_k_kwargs(ks=(2, 3)))

    def test_ks_empty(self):
        with pytest.raises(ValueError, match="ks must not be empty"):
            MultiKSubsetResult(**_multi_k_kwargs(ks=[]))

    def test_ks_element_not_int(self):
        with pytest.raises(TypeError, match="ks elements must be int"):
            MultiKSubsetResult(**_multi_k_kwargs(ks=[2.0, 3]))

    def test_ks_value_too_small(self):
        with pytest.raises(ValueError, match="ks values must be >= 2"):
            MultiKSubsetResult(**_multi_k_kwargs(ks=[1, 3]))

    def test_ks_not_sorted(self):
        with pytest.raises(ValueError, match="ks must be sorted"):
            MultiKSubsetResult(**_multi_k_kwargs(ks=[3, 2]))

    def test_ks_duplicates(self):
        with pytest.raises(ValueError, match="ks must not contain duplicates"):
            MultiKSubsetResult(**_multi_k_kwargs(ks=[2, 2]))

    def test_per_k_results_not_dict(self):
        with pytest.raises(TypeError, match="per_k_results must be dict"):
            MultiKSubsetResult(**_multi_k_kwargs(per_k_results=[]))

    def test_per_k_results_keys_mismatch(self):
        with pytest.raises(ValueError, match="per_k_results keys .* do not match"):
            # ks says [2, 3] but per_k_results only has 2
            kw = _multi_k_kwargs()
            del kw["per_k_results"][3]
            MultiKSubsetResult(**kw)

    def test_per_k_results_value_wrong_type(self):
        with pytest.raises(TypeError, match="must be SubsetPermutationResult"):
            kw = _multi_k_kwargs()
            kw["per_k_results"][2] = "not a result"
            MultiKSubsetResult(**kw)

    def test_cross_validate_coefficient_name(self):
        """Per-k results must share parent's coefficient_name."""
        bp = _bp()
        wrong_result = SubsetPermutationResult(
            coefficient_name="scotts_pi",  # Mismatch!
            ci_method="percentile",
            balance_profile=bp,
            subset_values=dict(SUBSETS_K2),
            subjects_per_subset=dict(SUBJECTS_K2),
            n_subsets=3,
            subset_k=2,
            mean_over_subsets=SUBSET_K2_MEAN,
            bootstrap_mean_method=None,
            bootstrap_ci_of_mean=None,
            bootstrap_distribution_of_mean=None,
            bootstrap_cis_per_subset=_cis_for(SUBSETS_K2),
            bootstrap_distributions_per_subset=_dists_for(SUBSETS_K2),
        )
        with pytest.raises(ValueError, match="coefficient_name="):
            MultiKSubsetResult(
                coefficient_name="fleiss_kappa",
                ci_method="percentile",
                balance_profile=bp,
                per_k_results={2: wrong_result},
                ks=[2],
            )

    def test_cross_validate_balance_profile_identity(self):
        """Per-k results must share the same balance_profile object."""
        bp_parent = _bp()
        bp_child = _bp()  # Different object, same content
        result_k2 = SubsetPermutationResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp_child,  # Different object!
            subset_values=dict(SUBSETS_K2),
            subjects_per_subset=dict(SUBJECTS_K2),
            n_subsets=3,
            subset_k=2,
            mean_over_subsets=SUBSET_K2_MEAN,
            bootstrap_mean_method=None,
            bootstrap_ci_of_mean=None,
            bootstrap_distribution_of_mean=None,
            bootstrap_cis_per_subset=_cis_for(SUBSETS_K2),
            bootstrap_distributions_per_subset=_dists_for(SUBSETS_K2),
        )
        with pytest.raises(ValueError, match="different balance_profile"):
            MultiKSubsetResult(
                coefficient_name="fleiss_kappa",
                ci_method="percentile",
                balance_profile=bp_parent,
                per_k_results={2: result_k2},
                ks=[2],
            )

    def test_cross_validate_subset_k_matches_key(self):
        """Each per-k result's subset_k must match its dict key."""
        bp = _bp()
        result_wrong_k = SubsetPermutationResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            subset_values=dict(SUBSETS_K2),
            subjects_per_subset=dict(SUBJECTS_K2),
            n_subsets=3,
            subset_k=2,  # subset_k=2 but stored under key 3
            mean_over_subsets=SUBSET_K2_MEAN,
            bootstrap_mean_method=None,
            bootstrap_ci_of_mean=None,
            bootstrap_distribution_of_mean=None,
            bootstrap_cis_per_subset=_cis_for(SUBSETS_K2),
            bootstrap_distributions_per_subset=_dists_for(SUBSETS_K2),
        )
        with pytest.raises(ValueError, match="subset_k=2, expected 3"):
            MultiKSubsetResult(
                coefficient_name="fleiss_kappa",
                ci_method="percentile",
                balance_profile=bp,
                per_k_results={3: result_wrong_k},
                ks=[3],
            )

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[MultiKSubsetResult\] Error"):
            MultiKSubsetResult(**_multi_k_kwargs(ks=[]))


class TestMultiKSubsetResultContainerProtocol:
    """__getitem__, __len__, __iter__."""

    def test_getitem(self):
        r = MultiKSubsetResult(**_multi_k_kwargs())
        assert r[2].subset_k == 2
        assert r[3].subset_k == 3

    def test_getitem_missing_key(self):
        r = MultiKSubsetResult(**_multi_k_kwargs())
        with pytest.raises(KeyError, match="k=5 not in results"):
            r[5]

    def test_len(self):
        r = MultiKSubsetResult(**_multi_k_kwargs())
        assert len(r) == 2

    def test_iter(self):
        r = MultiKSubsetResult(**_multi_k_kwargs())
        assert list(r) == [2, 3]


class TestMultiKSubsetResultStr:
    """__str__ output."""

    def test_basic(self):
        s = str(MultiKSubsetResult(**_multi_k_kwargs()))
        assert s.startswith("Subsets, Fleiss' κ:")
        assert "k=2   3 subsets   mean 0.500" in s
        assert "k=3   1 subset    mean 0.550" in s
        assert "no interval for the mean" in s

    def test_with_bootstrap(self):
        """Bootstrap CI shown per k when available."""
        bp = _bp()
        boot_dist = np.array([0.48, 0.50, 0.52])
        result_k2 = SubsetPermutationResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            subset_values=dict(SUBSETS_K2),
            subjects_per_subset=dict(SUBJECTS_K2),
            n_subsets=3,
            subset_k=2,
            mean_over_subsets=SUBSET_K2_MEAN,
            bootstrap_mean_method="correlated",
            bootstrap_ci_of_mean=_ci(0.45, 0.55, n_resamples=3),
            bootstrap_distribution_of_mean=boot_dist,
            bootstrap_cis_per_subset=(cis_k2 := {
                (0, 1): _ci(0.50, 0.70),
                (0, 2): _ci(0.30, 0.50),
                (1, 2): _ci(0.40, 0.60),
            }),
            bootstrap_distributions_per_subset=_dists_matching(cis_k2),
        )
        result_k3 = SubsetPermutationResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            subset_values=dict(SUBSETS_K3),
            subjects_per_subset=dict(SUBJECTS_K3),
            n_subsets=1,
            subset_k=3,
            mean_over_subsets=SUBSET_K3_MEAN,
            bootstrap_mean_method=None,
            bootstrap_ci_of_mean=None,
            bootstrap_distribution_of_mean=None,
            bootstrap_cis_per_subset=_cis_for(SUBSETS_K3),
            bootstrap_distributions_per_subset=_dists_for(SUBSETS_K3),
        )
        r = MultiKSubsetResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            per_k_results={2: result_k2, 3: result_k3},
            ks=[2, 3],
        )
        s = str(r)
        assert "mean 0.500   [0.450, 0.550]" in s
        assert "intervals: 95% CI, percentile, n=3, correlated bootstrap" in s
        # k=3 has no bootstrap of the mean, so its line ends at the mean
        assert any(line.rstrip().endswith("mean 0.550") for line in s.splitlines())

    def test_undefined_interval_is_named_not_omitted(self):
        """A k that bootstrapped and got nothing must not print like a k
        that was never bootstrapped at all."""
        bp = _bp()
        result_k2 = SubsetPermutationResult(**_subset_kwargs(
            balance_profile=bp,
            bootstrap_mean_method="correlated",
            bootstrap_ci_of_mean=_undefined_ci(),
            bootstrap_distribution_of_mean=np.full(200, np.nan),
        ))
        r = MultiKSubsetResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            per_k_results={2: result_k2},
            ks=[2],
        )
        assert "no CI: every draw failed" in str(r)

    def test_failed_draws_are_named(self):
        """An interval resting on a fraction of its draws says so."""
        bp = _bp()
        result_k2 = SubsetPermutationResult(**_subset_kwargs(
            balance_profile=bp,
            bootstrap_mean_method="correlated",
            bootstrap_ci_of_mean=_ci(0.45, 0.55, n_failed=312),
            bootstrap_distribution_of_mean=_draws_for(
                _ci(0.45, 0.55, n_failed=312)
            ),
        ))
        r = MultiKSubsetResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            per_k_results={2: result_k2},
            ks=[2],
        )
        assert "312 failed" in str(r)


class TestMultiKSubsetResultSummary:
    """summary() returns per-k summaries."""

    def test_basic(self):
        d = MultiKSubsetResult(**_multi_k_kwargs()).summary()
        assert d["coefficient_name"] == "fleiss_kappa"
        assert d["ci_method"] == "percentile"
        assert d["ks"] == [2, 3]
        assert 2 in d["per_k"]
        assert 3 in d["per_k"]
        assert d["per_k"][2]["subset_k"] == 2
        assert d["per_k"][3]["subset_k"] == 3

    def test_inherits_balance(self):
        d = MultiKSubsetResult(**_multi_k_kwargs()).summary()
        assert isinstance(d["balance"], dict)


# ============================================================================ #
#                           FleissKappaResult                                  #
# ============================================================================ #


class TestFleissKappaResult:
    """Fleiss' kappa: ci_method includes parametric_null, weights_type, is_exact."""

    def test_valid_construction(self):
        r = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_pooled(0.6),
        )
        assert r.point == 0.6
        assert r.ci_method == "bca"
        assert r.weights_type == "identity"

    def test_is_point_estimate_result(self):
        """FleissKappaResult inherits from PointEstimateResult."""
        r = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_pooled(0.6),
        )
        assert isinstance(r, PointEstimateResult)

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            FleissKappaResult(
                0.6, None, None, None,
                ci_method="bca", weights_type="identity",
                **_pooled(0.6),
            )

    def test_ci_method_parametric_null_allowed(self):
        r = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="parametric_null", weights_type="identity",
            **_pooled(0.6),
        )
        assert r.ci_method == "parametric_null"

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            FleissKappaResult(
                0.6, None, None, _bp(),
                ci_method=42, weights_type="identity",
                **_pooled(0.6),
            )

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            FleissKappaResult(
                0.6, None, None, _bp(),
                ci_method="invalid", weights_type="identity",
                **_pooled(0.6),
            )

    @pytest.mark.parametrize("wt", ["identity", "linear", "quadratic", "custom"])
    def test_valid_weights_types(self, wt):
        r = FleissKappaResult(
            0.6, None, None, _bp(), ci_method="bca", weights_type=wt,
            **_pooled(0.6, weighted=wt != "identity"),
        )
        assert r.weights_type == wt

    def test_weights_type_not_str(self):
        with pytest.raises(TypeError, match="weights_type must be str"):
            FleissKappaResult(
                0.6, None, None, _bp(), ci_method="bca", weights_type=42,
                **_pooled(0.6),
            )

    def test_weights_type_invalid(self):
        with pytest.raises(ValueError, match="weights_type must be"):
            FleissKappaResult(
                0.6, None, None, _bp(),
                ci_method="bca", weights_type="ordinal",
                **_pooled(0.6),
            )

    def test_is_exact_balanced(self):
        r = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="percentile", weights_type="identity",
            **_pooled(0.6),
        )
        assert r.is_exact is True

    def test_is_exact_unbalanced(self):
        r = FleissKappaResult(
            0.6, None, None, _unbalanced_bp(),
            ci_method="percentile", weights_type="identity",
            **_pooled(0.6),
        )
        assert r.is_exact is False

    def test_str_balanced_unweighted(self):
        s = str(FleissKappaResult(
            0.523, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_pooled(0.523),
        ))
        assert "Fleiss' κ" in s
        assert "unweighted" in s
        assert "exact formula" in s

    def test_str_unbalanced_weighted(self):
        s = str(FleissKappaResult(
            0.523, None, None, _unbalanced_bp(),
            ci_method="bca", weights_type="linear",
            **_pooled(0.523, weighted=True),
        ))
        assert "linear weights" in s
        assert "approximation (unbalanced data)" in s

    def test_summary_keys(self):
        d = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_pooled(0.6),
        ).summary()
        assert d["is_exact"] is True
        assert d["ci_method"] == "bca"
        assert d["weights_type"] == "identity"
        assert "balance" in d

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[FleissKappaResult\] Error"):
            FleissKappaResult(
                0.6, None, None, None,
                ci_method="bca", weights_type="identity",
                **_pooled(0.6),
            )

    def test_parent_error_chaining(self):
        """Parent TypeError wrapped with Fleiss prefix."""
        with pytest.raises(TypeError, match=r"\[FleissKappaResult\] Error"):
            FleissKappaResult(
                "bad", None, None, _bp(),
                ci_method="bca", weights_type="identity",
                **_pooled(0.6),
            )


class TestFleissKappaResultParts:
    """Po and Pe always there, weighted as the result is, Pe from the profile."""

    def _result(self, weights_type="identity", **fields):
        kwargs = dict(_pooled(0.6, weighted=weights_type != "identity"), **fields)
        return FleissKappaResult(
            0.6, None, None, _bp(), ci_method="bca", weights_type=weights_type,
            **kwargs,
        )

    def test_parts_not_computed_refused(self):
        with pytest.raises(TypeError, match=r"parts must be agreement parts \(Po and Pe\)"):
            self._result(parts=NotComputed("no reason"))

    def test_disagreement_parts_refused(self):
        with pytest.raises(TypeError, match="parts must be agreement parts"):
            self._result(parts=AgreementParts(0.2, 0.5, "disagreement", False))

    def test_weighted_parts_on_unweighted_result_refused(self):
        parts = _pooled(0.6, weighted=True)["parts"]
        with pytest.raises(ValueError, match="parts.weighted is True but weights_type is 'identity'"):
            self._result(parts=parts)

    def test_unweighted_parts_on_weighted_result_refused(self):
        parts = _pooled(0.6)["parts"]
        with pytest.raises(ValueError, match="parts.weighted is False but weights_type is 'linear'"):
            self._result(weights_type="linear", parts=parts)

    def test_computed_kappa_max_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be NotComputed for Fleiss' kappa"):
            self._result(kappa_max=1.0)

    def test_pe_from_other_data_refused(self):
        """Parts that rebuild the point, but not from the profile's shares."""
        other = AgreementParts(0.8, 0.5, "agreement", False)
        with pytest.raises(ValueError, match="describe different data"):
            self._result(parts=other)

    def test_weighted_pe_not_checked_against_shares(self):
        """Weighted Pe is p W p, which needs the weights; not checked here."""
        weighted = AgreementParts(0.8, 0.5, "agreement", True)
        assert self._result(weights_type="linear", parts=weighted).parts is weighted


# ============================================================================ #
#                            ScottsPiResult                                    #
# ============================================================================ #


class TestScottsPiResult:
    """Scott's pi: like Fleiss but no weights_type."""

    def test_valid_construction(self):
        r = ScottsPiResult(0.5, None, None, _bp(), ci_method="percentile", **_pooled(0.5))
        assert r.point == 0.5

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            ScottsPiResult(0.5, None, None, None, ci_method="bca", **_pooled(0.5))

    def test_ci_method_parametric_null_allowed(self):
        r = ScottsPiResult(0.5, None, None, _bp(), ci_method="parametric_null", **_pooled(0.5))
        assert r.ci_method == "parametric_null"

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            ScottsPiResult(0.5, None, None, _bp(), ci_method=42, **_pooled(0.5))

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            ScottsPiResult(0.5, None, None, _bp(), ci_method="wrong", **_pooled(0.5))

    def test_is_exact_balanced(self):
        assert ScottsPiResult(
            0.5, None, None, _bp(), ci_method="bca",
            **_pooled(0.5),
        ).is_exact is True

    def test_is_exact_unbalanced(self):
        assert ScottsPiResult(
            0.5, None, None, _unbalanced_bp(), ci_method="bca",
            **_pooled(0.5),
        ).is_exact is False

    def test_str_balanced(self):
        s = str(ScottsPiResult(0.5, None, None, _bp(), ci_method="bca", **_pooled(0.5)))
        assert "Scott's π" in s
        assert "exact formula" in s

    def test_str_unbalanced(self):
        s = str(ScottsPiResult(
            0.5, None, None, _unbalanced_bp(), ci_method="bca",
            **_pooled(0.5),
        ))
        assert "approximation (unbalanced data)" in s

    def test_summary_keys(self):
        d = ScottsPiResult(0.5, None, None, _bp(), ci_method="bca", **_pooled(0.5)).summary()
        assert d["is_exact"] is True
        assert d["ci_method"] == "bca"
        assert "weights_type" not in d

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[ScottsPiResult\] Error"):
            ScottsPiResult(0.5, None, None, None, ci_method="bca", **_pooled(0.5))


class TestScottsPiResultParts:
    """Po and Pe always there, unweighted, Pe from the profile."""

    def _result(self, **fields):
        kwargs = dict(_pooled(0.5), **fields)
        return ScottsPiResult(0.5, None, None, _bp(), ci_method="bca", **kwargs)

    def test_valid(self):
        assert self._result().parts.kind == "agreement"

    def test_parts_not_computed_refused(self):
        with pytest.raises(TypeError, match="parts must be agreement parts"):
            self._result(parts=NotComputed("no reason"))

    def test_weighted_parts_refused(self):
        with pytest.raises(ValueError, match="no weighted form"):
            self._result(parts=_pooled(0.5, weighted=True)["parts"])

    def test_computed_kappa_max_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be NotComputed for Scott's pi"):
            self._result(kappa_max=1.0)

    def test_pe_from_other_data_refused(self):
        other = AgreementParts(0.75, 0.5, "agreement", False)
        with pytest.raises(ValueError, match="describe different data"):
            self._result(parts=other)


# ============================================================================ #
#                            GwetsAC1Result                                    #
# ============================================================================ #


class TestGwetsAC1Result:
    """Gwet's AC1: no parametric_null (explicit 'not implemented' rejection)."""

    def test_valid_construction(self):
        r = GwetsAC1Result(0.7, None, None, _bp(), ci_method="bca", **_gwet(0.7))
        assert r.point == 0.7

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            GwetsAC1Result(0.7, None, None, None, ci_method="bca", **_gwet(0.7))

    def test_parametric_null_rejected(self):
        """Gwet's AC1 explicitly rejects parametric_null with 'not implemented'."""
        with pytest.raises(ValueError, match="parametric_null.*not implemented"):
            GwetsAC1Result(0.7, None, None, _bp(), ci_method="parametric_null", **_gwet(0.7))

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            GwetsAC1Result(0.7, None, None, _bp(), ci_method=42, **_gwet(0.7))

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            GwetsAC1Result(0.7, None, None, _bp(), ci_method="invalid", **_gwet(0.7))

    def test_is_exact_balanced(self):
        assert GwetsAC1Result(
            0.7, None, None, _bp(), ci_method="bca",
            **_gwet(0.7),
        ).is_exact is True

    def test_is_exact_unbalanced(self):
        assert GwetsAC1Result(
            0.7, None, None, _unbalanced_bp(), ci_method="bca",
            **_gwet(0.7),
        ).is_exact is False

    def test_str_balanced(self):
        s = str(GwetsAC1Result(0.7, None, None, _bp(), ci_method="bca", **_gwet(0.7)))
        assert "Gwet's AC1" in s
        assert "exact formula" in s

    def test_summary_keys(self):
        d = GwetsAC1Result(0.7, None, None, _bp(), ci_method="bca", **_gwet(0.7)).summary()
        assert d["is_exact"] is True
        assert d["ci_method"] == "bca"
        assert "weights_type" not in d

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[GwetsAC1Result\] Error"):
            GwetsAC1Result(0.7, None, None, None, ci_method="bca", **_gwet(0.7))


class TestGwetsAC1ResultParts:
    """Po and Pe always there, unweighted, Pe from the profile's shares."""

    def _result(self, **fields):
        kwargs = dict(_gwet(0.7), **fields)
        return GwetsAC1Result(0.7, None, None, _bp(), ci_method="bca", **kwargs)

    def test_valid(self):
        assert self._result().parts.reproduces(0.7)

    def test_parts_not_computed_refused(self):
        with pytest.raises(TypeError, match="parts must be agreement parts"):
            self._result(parts=NotComputed("no reason"))

    def test_weighted_parts_refused(self):
        weighted = _fields(0.7, _gwet(0.7)["parts"].chance, weighted=True)["parts"]
        with pytest.raises(ValueError, match="Gwet's AC2, a different coefficient"):
            self._result(parts=weighted)

    def test_computed_kappa_max_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be NotComputed for Gwet's AC1"):
            self._result(kappa_max=1.0)

    def test_fleiss_pe_refused(self):
        """Same data, the wrong chance term: Fleiss' Pe in an AC1 result."""
        with pytest.raises(ValueError, match="describe different data"):
            self._result(parts=_pooled(0.7)["parts"])


# ============================================================================ #
#                         BrennanPredigerResult                                #
# ============================================================================ #


class TestBrennanPredigerResult:
    """Brennan-Prediger: percentile/bca only, has weights_type."""

    def test_valid_construction(self):
        r = BrennanPredigerResult(
            0.8, None, None, _bp(),
            ci_method="percentile", weights_type="identity",
            **_uniform(0.8),
        )
        assert r.point == 0.8
        assert r.weights_type == "identity"

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            BrennanPredigerResult(
                0.8, None, None, None,
                ci_method="bca", weights_type="identity",
                **_uniform(0.8),
            )

    def test_parametric_null_not_allowed(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            BrennanPredigerResult(
                0.8, None, None, _bp(),
                ci_method="parametric_null", weights_type="identity",
                **_uniform(0.8),
            )

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            BrennanPredigerResult(
                0.8, None, None, _bp(),
                ci_method=42, weights_type="identity",
                **_uniform(0.8),
            )

    @pytest.mark.parametrize("wt", ["identity", "linear", "quadratic", "custom"])
    def test_valid_weights_types(self, wt):
        r = BrennanPredigerResult(
            0.8, None, None, _bp(), ci_method="bca", weights_type=wt,
            **_uniform(0.8, weighted=wt != "identity"),
        )
        assert r.weights_type == wt

    def test_weights_type_not_str(self):
        with pytest.raises(TypeError, match="weights_type must be str"):
            BrennanPredigerResult(
                0.8, None, None, _bp(), ci_method="bca", weights_type=42,
                **_uniform(0.8),
            )

    def test_weights_type_invalid(self):
        with pytest.raises(ValueError, match="weights_type must be"):
            BrennanPredigerResult(
                0.8, None, None, _bp(),
                ci_method="bca", weights_type="ordinal",
                **_uniform(0.8),
            )

    def test_is_exact_balanced(self):
        assert BrennanPredigerResult(
            0.8, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_uniform(0.8),
        ).is_exact is True

    def test_is_exact_unbalanced(self):
        assert BrennanPredigerResult(
            0.8, None, None, _unbalanced_bp(),
            ci_method="bca", weights_type="identity",
            **_uniform(0.8),
        ).is_exact is False

    def test_str_balanced_unweighted(self):
        s = str(BrennanPredigerResult(
            0.835, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_uniform(0.835),
        ))
        assert "Brennan-Prediger" in s
        assert "unweighted" in s
        assert "exact formula" in s

    def test_summary_keys(self):
        d = BrennanPredigerResult(
            0.8, None, None, _bp(),
            ci_method="bca", weights_type="linear",
            **_uniform(0.8, weighted=True),
        ).summary()
        assert d["is_exact"] is True
        assert d["ci_method"] == "bca"
        assert d["weights_type"] == "linear"

    def test_error_prefix(self):
        with pytest.raises(
            ValueError, match=r"\[BrennanPredigerResult\] Error"
        ):
            BrennanPredigerResult(
                0.8, None, None, None,
                ci_method="bca", weights_type="identity",
                **_uniform(0.8),
            )


class TestBrennanPredigerResultParts:
    """Po and Pe always there, weighted as the result is, Pe from K."""

    def _result(self, weights_type="identity", weight_matrix=None, **fields):
        kwargs = dict(_uniform(0.8, weighted=weights_type != "identity"), **fields)
        return BrennanPredigerResult(
            0.8, None, None, _bp(), ci_method="bca",
            weights_type=weights_type, weight_matrix=weight_matrix, **kwargs,
        )

    def test_valid(self):
        assert self._result().parts.chance == pytest.approx(1 / 3)

    def test_parts_not_computed_refused(self):
        with pytest.raises(TypeError, match="parts must be agreement parts"):
            self._result(parts=NotComputed("no reason"))

    def test_weighted_parts_on_unweighted_result_refused(self):
        with pytest.raises(ValueError, match="parts.weighted is True but weights_type is 'identity'"):
            self._result(parts=_uniform(0.8, weighted=True)["parts"])

    def test_computed_kappa_max_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be NotComputed for Brennan-Prediger"):
            self._result(kappa_max=1.0)

    def test_pe_for_other_k_refused(self):
        """1/2 is Pe for two categories; the profile has three."""
        with pytest.raises(ValueError, match=r"1/K with K = 3 is 0.333333"):
            self._result(parts=_fields(0.8, 0.5)["parts"])

    def test_weighted_pe_is_sum_of_weights_over_k_squared(self):
        w = [[1.0, 0.5, 0.0], [0.5, 1.0, 0.5], [0.0, 0.5, 1.0]]
        pe = sum(map(sum, w)) / 9
        r = self._result("linear", w, parts=_fields(0.8, pe, weighted=True)["parts"])
        assert r.parts.chance == pytest.approx(5 / 9)
        with pytest.raises(ValueError, match=r"sum\(W\)/K\^2 is 0.555556"):
            self._result("linear", w, parts=_fields(0.8, 1 / 3, weighted=True)["parts"])

    def test_weight_matrix_for_other_k_refused(self):
        w = [[1.0, 0.5], [0.5, 1.0]]
        with pytest.raises(ValueError, match="weight_matrix is 2x2 but the agreement profile has 3"):
            self._result("linear", w)


# ============================================================================ #
#                        KrippendorffAlphaResult                               #
# ============================================================================ #


class TestKrippendorffAlphaResult:
    """Krippendorff's alpha: is_exact always True, has weights_type."""

    def test_valid_construction(self):
        r = KrippendorffAlphaResult(
            0.65, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_alpha(0.65),
        )
        assert r.point == 0.65
        assert r.weights_type == "identity"

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            KrippendorffAlphaResult(
                0.65, None, None, None,
                ci_method="bca", weights_type="identity",
                **_alpha(0.65),
            )

    def test_parametric_null_not_allowed(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            KrippendorffAlphaResult(
                0.65, None, None, _bp(),
                ci_method="parametric_null", weights_type="identity",
                **_alpha(0.65),
            )

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            KrippendorffAlphaResult(
                0.65, None, None, _bp(),
                ci_method=42, weights_type="identity",
                **_alpha(0.65),
            )

    @pytest.mark.parametrize("wt", ["identity", "linear", "quadratic", "custom"])
    def test_valid_weights_types(self, wt):
        r = KrippendorffAlphaResult(
            0.65, None, None, _bp(), ci_method="bca", weights_type=wt,
            **_alpha(0.65, weighted=wt != "identity"),
        )
        assert r.weights_type == wt

    def test_weights_type_not_str(self):
        with pytest.raises(TypeError, match="weights_type must be str"):
            KrippendorffAlphaResult(
                0.65, None, None, _bp(), ci_method="bca", weights_type=42,
                **_alpha(0.65),
            )

    def test_weights_type_invalid(self):
        with pytest.raises(ValueError, match="weights_type must be"):
            KrippendorffAlphaResult(
                0.65, None, None, _bp(),
                ci_method="bca", weights_type="ordinal",
                **_alpha(0.65),
            )

    def test_is_exact_always_true_balanced(self):
        """Alpha is always exact, even when balanced."""
        assert KrippendorffAlphaResult(
            0.65, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_alpha(0.65),
        ).is_exact is True

    def test_is_exact_always_true_unbalanced(self):
        """Alpha handles unbalanced data natively -- still exact."""
        assert KrippendorffAlphaResult(
            0.65, None, None, _unbalanced_bp(),
            ci_method="bca", weights_type="identity",
            **_alpha(0.65),
        ).is_exact is True

    def test_str_balanced(self):
        s = str(KrippendorffAlphaResult(
            0.625, None, None, _bp(),
            ci_method="bca", weights_type="identity",
            **_alpha(0.625),
        ))
        assert s.startswith("Krippendorff's α = 0.625")
        assert "exact formula" in s
        # Balanced: no "(unbalanced data)" suffix
        assert "unbalanced" not in s

    def test_str_unbalanced(self):
        """Unbalanced shows 'exact formula (unbalanced data)' -- NOT 'approximation'."""
        s = str(KrippendorffAlphaResult(
            0.625, None, None, _unbalanced_bp(),
            ci_method="bca", weights_type="identity",
            **_alpha(0.625),
        ))
        assert "exact formula (unbalanced data)" in s

    def test_summary_keys(self):
        d = KrippendorffAlphaResult(
            0.65, None, None, _bp(),
            ci_method="bca", weights_type="quadratic",
            **_alpha(0.65, weighted=True),
        ).summary()
        assert d["is_exact"] is True
        assert d["ci_method"] == "bca"
        assert d["weights_type"] == "quadratic"

    def test_error_prefix(self):
        with pytest.raises(
            ValueError, match=r"\[KrippendorffAlphaResult\] Error"
        ):
            KrippendorffAlphaResult(
                0.65, None, None, None,
                ci_method="bca", weights_type="identity",
                **_alpha(0.65),
            )


class TestKrippendorffAlphaResultParts:
    """Do and De always there, weighted as the result is."""

    def _result(self, weights_type="identity", **fields):
        kwargs = dict(_alpha(0.65, weighted=weights_type != "identity"), **fields)
        return KrippendorffAlphaResult(
            0.65, None, None, _bp(), ci_method="bca", weights_type=weights_type,
            **kwargs,
        )

    def test_valid(self):
        assert self._result().parts.kind == "disagreement"

    def test_agreement_parts_refused(self):
        """Po and Pe that happen to rebuild alpha are still the wrong terms."""
        with pytest.raises(TypeError, match=r"disagreement parts \(Do and De\)"):
            self._result(parts=_fields(0.65, 0.5)["parts"])

    def test_parts_not_computed_refused(self):
        with pytest.raises(TypeError, match="disagreement parts"):
            self._result(parts=NotComputed("no reason"))

    def test_unweighted_parts_on_weighted_result_refused(self):
        with pytest.raises(ValueError, match="parts.weighted is False but weights_type is 'quadratic'"):
            self._result("quadratic", parts=_alpha(0.65)["parts"])

    def test_computed_kappa_max_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be NotComputed for Krippendorff's alpha"):
            self._result(kappa_max=1.0)

    def test_parts_that_do_not_rebuild_alpha_refused(self):
        with pytest.raises(ValueError, match="make 0.600000, but the point estimate is 0.650000"):
            self._result(parts=_alpha(0.6)["parts"])


# ============================================================================ #
#                       IndividualMarginalResult                               #
# ============================================================================ #


class TestIndividualMarginalResult:
    """Base for Cohen and Conger: balance_profile and n_raters must not be None."""

    def test_valid_construction(self):
        bp = _bp_with_n_raters(3)
        r = IndividualMarginalResult(0.6, None, None, bp, **_PANEL_NOT_UNDER_TEST)
        assert r.point == 0.6
        assert r.balance_profile is bp

    def test_is_point_estimate_result(self):
        r = IndividualMarginalResult(0.6, None, None, _bp_with_n_raters(3), **_PANEL_NOT_UNDER_TEST)
        assert isinstance(r, PointEstimateResult)

    def test_balance_profile_none_rejected(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            IndividualMarginalResult(0.6, None, None, None, **_PANEL_NOT_UNDER_TEST)

    def test_n_raters_none_rejected(self):
        """BalanceProfile from NijMatrix has n_raters=None -- rejected."""
        with pytest.raises(ValueError, match="n_raters must not be None"):
            IndividualMarginalResult(0.6, None, None, _bp(), **_PANEL_NOT_UNDER_TEST)

    def test_error_prefix(self):
        with pytest.raises(
            ValueError, match=r"\[IndividualMarginalResult\] Error"
        ):
            IndividualMarginalResult(0.6, None, None, None, **_PANEL_NOT_UNDER_TEST)

    def test_parent_validation_propagated(self):
        """Point type error wraps with IndividualMarginalResult prefix."""
        with pytest.raises(
            TypeError, match=r"\[IndividualMarginalResult\] Error"
        ):
            IndividualMarginalResult("bad", None, None, _bp_with_n_raters(3), **_PANEL_NOT_UNDER_TEST)


class TestIndividualMarginalResultProfile:
    """These coefficients need each rater's marginals, one row per rater."""

    def test_profile_from_counts_refused(self):
        with pytest.raises(TypeError, match="rater_marginals must be computed"):
            IndividualMarginalResult(0.6, None, None, _bp_with_n_raters(3), **_NOT_UNDER_TEST)

    def test_rater_count_must_match(self):
        fields = dict(_NOT_UNDER_TEST, agreement_profile=_PAIR)
        with pytest.raises(ValueError, match="marginals for 2 raters but balance_profile has 3"):
            IndividualMarginalResult(0.6, None, None, _bp_with_n_raters(3), **fields)


# ============================================================================ #
#                           CohenKappaResult                                   #
# ============================================================================ #


class TestCohenKappaResult:
    """Cohen's kappa: inherits IndividualMarginalResult, has weights_type."""

    def test_valid_construction(self):
        bp = _bp_with_n_raters(2)
        r = CohenKappaResult(
            0.7, None, None, bp,
            ci_method="percentile", weights_type="identity",
            **_cohen(0.7),
        )
        assert r.point == 0.7
        assert r.weights_type == "identity"

    def test_is_individual_marginal_result(self):
        """Cohen inherits from IndividualMarginalResult (not PointEstimateResult directly)."""
        bp = _bp_with_n_raters(2)
        r = CohenKappaResult(
            0.7, None, None, bp,
            ci_method="bca", weights_type="identity",
            **_cohen(0.7),
        )
        assert isinstance(r, IndividualMarginalResult)
        assert isinstance(r, PointEstimateResult)

    def test_balance_profile_none_rejected(self):
        """Inherited from IndividualMarginalResult."""
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            CohenKappaResult(
                0.7, None, None, None,
                ci_method="bca", weights_type="identity",
                **_cohen(0.7),
            )

    def test_n_raters_none_rejected(self):
        """Inherited from IndividualMarginalResult: n_raters=None rejected."""
        with pytest.raises(ValueError, match="n_raters must not be None"):
            CohenKappaResult(
                0.7, None, None, _bp(),
                ci_method="bca", weights_type="identity",
                **_cohen(0.7),
            )

    def test_parametric_null_allowed(self):
        bp = _bp_with_n_raters(2)
        r = CohenKappaResult(
            0.7, None, None, bp,
            ci_method="parametric_null", weights_type="identity",
            **_cohen(0.7),
        )
        assert r.ci_method == "parametric_null"

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(2),
                ci_method=42, weights_type="identity",
                **_cohen(0.7),
            )

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(2),
                ci_method="invalid", weights_type="identity",
                **_cohen(0.7),
            )

    @pytest.mark.parametrize("wt", ["identity", "linear", "quadratic", "custom"])
    def test_valid_weights_types(self, wt):
        bp = _bp_with_n_raters(2)
        r = CohenKappaResult(
            0.7, None, None, bp, ci_method="bca", weights_type=wt,
            **_cohen(0.7, weighted=wt != "identity"),
        )
        assert r.weights_type == wt

    def test_weights_type_not_str(self):
        with pytest.raises(TypeError, match="weights_type must be str"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(2),
                ci_method="bca", weights_type=42,
                **_cohen(0.7),
            )

    def test_weights_type_invalid(self):
        with pytest.raises(ValueError, match="weights_type must be"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(2),
                ci_method="bca", weights_type="ordinal",
                **_cohen(0.7),
            )

    def test_str_unweighted(self):
        bp = _bp_with_n_raters(2)
        s = str(CohenKappaResult(
            0.723, None, None, bp,
            ci_method="percentile", weights_type="identity",
            **_cohen(0.723),
        ))
        assert "Cohen's κ" in s
        assert "unweighted" in s

    def test_str_weighted(self):
        bp = _bp_with_n_raters(2)
        s = str(CohenKappaResult(
            0.723, None, None, bp,
            ci_method="bca", weights_type="quadratic",
            **_cohen(0.723, weighted=True),
        ))
        assert "quadratic weights" in s

    def test_summary_keys(self):
        bp = _bp_with_n_raters(2)
        d = CohenKappaResult(
            0.7, None, None, bp,
            ci_method="bca", weights_type="linear",
            **_cohen(0.7, weighted=True),
        ).summary()
        assert d["ci_method"] == "bca"
        assert d["weights_type"] == "linear"
        # Cohen has no is_exact
        assert "is_exact" not in d

    def test_error_prefix(self):
        """Cohen's own prefix wraps parent errors."""
        with pytest.raises(ValueError, match=r"\[CohenKappaResult\] Error"):
            CohenKappaResult(
                0.7, None, None, None,
                ci_method="bca", weights_type="identity",
                **_cohen(0.7),
            )


class TestCohenKappaResultParts:
    """Po, Pe and kappa_max checked against the two raters' marginals."""

    def _result(self, point=0.7, weights_type="identity", weight_matrix=None, **fields):
        kwargs = dict(_cohen(point, weighted=weights_type != "identity"), **fields)
        return CohenKappaResult(
            point, None, None, _bp_with_n_raters(2), ci_method="bca",
            weights_type=weights_type, weight_matrix=weight_matrix, **kwargs,
        )

    def test_valid(self):
        r = self._result()
        assert r.kappa_max == 1.0
        assert r.parts.chance == pytest.approx(0.375)

    def test_three_raters_refused(self):
        with pytest.raises(ValueError, match="one pair of raters, got 3"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(3), ci_method="bca",
                weights_type="identity", **dict(_cohen(0.7), agreement_profile=_PANEL),
            )

    def test_kappa_max_not_computed_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be computed for Cohen's kappa"):
            self._result(kappa_max=NotComputed("no reason"))

    def test_kappa_max_the_marginals_do_not_give_refused(self):
        with pytest.raises(ValueError, match="kappa_max is 0.900000, but the raters' marginals allow 1.000000"):
            self._result(kappa_max=0.9)

    def test_kappa_above_its_ceiling_refused(self):
        with pytest.raises(ValueError, match="exceeds kappa_max"):
            self._result(kappa_max=0.6)

    def test_pe_from_other_data_refused(self):
        with pytest.raises(ValueError, match="describe different data"):
            self._result(parts=_fields(0.7, 0.5)["parts"])

    def test_weighted_pe_is_p_a_w_p_b(self):
        w = [[1.0, 0.5, 0.0], [0.5, 1.0, 0.5], [0.0, 0.5, 1.0]]
        p_a, p_b = _marginals(_PAIR)
        pe = float(p_a @ np.array(w) @ p_b)
        r = self._result(
            weights_type="linear", weight_matrix=w,
            parts=_fields(0.7, pe, weighted=True)["parts"],
        )
        assert r.parts.chance == pytest.approx(0.625)
        with pytest.raises(ValueError, match="describe different data"):
            self._result(
                weights_type="linear", weight_matrix=w,
                parts=_fields(0.7, 0.375, weighted=True)["parts"],
            )

    def test_undefined_kappa_has_undefined_ceiling(self):
        """Fewer than two shared subjects: kappa, its parts and kappa_max all NaN."""
        nan = float("nan")
        r = self._result(
            point=nan, kappa_max=nan,
            parts=AgreementParts(nan, nan, "agreement", False),
        )
        assert math.isnan(r.kappa_max)


# ============================================================================ #
#                          CongerKappaResult                                   #
# ============================================================================ #


class TestCongerKappaResult:
    """Conger's kappa: inherits IndividualMarginalResult, always unweighted."""

    def test_valid_construction(self):
        bp = _complete_bp(3)
        r = CongerKappaResult(0.6, None, None, bp, ci_method="percentile", **_conger(0.6))
        assert r.point == 0.6

    def test_is_individual_marginal_result(self):
        """Conger inherits from IndividualMarginalResult."""
        bp = _complete_bp(3)
        r = CongerKappaResult(0.6, None, None, bp, ci_method="bca", **_conger(0.6))
        assert isinstance(r, IndividualMarginalResult)

    def test_balance_profile_none_rejected(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            CongerKappaResult(0.6, None, None, None, ci_method="bca", **_conger(0.6))

    def test_n_raters_none_rejected(self):
        """Inherited from IndividualMarginalResult: n_raters=None rejected."""
        with pytest.raises(ValueError, match="n_raters must not be None"):
            CongerKappaResult(0.6, None, None, _bp(), ci_method="bca", **_conger(0.6))

    def test_parametric_null_not_allowed(self):
        bp = _complete_bp(3)
        with pytest.raises(ValueError, match="ci_method must be"):
            CongerKappaResult(0.6, None, None, bp, ci_method="parametric_null", **_conger(0.6))

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            CongerKappaResult(
                0.6, None, None, _complete_bp(3), ci_method=42,
                **_conger(0.6),
            )

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            CongerKappaResult(
                0.6, None, None, _complete_bp(3), ci_method="invalid",
                **_conger(0.6),
            )

    def test_str_shows_n_raters(self):
        bp = _complete_bp(3)
        s = str(CongerKappaResult(0.612, None, None, bp, ci_method="percentile", **_conger(0.612)))
        assert "Conger's κ" in s
        assert "3 raters, every rater rated every subject" in s

    def test_summary_keys(self):
        bp = _complete_bp(3)
        d = CongerKappaResult(0.6, None, None, bp, ci_method="bca", **_conger(0.6)).summary()
        assert d["ci_method"] == "bca"
        # No weights_type or is_exact
        assert "weights_type" not in d
        assert "is_exact" not in d

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[CongerKappaResult\] Error"):
            CongerKappaResult(0.6, None, None, None, ci_method="bca", **_conger(0.6))


class TestCongerKappaResultParts:
    """Pe from the raters' marginals; kappa_max for a pair only."""

    @pytest.mark.parametrize("fill_rate", [0.9, None])
    def test_incomplete_data_refused(self, fill_rate):
        bp = dataclasses.replace(_bp_with_n_raters(3), fill_rate=fill_rate)
        with pytest.raises(ValueError, match="Conger's kappa is for complete data"):
            CongerKappaResult(0.6, None, None, bp, ci_method="bca", **_conger(0.6))

    def test_valid_panel(self):
        r = CongerKappaResult(0.6, None, None, _complete_bp(3), ci_method="bca", **_conger(0.6))
        assert isinstance(r.kappa_max, NotComputed)

    def test_panel_kappa_max_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be NotComputed with 3 raters"):
            CongerKappaResult(
                0.6, None, None, _complete_bp(3), ci_method="bca",
                **dict(_conger(0.6), kappa_max=1.0),
            )

    def test_pair_gets_cohens_ceiling(self):
        """Two raters: Conger's Pe is Cohen's, and so is the ceiling."""
        r = CongerKappaResult(0.7, None, None, _complete_bp(2), ci_method="bca", **_cohen(0.7))
        assert r.kappa_max == 1.0

    def test_pair_without_kappa_max_refused(self):
        with pytest.raises(TypeError, match="kappa_max must be computed with 2 raters"):
            CongerKappaResult(
                0.7, None, None, _complete_bp(2), ci_method="bca",
                **dict(_cohen(0.7), kappa_max=NotComputed("no reason")),
            )

    def test_fleiss_pe_refused(self):
        """Pooled Pe without Conger's variance term is Fleiss', not Conger's."""
        m = _marginals(_PANEL)
        pooled = float((m.mean(axis=0) ** 2).sum())
        with pytest.raises(ValueError, match="Conger's chance term"):
            CongerKappaResult(
                0.6, None, None, _complete_bp(3), ci_method="bca",
                **dict(_conger(0.6), parts=_fields(0.6, pooled)["parts"]),
            )

    def test_weighted_parts_refused(self):
        with pytest.raises(ValueError, match="weighted Conger's kappa is not implemented"):
            CongerKappaResult(
                0.6, None, None, _complete_bp(3), ci_method="bca",
                **dict(_conger(0.6), parts=_fields(0.6, 0.3, weighted=True)["parts"]),
            )


# ============================================================================ #
#                           LightKappaResult                                   #
# ============================================================================ #


class TestLightKappaResult:
    """Light's kappa: adds n_pairs, no parametric_null."""

    def test_valid_construction(self):
        bp = _complete_bp(3)
        expected_pairs = bp.n_raters * (bp.n_raters - 1) // 2
        r = LightKappaResult(
            0.55, None, None, bp,
            n_pairs=expected_pairs, ci_method="percentile",
            **_light(0.55),
        )
        assert r.point == 0.55
        assert r.n_pairs == expected_pairs

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            LightKappaResult(
                0.55, None, None, None, n_pairs=3, ci_method="bca",
                **_light(0.55),
            )

    def test_n_raters_none_rejected(self):
        """Light's kappa requires n_raters (rater identity needed)."""
        with pytest.raises(ValueError, match="n_raters must not be None"):
            LightKappaResult(
                0.55, None, None, _bp(), n_pairs=3, ci_method="bca",
                **_light(0.55),
            )

    def test_n_pairs_must_match(self):
        """n_pairs must equal C(n_raters, 2)."""
        bp = _complete_bp(3)
        with pytest.raises(ValueError, match="n_pairs must be"):
            LightKappaResult(
                0.55, None, None, bp, n_pairs=999, ci_method="bca",
                **_light(0.55),
            )

    def test_n_pairs_not_int(self):
        with pytest.raises(TypeError, match="n_pairs must be int"):
            LightKappaResult(
                0.55, None, None, _complete_bp(3),
                n_pairs=3.0, ci_method="bca",
                **_light(0.55),
            )

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            LightKappaResult(
                0.55, None, None, _complete_bp(3),
                n_pairs=3, ci_method=42,
                **_light(0.55),
            )

    def test_parametric_null_not_allowed(self):
        bp = _complete_bp(3)
        n = bp.n_raters * (bp.n_raters - 1) // 2
        with pytest.raises(ValueError, match="ci_method must be"):
            LightKappaResult(
                0.55, None, None, bp,
                n_pairs=n, ci_method="parametric_null",
                **_light(0.55),
            )

    def test_str_shows_n_pairs(self):
        bp = _complete_bp(3)
        n = bp.n_raters * (bp.n_raters - 1) // 2
        s = str(LightKappaResult(
            0.587, None, None, bp, n_pairs=n, ci_method="percentile",
            **_light(0.587),
        ))
        assert "Light's κ" in s
        assert "mean of 3 pair kappas" in s
        assert "3 raters, every rater rated every subject" in s

    def test_summary_keys(self):
        bp = _complete_bp(3)
        n = bp.n_raters * (bp.n_raters - 1) // 2
        d = LightKappaResult(
            0.55, None, None, bp, n_pairs=n, ci_method="bca",
            **_light(0.55),
        ).summary()
        assert d["n_pairs"] == n
        assert d["ci_method"] == "bca"

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[LightKappaResult\] Error"):
            LightKappaResult(
                0.55, None, None, None, n_pairs=3, ci_method="bca",
                **_light(0.55),
            )


class TestLightKappaResultPairs:
    """Each pair's parts and ceiling, and the mean of the pair kappas is the point."""

    @pytest.mark.parametrize("fill_rate", [0.9, None])
    def test_incomplete_data_refused(self, fill_rate):
        bp = dataclasses.replace(_bp_with_n_raters(3), fill_rate=fill_rate)
        with pytest.raises(ValueError, match="Light's kappa is for complete data"):
            LightKappaResult(0.55, None, None, bp, n_pairs=3, ci_method="bca", **_light(0.55))

    PAIRS = list(itertools.combinations(_PANEL.rater_marginals, 2))

    def _result(self, point=0.55, **fields):
        kwargs = dict(_light(point), **fields)
        return LightKappaResult(
            point, None, None, _complete_bp(3), n_pairs=3, ci_method="bca", **kwargs,
        )

    def _pairs_at(self, *kappas):
        return {p: _fields(k, 0.5)["parts"] for p, k in zip(self.PAIRS, kappas)}

    def test_valid(self):
        r = self._result()
        assert list(r.pair_parts) == self.PAIRS

    def test_mean_of_pairs_is_the_point(self):
        r = self._result(point=0.5, pair_parts=self._pairs_at(0.3, 0.5, 0.7))
        assert r.point == 0.5

    def test_mean_that_is_not_the_point_refused(self):
        with pytest.raises(ValueError, match="mean of the pair kappas is 0.500000, but the point estimate is 0.550000"):
            self._result(pair_parts=self._pairs_at(0.3, 0.5, 0.7))

    def test_undefined_pair_left_out_of_the_mean(self):
        """As light_kappa does: the mean is over the pairs that have a kappa."""
        nan = float("nan")
        pairs = self._pairs_at(0.5, 0.6, 0.0)
        pairs[self.PAIRS[2]] = AgreementParts(nan, nan, "agreement", False)
        r = self._result(point=0.55, pair_parts=pairs)
        assert r.point == 0.55

    def test_every_pair_undefined_means_nan(self):
        nan = float("nan")
        pairs = {p: AgreementParts(nan, nan, "agreement", False) for p in self.PAIRS}
        r = self._result(point=nan, pair_parts=pairs, pair_kappa_max={p: nan for p in self.PAIRS})
        assert math.isnan(r.point)

    def test_missing_pair_refused(self):
        pairs = dict(list(_light(0.55)["pair_parts"].items())[:2])
        with pytest.raises(ValueError, match="pair_parts must have one entry for each pair"):
            self._result(pair_parts=pairs)

    def test_pairs_out_of_order_refused(self):
        pairs = dict(reversed(list(_light(0.55)["pair_parts"].items())))
        with pytest.raises(ValueError, match="in order"):
            self._result(pair_parts=pairs)

    def test_weighted_pair_parts_refused(self):
        pairs = {p: _fields(0.55, 0.5, weighted=True)["parts"] for p in self.PAIRS}
        with pytest.raises(TypeError, match="must be unweighted agreement parts"):
            self._result(pair_parts=pairs)

    def test_pair_above_its_ceiling_refused(self):
        ceilings = {p: 0.5 for p in self.PAIRS}
        with pytest.raises(ValueError, match="above its ceiling 0.500000"):
            self._result(pair_kappa_max=ceilings)

    @pytest.mark.parametrize("name", ["parts", "kappa_max"])
    def test_single_parts_or_ceiling_refused(self, name):
        value = _fields(0.55, 0.5)["parts"] if name == "parts" else 1.0
        with pytest.raises(TypeError, match=f"{name} must be NotComputed for Light's kappa"):
            self._result(**{name: value})

    def test_summary_lists_the_pairs(self):
        d = self._result().summary()
        assert [(e["rater_a"], e["rater_b"]) for e in d["pairs"]] == self.PAIRS
        assert d["pairs"][0]["kappa_max"] == 1.0
        assert d["parts"]["not_computed"] == "a mean of pairs"


class TestPrintedForm:
    """Three lines: the number, what it is made of, the data it was computed on."""

    def test_fleiss_three_lines(self):
        ci = _ci(0.5, 0.7)
        r = FleissKappaResult(
            0.6, ci, _draws_for(ci), _bp(), ci_method="percentile",
            weights_type="identity", **_pooled(0.6),
        )
        lines = str(r).splitlines()
        assert lines[0] == "Fleiss' κ = 0.600 [0.500, 0.700] (95% CI, percentile, n=1000)"
        assert lines[1] == f"  {r.parts}; κ_max not computed (pooled shares)"
        assert lines[2] == "  4 subjects; at least 3 raters; unweighted; exact formula"

    def test_counts_only_say_at_least(self):
        """Built from counts, the rater count is a lower bound."""
        r = FleissKappaResult(
            0.6, None, None, _bp(), ci_method="bca", weights_type="identity", **_pooled(0.6),
        )
        assert _bp().n_raters is None
        assert "at least 3 raters" in str(r)

    def test_cohen_ceiling_and_both_raters(self):
        r = CohenKappaResult(
            0.7, None, None, _bp_with_n_raters(2), ci_method="bca",
            weights_type="identity", **_cohen(0.7),
        )
        assert "κ_max = 1.000 (these raters' marginals)" in str(r)
        assert "4 subjects; unweighted" in str(r)

    def test_every_subject_is_used(self):
        """A subject the coefficient left out is an error, not a printed count."""
        fewer = AgreementProfile.from_nij_matrix(NijMatrix(BALANCED_NIJ[:3].copy()))
        with pytest.raises(
            ValueError, match="the data has 4 subjects, but the agreement profile uses 3",
        ):
            FleissKappaResult(
                0.6, None, None, _bp(), ci_method="bca", weights_type="identity",
                **dict(_pooled(0.6), agreement_profile=fewer),
            )

    def test_conger_panel_has_no_ceiling(self):
        r = CongerKappaResult(0.6, None, None, _complete_bp(3), ci_method="bca", **_conger(0.6))
        assert "κ_max not computed (3 or more raters)" in str(r)

    def test_light_spread_and_undefined_pairs(self):
        nan = float("nan")
        pairs = list(itertools.combinations(_PANEL.rater_marginals, 2))
        parts = {
            pairs[0]: _fields(0.4, 0.5)["parts"],
            pairs[1]: _fields(0.6, 0.5)["parts"],
            pairs[2]: AgreementParts(nan, nan, "agreement", False),
        }
        r = LightKappaResult(
            0.5, None, None, _complete_bp(3), n_pairs=3, ci_method="bca",
            **dict(_light(0.5), pair_parts=parts),
        )
        assert "mean of 2 of 3 pair kappas (the rest undefined), 0.40 to 0.60" in str(r)

    def test_light_equal_pairs_print_one_value(self):
        r = LightKappaResult(
            0.5, None, None, _complete_bp(3), n_pairs=3, ci_method="bca", **_light(0.5),
        )
        assert "mean of 3 pair kappas, all 0.50" in str(r)

    def test_profile_not_printed_with_the_result(self):
        """The tables print on their own: print(result.agreement_profile)."""
        r = FleissKappaResult(
            0.6, None, None, _bp(), ci_method="bca", weights_type="identity", **_pooled(0.6),
        )
        assert len(str(r).splitlines()) == 3
        assert "How often" not in str(r)


class TestPermutationPrintedForm:
    """Every pair or subset listed up to the limit; to_frame() holds them all."""

    def test_rows_and_frame_agree(self):
        r = PairwisePermutationResult(**_pairwise_kwargs())
        frame = r.to_frame()
        assert list(frame["raters"]) == ["0, 1", "0, 2", "1, 2"]
        assert list(frame["estimate"]) == [0.6, 0.4, 0.5]
        assert list(frame["subjects"]) == [SUBJECTS_PER_PAIR[p] for p in PAIRS]
        rows = [line for line in str(r).splitlines() if line.startswith("    ")][1:]
        for line, (_, row) in zip(rows, frame.iterrows()):
            assert line.split()[:3] == [*row["raters"].split(), f"{row['estimate']:.3f}"]
            assert f"[{row['ci_lower']:.3f}, {row['ci_upper']:.3f}]" in line

    def test_frame_round_trips_through_csv(self, tmp_path):
        import pandas as pd

        frame = SubsetPermutationResult(**_subset_kwargs()).to_frame()
        path = tmp_path / "subsets.csv"
        frame.to_csv(path, index=False)
        back = pd.read_csv(path)
        pd.testing.assert_frame_equal(back, frame, check_dtype=False)

    def test_at_the_limit_every_row_is_listed(self, monkeypatch):
        import interrater.contracts as contracts

        monkeypatch.setattr(contracts, "MAX_ROWS_PRINTED", 3)
        s = str(PairwisePermutationResult(**_pairwise_kwargs()))
        assert "too many to list" not in s
        assert all(f"    {p[0]}, {p[1]}" in s for p in PAIRS)

    def test_above_the_limit_points_to_to_frame(self, monkeypatch):
        import interrater.contracts as contracts

        monkeypatch.setattr(contracts, "MAX_ROWS_PRINTED", 2)
        s = str(PairwisePermutationResult(**_pairwise_kwargs()))
        assert '3 pairs, too many to list here: result.to_frame() has them all, e.g. result.to_frame().to_csv("pairs.csv")' in s
        assert "    0, 1" not in s

    def test_limit_lives_in_contracts(self):
        from interrater.contracts import MAX_ROWS_PRINTED

        assert MAX_ROWS_PRINTED == 100

    def test_warning_only_for_the_independent_mean(self):
        correlated = str(PairwisePermutationResult(**_pairwise_with_bootstrap()))
        independent = str(PairwisePermutationResult(**_pairwise_with_bootstrap(
            bootstrap_mean_method="independent",
        )))
        assert "!" not in correlated
        lines = independent.splitlines()
        assert lines[0].endswith("independent bootstrap")
        assert lines[1].startswith("  ! independent bootstrap: ignores the correlation")

    def test_fallback_and_failures_noted_on_the_row(self):
        cis = _cis_for(PAIRS)
        cis[(0, 2)] = _ci(0.3, 0.5, method_intended="bca", method_used="percentile", n_failed=12, n_resamples=1012)
        s = str(PairwisePermutationResult(**_pairwise_kwargs(bootstrap_cis_per_pair=cis)))
        row = next(line for line in s.splitlines() if line.startswith("    0, 2"))
        assert row.endswith("percentile [bca requested]; 12 failed")

    def test_short_intervals_said_once_under_the_list(self):
        cis = {p: _ci(v - 0.1, v + 0.1, n_resamples=100) for p, v in PAIRS.items()}
        s = str(PairwisePermutationResult(**_pairwise_kwargs(bootstrap_cis_per_pair=cis)))
        assert s.splitlines()[-1] == "  3 of 3 intervals rest on fewer usable draws than the 250 recommended"

    def test_undefined_pair_says_so(self):
        values = {**PAIRS, (1, 2): float("nan")}
        cis = {**_cis_for(PAIRS), (1, 2): _undefined_ci(method_intended="percentile")}
        r = PairwisePermutationResult(**_pairwise_kwargs(
            pairwise_values=values, bootstrap_cis_per_pair=cis, mean_over_pairs=0.5,
        ))
        row = next(line for line in str(r).splitlines() if line.startswith("    1, 2"))
        assert row.split()[2:4] == ["undefined", "[undefined]"]
        assert row.endswith("no method ran [percentile requested]; all 200 failed")
