"""
Tests for the analysis_results sub-package.

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

import numpy as np
import pytest

from interrater.intermediates.analysis_results import (
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
from interrater.intermediates.analysis_results.individual_marginal_result import (
    IndividualMarginalResult,
)
from interrater.intermediates.analysis_results.light_kappa_result import (
    LightKappaResult,
)
from interrater.intermediates.analysis_results.permutation_result import (
    PermutationResult,
)
from interrater.intermediates.analysis_results.pairwise_permutation_result import (
    PairwisePermutationResult,
)
from interrater.intermediates.analysis_results.subset_permutation_result import (
    SubsetPermutationResult,
)
from interrater.intermediates.analysis_results.multi_k_subset_result import (
    MultiKSubsetResult,
)
from interrater.intermediates.estimator_parts.nij_matrix import NijMatrix
from interrater.intermediates.data_profiles.balance_profile import BalanceProfile

from .conftest import BALANCED_NIJ, UNBALANCED_NIJ


# -- Helpers ----------------------------------------------------------------- #


def _bp() -> BalanceProfile:
    """Balanced profile: 4 subjects, 3 categories, 3 raters each."""
    return BalanceProfile.from_nij_matrix(NijMatrix(BALANCED_NIJ.copy()))


def _unbalanced_bp() -> BalanceProfile:
    """Unbalanced profile: 4 subjects, 2 categories, variable raters."""
    return BalanceProfile.from_nij_matrix(NijMatrix(UNBALANCED_NIJ.copy()))


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


# -- Pairwise test data: 3 raters → C(3,2) = 3 pairs ----------------------- #

PAIRS = {(0, 1): 0.6, (0, 2): 0.4, (1, 2): 0.5}
SUBJECTS_PER_PAIR = {(0, 1): 100, (0, 2): 80, (1, 2): 90}
PAIRWISE_MEAN = 0.5  # (0.6 + 0.4 + 0.5) / 3


def _pairwise_kwargs(**overrides):
    """Defaults for PairwisePermutationResult without bootstrap."""
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
        bootstrap_cis_per_pair=None,
        bootstrap_distributions_per_pair=None,
    )
    base.update(overrides)
    return base


def _pairwise_with_bootstrap(**overrides):
    """Defaults for PairwisePermutationResult with correlated bootstrap."""
    boot_dist = np.array([0.48, 0.50, 0.52])
    base = _pairwise_kwargs(
        bootstrap_mean_method="correlated",
        bootstrap_ci_of_mean=(0.45, 0.55),
        bootstrap_distribution_of_mean=boot_dist,
        bootstrap_cis_per_pair={
            (0, 1): (0.50, 0.70),
            (0, 2): (0.30, 0.50),
            (1, 2): (0.40, 0.60),
        },
        bootstrap_distributions_per_pair=None,
    )
    base.update(overrides)
    return base


# -- Subset test data: k=2 (same structure as pairs) and k=3 --------------- #

SUBSETS_K2 = {(0, 1): 0.6, (0, 2): 0.4, (1, 2): 0.5}
SUBJECTS_K2 = {(0, 1): 100, (0, 2): 80, (1, 2): 90}
SUBSET_K2_MEAN = 0.5

SUBSETS_K3 = {(0, 1, 2): 0.55}
SUBJECTS_K3 = {(0, 1, 2): 70}
SUBSET_K3_MEAN = 0.55


def _subset_kwargs(**overrides):
    """Defaults for SubsetPermutationResult (k=2) without bootstrap."""
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
        bootstrap_cis_per_subset=None,
        bootstrap_distributions_per_subset=None,
    )
    base.update(overrides)
    return base


def _subset_with_bootstrap(**overrides):
    """Defaults for SubsetPermutationResult (k=2) with correlated bootstrap."""
    boot_dist = np.array([0.48, 0.50, 0.52])
    base = _subset_kwargs(
        bootstrap_mean_method="correlated",
        bootstrap_ci_of_mean=(0.45, 0.55),
        bootstrap_distribution_of_mean=boot_dist,
        bootstrap_cis_per_subset={
            (0, 1): (0.50, 0.70),
            (0, 2): (0.30, 0.50),
            (1, 2): (0.40, 0.60),
        },
        bootstrap_distributions_per_subset=None,
    )
    base.update(overrides)
    return base


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
        bootstrap_cis_per_subset=None,
        bootstrap_distributions_per_subset=None,
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
        bootstrap_cis_per_subset=None,
        bootstrap_distributions_per_subset=None,
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
        )
        assert r.point == 0.75
        assert r.confidence_interval is None
        assert r.bootstrap_distribution is None
        assert r.balance_profile is None

    def test_with_ci(self):
        """Point estimate with confidence interval."""
        r = PointEstimateResult(
            point=0.5,
            confidence_interval=(0.3, 0.7),
            bootstrap_distribution=None,
            balance_profile=None,
        )
        assert r.confidence_interval == (0.3, 0.7)

    def test_with_bootstrap(self):
        """Point estimate with bootstrap distribution."""
        boot = np.array([0.4, 0.5, 0.6])
        r = PointEstimateResult(
            point=0.5,
            confidence_interval=None,
            bootstrap_distribution=boot,
            balance_profile=None,
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
        )
        assert r.balance_profile is bp

    def test_all_fields(self):
        """All optional fields populated."""
        boot = np.array([0.3, 0.4, 0.5, 0.6, 0.7])
        bp = _bp()
        r = PointEstimateResult(
            point=0.5,
            confidence_interval=(0.3, 0.7),
            bootstrap_distribution=boot,
            balance_profile=bp,
        )
        assert r.point == 0.5
        assert r.confidence_interval == (0.3, 0.7)
        assert len(r.bootstrap_distribution) == 5
        assert r.balance_profile is bp

    def test_is_analysis_result(self):
        """PointEstimateResult is a subclass of AnalysisResult."""
        r = PointEstimateResult(0.5, None, None, None)
        assert isinstance(r, AnalysisResult)

    def test_is_frozen(self):
        """Frozen dataclass: cannot reassign attributes."""
        r = PointEstimateResult(0.5, None, None, None)
        with pytest.raises(AttributeError):
            r.point = 0.9

    def test_int_point_accepted(self):
        """Integer point values are accepted (int is numeric)."""
        r = PointEstimateResult(1, None, None, None)
        assert r.point == 1


class TestPointEstimateResultValidationRejects:
    """Each constraint in __post_init__ tested with one failing input."""

    def test_point_not_numeric(self):
        """Point must be int or float, not str."""
        with pytest.raises(TypeError, match="point must be float"):
            PointEstimateResult("0.5", None, None, None)

    def test_ci_not_tuple(self):
        """Confidence interval must be a tuple, not list."""
        with pytest.raises(ValueError, match="confidence_interval must be"):
            PointEstimateResult(0.5, [0.3, 0.7], None, None)

    def test_ci_wrong_length(self):
        """Confidence interval must have exactly 2 elements."""
        with pytest.raises(ValueError, match="confidence_interval must be"):
            PointEstimateResult(0.5, (0.3, 0.5, 0.7), None, None)

    def test_ci_bounds_not_numeric(self):
        """CI bounds must be numeric."""
        with pytest.raises(TypeError, match="CI bounds must be float"):
            PointEstimateResult(0.5, ("lo", "hi"), None, None)

    def test_ci_lower_exceeds_upper(self):
        """Lower bound must not exceed upper bound."""
        with pytest.raises(ValueError, match="lower bound .* > upper bound"):
            PointEstimateResult(0.5, (0.7, 0.3), None, None)

    def test_point_below_ci(self):
        """Point must fall within CI."""
        with pytest.raises(ValueError, match="point estimate .* outside CI"):
            PointEstimateResult(0.2, (0.3, 0.7), None, None)

    def test_point_above_ci(self):
        """Point must fall within CI."""
        with pytest.raises(ValueError, match="point estimate .* outside CI"):
            PointEstimateResult(0.8, (0.3, 0.7), None, None)

    def test_bootstrap_not_ndarray(self):
        """Bootstrap distribution must be np.ndarray."""
        with pytest.raises(TypeError, match="bootstrap_distribution must be"):
            PointEstimateResult(0.5, None, [0.4, 0.5, 0.6], None)

    def test_bootstrap_wrong_ndim(self):
        """Bootstrap distribution must be 1-dimensional."""
        with pytest.raises(ValueError, match="must be 1-dimensional"):
            PointEstimateResult(0.5, None, np.array([[0.4, 0.5]]), None)

    def test_bootstrap_empty(self):
        """Bootstrap distribution must not be empty."""
        with pytest.raises(ValueError, match="must have at least 1 value"):
            PointEstimateResult(0.5, None, np.array([]), None)

    def test_balance_profile_wrong_type(self):
        """Balance profile must be BalanceProfile, not dict."""
        with pytest.raises(TypeError, match="balance_profile must be"):
            PointEstimateResult(0.5, None, None, {"n_raters": 3})

    def test_error_prefix(self):
        """Error messages start with [PointEstimateResult] Error."""
        with pytest.raises(TypeError, match=r"\[PointEstimateResult\] Error"):
            PointEstimateResult("bad", None, None, None)


class TestPointEstimateResultStr:
    """The compact display string under different field combinations."""

    def test_point_only(self):
        """Just the point value, no CI or bootstrap."""
        s = str(PointEstimateResult(0.523, None, None, None))
        assert s == "= 0.523"

    def test_with_ci(self):
        """Point + CI, no bootstrap."""
        s = str(PointEstimateResult(0.523, (0.481, 0.565), None, None))
        assert "95% CI" in s
        assert "0.481" in s
        assert "0.565" in s

    def test_with_ci_and_bootstrap(self):
        """Point + CI + bootstrap count."""
        boot = np.ones(1000)
        s = str(PointEstimateResult(1.0, (1.0, 1.0), boot, None))
        assert "n_bootstrap=1000" in s

    def test_bootstrap_without_ci(self):
        """Bootstrap count shown even without CI."""
        boot = np.array([0.5])
        s = str(PointEstimateResult(0.5, None, boot, None))
        assert "n_bootstrap=1" in s


class TestPointEstimateResultSummary:
    """The plain-types dict under different field combinations."""

    def test_point_only(self):
        """Minimal summary has just 'point'."""
        d = PointEstimateResult(0.75, None, None, None).summary()
        assert d == {"point": 0.75}

    def test_with_ci(self):
        """CI adds ci_lower and ci_upper."""
        d = PointEstimateResult(0.5, (0.3, 0.7), None, None).summary()
        assert d["point"] == 0.5
        assert d["ci_lower"] == 0.3
        assert d["ci_upper"] == 0.7

    def test_with_bootstrap(self):
        """Bootstrap adds n, mean, std."""
        boot = np.array([1.0, 2.0, 3.0])
        d = PointEstimateResult(2.0, (1.0, 3.0), boot, None).summary()
        assert d["bootstrap_n"] == 3
        assert abs(d["bootstrap_mean"] - 2.0) < 1e-12
        assert abs(d["bootstrap_std"] - np.std(boot)) < 1e-12

    def test_with_balance_profile(self):
        """Balance profile adds nested 'balance' dict."""
        bp = _bp()
        d = PointEstimateResult(0.5, None, None, bp).summary()
        assert "balance" in d
        assert isinstance(d["balance"], dict)

    def test_all_fields(self):
        """Full summary has all keys."""
        boot = np.array([0.4, 0.5, 0.6])
        bp = _bp()
        d = PointEstimateResult(0.5, (0.4, 0.6), boot, bp).summary()
        expected_keys = {
            "point", "ci_lower", "ci_upper",
            "bootstrap_n", "bootstrap_mean", "bootstrap_std",
            "balance",
        }
        assert set(d.keys()) == expected_keys


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
        )
        assert r.null_distribution is None
        assert r.z_statistic is None
        assert r.null_variance is None

    def test_is_frozen(self):
        """Frozen dataclass: cannot reassign attributes."""
        r = HypothesisTestResult(
            0.001, 0.05, True, "parametric_null", 0.65,
            z_statistic=3.42, null_variance=0.036,
        )
        with pytest.raises(AttributeError):
            r.p_value = 0.01

    def test_boundary_p_equals_alpha_reject(self):
        """At p == alpha, either decision is valid -- reject accepted."""
        r = HypothesisTestResult(
            p_value=0.05, alpha=0.05, reject_null=True,
            method="permutation", point=0.3, n_permutations=100,
        )
        assert r.reject_null is True

    def test_boundary_p_equals_alpha_no_reject(self):
        """At p == alpha, either decision is valid -- fail to reject accepted."""
        r = HypothesisTestResult(
            p_value=0.05, alpha=0.05, reject_null=False,
            method="permutation", point=0.3, n_permutations=100,
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
            )

    def test_p_value_negative(self):
        with pytest.raises(ValueError, match=r"p_value must be in \[0, 1\]"):
            HypothesisTestResult(
                p_value=-0.1, alpha=0.05, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
            )

    def test_p_value_above_one(self):
        with pytest.raises(ValueError, match=r"p_value must be in \[0, 1\]"):
            HypothesisTestResult(
                p_value=1.5, alpha=0.05, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
            )

    # --- alpha -------------------------------------------------------------- #

    def test_alpha_not_numeric(self):
        with pytest.raises(TypeError, match="alpha must be float"):
            HypothesisTestResult(
                p_value=0.5, alpha="0.05", reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
            )

    def test_alpha_zero(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.0, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
            )

    def test_alpha_one(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            HypothesisTestResult(
                p_value=0.5, alpha=1.0, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
            )

    # --- reject_null -------------------------------------------------------- #

    def test_reject_null_not_bool(self):
        with pytest.raises(TypeError, match="reject_null must be bool"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=0,
                method="permutation", point=0.1, n_permutations=100,
            )

    def test_reject_null_false_but_p_below_alpha(self):
        """p < alpha but reject_null=False is inconsistent."""
        with pytest.raises(ValueError, match="reject_null is False but"):
            HypothesisTestResult(
                p_value=0.01, alpha=0.05, reject_null=False,
                method="permutation", point=0.5, n_permutations=100,
            )

    def test_reject_null_true_but_p_above_alpha(self):
        """p > alpha but reject_null=True is inconsistent."""
        with pytest.raises(ValueError, match="reject_null is True but"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=True,
                method="permutation", point=0.1, n_permutations=100,
            )

    # --- method ------------------------------------------------------------- #

    def test_method_not_str(self):
        with pytest.raises(TypeError, match="method must be str"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=False,
                method=42, point=0.1,
            )

    def test_method_invalid(self):
        with pytest.raises(ValueError, match="method must be"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=False,
                method="custom", point=0.1,
            )

    # --- point -------------------------------------------------------------- #

    def test_point_not_numeric(self):
        with pytest.raises(TypeError, match="point must be float"):
            HypothesisTestResult(
                p_value=0.5, alpha=0.05, reject_null=False,
                method="permutation", point="0.1", n_permutations=100,
            )

    # --- Cross-validation: method-specific required fields ------------------ #

    def test_parametric_null_requires_z_statistic(self):
        with pytest.raises(ValueError, match="z_statistic must not be None"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                null_variance=0.036,
            )

    def test_parametric_null_requires_null_variance(self):
        with pytest.raises(ValueError, match="null_variance must not be None"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                z_statistic=3.42,
            )

    def test_permutation_requires_n_permutations(self):
        with pytest.raises(ValueError, match="n_permutations must not be None"):
            HypothesisTestResult(
                p_value=0.04, alpha=0.05, reject_null=True,
                method="permutation", point=0.5,
            )

    # --- Optional field type/value checks ----------------------------------- #

    def test_z_statistic_not_numeric(self):
        with pytest.raises(TypeError, match="z_statistic must be float"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                z_statistic="3.42", null_variance=0.036,
            )

    def test_null_variance_not_numeric(self):
        with pytest.raises(TypeError, match="null_variance must be float"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                z_statistic=3.42, null_variance="0.036",
            )

    def test_null_variance_negative(self):
        with pytest.raises(ValueError, match="null_variance must be >= 0"):
            HypothesisTestResult(
                p_value=0.001, alpha=0.05, reject_null=True,
                method="parametric_null", point=0.65,
                z_statistic=3.42, null_variance=-0.01,
            )

    def test_null_distribution_not_ndarray(self):
        with pytest.raises(TypeError, match="null_distribution must be"):
            HypothesisTestResult(
                p_value=0.04, alpha=0.05, reject_null=True,
                method="permutation", point=0.5,
                null_distribution=[0.1, 0.2], n_permutations=100,
            )

    def test_n_permutations_not_int(self):
        with pytest.raises(TypeError, match="n_permutations must be int"):
            HypothesisTestResult(
                p_value=0.04, alpha=0.05, reject_null=True,
                method="permutation", point=0.5,
                n_permutations=100.0,
            )

    def test_n_permutations_zero(self):
        with pytest.raises(ValueError, match="n_permutations must be >= 1"):
            HypothesisTestResult(
                p_value=0.04, alpha=0.05, reject_null=True,
                method="permutation", point=0.5,
                n_permutations=0,
            )

    # --- Error prefix ------------------------------------------------------- #

    def test_error_prefix(self):
        with pytest.raises(TypeError, match=r"\[HypothesisTestResult\] Error"):
            HypothesisTestResult(
                p_value="bad", alpha=0.05, reject_null=False,
                method="permutation", point=0.1, n_permutations=100,
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
        )
        s = str(r)
        assert "z = 3.42" in s
        assert "reject" in s
        assert "α = 0.05" in s
        assert "point estimate" in s

    def test_permutation_format(self):
        """Permutation: shows n_permutations."""
        r = HypothesisTestResult(
            p_value=0.003,
            alpha=0.05,
            reject_null=True,
            method="permutation",
            point=0.5,
            n_permutations=5000,
        )
        s = str(r)
        assert "5000 permutations" in s
        assert "reject" in s
        assert "point estimate" in s

    def test_very_small_p_value(self):
        """p < 0.000001 uses 'p < 0.000001' format."""
        r = HypothesisTestResult(
            1e-10, 0.05, True, "parametric_null", 0.9,
            z_statistic=8.0, null_variance=0.012,
        )
        s = str(r)
        assert "p < 0.000001" in s

    def test_normal_p_value_format(self):
        """Normal p-values show 6 decimal places."""
        r = HypothesisTestResult(
            0.042, 0.05, True, "permutation", 0.5, n_permutations=1000,
        )
        s = str(r)
        assert "p = 0.042000" in s


class TestHypothesisTestResultSummary:
    """summary() returns method-specific keys."""

    def test_required_keys_always_present(self):
        """Permutation with only required extras has base keys + n_permutations."""
        d = HypothesisTestResult(
            0.5, 0.05, False, "permutation", 0.1, n_permutations=100,
        ).summary()
        assert set(d.keys()) == {
            "p_value", "alpha", "reject_null", "method", "point",
            "n_permutations",
        }

    def test_parametric_keys(self):
        """z_statistic and null_variance included when present."""
        d = HypothesisTestResult(
            0.001, 0.05, True, "parametric_null", 0.65,
            z_statistic=3.42, null_variance=0.036,
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
        ).summary()
        assert d["n_permutations"] == 5000
        assert "null_distribution" not in d
        assert "z_statistic" not in d


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

    def test_no_bootstrap(self):
        r = PairwisePermutationResult(**_pairwise_kwargs())
        assert r.n_pairs == 3
        assert r.mean_over_pairs == PAIRWISE_MEAN
        assert r.bootstrap_mean_method is None

    def test_with_bootstrap(self):
        r = PairwisePermutationResult(**_pairwise_with_bootstrap())
        assert r.bootstrap_mean_method == "correlated"
        assert r.bootstrap_ci_of_mean == (0.45, 0.55)
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

    def test_bootstrap_ci_lower_exceeds_upper(self):
        with pytest.raises(ValueError, match="lower bound .* > upper bound"):
            PairwisePermutationResult(**_pairwise_with_bootstrap(
                bootstrap_ci_of_mean=(0.7, 0.3),
            ))

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
        dists = {
            (0, 1): np.array([0.5, 0.6]),
            (0, 2): np.array([0.3, 0.4]),
            (1, 2): np.array([0.4, 0.5]),
        }
        r = PairwisePermutationResult(**_pairwise_with_bootstrap(
            bootstrap_distributions_per_pair=dists,
        ))
        result = r.get_bootstrap_distributions_per_pair()
        assert result is dists

    def test_raises_when_no_bootstrap(self):
        """RuntimeError when bootstrap was not performed."""
        r = PairwisePermutationResult(**_pairwise_kwargs())
        with pytest.raises(RuntimeError, match="bootstrap was not performed"):
            r.get_bootstrap_distributions_per_pair()

    def test_raises_when_not_stored(self):
        """RuntimeError with threshold info when distributions not stored."""
        r = PairwisePermutationResult(**_pairwise_with_bootstrap())
        # bootstrap was performed but distributions_per_pair is None
        with pytest.raises(RuntimeError, match="MAX_PAIRS_FOR_DISTRIBUTIONS"):
            r.get_bootstrap_distributions_per_pair()


class TestPairwisePermutationResultStr:
    """__str__ output."""

    def test_no_bootstrap(self):
        s = str(PairwisePermutationResult(**_pairwise_kwargs()))
        assert "Pairwise (cohen_kappa)" in s
        assert "mean = 0.500" in s
        assert "n_pairs=3" in s

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
        assert "bootstrap_n" in d
        assert d["pair_distributions_stored"] is False

    def test_pair_distributions_stored_flag(self):
        """When bootstrap_distributions_per_pair is populated, flag is True."""
        dists = {
            (0, 1): np.array([0.5, 0.6]),
            (0, 2): np.array([0.3, 0.4]),
            (1, 2): np.array([0.4, 0.5]),
        }
        d = PairwisePermutationResult(**_pairwise_with_bootstrap(
            bootstrap_distributions_per_pair=dists,
        )).summary()
        assert d["pair_distributions_stored"] is True


# ============================================================================ #
#                      SubsetPermutationResult                                 #
# ============================================================================ #


class TestSubsetPermutationResultConstruction:
    """Valid construction with and without bootstrap."""

    def test_no_bootstrap(self):
        r = SubsetPermutationResult(**_subset_kwargs())
        assert r.subset_k == 2
        assert r.n_subsets == 3
        assert r.mean_over_subsets == SUBSET_K2_MEAN
        assert r.bootstrap_mean_method is None

    def test_with_bootstrap(self):
        r = SubsetPermutationResult(**_subset_with_bootstrap())
        assert r.bootstrap_mean_method == "correlated"
        assert r.bootstrap_ci_of_mean == (0.45, 0.55)
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

    def test_bootstrap_ci_lower_exceeds_upper(self):
        with pytest.raises(ValueError, match="lower bound .* > upper bound"):
            SubsetPermutationResult(**_subset_with_bootstrap(
                bootstrap_ci_of_mean=(0.7, 0.3),
            ))

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
        dists = {
            (0, 1): np.array([0.5, 0.6]),
            (0, 2): np.array([0.3, 0.4]),
            (1, 2): np.array([0.4, 0.5]),
        }
        r = SubsetPermutationResult(**_subset_with_bootstrap(
            bootstrap_distributions_per_subset=dists,
        ))
        result = r.get_bootstrap_distributions_per_subset()
        assert result is dists

    def test_raises_when_no_bootstrap(self):
        """RuntimeError when bootstrap was not performed."""
        r = SubsetPermutationResult(**_subset_kwargs())
        with pytest.raises(RuntimeError, match="bootstrap was not performed"):
            r.get_bootstrap_distributions_per_subset()

    def test_raises_when_not_stored(self):
        """RuntimeError with threshold info when distributions not stored."""
        r = SubsetPermutationResult(**_subset_with_bootstrap())
        with pytest.raises(
            RuntimeError, match="MAX_SUBSETS_FOR_DISTRIBUTIONS"
        ):
            r.get_bootstrap_distributions_per_subset()


class TestSubsetPermutationResultStr:
    """__str__ output."""

    def test_no_bootstrap(self):
        s = str(SubsetPermutationResult(**_subset_kwargs()))
        assert "Subset (k=2, fleiss_kappa)" in s
        assert "mean = 0.500" in s
        assert "n_subsets=3" in s

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
        assert "bootstrap_n" in d
        assert d["subset_distributions_stored"] is False

    def test_subset_distributions_stored_flag(self):
        """When distributions_per_subset is populated, flag is True."""
        dists = {
            (0, 1): np.array([0.5, 0.6]),
            (0, 2): np.array([0.3, 0.4]),
            (1, 2): np.array([0.4, 0.5]),
        }
        d = SubsetPermutationResult(**_subset_with_bootstrap(
            bootstrap_distributions_per_subset=dists,
        )).summary()
        assert d["subset_distributions_stored"] is True


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
            bootstrap_cis_per_subset=None,
            bootstrap_distributions_per_subset=None,
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
            bootstrap_cis_per_subset=None,
            bootstrap_distributions_per_subset=None,
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
            bootstrap_cis_per_subset=None,
            bootstrap_distributions_per_subset=None,
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
            bootstrap_cis_per_subset=None,
            bootstrap_distributions_per_subset=None,
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
        assert "MultiK Subset (fleiss_kappa, ks=[2, 3]):" in s
        assert "k=2: mean=0.500 (3 subsets)" in s
        assert "k=3: mean=0.550 (1 subsets)" in s

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
            bootstrap_ci_of_mean=(0.45, 0.55),
            bootstrap_distribution_of_mean=boot_dist,
            bootstrap_cis_per_subset={
                (0, 1): (0.50, 0.70),
                (0, 2): (0.30, 0.50),
                (1, 2): (0.40, 0.60),
            },
            bootstrap_distributions_per_subset=None,
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
            bootstrap_cis_per_subset=None,
            bootstrap_distributions_per_subset=None,
        )
        r = MultiKSubsetResult(
            coefficient_name="fleiss_kappa",
            ci_method="percentile",
            balance_profile=bp,
            per_k_results={2: result_k2, 3: result_k3},
            ks=[2, 3],
        )
        s = str(r)
        assert "95% CI: [0.450, 0.550]" in s
        assert "correlated" in s


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
        )
        assert r.point == 0.6
        assert r.ci_method == "bca"
        assert r.weights_type == "identity"

    def test_is_point_estimate_result(self):
        """FleissKappaResult inherits from PointEstimateResult."""
        r = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="bca", weights_type="identity",
        )
        assert isinstance(r, PointEstimateResult)

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            FleissKappaResult(
                0.6, None, None, None,
                ci_method="bca", weights_type="identity",
            )

    def test_ci_method_parametric_null_allowed(self):
        r = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="parametric_null", weights_type="identity",
        )
        assert r.ci_method == "parametric_null"

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            FleissKappaResult(
                0.6, None, None, _bp(),
                ci_method=42, weights_type="identity",
            )

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            FleissKappaResult(
                0.6, None, None, _bp(),
                ci_method="invalid", weights_type="identity",
            )

    @pytest.mark.parametrize("wt", ["identity", "linear", "quadratic", "custom"])
    def test_valid_weights_types(self, wt):
        r = FleissKappaResult(
            0.6, None, None, _bp(), ci_method="bca", weights_type=wt,
        )
        assert r.weights_type == wt

    def test_weights_type_not_str(self):
        with pytest.raises(TypeError, match="weights_type must be str"):
            FleissKappaResult(
                0.6, None, None, _bp(), ci_method="bca", weights_type=42,
            )

    def test_weights_type_invalid(self):
        with pytest.raises(ValueError, match="weights_type must be"):
            FleissKappaResult(
                0.6, None, None, _bp(),
                ci_method="bca", weights_type="ordinal",
            )

    def test_is_exact_balanced(self):
        r = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="percentile", weights_type="identity",
        )
        assert r.is_exact is True

    def test_is_exact_unbalanced(self):
        r = FleissKappaResult(
            0.6, None, None, _unbalanced_bp(),
            ci_method="percentile", weights_type="identity",
        )
        assert r.is_exact is False

    def test_str_balanced_unweighted(self):
        s = str(FleissKappaResult(
            0.523, None, None, _bp(),
            ci_method="bca", weights_type="identity",
        ))
        assert "Fleiss' κ" in s
        assert "unweighted" in s
        assert "exact formula" in s

    def test_str_unbalanced_weighted(self):
        s = str(FleissKappaResult(
            0.523, None, None, _unbalanced_bp(),
            ci_method="bca", weights_type="linear",
        ))
        assert "weighted (linear)" in s
        assert "approximation (unbalanced data)" in s

    def test_summary_keys(self):
        d = FleissKappaResult(
            0.6, None, None, _bp(),
            ci_method="bca", weights_type="identity",
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
            )

    def test_parent_error_chaining(self):
        """Parent TypeError wrapped with Fleiss prefix."""
        with pytest.raises(TypeError, match=r"\[FleissKappaResult\] Error"):
            FleissKappaResult(
                "bad", None, None, _bp(),
                ci_method="bca", weights_type="identity",
            )


# ============================================================================ #
#                            ScottsPiResult                                    #
# ============================================================================ #


class TestScottsPiResult:
    """Scott's pi: like Fleiss but no weights_type."""

    def test_valid_construction(self):
        r = ScottsPiResult(0.5, None, None, _bp(), ci_method="percentile")
        assert r.point == 0.5

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            ScottsPiResult(0.5, None, None, None, ci_method="bca")

    def test_ci_method_parametric_null_allowed(self):
        r = ScottsPiResult(0.5, None, None, _bp(), ci_method="parametric_null")
        assert r.ci_method == "parametric_null"

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            ScottsPiResult(0.5, None, None, _bp(), ci_method=42)

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            ScottsPiResult(0.5, None, None, _bp(), ci_method="wrong")

    def test_is_exact_balanced(self):
        assert ScottsPiResult(
            0.5, None, None, _bp(), ci_method="bca",
        ).is_exact is True

    def test_is_exact_unbalanced(self):
        assert ScottsPiResult(
            0.5, None, None, _unbalanced_bp(), ci_method="bca",
        ).is_exact is False

    def test_str_balanced(self):
        s = str(ScottsPiResult(0.5, None, None, _bp(), ci_method="bca"))
        assert "Scott's π" in s
        assert "exact formula" in s

    def test_str_unbalanced(self):
        s = str(ScottsPiResult(
            0.5, None, None, _unbalanced_bp(), ci_method="bca",
        ))
        assert "approximation (unbalanced data)" in s

    def test_summary_keys(self):
        d = ScottsPiResult(0.5, None, None, _bp(), ci_method="bca").summary()
        assert d["is_exact"] is True
        assert d["ci_method"] == "bca"
        assert "weights_type" not in d

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[ScottsPiResult\] Error"):
            ScottsPiResult(0.5, None, None, None, ci_method="bca")


# ============================================================================ #
#                            GwetsAC1Result                                    #
# ============================================================================ #


class TestGwetsAC1Result:
    """Gwet's AC1: no parametric_null (explicit 'not implemented' rejection)."""

    def test_valid_construction(self):
        r = GwetsAC1Result(0.7, None, None, _bp(), ci_method="bca")
        assert r.point == 0.7

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            GwetsAC1Result(0.7, None, None, None, ci_method="bca")

    def test_parametric_null_rejected(self):
        """Gwet's AC1 explicitly rejects parametric_null with 'not implemented'."""
        with pytest.raises(ValueError, match="parametric_null.*not implemented"):
            GwetsAC1Result(0.7, None, None, _bp(), ci_method="parametric_null")

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            GwetsAC1Result(0.7, None, None, _bp(), ci_method=42)

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            GwetsAC1Result(0.7, None, None, _bp(), ci_method="invalid")

    def test_is_exact_balanced(self):
        assert GwetsAC1Result(
            0.7, None, None, _bp(), ci_method="bca",
        ).is_exact is True

    def test_is_exact_unbalanced(self):
        assert GwetsAC1Result(
            0.7, None, None, _unbalanced_bp(), ci_method="bca",
        ).is_exact is False

    def test_str_balanced(self):
        s = str(GwetsAC1Result(0.7, None, None, _bp(), ci_method="bca"))
        assert "Gwet's AC1" in s
        assert "exact formula" in s

    def test_summary_keys(self):
        d = GwetsAC1Result(0.7, None, None, _bp(), ci_method="bca").summary()
        assert d["is_exact"] is True
        assert d["ci_method"] == "bca"
        assert "weights_type" not in d

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[GwetsAC1Result\] Error"):
            GwetsAC1Result(0.7, None, None, None, ci_method="bca")


# ============================================================================ #
#                         BrennanPredigerResult                                #
# ============================================================================ #


class TestBrennanPredigerResult:
    """Brennan-Prediger: percentile/bca only, has weights_type."""

    def test_valid_construction(self):
        r = BrennanPredigerResult(
            0.8, None, None, _bp(),
            ci_method="percentile", weights_type="identity",
        )
        assert r.point == 0.8
        assert r.weights_type == "identity"

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            BrennanPredigerResult(
                0.8, None, None, None,
                ci_method="bca", weights_type="identity",
            )

    def test_parametric_null_not_allowed(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            BrennanPredigerResult(
                0.8, None, None, _bp(),
                ci_method="parametric_null", weights_type="identity",
            )

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            BrennanPredigerResult(
                0.8, None, None, _bp(),
                ci_method=42, weights_type="identity",
            )

    @pytest.mark.parametrize("wt", ["identity", "linear", "quadratic", "custom"])
    def test_valid_weights_types(self, wt):
        r = BrennanPredigerResult(
            0.8, None, None, _bp(), ci_method="bca", weights_type=wt,
        )
        assert r.weights_type == wt

    def test_weights_type_not_str(self):
        with pytest.raises(TypeError, match="weights_type must be str"):
            BrennanPredigerResult(
                0.8, None, None, _bp(), ci_method="bca", weights_type=42,
            )

    def test_weights_type_invalid(self):
        with pytest.raises(ValueError, match="weights_type must be"):
            BrennanPredigerResult(
                0.8, None, None, _bp(),
                ci_method="bca", weights_type="ordinal",
            )

    def test_is_exact_balanced(self):
        assert BrennanPredigerResult(
            0.8, None, None, _bp(),
            ci_method="bca", weights_type="identity",
        ).is_exact is True

    def test_is_exact_unbalanced(self):
        assert BrennanPredigerResult(
            0.8, None, None, _unbalanced_bp(),
            ci_method="bca", weights_type="identity",
        ).is_exact is False

    def test_str_balanced_unweighted(self):
        s = str(BrennanPredigerResult(
            0.835, None, None, _bp(),
            ci_method="bca", weights_type="identity",
        ))
        assert "Brennan-Prediger" in s
        assert "unweighted" in s
        assert "exact formula" in s

    def test_summary_keys(self):
        d = BrennanPredigerResult(
            0.8, None, None, _bp(),
            ci_method="bca", weights_type="linear",
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
            )


# ============================================================================ #
#                        KrippendorffAlphaResult                               #
# ============================================================================ #


class TestKrippendorffAlphaResult:
    """Krippendorff's alpha: is_exact always True, has weights_type."""

    def test_valid_construction(self):
        r = KrippendorffAlphaResult(
            0.65, None, None, _bp(),
            ci_method="bca", weights_type="identity",
        )
        assert r.point == 0.65
        assert r.weights_type == "identity"

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            KrippendorffAlphaResult(
                0.65, None, None, None,
                ci_method="bca", weights_type="identity",
            )

    def test_parametric_null_not_allowed(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            KrippendorffAlphaResult(
                0.65, None, None, _bp(),
                ci_method="parametric_null", weights_type="identity",
            )

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            KrippendorffAlphaResult(
                0.65, None, None, _bp(),
                ci_method=42, weights_type="identity",
            )

    @pytest.mark.parametrize("wt", ["identity", "linear", "quadratic", "custom"])
    def test_valid_weights_types(self, wt):
        r = KrippendorffAlphaResult(
            0.65, None, None, _bp(), ci_method="bca", weights_type=wt,
        )
        assert r.weights_type == wt

    def test_weights_type_not_str(self):
        with pytest.raises(TypeError, match="weights_type must be str"):
            KrippendorffAlphaResult(
                0.65, None, None, _bp(), ci_method="bca", weights_type=42,
            )

    def test_weights_type_invalid(self):
        with pytest.raises(ValueError, match="weights_type must be"):
            KrippendorffAlphaResult(
                0.65, None, None, _bp(),
                ci_method="bca", weights_type="ordinal",
            )

    def test_is_exact_always_true_balanced(self):
        """Alpha is always exact, even when balanced."""
        assert KrippendorffAlphaResult(
            0.65, None, None, _bp(),
            ci_method="bca", weights_type="identity",
        ).is_exact is True

    def test_is_exact_always_true_unbalanced(self):
        """Alpha handles unbalanced data natively -- still exact."""
        assert KrippendorffAlphaResult(
            0.65, None, None, _unbalanced_bp(),
            ci_method="bca", weights_type="identity",
        ).is_exact is True

    def test_str_balanced(self):
        s = str(KrippendorffAlphaResult(
            0.625, None, None, _bp(),
            ci_method="bca", weights_type="identity",
        ))
        assert "Krippendorff's alpha" in s
        assert "exact formula" in s
        # Balanced: no "(unbalanced data)" suffix
        assert "unbalanced" not in s

    def test_str_unbalanced(self):
        """Unbalanced shows 'exact formula (unbalanced data)' -- NOT 'approximation'."""
        s = str(KrippendorffAlphaResult(
            0.625, None, None, _unbalanced_bp(),
            ci_method="bca", weights_type="identity",
        ))
        assert "exact formula (unbalanced data)" in s

    def test_summary_keys(self):
        d = KrippendorffAlphaResult(
            0.65, None, None, _bp(),
            ci_method="bca", weights_type="quadratic",
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
            )


# ============================================================================ #
#                       IndividualMarginalResult                               #
# ============================================================================ #


class TestIndividualMarginalResult:
    """Base for Cohen and Conger: balance_profile and n_raters must not be None."""

    def test_valid_construction(self):
        bp = _bp_with_n_raters(3)
        r = IndividualMarginalResult(0.6, None, None, bp)
        assert r.point == 0.6
        assert r.balance_profile is bp

    def test_is_point_estimate_result(self):
        r = IndividualMarginalResult(0.6, None, None, _bp_with_n_raters(3))
        assert isinstance(r, PointEstimateResult)

    def test_balance_profile_none_rejected(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            IndividualMarginalResult(0.6, None, None, None)

    def test_n_raters_none_rejected(self):
        """BalanceProfile from NijMatrix has n_raters=None -- rejected."""
        with pytest.raises(ValueError, match="n_raters must not be None"):
            IndividualMarginalResult(0.6, None, None, _bp())

    def test_error_prefix(self):
        with pytest.raises(
            ValueError, match=r"\[IndividualMarginalResult\] Error"
        ):
            IndividualMarginalResult(0.6, None, None, None)

    def test_parent_validation_propagated(self):
        """Point type error wraps with IndividualMarginalResult prefix."""
        with pytest.raises(
            TypeError, match=r"\[IndividualMarginalResult\] Error"
        ):
            IndividualMarginalResult("bad", None, None, _bp_with_n_raters(3))


# ============================================================================ #
#                           CohenKappaResult                                   #
# ============================================================================ #


class TestCohenKappaResult:
    """Cohen's kappa: inherits IndividualMarginalResult, has weights_type."""

    def test_valid_construction(self):
        bp = _bp_with_n_raters(3)
        r = CohenKappaResult(
            0.7, None, None, bp,
            ci_method="percentile", weights_type="identity",
        )
        assert r.point == 0.7
        assert r.weights_type == "identity"

    def test_is_individual_marginal_result(self):
        """Cohen inherits from IndividualMarginalResult (not PointEstimateResult directly)."""
        bp = _bp_with_n_raters(3)
        r = CohenKappaResult(
            0.7, None, None, bp,
            ci_method="bca", weights_type="identity",
        )
        assert isinstance(r, IndividualMarginalResult)
        assert isinstance(r, PointEstimateResult)

    def test_balance_profile_none_rejected(self):
        """Inherited from IndividualMarginalResult."""
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            CohenKappaResult(
                0.7, None, None, None,
                ci_method="bca", weights_type="identity",
            )

    def test_n_raters_none_rejected(self):
        """Inherited from IndividualMarginalResult: n_raters=None rejected."""
        with pytest.raises(ValueError, match="n_raters must not be None"):
            CohenKappaResult(
                0.7, None, None, _bp(),
                ci_method="bca", weights_type="identity",
            )

    def test_parametric_null_allowed(self):
        bp = _bp_with_n_raters(3)
        r = CohenKappaResult(
            0.7, None, None, bp,
            ci_method="parametric_null", weights_type="identity",
        )
        assert r.ci_method == "parametric_null"

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(3),
                ci_method=42, weights_type="identity",
            )

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(3),
                ci_method="invalid", weights_type="identity",
            )

    @pytest.mark.parametrize("wt", ["identity", "linear", "quadratic", "custom"])
    def test_valid_weights_types(self, wt):
        bp = _bp_with_n_raters(3)
        r = CohenKappaResult(
            0.7, None, None, bp, ci_method="bca", weights_type=wt,
        )
        assert r.weights_type == wt

    def test_weights_type_not_str(self):
        with pytest.raises(TypeError, match="weights_type must be str"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(3),
                ci_method="bca", weights_type=42,
            )

    def test_weights_type_invalid(self):
        with pytest.raises(ValueError, match="weights_type must be"):
            CohenKappaResult(
                0.7, None, None, _bp_with_n_raters(3),
                ci_method="bca", weights_type="ordinal",
            )

    def test_str_unweighted(self):
        bp = _bp_with_n_raters(3)
        s = str(CohenKappaResult(
            0.723, None, None, bp,
            ci_method="percentile", weights_type="identity",
        ))
        assert "Cohen's κ" in s
        assert "unweighted" in s

    def test_str_weighted(self):
        bp = _bp_with_n_raters(3)
        s = str(CohenKappaResult(
            0.723, None, None, bp,
            ci_method="bca", weights_type="quadratic",
        ))
        assert "weighted (quadratic)" in s

    def test_summary_keys(self):
        bp = _bp_with_n_raters(3)
        d = CohenKappaResult(
            0.7, None, None, bp,
            ci_method="bca", weights_type="linear",
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
            )


# ============================================================================ #
#                          CongerKappaResult                                   #
# ============================================================================ #


class TestCongerKappaResult:
    """Conger's kappa: inherits IndividualMarginalResult, always unweighted."""

    def test_valid_construction(self):
        bp = _bp_with_n_raters(3)
        r = CongerKappaResult(0.6, None, None, bp, ci_method="percentile")
        assert r.point == 0.6

    def test_is_individual_marginal_result(self):
        """Conger inherits from IndividualMarginalResult."""
        bp = _bp_with_n_raters(3)
        r = CongerKappaResult(0.6, None, None, bp, ci_method="bca")
        assert isinstance(r, IndividualMarginalResult)

    def test_balance_profile_none_rejected(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            CongerKappaResult(0.6, None, None, None, ci_method="bca")

    def test_n_raters_none_rejected(self):
        """Inherited from IndividualMarginalResult: n_raters=None rejected."""
        with pytest.raises(ValueError, match="n_raters must not be None"):
            CongerKappaResult(0.6, None, None, _bp(), ci_method="bca")

    def test_parametric_null_not_allowed(self):
        bp = _bp_with_n_raters(3)
        with pytest.raises(ValueError, match="ci_method must be"):
            CongerKappaResult(0.6, None, None, bp, ci_method="parametric_null")

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            CongerKappaResult(
                0.6, None, None, _bp_with_n_raters(3), ci_method=42,
            )

    def test_ci_method_invalid(self):
        with pytest.raises(ValueError, match="ci_method must be"):
            CongerKappaResult(
                0.6, None, None, _bp_with_n_raters(3), ci_method="invalid",
            )

    def test_str_shows_n_raters(self):
        bp = _bp_with_n_raters(3)
        s = str(CongerKappaResult(0.612, None, None, bp, ci_method="percentile"))
        assert "Conger's κ" in s
        assert "n_raters=3" in s

    def test_summary_keys(self):
        bp = _bp_with_n_raters(3)
        d = CongerKappaResult(0.6, None, None, bp, ci_method="bca").summary()
        assert d["ci_method"] == "bca"
        # No weights_type or is_exact
        assert "weights_type" not in d
        assert "is_exact" not in d

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[CongerKappaResult\] Error"):
            CongerKappaResult(0.6, None, None, None, ci_method="bca")


# ============================================================================ #
#                           LightKappaResult                                   #
# ============================================================================ #


class TestLightKappaResult:
    """Light's kappa: adds n_pairs, no parametric_null."""

    def test_valid_construction(self):
        bp = _bp_with_n_raters(3)
        expected_pairs = bp.n_raters * (bp.n_raters - 1) // 2
        r = LightKappaResult(
            0.55, None, None, bp,
            n_pairs=expected_pairs, ci_method="percentile",
        )
        assert r.point == 0.55
        assert r.n_pairs == expected_pairs

    def test_balance_profile_required(self):
        with pytest.raises(ValueError, match="balance_profile must not be None"):
            LightKappaResult(
                0.55, None, None, None, n_pairs=3, ci_method="bca",
            )

    def test_n_raters_none_rejected(self):
        """Light's kappa requires n_raters (rater identity needed)."""
        with pytest.raises(ValueError, match="n_raters must not be None"):
            LightKappaResult(
                0.55, None, None, _bp(), n_pairs=3, ci_method="bca",
            )

    def test_n_pairs_must_match(self):
        """n_pairs must equal C(n_raters, 2)."""
        bp = _bp_with_n_raters(3)
        with pytest.raises(ValueError, match="n_pairs must be"):
            LightKappaResult(
                0.55, None, None, bp, n_pairs=999, ci_method="bca",
            )

    def test_n_pairs_not_int(self):
        with pytest.raises(TypeError, match="n_pairs must be int"):
            LightKappaResult(
                0.55, None, None, _bp_with_n_raters(3),
                n_pairs=3.0, ci_method="bca",
            )

    def test_ci_method_not_str(self):
        with pytest.raises(TypeError, match="ci_method must be str"):
            LightKappaResult(
                0.55, None, None, _bp_with_n_raters(3),
                n_pairs=3, ci_method=42,
            )

    def test_parametric_null_not_allowed(self):
        bp = _bp_with_n_raters(3)
        n = bp.n_raters * (bp.n_raters - 1) // 2
        with pytest.raises(ValueError, match="ci_method must be"):
            LightKappaResult(
                0.55, None, None, bp,
                n_pairs=n, ci_method="parametric_null",
            )

    def test_str_shows_n_pairs(self):
        bp = _bp_with_n_raters(3)
        n = bp.n_raters * (bp.n_raters - 1) // 2
        s = str(LightKappaResult(
            0.587, None, None, bp, n_pairs=n, ci_method="percentile",
        ))
        assert "Light's κ" in s
        assert "n_pairs=3" in s
        assert "n_raters=3" in s

    def test_summary_keys(self):
        bp = _bp_with_n_raters(3)
        n = bp.n_raters * (bp.n_raters - 1) // 2
        d = LightKappaResult(
            0.55, None, None, bp, n_pairs=n, ci_method="bca",
        ).summary()
        assert d["n_pairs"] == n
        assert d["ci_method"] == "bca"

    def test_error_prefix(self):
        with pytest.raises(ValueError, match=r"\[LightKappaResult\] Error"):
            LightKappaResult(
                0.55, None, None, None, n_pairs=3, ci_method="bca",
            )
