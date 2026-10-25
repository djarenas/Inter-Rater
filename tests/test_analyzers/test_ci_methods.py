"""Tests for ci_methods.py -- percentile_ci, bca_ci, _acceleration.

Deterministic tests only: no bootstrap loops. Tolerance is 1e-12 for
exact computations, relaxed where floating-point order matters.

Both public functions return a BootstrapCIResult rather than a tuple, so
these tests read ``.lower`` and ``.upper`` and also check what the result
says about itself -- which method was asked for, which one ran, and how
many draws failed. A test that only checks the bounds would pass just as
happily on an interval that silently fell back.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import yaml

from interrater.analyzers.ci_methods import (
    percentile_ci,
    bca_ci,
    _acceleration,
)
from interrater.intermediates.results.ci_results import BootstrapCIResult


# A jackknife with some spread in it, for the BCa tests that are not about
# acceleration itself. Any set of distinct values gives a finite ``a``.
_JACKKNIFE = np.linspace(0.40, 0.60, 30)


_FIXTURES = Path(__file__).resolve().parents[1] / "fixtures"


def _law_school_data() -> tuple[np.ndarray, np.ndarray]:
    """LSAT and GPA for Efron & Tibshirani's 15 law schools."""
    with open(_FIXTURES / "efron1986_lawschool_data.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    lsat = np.array([float(row["LSAT"]) for row in rows])
    gpa = np.array([float(row["GPA"]) for row in rows])
    return lsat, gpa


def _law_school_expected() -> dict:
    """The constants the paper publishes, with their citation."""
    with open(_FIXTURES / "efron1986_lawschool_expected.yaml") as handle:
        return yaml.safe_load(handle)


# ------------------------------------------------------------------ #
# Module-level _ERROR_PREFIX                                           #
# ------------------------------------------------------------------ #


class TestErrorPrefix:
    """_ERROR_PREFIX lives at module level, not inside functions."""

    def test_error_prefix_at_module_level(self):
        import interrater.analyzers.ci_methods as mod

        assert hasattr(mod, "_ERROR_PREFIX")
        assert mod._ERROR_PREFIX == "[ci_methods] Error"

    def test_error_prefix_appears_in_errors(self):
        with pytest.raises(ValueError, match=r"\[ci_methods\] Error"):
            percentile_ci(np.array([]), 0.05)


# ------------------------------------------------------------------ #
# percentile_ci                                                        #
# ------------------------------------------------------------------ #


class TestPercentileCi:
    """Tests for percentile_ci()."""

    def test_known_symmetric_distribution(self):
        """np.arange(1, 101): 2.5th percentile=3.475, 97.5th=97.525."""
        dist = np.arange(1, 101, dtype=float)
        ci = percentile_ci(dist, 0.05)
        assert ci.lower == pytest.approx(3.475, abs=1e-12)
        assert ci.upper == pytest.approx(97.525, abs=1e-12)

    def test_all_same_values(self):
        """Constant distribution -> CI collapses to the single value."""
        dist = np.full(1000, 0.75)
        ci = percentile_ci(dist, 0.05)
        assert ci.lower == pytest.approx(0.75, abs=1e-12)
        assert ci.upper == pytest.approx(0.75, abs=1e-12)

    def test_alpha_10_percent(self):
        """alpha=0.10 -> 5th and 95th percentile."""
        dist = np.arange(1, 101, dtype=float)
        ci = percentile_ci(dist, 0.10)
        assert ci.lower == pytest.approx(5.95, abs=1e-12)
        assert ci.upper == pytest.approx(95.05, abs=1e-12)

    def test_small_alpha_gives_wider_ci(self):
        """Smaller alpha -> wider CI."""
        dist = np.arange(1, 1001, dtype=float)
        ci_05 = percentile_ci(dist, 0.05)
        ci_01 = percentile_ci(dist, 0.01)
        assert ci_01.lower < ci_05.lower
        assert ci_01.upper > ci_05.upper

    def test_single_element(self):
        """Single element: CI collapses."""
        dist = np.array([0.5])
        ci = percentile_ci(dist, 0.05)
        assert ci.lower == pytest.approx(0.5, abs=1e-12)
        assert ci.upper == pytest.approx(0.5, abs=1e-12)

    def test_returns_bootstrap_ci_result(self):
        dist = np.arange(1, 101, dtype=float)
        ci = percentile_ci(dist, 0.05)
        assert isinstance(ci, BootstrapCIResult)
        assert isinstance(ci.lower, float)
        assert isinstance(ci.upper, float)
        assert ci.as_tuple() == (ci.lower, ci.upper)

    def test_reports_its_own_method(self):
        ci = percentile_ci(np.arange(1, 101, dtype=float), 0.05)
        assert ci.method_intended == "percentile"
        assert ci.method_used == "percentile"
        assert ci.degraded is False
        assert ci.acceleration is None
        assert ci.bias_correction is None

    def test_carries_alpha_and_level(self):
        ci = percentile_ci(np.arange(1, 101, dtype=float), 0.10)
        assert ci.alpha == pytest.approx(0.10)
        assert ci.confidence_level == pytest.approx(0.90)

    # --- Failed resamples ---

    def test_counts_failed_draws(self):
        """NaN draws are excluded from the bounds and counted."""
        dist = np.arange(1, 101, dtype=float)
        dist[:10] = np.nan
        ci = percentile_ci(dist, 0.05)
        assert ci.n_resamples == 100
        assert ci.n_failed == 10
        assert ci.n_usable == 90
        assert ci.is_defined

    def test_bounds_computed_from_usable_draws_only(self):
        """Dropping the NaNs must give the same answer as never having them."""
        clean = np.arange(11, 101, dtype=float)
        holed = np.arange(1, 101, dtype=float)
        holed[:10] = np.nan

        assert percentile_ci(holed, 0.05).lower == pytest.approx(
            percentile_ci(clean, 0.05).lower, abs=1e-12
        )
        assert percentile_ci(holed, 0.05).upper == pytest.approx(
            percentile_ci(clean, 0.05).upper, abs=1e-12
        )

    def test_infinite_draws_count_as_failed(self):
        dist = np.arange(1, 101, dtype=float)
        dist[0] = np.inf
        ci = percentile_ci(dist, 0.05)
        assert ci.n_failed == 1
        assert np.isfinite(ci.lower) and np.isfinite(ci.upper)

    def test_all_draws_failed_means_no_method_ran(self):
        """An all-NaN distribution is a result, not an error."""
        ci = percentile_ci(np.full(200, np.nan), 0.05)
        assert ci.method_intended == "percentile"
        assert ci.method_used is None
        assert ci.is_defined is False
        assert ci.n_resamples == 200
        assert ci.n_failed == 200
        assert ci.degraded is False  # an absence, not a substitution

    # --- Guards ---

    def test_empty_distribution_raises(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            percentile_ci(np.array([]), 0.05)

    def test_alpha_zero_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            percentile_ci(np.arange(100.0), 0.0)

    def test_alpha_one_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            percentile_ci(np.arange(100.0), 1.0)

    def test_alpha_negative_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            percentile_ci(np.arange(100.0), -0.1)

    def test_alpha_greater_than_one_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            percentile_ci(np.arange(100.0), 1.5)


# ------------------------------------------------------------------ #
# _acceleration                                                        #
# ------------------------------------------------------------------ #


class TestAcceleration:
    """Tests for _acceleration() -- the jackknife skewness of eq. (8.5)."""

    def test_symmetric_jackknife_gives_zero(self):
        """Symmetric leave-one-out values have no skewness to correct."""
        values = np.array([0.4, 0.45, 0.5, 0.55, 0.6])
        assert _acceleration(values) == pytest.approx(0.0, abs=1e-12)

    def test_matches_hand_computed_value(self):
        """a = sum(d^3) / (6 * [sum(d^2)]^1.5), d = mean - value."""
        values = np.array([0.1, 0.2, 0.9])
        deviations = values.mean() - values
        expected = float(
            np.sum(deviations ** 3) / (6.0 * np.sum(deviations ** 2) ** 1.5)
        )
        assert _acceleration(values) == pytest.approx(expected, abs=1e-12)

    def test_sign_follows_skewness(self):
        """Left-skewed values and their mirror image give opposite signs."""
        values = np.array([0.1, 0.2, 0.9])
        mirrored = 1.0 - values
        assert _acceleration(values) == pytest.approx(
            -_acceleration(mirrored), abs=1e-12
        )

    def test_constant_jackknife_returns_none(self):
        """No spread means no skewness -- not zero, undefined."""
        assert _acceleration(np.full(20, 0.6)) is None

    def test_single_value_returns_none(self):
        assert _acceleration(np.array([0.5])) is None

    def test_empty_returns_none(self):
        assert _acceleration(np.array([])) is None

    def test_nan_values_are_dropped_not_propagated(self):
        """A degenerate leave-one-out dataset must not poison the rest."""
        clean = np.array([0.1, 0.2, 0.9])
        holed = np.array([0.1, 0.2, 0.9, np.nan, np.nan])
        assert _acceleration(holed) == pytest.approx(
            _acceleration(clean), abs=1e-12
        )

    def test_all_nan_returns_none(self):
        assert _acceleration(np.full(10, np.nan)) is None


# ------------------------------------------------------------------ #
# bca_ci                                                               #
# ------------------------------------------------------------------ #


class TestBcaCi:
    """Tests for bca_ci()."""

    def test_known_symmetric_close_to_percentile(self):
        """Symmetric dist -> BCa ≈ percentile (bias & skewness ≈ 0)."""
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.6, scale=0.1, size=5000)
        point = 0.6

        bca = bca_ci(dist, point, 0.05, _JACKKNIFE)
        pct = percentile_ci(dist, 0.05)

        assert bca.lower == pytest.approx(pct.lower, abs=0.02)
        assert bca.upper == pytest.approx(pct.upper, abs=0.02)

    def test_skewed_distribution_shifts_ci(self):
        """Skewed dist -> BCa adjusts bounds relative to percentile."""
        rng = np.random.default_rng(42)
        dist = rng.exponential(scale=0.1, size=5000) + 0.3
        point = 0.3 + 0.1

        bca = bca_ci(dist, point, 0.05, _JACKKNIFE)
        pct = percentile_ci(dist, 0.05)

        assert bca.lower < bca.upper
        assert pct.lower < pct.upper
        assert not (
            bca.lower == pytest.approx(pct.lower, abs=1e-12)
            and bca.upper == pytest.approx(pct.upper, abs=1e-12)
        )

    def test_returns_bootstrap_ci_result(self):
        rng = np.random.default_rng(42)
        dist = rng.normal(size=200)
        ci = bca_ci(dist, 0.0, 0.05, _JACKKNIFE)
        assert isinstance(ci, BootstrapCIResult)
        assert isinstance(ci.lower, float)
        assert isinstance(ci.upper, float)

    def test_reports_diagnostics_when_it_runs(self):
        """A real BCa interval carries the two constants it was built from."""
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=1000)
        ci = bca_ci(dist, 0.5, 0.05, _JACKKNIFE)

        assert ci.method_intended == "bca"
        assert ci.method_used == "bca"
        assert ci.degraded is False
        assert ci.acceleration is not None
        assert ci.bias_correction is not None
        assert np.isfinite(ci.acceleration)
        assert np.isfinite(ci.bias_correction)

    def test_acceleration_matches_the_jackknife_it_was_given(self):
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=1000)
        jackknife = np.array([0.1, 0.2, 0.9, 0.45, 0.5])
        ci = bca_ci(dist, 0.5, 0.05, jackknife)
        assert ci.acceleration == pytest.approx(
            _acceleration(jackknife), abs=1e-12
        )

    def test_ties_at_the_point_count_half(self):
        """
        Draws centred on the point, with 40% tied at it: no bias to correct.

        Counting the ties as below (the old `<=`) read this as 70% of draws
        below the point and shifted the interval up (z0 = +0.52). Half
        below and half above, as scipy.stats.bootstrap does, gives 50%.
        """
        dist = np.repeat([0.4, 0.5, 0.6], [300, 400, 300])
        ci = bca_ci(dist, 0.5, 0.05, _JACKKNIFE)
        assert ci.method_used == "bca"
        assert ci.bias_correction == pytest.approx(0.0, abs=1e-12)

    def test_ties_at_the_point_hand_worked(self):
        """100 below, 200 tied, 700 above: (100 + 100) / 1000 = 20% below."""
        from scipy.stats import norm

        dist = np.repeat([0.4, 0.5, 0.6], [100, 200, 700])
        ci = bca_ci(dist, 0.5, 0.05, _JACKKNIFE)
        assert ci.bias_correction == pytest.approx(norm.ppf(0.2), abs=1e-12)

    def test_extreme_inputs_do_not_raise(self):
        """A point far from the distribution mass must not crash BCa."""
        dist = np.full(100, 0.99)
        dist[0] = 0.98
        ci = bca_ci(dist, 0.01, 0.05, _JACKKNIFE)
        assert isinstance(ci.lower, float)
        assert isinstance(ci.upper, float)

    # --- Falling back ---

    def test_falls_back_when_acceleration_undefined(self):
        """A constant jackknife gives no acceleration, so BCa steps down."""
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=1000)
        ci = bca_ci(dist, 0.5, 0.05, np.full(20, 0.6))

        assert ci.method_intended == "bca"
        assert ci.method_used == "percentile"
        assert ci.degraded is True

    def test_fallback_discards_both_corrections(self):
        """'percentile' means percentile -- z0 is dropped along with a."""
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=1000)
        fallen = bca_ci(dist, 0.9, 0.05, np.full(20, 0.6))
        plain = percentile_ci(dist, 0.05)

        assert fallen.acceleration is None
        assert fallen.bias_correction is None
        assert fallen.lower == pytest.approx(plain.lower, abs=1e-12)
        assert fallen.upper == pytest.approx(plain.upper, abs=1e-12)

    def test_all_draws_failed_means_no_method_ran(self):
        ci = bca_ci(np.full(200, np.nan), 0.5, 0.05, _JACKKNIFE)
        assert ci.method_intended == "bca"
        assert ci.method_used is None
        assert ci.is_defined is False
        assert ci.n_failed == 200
        assert ci.acceleration is None

    def test_counts_failed_draws(self):
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=1000)
        dist[:25] = np.nan
        ci = bca_ci(dist, 0.5, 0.05, _JACKKNIFE)
        assert ci.n_resamples == 1000
        assert ci.n_failed == 25
        assert ci.n_usable == 975

    # --- Guards ---

    def test_missing_jackknife_raises(self):
        """Asking for BCa without the means to compute it is a caller bug."""
        rng = np.random.default_rng(42)
        dist = rng.normal(size=200)
        with pytest.raises(ValueError, match="requires jackknife_estimates"):
            bca_ci(dist, 0.0, 0.05, None)

    def test_empty_jackknife_raises(self):
        rng = np.random.default_rng(42)
        dist = rng.normal(size=200)
        with pytest.raises(ValueError, match="requires jackknife_estimates"):
            bca_ci(dist, 0.0, 0.05, np.array([]))

    def test_empty_distribution_raises(self):
        with pytest.raises(ValueError, match="cannot be empty"):
            bca_ci(np.array([]), 0.5, 0.05, _JACKKNIFE)

    def test_alpha_zero_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            bca_ci(np.arange(100.0), 0.5, 0.0, _JACKKNIFE)

    def test_alpha_one_raises(self):
        with pytest.raises(ValueError, match=r"alpha must be in \(0, 1\)"):
            bca_ci(np.arange(100.0), 0.5, 1.0, _JACKKNIFE)


# ------------------------------------------------------------------ #
# Validation against Efron & Tibshirani (1986), Example 3              #
# ------------------------------------------------------------------ #


class TestEfronLawSchool:
    """
    Pin the acceleration to the published literature.

    Efron & Tibshirani (1986) Example 3, p. 71, works the law school data:
    15 schools, correlation between LSAT score and GPA, theta_hat = .776.
    The paper reports a = -.0817, z0 = -.0927, and a central 90% interval
    of [.43, .92] from B = 100,000 draws.

    The paper also prints the fifteen empirical influence values its
    acceleration was computed from, which makes the first test below
    exact: no bootstrap, no RNG, no tolerance beyond the paper's own four
    decimal places. It tests the formula separately from its input, which
    matters because the two differ -- Efron computes influence values from
    the infinitesimal jackknife of eq. (10.11), and this package uses the
    delete-one jackknife of the 1993 book.

    Data and published values live in ``tests/fixtures/`` rather than
    inline, so that the numbers and their citation sit in one place. They
    are not in ``example_data/``: that holds ratings datasets a user can
    load and run an analyzer on, and this is 15 pairs of test scores with
    no raters, no categories, and no agreement coefficient. It validates
    ci_methods directly, underneath the analyzer layer.
    """

    @staticmethod
    def _correlation(x: np.ndarray, y: np.ndarray) -> float:
        return float(np.corrcoef(x, y)[0, 1])

    def test_formula_reproduces_published_acceleration(self):
        """Given the paper's own influence values, eq. (8.5) returns its a.

        Exact to the four decimals the paper prints. This is the test that
        says our arithmetic is the published arithmetic; everything else
        about acceleration is a question of what we feed it.

        The influence values cannot be passed in as they stand. They are
        already deviations, and _acceleration takes *estimates* and forms
        its own deviations as ``mean - value`` -- hand it U_i directly and
        every deviation comes back negated, which cubes to the right
        magnitude with the wrong sign. The two are related by

            U_i = (n - 1) * (theta_(.) - theta_(i))

        so inverting that gives leave-one-out estimates whose deviations
        are the published influence values again. The centre is arbitrary,
        since the function re-centres on the mean; theta_hat is used for
        readability.
        """
        expected = _law_school_expected()
        influence = np.array(expected["empirical_influence_values"])
        n = len(influence)
        pseudo_estimates = expected["point"] - influence / (n - 1)

        assert _acceleration(pseudo_estimates) == pytest.approx(
            expected["acceleration"], abs=5e-5
        )

    def test_influence_values_are_not_estimates(self):
        """Passing influence values straight in negates the result.

        Pinned deliberately. The sign of ``a`` leans the interval, so a
        future refactor that quietly changes which direction _acceleration
        subtracts would otherwise pass every other test in this file.
        """
        expected = _law_school_expected()
        influence = np.array(expected["empirical_influence_values"])
        n = len(influence)
        pseudo_estimates = expected["point"] - influence / (n - 1)

        assert _acceleration(influence) == pytest.approx(
            -_acceleration(pseudo_estimates), abs=1e-12
        )

    def test_point_estimate_matches_published(self):
        """theta_hat = .776 (p. 56). A wrong dataset fails here first."""
        lsat, gpa = _law_school_data()
        expected = _law_school_expected()
        assert self._correlation(lsat, gpa) == pytest.approx(
            expected["point"], abs=5e-4
        )

    def test_delete_one_jackknife_lands_near_published_acceleration(self):
        """Our jackknife is not Efron's, and the gap is small and stable.

        Efron's -.0817 comes from the infinitesimal jackknife; the
        delete-one jackknife this package uses gives -.0757. Both are
        standard. This pins the difference so that a future change to
        subject_jackknife cannot widen it unnoticed.
        """
        lsat, gpa = _law_school_data()
        expected = _law_school_expected()
        jackknife = np.array([
            self._correlation(np.delete(lsat, i), np.delete(gpa, i))
            for i in range(len(lsat))
        ])
        acceleration = _acceleration(jackknife)

        assert acceleration == pytest.approx(-0.0757, abs=5e-4)
        assert acceleration == pytest.approx(
            expected["acceleration"], abs=0.01
        )
        # Same sign as the paper -- the thing we got wrong before 2.0.1.
        assert acceleration < 0

    def test_bca_interval_reproduces_published_bounds(self):
        """The published interval is CENTRAL 90%, so alpha = 0.10.

        Comparing a 95% interval against [.43, .92] looks like a failure
        and is not one. B here is far below the paper's 100,000, which it
        calls "about 100 times more than actually necessary", so the
        tolerance is Monte Carlo slack rather than a claim about accuracy.
        """
        lsat, gpa = _law_school_data()
        expected = _law_school_expected()
        n = len(lsat)
        point = self._correlation(lsat, gpa)

        rng = np.random.default_rng(1986)
        indices = rng.integers(0, n, size=(10_000, n))
        draws = np.array([
            self._correlation(lsat[i], gpa[i]) for i in indices
        ])
        jackknife = np.array([
            self._correlation(np.delete(lsat, i), np.delete(gpa, i))
            for i in range(n)
        ])

        ci = bca_ci(draws, point, expected["alpha"], jackknife)
        published_lower, published_upper = expected["ci_bca"]

        assert ci.method_used == "bca"
        assert ci.lower == pytest.approx(published_lower, abs=0.03)
        assert ci.upper == pytest.approx(published_upper, abs=0.03)
        assert ci.bias_correction == pytest.approx(
            expected["bias_correction_z0"], abs=0.02
        )


# ------------------------------------------------------------------ #
# Cross-method consistency                                             #
# ------------------------------------------------------------------ #


class TestCrossMethodConsistency:
    """Sanity checks across both CI methods."""

    def test_both_methods_return_lower_le_upper(self):
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=1000)

        pct = percentile_ci(dist, 0.05)
        bca = bca_ci(dist, 0.5, 0.05, _JACKKNIFE)

        assert pct.lower <= pct.upper
        assert bca.lower <= bca.upper

    def test_wider_alpha_gives_narrower_ci_percentile(self):
        dist = np.arange(1, 1001, dtype=float)
        ci_05 = percentile_ci(dist, 0.05)
        ci_10 = percentile_ci(dist, 0.10)
        assert (ci_05.upper - ci_05.lower) > (ci_10.upper - ci_10.lower)

    def test_wider_alpha_gives_narrower_ci_bca(self):
        rng = np.random.default_rng(42)
        dist = rng.normal(loc=0.5, scale=0.1, size=2000)
        ci_05 = bca_ci(dist, 0.5, 0.05, _JACKKNIFE)
        ci_10 = bca_ci(dist, 0.5, 0.10, _JACKKNIFE)
        assert (ci_05.upper - ci_05.lower) > (ci_10.upper - ci_10.lower)

    def test_percentile_never_falls_back(self):
        """Only BCa can degrade; percentile either runs or nothing does."""
        rng = np.random.default_rng(42)
        for dist in (
            rng.normal(loc=0.5, scale=0.1, size=500),
            np.full(500, 0.75),
            np.concatenate([np.full(400, np.nan), np.full(100, 0.5)]),
        ):
            ci = percentile_ci(dist, 0.05)
            assert ci.method_used in ("percentile", None)
            assert ci.degraded is False
