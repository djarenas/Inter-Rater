"""
Integration tests for pairwise permutation analysis.

Pairwise permutation analysis on SR1 prelim (4 raters, balanced) with
both bootstrap methods to verify the CI width relationship.

With 4 raters, there are C(4, 2) = 6 rater pairs.  The mean of all
pairwise Cohen kappas is Light's kappa (0.496, validated in
test_research_data.py).

Two bootstrap methods:
  - correlated:   resample subjects globally, recompute all pairs on the
                  same resampled data, then average → captures shared-subject
                  dependence → wider (correct) CI
  - independent:  sample from per-pair bootstrap distributions independently
                  → ignores cross-pair correlation → narrower (anti-conservative) CI

Point estimates (mean_over_pairs) are deterministic and asserted to 3
decimal places (abs=5e-4).  CI width comparison is structural (no specific
values asserted -- just correlated > independent).
"""

import pytest

from interrater import PairwisePermutationAnalysis
from interrater.analyzers.configs import (
    PairwisePermutationConfig,
    CohenKappaConfig,
    PermutationBootstrapMeanConfig,
)

# Tolerance: 3 decimal places.
TOL = 5e-4


# ── SR1 prelim pairwise (C(4,2) = 6 pairs) ──────────────────────────── #
#
# Source: Systematic Review 1, preliminary screening (199 subjects).
# Mean pairwise Cohen kappa = Light's kappa = 0.496 (V2 validated).


class TestSR1PrelimPairwise:
    """Pairwise permutation analysis on SR1 prelim (199 subjects, 4 raters, K=3)."""

    def test_cohen_correlated(self, sr1_prelim_ratings):
        """Mean pairwise Cohen = 0.496, correlated bootstrap for correct CI.

        Correlated bootstrap resamples subjects globally so all 6 pair
        statistics are computed on the same resampled subjects.  This
        preserves the cross-pair dependence and produces a correctly-
        sized confidence interval.
        """
        config = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(
                n_resamples=50, random_state=42,
            ),
            bootstrap_mean_config=PermutationBootstrapMeanConfig(
                method="correlated",
                n_bootstraps=50,
                random_state=99,
            ),
        )
        result = PairwisePermutationAnalysis(sr1_prelim_ratings, config).run()

        # Point estimate (deterministic)
        assert result.mean_over_pairs == pytest.approx(0.496, abs=TOL)

        # CI structure
        ci = result.bootstrap_ci_of_mean
        assert ci is not None
        assert ci[0] < result.mean_over_pairs < ci[1]

        # Bootstrap method recorded correctly
        assert result.bootstrap_mean_method == "correlated"

    def test_cohen_independent(self, sr1_prelim_ratings):
        """Mean pairwise Cohen = 0.496, independent bootstrap (narrower CI).

        Independent bootstrap samples from each pair's bootstrap distribution
        independently, ignoring the shared-subject dependence.  The mean is
        the same but the CI is narrower (anti-conservative).
        """
        config = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(
                n_resamples=50, random_state=42,
            ),
            bootstrap_mean_config=PermutationBootstrapMeanConfig(
                method="independent",
                n_bootstraps=50,
                random_state=99,
            ),
        )
        result = PairwisePermutationAnalysis(sr1_prelim_ratings, config).run()

        # Point estimate (deterministic, same as correlated)
        assert result.mean_over_pairs == pytest.approx(0.496, abs=TOL)

        # CI structure
        ci = result.bootstrap_ci_of_mean
        assert ci is not None
        assert ci[0] < result.mean_over_pairs < ci[1]

        # Bootstrap method recorded correctly
        assert result.bootstrap_mean_method == "independent"

    def test_correlated_wider_than_independent(self, sr1_prelim_ratings):
        """Correlated bootstrap CI is wider than independent.

        This is a structural property: the correlated method preserves
        cross-pair dependence from shared subjects, which inflates
        variance of the mean.  The independent method treats each pair's
        distribution as if it were drawn from a separate experiment.

        With 50 bootstraps the ratio is approximately 1.8x.
        """
        cohen_config = CohenKappaConfig(n_resamples=50, random_state=42)

        config_corr = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=cohen_config,
            bootstrap_mean_config=PermutationBootstrapMeanConfig(
                method="correlated",
                n_bootstraps=50,
                random_state=99,
            ),
        )
        result_corr = PairwisePermutationAnalysis(
            sr1_prelim_ratings, config_corr
        ).run()

        config_ind = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=cohen_config,
            bootstrap_mean_config=PermutationBootstrapMeanConfig(
                method="independent",
                n_bootstraps=50,
                random_state=99,
            ),
        )
        result_ind = PairwisePermutationAnalysis(
            sr1_prelim_ratings, config_ind
        ).run()

        ci_corr = result_corr.bootstrap_ci_of_mean
        ci_ind = result_ind.bootstrap_ci_of_mean

        width_corr = ci_corr[1] - ci_corr[0]
        width_ind = ci_ind[1] - ci_ind[0]

        assert width_corr > width_ind

    def test_n_pairs(self, sr1_prelim_ratings):
        """Result reports n_pairs = C(4, 2) = 6."""
        config = PairwisePermutationConfig(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(
                n_resamples=50, random_state=42,
            ),
        )
        result = PairwisePermutationAnalysis(sr1_prelim_ratings, config).run()

        assert result.n_pairs == 6
        assert len(result.pairwise_values) == 6
