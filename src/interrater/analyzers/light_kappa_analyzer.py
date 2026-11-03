"""Analyzer for Light's kappa coefficient."""

from __future__ import annotations

from math import comb

import numpy as np

from .individual_marginal_analyzer import IndividualMarginalAnalyzer, check_complete

from ..data_objects.ratings import Ratings
from ..intermediates.results import LightKappaResult, HypothesisTestResult
from ..intermediates.data_profiles.balance_profile import BalanceProfile
from .configs import LightKappaConfig

from ..intermediates.estimator_parts.estimator_input import EstimatorInput
from .agreement_functions.individual_marginals import (
    cohen_kappa_max,
    light_kappa,
    light_pair_parts,
    rater_marginals,
    shared_rows,
)
from ..intermediates.data_profiles import AgreementProfile
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.light_kappa_result import (
    KAPPA_MAX_NOT_COMPUTED,
    PARTS_NOT_COMPUTED,
)

_ERROR_PREFIX = "[LightKappaAnalyzer] Error"


class LightKappaAnalyzer(IndividualMarginalAnalyzer):
    """
    I am the analyzer for Light's kappa, the mean of all C(n, 2) pairwise
    Cohen's kappas.

    Light's kappa (Light 1971) keeps rater identity through Cohen's kappa:
    each pair uses its own marginals. The mean then compresses the C(n, 2)
    values into a single summary. For per-pair detail, spread, and per-rater
    diagnostics, use the pairwise permutation analyzer -- Light's kappa is
    the headline number, not the diagnostic.

    At exactly two raters (n_pairs=1), Light's kappa reduces to Cohen's
    unweighted kappa (there is only one pair).

    I validate that:
    - Config is LightKappaConfig

    I compute:
    - Point estimate of Light's kappa (always unweighted)
    - Bootstrap confidence interval (percentile or BCa -- no parametric null)
    - Return LightKappaResult with all metadata including n_pairs

    There is no NijMatrix path. Light's kappa reads rater identity -- who said
    what -- and counts discard that.

    There is no weighted form in this package. Light (1971) defined the coefficient as the
    mean of unweighted pairwise kappas.

    There is no parametric null variance. The C(n, 2) pairwise kappas are
    correlated (pairs sharing a rater are not independent), so confidence
    intervals and hypothesis tests use bootstrap and permutation methods
    respectively.

    Attributes
    ----------
    ratings : Ratings
        The ratings data (from Analyzer, always required).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    config : LightKappaConfig
        The analysis configuration for Light's kappa.
    estimator_input : EstimatorInput
        Built from ratings, carrying category_index, n_categories, and
        weights=None (from IndividualMarginalAnalyzer).

    References
    ----------
    Light, R. J. (1971). Measures of response agreement for qualitative
        data: Some generalizations and alternatives. Psychological
        Bulletin, 76(5), 365--377.
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    # Populated at parent class IndividualMarginalAnalyzer
    estimator_input: EstimatorInput
    # Checked at this level
    config: LightKappaConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(self, ratings: Ratings, config: LightKappaConfig) -> None:
        """
        Initialize Light's kappa analyzer: validate config.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze. Must contain at least two raters.
        config : LightKappaConfig
            The analysis configuration for Light's kappa.

        Raises
        ------
        TypeError
            If config is not LightKappaConfig.
        ValueError
            If ratings is None, or any rating is missing (after injection):
            Light (1971) defines the coefficient for complete data.
        """
        # Call parent __init__ (validates ratings not None, validates config is
        # IndividualMarginalConfig, computes BalanceProfile, builds
        # EstimatorInput with weights=None)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate config type
        self._validate_light_config()

        # Complete data only, checked after any injection
        try:
            check_complete(self.ratings.category_index, "Light's kappa")
        except ValueError as exc:
            raise ValueError(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: kappa = 0 using permutation.

        Only "permutation" is supported -- there is no closed-form null
        variance for Light's kappa.

        Returns
        -------
        HypothesisTestResult
            Test statistic, p-value, and metadata.

        Raises
        ------
        ValueError
            If config.hypothesis_test_config is None.
            If method is not "permutation".
        """
        # Guard: hypothesis_test_config must be set
        if self.config.hypothesis_test_config is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis_test_config is None. "
                f"Set hypothesis_test_config in your config to use "
                f"hypothesis testing."
            )

        ht_config = self.config.hypothesis_test_config
        n_categories = self.ratings.semantics.n_categories

        # Compute point estimate
        point = light_kappa(self.estimator_input)

        if ht_config.method == "permutation":
            # Permutation test: shuffle ratings to build null distribution
            from .resampling_functions.permutation_resampling import (
                within_rater_shuffle,
                pool_and_redistribute,
            )

            scheme_map = {
                "within_rater_shuffle": within_rater_shuffle,
                "pool_and_redistribute": pool_and_redistribute,
            }

            scheme_fn = scheme_map.get(ht_config.scheme)
            if scheme_fn is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: unknown permutation scheme "
                    f"'{ht_config.scheme}'. Expected one of: "
                    f"{set(scheme_map.keys())}."
                )

            # Build null distribution
            null_values = []

            for permuted_array in scheme_fn(
                self.ratings.category_index,
                n_permutations=ht_config.n_permutations,
                random_state=ht_config.random_state,
            ):
                ei = EstimatorInput(
                    category_index=permuted_array,
                    n_categories=n_categories,
                )
                kappa_perm = light_kappa(ei)
                null_values.append(kappa_perm)

            null_distribution = np.array(null_values)

            counts = self._permutation_p_value(
                null_distribution, point, ht_config.alternative
            )

            return HypothesisTestResult(
                alpha=ht_config.alpha,
                reject_null=self._reject_null(counts["p_value"], ht_config.alpha),
                method="permutation",
                coefficient="light_kappa",
                point=point,
                null_distribution=null_distribution,
                n_permutations=ht_config.n_permutations,
                **counts,
            )

        else:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown hypothesis test method "
                f"'{ht_config.method}'. Light's kappa has no closed-form null "
                f"variance, so only 'permutation' is available. "
                f"Got '{ht_config.method}'."
            )

    # Private helpers -------------------------------------------------------- #

    def _validate_light_config(self) -> None:
        """
        Validate that self.config is a LightKappaConfig.

        Raises
        ------
        TypeError
            If config is not LightKappaConfig.
        """
        if not isinstance(self.config, LightKappaConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be LightKappaConfig, "
                f"got {type(self.config).__name__}"
            )

    def _pairs(self) -> tuple[dict, dict]:
        """
        Each pair's Po and Pe, from the function light_kappa is built from,
        and each pair's Cohen ceiling on the subjects the pair shares (NaN
        where they share fewer than two). Keyed by the pair's rater names.
        """
        names = [str(r) for r in self.ratings.semantics.raters]
        pair_parts, pair_kappa_max = {}, {}
        for (i, j), (po, pe) in light_pair_parts(self.estimator_input).items():
            key = (names[i], names[j])
            pair_parts[key] = AgreementParts(
                observed=po, chance=pe, kind="agreement", weighted=False,
            )
            pair_data = self.estimator_input.select_raters(np.array([i, j]))
            rows = shared_rows(pair_data, 0, 1)
            if rows.size < 2:
                pair_kappa_max[key] = float("nan")
            else:
                marginals = rater_marginals(pair_data, rows)
                pair_kappa_max[key] = cohen_kappa_max(marginals[0], marginals[1])
        return pair_parts, pair_kappa_max

    def _estimate_from_array(self, category_index: np.ndarray) -> float:
        """
        Light's kappa from a ratings array.

        The only coefficient-specific step in my resampling. Everything that
        draws samples is written once on Analyzer in terms of this.

        I build my own EstimatorInput because the array I am handed is a
        resample, not the original data. Light's kappa is the mean of the
        C(m,2) pairwise Cohen's kappas, so every rater pair is re-formed from
        the column positions of whatever array I am given.
        """
        ei = EstimatorInput(
            category_index=category_index,
            n_categories=self.ratings.semantics.n_categories,
        )
        return light_kappa(ei)

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> LightKappaResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.

        Only percentile and bca intervals: there is no closed-form null
        variance for Light's kappa.
        """
        n_pairs = comb(self.balance_profile.n_raters, 2)

        pair_parts, pair_kappa_max = self._pairs()
        return LightKappaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            parts=NotComputed(PARTS_NOT_COMPUTED),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=AgreementProfile.from_ratings(self.ratings),
            n_pairs=n_pairs,
            ci_method=self.config.ci_method,
            pair_parts=pair_parts,
            pair_kappa_max=pair_kappa_max,
        )
