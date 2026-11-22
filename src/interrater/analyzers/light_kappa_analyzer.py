"""Analyzer for Light's kappa coefficient."""

from __future__ import annotations

from math import comb

import numpy as np

from .individual_marginal_analyzer import IndividualMarginalAnalyzer, check_complete

from ..data_objects.ratings import Ratings
from ..intermediates.results import LightKappaResult, HypothesisTestResult
from ..intermediates.data_profiles.balance_profile import BalanceProfile
from ..intermediates.data_profiles.rater_profile import RaterProfile
from .configs import LightKappaConfig

from ..intermediates.estimator_parts.estimator_input import EstimatorInput
from .agreement_functions.individual_marginals import (
    cohen_kappa_max,
    light_kappa,
    light_pair_parts,
    rater_marginals,
    shared_rows,
)
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

    Light's kappa reads rater identity -- who said what -- and counts discard
    that.

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
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped. Always present here:
        this coefficient needs rater identity.
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

    # Label in hypothesis test results
    _COEFFICIENT_KEY = "light_kappa"

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
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
            If any rating is missing (after injection):
            Light (1971) defines the coefficient for complete data.
        """
        # Call parent __init__ (validates ratings and checks that config is
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
        Test H_0: kappa = 0 by permutation (``_permutation_test``).

        There is no closed-form null variance for Light's kappa.
        """
        return self._permutation_test()

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

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Light's kappa from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples is written once on Analyzer in terms of this. Light's
        kappa is the mean of the C(m,2) pairwise Cohen's kappas, so every
        rater pair is re-formed from the column positions of data.
        """
        return light_kappa(data)

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> LightKappaResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.

        Only percentile and bca intervals: there is no closed-form null
        variance for Light's kappa.
        """
        n_pairs = comb(self.rater_profile.n_raters, 2)

        pair_parts, pair_kappa_max = self._pairs()
        return LightKappaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=NotComputed(PARTS_NOT_COMPUTED),
            kappa_max=NotComputed(KAPPA_MAX_NOT_COMPUTED),
            agreement_profile=self._agreement_profile(),
            n_pairs=n_pairs,
            ci_method=self.config.ci_method,
            pair_parts=pair_parts,
            pair_kappa_max=pair_kappa_max,
        )
