"""Analyzer for Conger's kappa coefficient."""

from __future__ import annotations

import numpy as np

from .individual_marginal_analyzer import IndividualMarginalAnalyzer, check_complete

from ..data_objects.ratings import Ratings
from ..intermediates.results import CongerKappaResult, HypothesisTestResult
from ..intermediates.data_profiles.balance_profile import BalanceProfile
from ..intermediates.data_profiles.rater_profile import RaterProfile
from .configs import CongerKappaConfig

from ..intermediates.estimator_parts.estimator_input import EstimatorInput
from .agreement_functions.individual_marginals import (
    cohen_kappa_max,
    conger_kappa,
    conger_parts,
    rater_marginals,
)
from ..intermediates.not_computed import NotComputed
from ..intermediates.results import AgreementParts
from ..intermediates.results.analyzer_results.conger_kappa_result import (
    KAPPA_MAX_NOT_COMPUTED,
)

_ERROR_PREFIX = "[CongerKappaAnalyzer] Error"


class CongerKappaAnalyzer(IndividualMarginalAnalyzer):
    """
    I am the analyzer for Conger's kappa, the multi-rater agreement coefficient
    with individual marginals.

    Conger's kappa (Conger 1980) generalises Cohen's kappa to two or more
    raters while preserving each rater's own category distribution in
    the chance term.

    At exactly two raters, Conger's kappa reduces to Cohen's unweighted kappa.

    I validate that:
    - Ratings and config pass my parent's checks, which refuse ordered or
      weighted semantics
    - Config is CongerKappaConfig
    - Every rater rated every subject, after any injection

    I compute:
    - Point estimate of Conger's kappa (always unweighted)
    - Bootstrap confidence interval (percentile or BCa -- no parametric null)
    - Return a CongerKappaResult

    Conger's kappa reads rater identity -- who said what -- and counts discard
    that. If you have only counts, you want Fleiss'
    kappa, not Conger's.

    The package implements no weighted form for Conger's kappa. Conger
    (1980) defined the coefficient for nominal data only. I refuse ordered
    or weighted semantics (``_refuse_weights_if_unsupported``), and
    ``conger_kappa`` raises NotImplementedError on non-identity weights.

    The package implements no parametric null for Conger's kappa, so
    confidence intervals use bootstrap (percentile or BCa) and hypothesis
    tests use permutation only.

    Attributes
    ----------
    ratings : Ratings
        The validated ratings data (from Analyzer).
    balance_profile : BalanceProfile
        Describes the balance structure of the data (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped. Always present here:
        this coefficient needs rater identity.
    weights : None
        Always None: Conger's kappa has no weighted form here (from
        Analyzer).
    estimator_input : EstimatorInput
        What my coefficient is computed from: category_index and
        n_categories, with weights=None (from Analyzer).
    config : CongerKappaConfig
        The analysis configuration for Conger's kappa.

    References
    ----------
    Conger, A. J. (1980). Integration and generalization of kappas for
        multiple raters. Psychological Bulletin, 88(2), 322--328.
    """

    # Label in hypothesis test results
    _COEFFICIENT_KEY = "conger_kappa"

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    weights: None
    estimator_input: EstimatorInput
    # Checked at this level
    config: CongerKappaConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(self, ratings: Ratings, config: CongerKappaConfig) -> None:
        """
        Run the parent's checks, then check the config type and that the
        data are complete.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze. Must contain at least two raters,
            and every rater must have rated every subject.
        config : CongerKappaConfig
            The analysis configuration for Conger's kappa.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object, or config is not
            CongerKappaConfig.
        ValueError
            If the semantics are ordered or weighted, or any rating is
            missing (after injection).
        """
        # Call parent __init__ (validates ratings and checks that config is
        # IndividualMarginalConfig, computes BalanceProfile, builds
        # EstimatorInput with weights=None)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate config type
        self._validate_conger_config()

        # Complete data only, checked after any injection
        try:
            check_complete(self.ratings.category_index, "Conger's kappa")
        except ValueError as exc:
            raise ValueError(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

    def hypothesis_test(self) -> HypothesisTestResult:
        """
        Test H_0: kappa = 0 by permutation (``_permutation_test``).

        The package implements no parametric null for Conger's kappa.
        """
        return self._permutation_test()

    # Private helpers -------------------------------------------------------- #

    def _validate_conger_config(self) -> None:
        """
        Validate that self.config is a CongerKappaConfig.

        Raises
        ------
        TypeError
            If config is not CongerKappaConfig.
        """
        if not isinstance(self.config, CongerKappaConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be CongerKappaConfig, "
                f"got {type(self.config).__name__}"
            )

    def _agreement_parts(self) -> AgreementParts:
        """
        Po and Pe from the same function conger_kappa is built from, so the
        result can check that they rebuild the point estimate.
        """
        po, pe = conger_parts(self.estimator_input)
        return AgreementParts(
            observed=po, chance=pe, kind="agreement", weighted=False
        )

    def _kappa_max(self) -> float | NotComputed:
        """
        Cohen's ceiling with two raters, where Conger's kappa is Cohen's;
        NotComputed with more.
        """
        if self.estimator_input.category_index.shape[1] != 2:
            return NotComputed(KAPPA_MAX_NOT_COMPUTED)
        marginals = rater_marginals(self.estimator_input)
        return cohen_kappa_max(marginals[0], marginals[1])

    def _estimate_from_input(self, data: EstimatorInput) -> float:
        """
        Conger's kappa from an EstimatorInput.

        The only coefficient-specific step in my resampling. Everything that
        draws samples is written once on Analyzer in terms of this. Conger's
        chance term is computed from each rater's own marginals, so the
        column positions of data matter.
        """
        return conger_kappa(data)

    def _result(
        self, point: float, ci, bootstrap_distribution: np.ndarray | None
    ) -> CongerKappaResult:
        """
        My result: the point, interval and draws ``estimate()`` (or a
        permutation analysis) computed, with my parts, ceiling and profiles.

        Only percentile and bca intervals: the package implements no
        parametric null for Conger's kappa.
        """
        return CongerKappaResult(
            point=point,
            confidence_interval=ci,
            bootstrap_distribution=bootstrap_distribution,
            balance_profile=self.balance_profile,
            rater_profile=self.rater_profile,
            parts=self._agreement_parts(),
            kappa_max=self._kappa_max(),
            agreement_profile=self._agreement_profile(),
            ci_method=self.config.ci_method,
        )
