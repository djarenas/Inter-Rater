"""Result for Conger's kappa analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from ...not_computed import NotComputed
from ..agreement_parts import AgreementParts
from .individual_marginal_result import (
    IndividualMarginalResult,
    _same,
    _unweighted_pair_kappa_max,
)

if TYPE_CHECKING:
    from ...data_profiles.agreement_profile import AgreementProfile
    from ...data_profiles.balance_profile import BalanceProfile
    from ..ci_results.ci_result import CIResult

_ERROR_PREFIX = "[CongerKappaResult] Error"

# Why kappa_max is not computed for 3 or more raters. The analyzer passes it in.
KAPPA_MAX_NOT_COMPUTED = (
    "with 3 or more raters the most agreement depends on fitting every "
    "rater's marginals at once, which has no simple formula; with 2 raters "
    "it is Cohen's"
)


@dataclass(frozen=True)
class CongerKappaResult(IndividualMarginalResult):
    """
    I am Conger's kappa for any number of raters (n >= 2).

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult.

    Conger's kappa generalises Cohen's kappa to many raters while keeping
    individual marginals. At exactly two raters I reduce to Cohen's unweighted
    kappa. I am always unweighted -- Conger (1980) defined the coefficient for
    nominal data only -- so I carry no weights_type field.

    I exist only for complete data, where every rater rated every subject
    (Conger 1980), so my balance profile must have a fill rate of 1.

    My parts are Po and Pe, unweighted, and Pe must be Conger's chance term
    on the raters' marginals in my agreement profile:
    Sum_j p_bar_j^2 - Sum_j var_r(p_rj) / n_raters. kappa_max is computed
    for two raters, where I am Cohen's kappa and his ceiling applies, and
    NotComputed for more.

    References
    ----------
    Conger, A. J. (1980). Integration and generalization of kappas for
        multiple raters. Psychological Bulletin, 88(2), 322--328.

    Attributes
    ----------
    point : float
        Conger's kappa (inherited from PointEstimateResult).
    confidence_interval : CIResult | None
        Bootstrap confidence interval (inherited).
    bootstrap_distribution : np.ndarray | NotComputed | None
        Bootstrap distribution (inherited); NotComputed when the draws were
        not kept.
    balance_profile : BalanceProfile
        Balance structure of the data (inherited).
    parts : AgreementParts
        Po and Pe, unweighted (inherited). Always computed.
    kappa_max : float | NotComputed
        Cohen's ceiling with two raters; NotComputed with more, with
        KAPPA_MAX_NOT_COMPUTED as the reason (inherited).
    agreement_profile : AgreementProfile
        Prevalence, each rater's marginals, specific agreement (inherited).
    ci_method : str
        Method used: "percentile" or "bca".
    """

    # Attribute declaration

    # What I am called: the config name, my printed name and my symbol.
    # Every other place that names this coefficient reads these.
    COEFFICIENT: ClassVar[str] = "conger_kappa"
    NAME:        ClassVar[str] = "Conger's κ"
    SYMBOL:      ClassVar[str] = "κ"

    # Inherited from IndividualMarginalResult (which validates
    # balance_profile is not None and n_raters is not None)
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    balance_profile: BalanceProfile | None
    parts: AgreementParts
    kappa_max: float | NotComputed
    agreement_profile: AgreementProfile

    # Populated at this level
    ci_method: str

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Conger's kappa result attributes."""

        # Validate ci_method (no parametric_null for Conger)
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or 'bca', "
                f"got '{self.ci_method}'. Conger's kappa has no closed-form "
                f"null variance, so 'parametric_null' is not available."
            )

        # Call parent validation (IndividualMarginalResult validates
        # balance_profile and n_raters; PointEstimateResult validates
        # point, confidence_interval, bootstrap_distribution)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._check_parts_fit_conger()

    def _check_parts_fit_conger(self) -> None:
        """Complete data; Po and Pe fit the raters' marginals; kappa_max only for a pair."""
        if self.balance_profile.fill_rate != 1.0:
            raise ValueError(
                f"{_ERROR_PREFIX}: Conger's kappa is for complete data, where "
                f"every rater rated every subject; balance_profile has a fill "
                f"rate of {self.balance_profile.fill_rate}"
            )
        if not isinstance(self.parts, AgreementParts) or self.parts.kind != "agreement":
            raise TypeError(
                f"{_ERROR_PREFIX}: parts must be agreement parts (Po and Pe), "
                f"got {self.parts!r}"
            )
        if self.parts.weighted:
            raise ValueError(
                f"{_ERROR_PREFIX}: parts must be unweighted; weighted Conger's "
                f"kappa is not implemented"
            )
        pair = self.balance_profile.n_raters == 2
        if pair and isinstance(self.kappa_max, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: kappa_max must be computed with 2 raters, "
                f"where Conger's kappa is Cohen's and his ceiling applies"
            )
        if not pair and not isinstance(self.kappa_max, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: kappa_max must be NotComputed with "
                f"{self.balance_profile.n_raters} raters "
                f"({KAPPA_MAX_NOT_COMPUTED}), got {self.kappa_max!r}"
            )

        marginals = self._marginal_rows()
        if not math.isnan(self.parts.chance):
            pe = float((marginals.mean(axis=0) ** 2).sum()) - float(
                marginals.var(axis=0, ddof=1).sum()
            ) / marginals.shape[0]
            if not _same(pe, self.parts.chance):
                raise ValueError(
                    f"{_ERROR_PREFIX}: Pe is {self.parts.chance:.6f}, but "
                    f"Conger's chance term on the raters' marginals in the "
                    f"agreement profile is {pe:.6f}; the parts and the "
                    f"profile describe different data"
                )
        if pair and not math.isnan(self.point):
            ceiling = _unweighted_pair_kappa_max(marginals[0], marginals[1])
            if not _same(ceiling, self.kappa_max):
                raise ValueError(
                    f"{_ERROR_PREFIX}: kappa_max is {self.kappa_max:.6f}, but "
                    f"the raters' marginals allow {ceiling:.6f}"
                )

    # Public API ------------------------------------------------------------- #

    # Printed form: see PointEstimateResult.__str__ ----------------------- #

    def _kappa_max_note(self) -> str:
        if isinstance(self.kappa_max, NotComputed):
            return "κ_max not computed (3 or more raters)"
        note = super()._kappa_max_note()
        return note if note == "κ_max undefined" else f"{note} (these raters' marginals)"

    def _data_line(self) -> str:
        return (
            f"{self.balance_profile.n_subjects:,} subjects; "
            f"{self._raters_phrase()}, every rater rated every subject"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and Conger-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'ci_method'.
        """
        result = super().summary()
        result.update({
            "ci_method": self.ci_method,
        })
        return result
