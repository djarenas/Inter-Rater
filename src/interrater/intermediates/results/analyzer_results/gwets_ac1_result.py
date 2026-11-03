"""Result for Gwet's AC1 analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ...not_computed import NotComputed
from ..agreement_parts import AgreementParts
from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ...data_profiles.agreement_profile import AgreementProfile
    from ...data_profiles.balance_profile import BalanceProfile
    from ..ci_results.ci_result import CIResult

_ERROR_PREFIX = "[GwetsAC1Result] Error"
_TOL = 1e-9

# Why kappa_max is not computed for me. The analyzer passes it in.
KAPPA_MAX_NOT_COMPUTED = (
    "Gwet's AC1 takes chance from category shares pooled over all raters, "
    "which perfect agreement can always match, so the maximum is always 1"
)


@dataclass(frozen=True)
class GwetsAC1Result(PointEstimateResult):
    """
    I am Gwet's AC1, exact or approximate, with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers whether this is the exact coefficient or a generalization
    applied to unbalanced data.

    When is_exact is True, I am the exact formula (balanced data, constant
    raters per subject). When False, I am a valid generalization using
    per-subject averaging of p_j, but not the original coefficient.

    I exist as a separate class (rather than just using PointEstimateResult)
    for type clarity. This is especially useful for error chaining.

    My parts are Po and Pe, unweighted, and Pe must be
    Sum_j p_j (1 - p_j) / (K - 1) over the category shares in my agreement
    profile -- so the parts and the profile must come from the same data.
    kappa_max is always NotComputed: this chance term puts no limit below 1
    on agreement.

    Attributes
    ----------
    point : float
        Gwet's AC1 (inherited from PointEstimateResult).
    confidence_interval : CIResult | None
        The bootstrap interval (inherited).
    bootstrap_distribution : np.ndarray | NotComputed | None
        Bootstrap distribution (inherited); NotComputed when the draws were
        not kept.
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Its is_balanced property
        determines whether this result is exact or approximate.
    parts : AgreementParts
        Po and Pe, unweighted (inherited). Always computed for AC1.
    kappa_max : NotComputed
        Always NotComputed, with KAPPA_MAX_NOT_COMPUTED as the reason
        (inherited).
    agreement_profile : AgreementProfile
        Prevalence, rater marginals, specific agreement (inherited).
    ci_method : str
        Method used: "percentile" or "bca".
    """

    # Attribute declaration

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    balance_profile: BalanceProfile | None
    parts: AgreementParts
    kappa_max: NotComputed
    agreement_profile: AgreementProfile

    # Populated at this level
    ci_method: str

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Gwet's AC1 result attributes."""

        # Gwet's AC1 always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for Gwet's AC1"
            )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method 'parametric_null' is not implemented "
                f"for Gwet's AC1"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or 'bca', "
                f"got '{self.ci_method}'"
            )

        # Call parent validation (validates point, confidence_interval,
        # bootstrap_distribution, balance_profile type)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._check_parts_fit_gwet()

    def _check_parts_fit_gwet(self) -> None:
        """Po and Pe are always there, unweighted, Pe from my profile."""
        if not isinstance(self.parts, AgreementParts) or self.parts.kind != "agreement":
            raise TypeError(
                f"{_ERROR_PREFIX}: parts must be agreement parts (Po and Pe), "
                f"got {self.parts!r}"
            )
        if self.parts.weighted:
            raise ValueError(
                f"{_ERROR_PREFIX}: parts must be unweighted; the weighted form "
                f"is Gwet's AC2, a different coefficient"
            )
        if not isinstance(self.kappa_max, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: kappa_max must be NotComputed for Gwet's "
                f"AC1 ({KAPPA_MAX_NOT_COMPUTED}), got {self.kappa_max!r}"
            )
        if math.isnan(self.parts.chance):
            return
        shares = list(self.agreement_profile.prevalence.values())
        pe = sum(p * (1 - p) for p in shares) / (len(shares) - 1)
        if abs(pe - self.parts.chance) > _TOL:
            raise ValueError(
                f"{_ERROR_PREFIX}: Pe is {self.parts.chance:.6f}, but the "
                f"agreement profile's category shares give {pe:.6f}; the parts "
                f"and the profile describe different data"
            )

    # Public API ------------------------------------------------------------- #

    # Printed form: see PointEstimateResult.__str__ ----------------------- #

    def _name(self) -> str:
        return "Gwet's AC1"

    def _kappa_max_note(self) -> str:
        return "κ_max not computed (pooled shares)"

    def _data_line(self) -> str:
        exact = "exact formula" if self.is_exact else "approximation (unbalanced data)"
        return f"{self._subjects_phrase()}; {self._raters_phrase()}; {exact}"

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and AC1-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'is_exact', 'ci_method'.
        """
        result = super().summary()
        result.update({
            "is_exact": self.is_exact,
            "ci_method": self.ci_method,
        })
        return result

    @property
    def is_exact(self) -> bool:
        """
        Whether this is the exact AC1 formula or an approximation.

        Returns True if data was balanced (constant raters per subject).
        Returns False if an approximation was used for unbalanced data.
        """
        return self.balance_profile.is_balanced
