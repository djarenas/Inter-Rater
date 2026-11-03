"""Result for Scott's pi analysis."""

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

_ERROR_PREFIX = "[ScottsPiResult] Error"
_TOL = 1e-9

# Why kappa_max is not computed for me. The analyzer passes it in.
KAPPA_MAX_NOT_COMPUTED = (
    "Scott's pi corrects with category shares pooled over both raters, "
    "so perfect agreement is always possible and the maximum is always 1"
)


@dataclass(frozen=True)
class ScottsPiResult(PointEstimateResult):
    """
    I am Scott's pi, exact or approximate, with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers whether this is the exact coefficient from Scott 1955
    or a generalization applied to unbalanced data.

    When is_exact is True, I am the exact formula (balanced data, constant
    raters per subject). When False, I am a valid generalization using
    per-subject averaging of p_j, but not the original coefficient.

    I exist as a separate class (rather than just using PointEstimateResult)
    for type clarity. Downstream code can see `result: ScottsPiResult` and
    know the coefficient and that metadata about balance is available.

    My parts are Po and Pe, unweighted, and Pe must be the sum of squared
    category shares in my agreement profile -- so the parts and the profile
    must come from the same data. kappa_max is always NotComputed: pooled
    shares put no limit below 1 on agreement.

    Attributes
    ----------
    point : float
        Scott's pi (inherited from PointEstimateResult).
    confidence_interval : CIResult | None
        The interval, bootstrap or parametric (inherited).
    bootstrap_distribution : np.ndarray | NotComputed | None
        Bootstrap distribution (inherited); NotComputed when the draws were
        not kept.
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Its is_balanced property
        determines whether this result is exact or approximate.
    parts : AgreementParts
        Po and Pe, unweighted (inherited). Always computed for Scott's pi.
    kappa_max : NotComputed
        Always NotComputed, with KAPPA_MAX_NOT_COMPUTED as the reason
        (inherited).
    agreement_profile : AgreementProfile
        Prevalence, rater marginals, specific agreement (inherited).
    ci_method : str
        Method used: "percentile", "bca", or "parametric_null".
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
        """Validate Scott's pi result attributes."""

        # Scott's pi always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for Scott's pi"
            )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca", "parametric_null"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile', 'bca', "
                f"or 'parametric_null', got '{self.ci_method}'"
            )

        # Call parent validation (validates point, confidence_interval,
        # bootstrap_distribution, balance_profile type)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._check_parts_fit_scott()

    def _check_parts_fit_scott(self) -> None:
        """Po and Pe are always there, unweighted, Pe from my profile."""
        if not isinstance(self.parts, AgreementParts) or self.parts.kind != "agreement":
            raise TypeError(
                f"{_ERROR_PREFIX}: parts must be agreement parts (Po and Pe), "
                f"got {self.parts!r}"
            )
        if self.parts.weighted:
            raise ValueError(
                f"{_ERROR_PREFIX}: parts must be unweighted; Scott's pi has "
                f"no weighted form"
            )
        if not isinstance(self.kappa_max, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: kappa_max must be NotComputed for Scott's "
                f"pi ({KAPPA_MAX_NOT_COMPUTED}), got {self.kappa_max!r}"
            )
        if math.isnan(self.parts.chance):
            return
        pe = sum(p * p for p in self.agreement_profile.prevalence.values())
        if abs(pe - self.parts.chance) > _TOL:
            raise ValueError(
                f"{_ERROR_PREFIX}: Pe is {self.parts.chance:.6f}, but the "
                f"agreement profile's category shares give {pe:.6f}; the parts "
                f"and the profile describe different data"
            )

    # Public API ------------------------------------------------------------- #

    # Printed form: see PointEstimateResult.__str__ ----------------------- #

    def _name(self) -> str:
        return "Scott's π"

    def _kappa_max_note(self) -> str:
        return "κ_max not computed (pooled shares)"

    def _data_line(self) -> str:
        exact = "exact formula" if self.is_exact else "approximation (unbalanced data)"
        return f"{self._subjects_phrase()}; {self._raters_phrase()}; {exact}"

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and Scott-specific metadata.

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
        Whether this is the exact Scott formula or an approximation.

        Returns True if data was balanced (constant raters per subject).
        Returns False if an approximation was used for unbalanced data.
        """
        return self.balance_profile.is_balanced
