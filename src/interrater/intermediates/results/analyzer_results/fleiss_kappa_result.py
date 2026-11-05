"""Result for Fleiss' kappa analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from ...not_computed import NotComputed
from ..agreement_parts import AgreementParts
from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ...data_profiles.agreement_profile import AgreementProfile
    from ...data_profiles.balance_profile import BalanceProfile
    from ...data_profiles.rater_profile import RaterProfile
    from ..ci_results.ci_result import CIResult

_ERROR_PREFIX = "[FleissKappaResult] Error"
_TOL = 1e-9

# Why kappa_max is not computed for me. The analyzer passes it in.
KAPPA_MAX_NOT_COMPUTED = (
    "Fleiss' kappa corrects with category shares pooled over all raters, "
    "so perfect agreement is always possible and the maximum is always 1"
)


@dataclass(frozen=True)
class FleissKappaResult(PointEstimateResult):
    """
    I am Fleiss' kappa, exact or approximate, weighted or unweighted,
    with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers whether this is the exact coefficient from Fleiss 1971
    or a generalization applied to unbalanced data.

    When is_exact is True, I am the exact formula (balanced data, constant
    raters per subject). When False, I am a valid generalization using
    per-subject averaging of p_j, but not the original coefficient.

    weights_type is on me (not just on the config) because weighted and
    unweighted Fleiss' kappa are different statistics with different
    interpretations. A result must be self-describing.

    My parts are Po and Pe, and I check them against the rest of me: they
    are weighted exactly when I am, and when unweighted, Pe is the sum of
    squared category shares in my agreement profile -- so the parts and the
    profile must come from the same data. kappa_max is always NotComputed:
    pooled shares put no limit below 1 on agreement.

    Attributes
    ----------
    point : float
        Fleiss' kappa (inherited from PointEstimateResult).
    confidence_interval : CIResult | None
        The interval, bootstrap or parametric (inherited).
    bootstrap_distribution : np.ndarray | NotComputed | None
        Bootstrap distribution (inherited); NotComputed when the draws were
        not kept.
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Its is_balanced property
        determines whether this result is exact or approximate.
    rater_profile : RaterProfile | NotComputed
        The raters (inherited). NotComputed when built from counts.
    parts : AgreementParts
        Po and Pe (inherited). Always computed for Fleiss' kappa.
    kappa_max : NotComputed
        Always NotComputed, with KAPPA_MAX_NOT_COMPUTED as the reason
        (inherited).
    agreement_profile : AgreementProfile
        Prevalence, rater marginals, specific agreement (inherited).
    ci_method : str
        Method used: "percentile", "bca", or "parametric_null".
    weights_type : str
        Weighting scheme: "identity" (unweighted), "linear", "quadratic",
        or "custom". Determines what statistic this result represents.
    weight_matrix : list[list[float]] | None, default None
        The weight matrix that was used, as nested lists (row i, column j is
        the credit for categories i and j). None when unweighted. Recorded so
        a report shows the numbers, not only the name of the scheme.
    """

    # Attribute declaration

    # What I am called: the config name, my printed name and my symbol.
    # Every other place that names this coefficient reads these.
    COEFFICIENT: ClassVar[str] = "fleiss_kappa"
    NAME:        ClassVar[str] = "Fleiss' κ"
    SYMBOL:      ClassVar[str] = "κ"

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    balance_profile: BalanceProfile | None
    rater_profile: RaterProfile | NotComputed
    parts: AgreementParts
    kappa_max: NotComputed
    agreement_profile: AgreementProfile

    # Populated at this level
    ci_method: str
    weights_type: str
    weight_matrix: list[list[float]] | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Fleiss' kappa result attributes."""

        # Fleiss always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for Fleiss' kappa"
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

        # Validate weights_type
        if not isinstance(self.weights_type, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: weights_type must be str, "
                f"got {type(self.weights_type).__name__}"
            )
        if self.weights_type not in {"identity", "linear", "quadratic", "custom"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: weights_type must be 'identity', 'linear', "
                f"'quadratic', or 'custom', got '{self.weights_type}'"
            )

        # Validate weight_matrix: none when unweighted, square when given
        if self.weight_matrix is not None:
            if self.weights_type == "identity":
                raise ValueError(
                    f"{_ERROR_PREFIX}: weight_matrix must be None when "
                    f"weights_type is 'identity' (unweighted)"
                )
            matrix = np.asarray(self.weight_matrix, dtype=float)
            if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
                raise ValueError(
                    f"{_ERROR_PREFIX}: weight_matrix must be square, "
                    f"got shape {matrix.shape}"
                )

        # Call parent validation (validates point, confidence_interval,
        # bootstrap_distribution, balance_profile type)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._check_parts_fit_fleiss()

    def _check_parts_fit_fleiss(self) -> None:
        """Po and Pe are always there, weighted as I am, Pe from my profile."""
        if not isinstance(self.parts, AgreementParts) or self.parts.kind != "agreement":
            raise TypeError(
                f"{_ERROR_PREFIX}: parts must be agreement parts (Po and Pe), "
                f"got {self.parts!r}"
            )
        weighted = self.weights_type != "identity"
        if self.parts.weighted != weighted:
            raise ValueError(
                f"{_ERROR_PREFIX}: parts.weighted is {self.parts.weighted} but "
                f"weights_type is '{self.weights_type}'"
            )
        if not isinstance(self.kappa_max, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: kappa_max must be NotComputed for Fleiss' "
                f"kappa ({KAPPA_MAX_NOT_COMPUTED}), got {self.kappa_max!r}"
            )
        if weighted or math.isnan(self.parts.chance):
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

    def _kappa_max_note(self) -> str:
        return "κ_max not computed (pooled shares)"

    def _data_line(self) -> str:
        exact = "exact formula" if self.is_exact else "approximation (unbalanced data)"
        return (
            f"{self._subjects_phrase()}; {self._raters_phrase()}; "
            f"{self._weights_phrase(self.weights_type)}; {exact}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and Fleiss-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'is_exact', 'ci_method',
            'weights_type', 'weight_matrix'.
        """
        result = super().summary()
        result.update({
            "is_exact": self.is_exact,
            "ci_method": self.ci_method,
            "weights_type": self.weights_type,
            "weight_matrix": self.weight_matrix,
        })
        return result

    @property
    def is_exact(self) -> bool:
        """
        Whether this is the exact Fleiss formula or an approximation.

        Returns True if data was balanced (constant raters per subject).
        Returns False if an approximation was used for unbalanced data.
        """
        return self.balance_profile.is_balanced
