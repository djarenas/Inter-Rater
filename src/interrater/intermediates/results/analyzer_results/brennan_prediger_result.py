"""Result for Brennan-Prediger analysis."""

from __future__ import annotations

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

_ERROR_PREFIX = "[BrennanPredigerResult] Error"
_TOL = 1e-9
# Printed on unbalanced data: the formula is Gwet's, not an approximation
UNBALANCED_DATA_LABEL = "unbalanced data: Gwet (2014) formula for missing ratings"

# Why kappa_max is not computed for me. The analyzer passes it in.
KAPPA_MAX_NOT_COMPUTED = (
    "Brennan-Prediger's chance term is fixed by the number of categories, "
    "so perfect agreement is always possible and the maximum is always 1"
)


@dataclass(frozen=True)
class BrennanPredigerResult(PointEstimateResult):
    """
    I am Brennan-Prediger's coefficient, exact or approximate, weighted or
    unweighted, with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers whether the data was balanced.

    When is_exact is True, the data was balanced and the formula is the
    standard Brennan-Prediger with Pe = 1/K. When False, I am Gwet's formula
    for missing ratings: Po averaged per subject, each subject divided by its
    own rater count, over the same Pe (Gwet, K. L. (2014). Handbook of
    Inter-Rater Reliability, 4th ed., as implemented in irrCAC's
    ``bp.coeff.raw``). Not an approximation; my printed form names the source.

    weights_type is on me (not just on the config) because weighted and
    unweighted Brennan-Prediger are different statistics with different
    interpretations. A result must be self-describing.

    My parts are Po and Pe, weighted exactly when I am. Pe depends only on
    the categories: 1/K unweighted, sum(W)/K^2 with my weight matrix -- so I
    check it against the K in my agreement profile. kappa_max is always
    NotComputed: a fixed chance term puts no limit below 1 on agreement.

    Attributes
    ----------
    point : float
        Brennan-Prediger coefficient (inherited from PointEstimateResult).
    confidence_interval : CIResult | None
        The bootstrap interval (inherited).
    bootstrap_distribution : np.ndarray | NotComputed | None
        Bootstrap distribution (inherited); NotComputed when the draws were
        not kept.
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Its is_balanced property
        determines whether this result is exact or approximate.
    rater_profile : RaterProfile | NotComputed
        The raters (inherited). NotComputed when built from counts.
    parts : AgreementParts
        Po and Pe (inherited). Always computed for Brennan-Prediger.
    kappa_max : NotComputed
        Always NotComputed, with KAPPA_MAX_NOT_COMPUTED as the reason
        (inherited).
    agreement_profile : AgreementProfile
        Prevalence, rater marginals, specific agreement (inherited).
    ci_method : str
        Method used: "percentile" or "bca".
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
    COEFFICIENT: ClassVar[str] = "brennan_prediger"
    NAME:        ClassVar[str] = "Brennan-Prediger"
    SYMBOL:      ClassVar[str] = "B-P"

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    balance_profile: BalanceProfile
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
        """Validate Brennan-Prediger result attributes."""

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or 'bca', "
                f"got '{self.ci_method}'"
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

        self._check_parts_fit_brennan_prediger()

    def _check_parts_fit_brennan_prediger(self) -> None:
        """Po and Pe are always there, weighted as I am, Pe from K."""
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
                f"{_ERROR_PREFIX}: kappa_max must be NotComputed for "
                f"Brennan-Prediger ({KAPPA_MAX_NOT_COMPUTED}), "
                f"got {self.kappa_max!r}"
            )
        k = len(self.agreement_profile.category_names)
        if not weighted:
            pe, formula = 1.0 / k, f"1/K with K = {k}"
        elif self.weight_matrix is not None:
            matrix = np.asarray(self.weight_matrix, dtype=float)
            if matrix.shape != (k, k):
                raise ValueError(
                    f"{_ERROR_PREFIX}: weight_matrix is {matrix.shape[0]}x"
                    f"{matrix.shape[1]} but the agreement profile has {k} "
                    f"categories"
                )
            pe, formula = float(matrix.sum()) / (k * k), "sum(W)/K^2"
        else:
            return
        if abs(pe - self.parts.chance) > _TOL:
            raise ValueError(
                f"{_ERROR_PREFIX}: Pe is {self.parts.chance:.6f}, but "
                f"{formula} is {pe:.6f}; the parts and the profile describe "
                f"different categories"
            )

    # Public API -------------------------------------------------------------  #

    # Printed form: see PointEstimateResult.__str__ ----------------------- #

    def _kappa_max_note(self) -> str:
        return "κ_max not computed (fixed chance term)"

    def _data_line(self) -> str:
        exact = "exact formula" if self.is_exact else UNBALANCED_DATA_LABEL
        return (
            f"{self._subjects_phrase()}; {self._raters_phrase()}; "
            f"{self._weights_phrase(self.weights_type)}; {exact}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and BP-specific metadata.

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
        Whether the data was balanced, where every reference agrees on the
        formula.

        Returns True if data was balanced (constant raters per subject).
        Returns False on unbalanced data, where I am Gwet's (2014) formula
        for missing ratings.
        """
        return self.balance_profile.is_balanced
