"""Result for Brennan-Prediger analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from ...not_computed import NotComputed
from ..agreement_parts import AgreementParts
from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile
    from ...data_profiles.rater_profile import RaterProfile
    from ...data_profiles.agreement_profile import AgreementProfile
    from ..ci_results.ci_result import CIResult

_ERROR_PREFIX = "[BrennanPredigerResult] Error"
_TOL = 1e-9
# Printed on unbalanced data: the formula is Gwet's, not an approximation
UNBALANCED_DATA_LABEL = "unbalanced data: Gwet (2016) formula for missing ratings"

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
    own rater count, over the same Pe (Gwet, K. L. (2016). Educational and
    Psychological Measurement, 76(4), 609-637, Appendix A). Not an
    approximation; my printed form names the source.

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
    rater_profile : RaterProfile
        The raters (inherited).
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
    _ACCEPTS_WEIGHTS: ClassVar[bool] = True

    # Inherited from parent class PanelResult (keyword-only)
    balance_profile: BalanceProfile = field(kw_only=True)
    rater_profile: RaterProfile = field(kw_only=True)
    ci_method: str = field(kw_only=True)
    weights_type: str = field(default="identity", kw_only=True)
    weight_matrix: list[list[float]] | None = field(default=None, kw_only=True)

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    parts: AgreementParts
    kappa_max: NotComputed
    agreement_profile: AgreementProfile

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Brennan-Prediger result attributes."""

        # Call parent validation (profiles, interval method, weighting,
        # point, interval and draws)
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
        n_categories = len(self.agreement_profile.category_names)
        if not weighted:
            pe, formula = 1.0 / n_categories, f"1/K with K = {n_categories}"
        elif self.weight_matrix is not None:
            matrix = np.asarray(self.weight_matrix, dtype=float)
            if matrix.shape != (n_categories, n_categories):
                raise ValueError(
                    f"{_ERROR_PREFIX}: weight_matrix is {matrix.shape[0]}x"
                    f"{matrix.shape[1]} but the agreement profile has "
                    f"{n_categories} categories"
                )
            pe = float(matrix.sum()) / (n_categories * n_categories)
            formula = "sum(W)/K^2"
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
            Keys from parent (including balance, ci_method, weights_type,
            weight_matrix) plus 'is_exact'.
        """
        result = super().summary()
        result.update({
            "is_exact": self.is_exact,
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
