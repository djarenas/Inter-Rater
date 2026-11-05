"""Result for Cohen's kappa analysis."""

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
    from ...data_profiles.rater_profile import RaterProfile
    from ..ci_results.ci_result import CIResult

_ERROR_PREFIX = "[CohenKappaResult] Error"


@dataclass(frozen=True)
class CohenKappaResult(IndividualMarginalResult):
    """
    I am Cohen's kappa, weighted or unweighted, for a single pair of raters.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. I add CI method and weights_type
    so anyone reading me knows what statistic I am.

    weights_type is on me (not just on the config) because weighted Cohen's
    kappa and unweighted Cohen's kappa are different statistics with different
    interpretations. A result must be self-describing.

    Cohen's kappa operates on paired data (two raters, fully crossed), so
    balance is guaranteed by construction. My balance_profile is always
    populated and always balanced. n_subjects is available via
    balance_profile.n_subjects.

    My rater_profile holds the two raters: their names, how many subjects
    each rated, and how many they share.

    My parts are Po and Pe, weighted exactly when I am, and kappa_max is
    always computed: two raters' marginals cap their agreement (Cohen 1960).
    I check both against the two raters' marginals in my agreement profile:
    Pe is Sum_j p_aj p_bj (p_a W p_b when weighted), and unweighted
    kappa_max is built from Sum_j min(p_aj, p_bj). The weighted maximum
    needs a linear program, so it is checked only through kappa <= kappa_max.

    References
    ----------
    Cohen, J. (1960). A coefficient of agreement for nominal scales.
        Educational and Psychological Measurement, 20(1), 37--46.
    Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision
        for scaled disagreement or partial credit. Psychological Bulletin,
        70(4), 213--220.

    Attributes
    ----------
    point : float
        Cohen's kappa (inherited from PointEstimateResult).
    confidence_interval : CIResult | None
        Bootstrap or parametric confidence interval (inherited).
    bootstrap_distribution : np.ndarray | NotComputed | None
        Bootstrap distribution (inherited); NotComputed when the draws were
        not kept.
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Always populated for
        Cohen's kappa. n_subjects available via balance_profile.n_subjects.
    rater_profile : RaterProfile
        The raters (inherited). Always a RaterProfile: Cohen's kappa needs
        rater identity.
    parts : AgreementParts
        Po and Pe (inherited). Always computed for Cohen's kappa.
    kappa_max : float
        The largest kappa these two raters' marginals allow (inherited).
        NaN when kappa itself is undefined.
    agreement_profile : AgreementProfile
        Prevalence, the two raters' marginals, specific agreement, and the
        prevalence and bias indices for two categories (inherited).
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
    COEFFICIENT: ClassVar[str] = "cohen_kappa"
    NAME:        ClassVar[str] = "Cohen's κ"
    SYMBOL:      ClassVar[str] = "κ"

    # Inherited from IndividualMarginalResult (which validates
    # balance_profile is not None and rater_profile is a RaterProfile)
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    balance_profile: BalanceProfile | None
    rater_profile: RaterProfile
    parts: AgreementParts
    kappa_max: float
    agreement_profile: AgreementProfile

    # Populated at this level
    ci_method: str
    weights_type: str
    weight_matrix: list[list[float]] | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Cohen's kappa result attributes."""

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

        # Call parent validation (IndividualMarginalResult validates
        # balance_profile and rater_profile; PointEstimateResult validates
        # point, confidence_interval, bootstrap_distribution)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        self._check_parts_fit_cohen()

    def _check_parts_fit_cohen(self) -> None:
        """Po, Pe and kappa_max fit the two raters' marginals."""
        if self.rater_profile.n_raters != 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: Cohen's kappa is for one pair of raters, "
                f"got {self.rater_profile.n_raters}"
            )
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
        if isinstance(self.kappa_max, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: kappa_max must be computed for Cohen's "
                f"kappa: two raters' marginals always cap it (Cohen 1960)"
            )
        if math.isnan(self.point):
            return

        p_a, p_b = self._marginal_rows()
        if not weighted:
            pe = float((p_a * p_b).sum())
        elif self.weight_matrix is not None:
            pe = float(p_a @ np.asarray(self.weight_matrix, dtype=float) @ p_b)
        else:
            return
        if not _same(pe, self.parts.chance):
            raise ValueError(
                f"{_ERROR_PREFIX}: Pe is {self.parts.chance:.6f}, but the "
                f"raters' marginals in the agreement profile give {pe:.6f}; "
                f"the parts and the profile describe different data"
            )
        if not weighted:
            ceiling = _unweighted_pair_kappa_max(p_a, p_b)
            if not _same(ceiling, self.kappa_max):
                raise ValueError(
                    f"{_ERROR_PREFIX}: kappa_max is {self.kappa_max:.6f}, but "
                    f"the raters' marginals allow {ceiling:.6f}"
                )

    # Public API ------------------------------------------------------------- #

    # Printed form: see PointEstimateResult.__str__ ----------------------- #

    def _kappa_max_note(self) -> str:
        note = super()._kappa_max_note()
        return note if note == "κ_max undefined" else f"{note} (these raters' marginals)"

    def _data_line(self) -> str:
        return (
            f"{self._subjects_phrase()}; "
            f"{self._weights_phrase(self.weights_type)}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and Cohen-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'ci_method', 'weights_type', 'weight_matrix'.
        """
        result = super().summary()
        result.update({
            "ci_method": self.ci_method,
            "weights_type": self.weights_type,
            "weight_matrix": self.weight_matrix,
        })
        return result
