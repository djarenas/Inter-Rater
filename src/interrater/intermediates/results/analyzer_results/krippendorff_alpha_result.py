"""Result for Krippendorff's alpha analysis."""

from __future__ import annotations

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

_ERROR_PREFIX = "[KrippendorffAlphaResult] Error"

# Why kappa_max is not computed for me. The analyzer passes it in.
KAPPA_MAX_NOT_COMPUTED = (
    "Krippendorff's alpha takes expected disagreement from all ratings "
    "pooled, so perfect agreement is always possible and the maximum is "
    "always 1"
)


@dataclass(frozen=True)
class KrippendorffAlphaResult(PointEstimateResult):
    """
    I am Krippendorff's alpha, weighted or unweighted, with data context.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult. balance_profile tells users
    and visualizers about the balance structure of the data.

    Unlike the kappa family, alpha handles unbalanced data natively via
    pair-pooling -- a subject rated by six raters contributes fifteen pairs
    and one rated by two contributes a single pair. There is no
    approximation involved: the unbalanced formula *is* the formula. So
    is_exact is always True, and the balance profile is informational
    context rather than a correctness flag.

    Alpha measures disagreement rather than agreement:

        alpha = 1 - Do / De

    This is not the same expression as (Po - Pe) / (1 - Pe) rearranged.
    The two use different units (pairs pooled vs. subjects averaged) and
    different chance models (finite-sample permutation vs. population
    rates).

    weights_type is on me (not just on the config) because weighted and
    unweighted alpha are different statistics with different
    interpretations. A result must be self-describing.

    My parts are Do and De, the observed and expected disagreement, and
    they are weighted exactly when I am. Unlike the kappa family, De is not
    checked against my agreement profile: De pools every pairable rating,
    while the profile's shares are averaged per subject, and the two differ
    on unbalanced data. kappa_max is always NotComputed.

    Attributes
    ----------
    point : float
        Krippendorff's alpha (inherited from PointEstimateResult).
    confidence_interval : CIResult | None
        The bootstrap interval (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Informational only --
        alpha handles both balanced and unbalanced data exactly.
    parts : AgreementParts
        Do and De, kind "disagreement" (inherited). Always computed for
        alpha.
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

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | None
    balance_profile: BalanceProfile | None
    parts: AgreementParts
    kappa_max: NotComputed
    agreement_profile: AgreementProfile

    # Populated at this level
    ci_method: str
    weights_type: str
    weight_matrix: list[list[float]] | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Krippendorff's alpha result attributes."""

        # Krippendorff's alpha always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for "
                f"Krippendorff's alpha"
            )

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

        self._check_parts_fit_alpha()

    def _check_parts_fit_alpha(self) -> None:
        """Do and De are always there, weighted as I am."""
        if not isinstance(self.parts, AgreementParts) or self.parts.kind != "disagreement":
            raise TypeError(
                f"{_ERROR_PREFIX}: parts must be disagreement parts (Do and "
                f"De), got {self.parts!r}"
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
                f"Krippendorff's alpha ({KAPPA_MAX_NOT_COMPUTED}), "
                f"got {self.kappa_max!r}"
            )

    # Public API -------------------------------------------------------------  #

    # Printed form: see PointEstimateResult.__str__ ----------------------- #

    def _name(self) -> str:
        return "Krippendorff's α"

    def _kappa_max_note(self) -> str:
        return "κ_max not computed (pooled ratings)"

    def _data_line(self) -> str:
        # Alpha is exact on unbalanced data too. Said so explicitly, so that
        # a reader used to Fleiss' "approximation (unbalanced data)" does not
        # assume the same here.
        exact = (
            "exact formula" if self.balance_profile.is_balanced
            else "exact formula (unbalanced data)"
        )
        return (
            f"{self._subjects_phrase()}; {self._raters_phrase()}; "
            f"{self._weights_phrase(self.weights_type)}; {exact}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and alpha-specific metadata.

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
        Whether this is the exact formula or an approximation.

        Always returns True for Krippendorff's alpha. Alpha handles
        unbalanced data natively via pair-pooling -- the unbalanced
        formula is the formula, not an approximation. The balance
        profile is available for informational context but does not
        affect the correctness of the result.
        """
        return True
