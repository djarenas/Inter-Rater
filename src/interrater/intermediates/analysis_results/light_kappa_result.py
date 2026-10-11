"""Result for Light's kappa analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import numpy as np

from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[LightKappaResult] Error"


@dataclass(frozen=True)
class LightKappaResult(PointEstimateResult):
    """
    I am Light's kappa -- the mean of all C(n, 2) pairwise Cohen's kappas.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult.

    I am the headline number: one summary of average pairwise agreement. I do
    not carry per-pair kappas, spread, or per-rater diagnostics -- those belong
    to the pairwise permutation analyzer, which shows which raters are the
    problem. I am what you report when you want a single number.

    At exactly two raters (n_pairs=1) I reduce to Cohen's unweighted kappa.

    I am always unweighted. Light (1971) defined the coefficient as the mean
    of unweighted pairwise kappas. I carry no weights_type field.

    References
    ----------
    Light, R. J. (1971). Measures of response agreement for qualitative
        data: Some generalizations and alternatives. Psychological
        Bulletin, 76(5), 365--377.

    Attributes
    ----------
    point : float
        Light's kappa (inherited from PointEstimateResult).
    confidence_interval : tuple[float, float] | None
        Bootstrap confidence interval (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited).
    n_pairs : int
        Number of rater pairs: C(n_raters, 2).
    ci_method : str
        Method used: "percentile" or "bca".
    """

    # Attribute declaration

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: Optional[tuple[float, float]]
    bootstrap_distribution: Optional[np.ndarray]
    balance_profile: Optional[BalanceProfile]

    # Populated at this level
    n_pairs: int
    ci_method: str

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Light's kappa result attributes."""

        # Light's kappa always requires a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for "
                f"Light's kappa"
            )

        # Validate n_pairs
        if not isinstance(self.n_pairs, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_pairs must be int, "
                f"got {type(self.n_pairs).__name__}"
            )
        if self.balance_profile.n_raters is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile.n_raters must not be "
                f"None for Light's kappa (rater identity is needed to "
                f"compute pairwise kappas)"
            )
        expected_pairs = self.balance_profile.n_raters * (self.balance_profile.n_raters - 1) // 2
        if self.n_pairs != expected_pairs:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_pairs must be C(n_raters, 2) = "
                f"{expected_pairs}, got {self.n_pairs}"
            )

        # Validate ci_method (no parametric_null for Light)
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or 'bca', "
                f"got '{self.ci_method}'. Light's kappa has no closed-form "
                f"null variance (pairwise kappas are correlated), so "
                f"'parametric_null' is not available."
            )

        # Call parent validation (validates point, confidence_interval,
        # bootstrap_distribution, balance_profile type)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Return a compact string representation with coefficient name.

        Format:
            Light's κ = 0.587 (95% CI: [0.521, 0.653], n_bootstrap=1000)
            n_subjects=80; n_raters=4; n_pairs=6; ci_method=percentile
        """
        parent_str = super().__str__()
        return (
            f"Light's κ {parent_str}\n"
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"n_raters={self.balance_profile.n_raters}; "
            f"n_pairs={self.n_pairs}; ci_method={self.ci_method}"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and Light-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'n_pairs', 'ci_method'.
        """
        result = super().summary()
        result.update({
            "n_pairs": self.n_pairs,
            "ci_method": self.ci_method,
        })
        return result
