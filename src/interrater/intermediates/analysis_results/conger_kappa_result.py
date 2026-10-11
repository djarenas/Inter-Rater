"""Result for Conger's kappa analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .individual_marginal_result import IndividualMarginalResult

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[CongerKappaResult] Error"


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

    References
    ----------
    Conger, A. J. (1980). Integration and generalization of kappas for
        multiple raters. Psychological Bulletin, 88(2), 322--328.

    Attributes
    ----------
    point : float
        Conger's kappa (inherited from PointEstimateResult).
    confidence_interval : tuple[float, float] | None
        Bootstrap confidence interval (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited).
    ci_method : str
        Method used: "percentile" or "bca".
    """

    # Attribute declaration

    # Inherited from IndividualMarginalResult (which validates
    # balance_profile is not None and n_raters is not None)
    point: float
    confidence_interval: tuple[float, float] | None
    bootstrap_distribution: np.ndarray | None
    balance_profile: BalanceProfile | None

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

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Return a compact string representation with coefficient name.

        Format:
            Conger's κ = 0.612 (95% CI: [0.551, 0.673], n_bootstrap=1000)
            n_subjects=80; n_raters=4; ci_method=percentile
        """
        parent_str = super().__str__()
        return (
            f"Conger's κ {parent_str}\n"
            f"n_subjects={self.balance_profile.n_subjects}; "
            f"n_raters={self.balance_profile.n_raters}; "
            f"ci_method={self.ci_method}"
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
