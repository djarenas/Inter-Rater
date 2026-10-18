"""Result base for individual marginal coefficients."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[IndividualMarginalResult] Error"


@dataclass(frozen=True)
class IndividualMarginalResult(PointEstimateResult):
    """
    I am the base for results of individual marginal coefficients --
    coefficients that use each rater's own category distribution rather than
    pooling marginals across all raters.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult.

    Individual marginal coefficients (Cohen, Conger, Light) require rater
    identity: you must know which rater gave which rating to compute
    per-rater marginals. This means they can only be built from Ratings,
    never from an NijMatrix (which discards rater identity). I enforce
    that constraint here: balance_profile must not be None, and
    balance_profile.n_raters must not be None.

    I exist so that these two guards are written once instead of repeated
    on every leaf. Downstream code can also check
    isinstance(result, IndividualMarginalResult) to know the coefficient
    came from rater-identity-aware analysis.

    Attributes
    ----------
    point : float
        The coefficient value (inherited from PointEstimateResult).
    confidence_interval : tuple[float, float] | None
        Bootstrap confidence interval (lower, upper) (inherited).
    bootstrap_distribution : np.ndarray | None
        Bootstrap distribution (inherited).
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Must not be None for
        individual marginal coefficients. Its n_raters must also not be
        None (rater identity is always known).
    """

    # Attribute declaration

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: tuple[float, float] | None
    bootstrap_distribution: np.ndarray | None
    balance_profile: BalanceProfile | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate individual marginal result constraints.

        Enforces that balance_profile and balance_profile.n_raters are
        not None -- individual marginal coefficients always require
        rater identity.
        """

        # Individual marginal coefficients always require a balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None for "
                f"individual marginal coefficients (rater identity required)."
            )

        # Individual marginal coefficients always require known rater count
        if self.balance_profile.n_raters is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile.n_raters must not be None "
                f"for individual marginal coefficients (rater identity required)"
            )

        # Call parent validation (validates point, confidence_interval,
        # bootstrap_distribution, balance_profile type)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc
