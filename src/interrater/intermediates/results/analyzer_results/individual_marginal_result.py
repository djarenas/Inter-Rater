"""Result base for individual marginal coefficients."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ...not_computed import NotComputed
from .point_estimate_result import PointEstimateResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile
    from ...data_profiles.rater_profile import RaterProfile
    from ..ci_results.ci_result import CIResult

_ERROR_PREFIX = "[IndividualMarginalResult] Error"
_TOL = 1e-9


def _same(a: float, b: float) -> bool:
    """Equal within tolerance, NaN equal only to NaN."""
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    return abs(a - b) <= _TOL


def _unweighted_pair_kappa_max(p_a: np.ndarray, p_b: np.ndarray) -> float:
    """
    Cohen's (1960) ceiling for one pair, from their marginals: Po can be at
    most Sum_j min(p_aj, p_bj), and Pe is Sum_j p_aj p_bj.
    """
    pe = float((p_a * p_b).sum())
    if np.isclose(pe, 1.0):
        return float("nan")
    return (float(np.minimum(p_a, p_b).sum()) - pe) / (1.0 - pe)


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
    rater_profile must be a RaterProfile, not NotComputed.

    For the same reason my agreement profile must carry each rater's
    marginals, one row per rater in my rater profile: they are what my
    children's chance terms are built from, and they check their parts
    against them.

    I exist so that these guards are written once instead of repeated
    on every leaf. Downstream code can also check
    isinstance(result, IndividualMarginalResult) to know the coefficient
    came from rater-identity-aware analysis.

    Attributes
    ----------
    point : float
        The coefficient value (inherited from PointEstimateResult).
    confidence_interval : CIResult | None
        The interval, bootstrap or parametric (inherited).
    bootstrap_distribution : np.ndarray | NotComputed | None
        Bootstrap distribution (inherited); NotComputed when the draws were
        not kept.
    balance_profile : BalanceProfile
        Balance structure of the data (inherited). Must not be None for
        individual marginal coefficients.
    rater_profile : RaterProfile
        The raters (inherited). Must be a RaterProfile, not NotComputed:
        rater identity is always known for these coefficients.
    """

    # Attribute declaration

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    balance_profile: BalanceProfile
    rater_profile: RaterProfile

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate individual marginal result constraints.

        Enforces that balance_profile is not None and rater_profile is a
        RaterProfile -- individual marginal coefficients always require
        rater identity.
        """

        # Individual marginal coefficients always require the raters
        if isinstance(self.rater_profile, NotComputed):
            raise ValueError(
                f"{_ERROR_PREFIX}: rater_profile must be a RaterProfile "
                f"for individual marginal coefficients (rater identity "
                f"required), got {self.rater_profile}"
            )

        # Call parent validation (validates point, confidence_interval,
        # bootstrap_distribution, balance_profile type)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Each rater's marginals, one row per rater
        marginals = self.agreement_profile.rater_marginals
        if isinstance(marginals, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: agreement_profile.rater_marginals must be "
                f"computed for individual marginal coefficients (rater "
                f"identity required), got {marginals}"
            )
        if len(marginals) != self.rater_profile.n_raters:
            raise ValueError(
                f"{_ERROR_PREFIX}: agreement_profile has marginals for "
                f"{len(marginals)} raters but rater_profile has "
                f"{self.rater_profile.n_raters}"
            )

    def _marginal_rows(self) -> np.ndarray:
        """Each rater's marginals from my profile, shape (n_raters, K)."""
        return np.array([
            list(row.values())
            for row in self.agreement_profile.rater_marginals.values()
        ], dtype=float)
