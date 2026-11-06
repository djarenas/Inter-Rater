"""Result for Light's kappa analysis."""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations
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

_ERROR_PREFIX = "[LightKappaResult] Error"
_TOL = 1e-9

# Why parts and kappa_max are not computed for me. The analyzer passes them in.
PARTS_NOT_COMPUTED = (
    "Light's kappa is the mean of the pair kappas, not one Po and Pe; each "
    "pair's are in pair_parts"
)
KAPPA_MAX_NOT_COMPUTED = (
    "the pairs may not reach their ceilings together, so the mean of the "
    "pair ceilings is an upper bound, not the maximum; each pair's ceiling "
    "is in pair_kappa_max"
)


@dataclass(frozen=True)
class LightKappaResult(PointEstimateResult):
    """
    I am Light's kappa -- the mean of all C(n, 2) pairwise Cohen's kappas.

    I inherit point, confidence_interval, bootstrap_distribution, and
    balance_profile from PointEstimateResult.

    I am the headline number: one summary of average pairwise agreement.
    Spread and per-rater diagnostics belong to the pairwise permutation
    analyzer, which shows which raters are the problem. I am what you report
    when you want a single number.

    I exist only for complete data, where every rater rated every subject
    (Light 1971), so my rater profile must have a fill rate of 1.

    I am not one Po and Pe, so my parts and kappa_max are NotComputed.
    Instead I carry each pair's: pair_parts holds every pair's Po and Pe,
    and pair_kappa_max its Cohen ceiling, both keyed by the pair's rater
    names. I check that they cover every pair of the raters in my agreement
    profile, that each pair's kappa is within its ceiling, and that the mean
    of the pair kappas (skipping undefined ones, as Light's kappa does) is
    my point estimate.

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
    confidence_interval : CIResult | None
        Bootstrap confidence interval (inherited).
    bootstrap_distribution : np.ndarray | NotComputed | None
        Bootstrap distribution (inherited); NotComputed when the draws were
        not kept.
    balance_profile : BalanceProfile
        Balance structure of the data (inherited).
    rater_profile : RaterProfile
        The raters (inherited). Always a RaterProfile: Light's kappa needs
        rater identity.
    parts : NotComputed
        Always NotComputed, with PARTS_NOT_COMPUTED as the reason (inherited).
    kappa_max : NotComputed
        Always NotComputed, with KAPPA_MAX_NOT_COMPUTED as the reason
        (inherited).
    agreement_profile : AgreementProfile
        Prevalence, each rater's marginals, specific agreement (inherited).
    n_pairs : int
        Number of rater pairs: C(n_raters, 2).
    ci_method : str
        Method used: "percentile" or "bca".
    pair_parts : dict[tuple[str, str], AgreementParts]
        Each pair's Po and Pe, unweighted, keyed by (rater_a, rater_b) in
        declared order. NaN terms for a pair sharing fewer than two subjects.
    pair_kappa_max : dict[tuple[str, str], float]
        Each pair's Cohen ceiling, same keys. NaN where the pair's kappa is.
    """

    # Attribute declaration

    # What I am called: the config name, my printed name and my symbol.
    # Every other place that names this coefficient reads these.
    COEFFICIENT: ClassVar[str] = "light_kappa"
    NAME:        ClassVar[str] = "Light's κ"
    SYMBOL:      ClassVar[str] = "κ"

    # Inherited from parent class PointEstimateResult
    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    parts: NotComputed
    kappa_max: NotComputed
    agreement_profile: AgreementProfile

    # Populated at this level
    n_pairs: int
    ci_method: str
    pair_parts: dict
    pair_kappa_max: dict

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Light's kappa result attributes."""

        # Validate n_pairs
        if not isinstance(self.n_pairs, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_pairs must be int, "
                f"got {type(self.n_pairs).__name__}"
            )
        from ...data_profiles.rater_profile import RaterProfile

        if not isinstance(self.rater_profile, RaterProfile):
            raise ValueError(
                f"{_ERROR_PREFIX}: rater_profile must be a RaterProfile for "
                f"Light's kappa (rater identity is needed to compute "
                f"pairwise kappas), got {self.rater_profile}"
            )
        expected_pairs = self.rater_profile.n_raters * (self.rater_profile.n_raters - 1) // 2
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

        self._check_pairs()

    def _check_pairs(self) -> None:
        """Complete data; every pair's parts and ceiling; their mean is my point."""
        if self.rater_profile.fill_rate != 1.0:
            raise ValueError(
                f"{_ERROR_PREFIX}: Light's kappa is for complete data, where "
                f"every rater rated every subject; rater_profile has a fill "
                f"rate of {self.rater_profile.fill_rate}"
            )
        for name in ("parts", "kappa_max"):
            if not isinstance(getattr(self, name), NotComputed):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must be NotComputed for Light's "
                    f"kappa; each pair's is in pair_{name}"
                )
        marginals = self.agreement_profile.rater_marginals
        if isinstance(marginals, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: agreement_profile.rater_marginals must be "
                f"computed for Light's kappa (rater identity required)"
            )
        pairs = list(combinations(marginals, 2))
        for name in ("pair_parts", "pair_kappa_max"):
            value = getattr(self, name)
            if not isinstance(value, dict):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must be a dict, "
                    f"got {type(value).__name__}"
                )
            if list(value) != pairs:
                raise ValueError(
                    f"{_ERROR_PREFIX}: {name} must have one entry for each "
                    f"pair of the raters in the agreement profile, in order: "
                    f"{pairs}, got {list(value)}"
                )

        kappas = []
        for pair in pairs:
            parts = self.pair_parts[pair]
            if not isinstance(parts, AgreementParts) or parts.kind != "agreement" or parts.weighted:
                raise TypeError(
                    f"{_ERROR_PREFIX}: pair_parts[{pair}] must be unweighted "
                    f"agreement parts (Po and Pe), got {parts!r}"
                )
            ceiling = self.pair_kappa_max[pair]
            if isinstance(ceiling, bool) or not isinstance(ceiling, (int, float)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: pair_kappa_max[{pair}] must be a "
                    f"number, got {type(ceiling).__name__}"
                )
            kappa = parts.coefficient()
            if math.isnan(kappa):
                continue
            if not math.isnan(ceiling) and kappa > ceiling + _TOL:
                raise ValueError(
                    f"{_ERROR_PREFIX}: pair {pair} has kappa {kappa:.6f} "
                    f"above its ceiling {ceiling:.6f}"
                )
            kappas.append(kappa)

        mean = float(np.mean(kappas)) if kappas else float("nan")
        if math.isnan(mean) or math.isnan(self.point):
            same = math.isnan(mean) and math.isnan(self.point)
        else:
            same = abs(mean - self.point) <= _TOL
        if not same:
            raise ValueError(
                f"{_ERROR_PREFIX}: the mean of the pair kappas is "
                f"{mean:.6f}, but the point estimate is {self.point:.6f}"
            )

    # Public API ------------------------------------------------------------- #

    # Printed form: see PointEstimateResult.__str__ ----------------------- #

    def _made_of_line(self) -> str:
        """The pair kappas the mean is made of, and where each pair's terms are."""
        kappas = [
            k for k in (p.coefficient() for p in self.pair_parts.values())
            if not math.isnan(k)
        ]
        where = "each pair's Po and Pe in pair_parts, κ_max in pair_kappa_max"
        if not kappas:
            return f"no pair kappa is defined; {where}"
        of = (
            f"{len(kappas)} pair kappas" if len(kappas) == self.n_pairs
            else f"{len(kappas)} of {self.n_pairs} pair kappas (the rest undefined)"
        )
        low, high = f"{min(kappas):.2f}", f"{max(kappas):.2f}"
        spread = f"all {low}" if low == high else f"{low} to {high}"
        return f"mean of {of}, {spread}; {where}"

    def _data_line(self) -> str:
        return (
            f"{self.balance_profile.n_subjects:,} subjects; "
            f"{self._raters_phrase()}, every rater rated every subject"
        )

    def summary(self) -> dict:
        """
        Return point, CI, bootstrap diagnostics, and Light-specific metadata.

        Returns
        -------
        dict
            Keys from parent (including balance) plus 'n_pairs',
            'ci_method', and 'pairs': one dict per pair with 'rater_a',
            'rater_b', its parts and its 'kappa_max'.
        """
        result = super().summary()
        result.update({
            "n_pairs": self.n_pairs,
            "ci_method": self.ci_method,
            "pairs": [
                {
                    "rater_a": a,
                    "rater_b": b,
                    **self.pair_parts[(a, b)].summary(),
                    "kappa_max": float(self.pair_kappa_max[(a, b)]),
                }
                for a, b in self.pair_parts
            ],
        })
        return result
