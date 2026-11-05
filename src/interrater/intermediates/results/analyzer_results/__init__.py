"""Single-coefficient results: one point estimate, one interval, one dataset.

Every class here descends from ``PointEstimateResult``, which is the shape an
analyzer returns when it has estimated one coefficient. ``IndividualMarginalResult``
is deliberately absent from ``__all__``: it is an internal base carrying the
invariant that Cohen and Conger share, not a result any analyzer produces.
"""

from .point_estimate_result import PointEstimateResult
from .fleiss_kappa_result import FleissKappaResult
from .scotts_pi_result import ScottsPiResult
from .gwets_ac1_result import GwetsAC1Result
from .brennan_prediger_result import BrennanPredigerResult
from .krippendorff_alpha_result import KrippendorffAlphaResult
from .cohen_kappa_result import CohenKappaResult
from .conger_kappa_result import CongerKappaResult
from .light_kappa_result import LightKappaResult

# Config name -> the coefficient's result class, built from each class's own
# COEFFICIENT. The one place a name like "fleiss_kappa" is looked up, for the
# result class and for what the coefficient is called (NAME, SYMBOL).
COEFFICIENT_RESULTS: dict[str, type[PointEstimateResult]] = {
    cls.COEFFICIENT: cls
    for cls in (
        FleissKappaResult,
        ScottsPiResult,
        GwetsAC1Result,
        BrennanPredigerResult,
        KrippendorffAlphaResult,
        CohenKappaResult,
        CongerKappaResult,
        LightKappaResult,
    )
}


def result_class_for(coefficient: str) -> type[PointEstimateResult]:
    """The result class for a config name such as "fleiss_kappa"; raises if unknown."""
    try:
        return COEFFICIENT_RESULTS[coefficient]
    except KeyError:
        raise ValueError(
            f"[analyzer_results] Error: unknown coefficient {coefficient!r}. "
            f"Known: {sorted(COEFFICIENT_RESULTS)}."
        ) from None


__all__ = [
    "COEFFICIENT_RESULTS",
    "result_class_for",
    "PointEstimateResult",
    "FleissKappaResult",
    "ScottsPiResult",
    "GwetsAC1Result",
    "BrennanPredigerResult",
    "KrippendorffAlphaResult",
    "CohenKappaResult",
    "CongerKappaResult",
    "LightKappaResult",
]
