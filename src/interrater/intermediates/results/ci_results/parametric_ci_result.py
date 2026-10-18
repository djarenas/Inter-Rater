"""Result for a confidence interval computed from a variance formula."""

from __future__ import annotations

import warnings
from dataclasses import dataclass

from .ci_result import CIResult

_ERROR_PREFIX = "[ParametricCIResult] Error"
_WARNING_PREFIX = "[ParametricCIResult] Warning"

# What the half-width was built from. Only 'null' exists today: the variance of
# the coefficient under H0, which is what the published formulas give. An
# unconstrained variance -- estimated at the observed value rather than at zero
# -- is the natural next entry, and is not implemented.
_VALID_VARIANCE_SOURCES = ("null",)


@dataclass(frozen=True)
class ParametricCIResult(CIResult):
    """
    I am an interval from a closed-form variance rather than a resampling.

    I have no resamples, so nothing about me can fail partway: either the
    formula applied or it did not. That is why I carry no counts, no method
    pair and no fallback -- there is no second method for me to fall back to,
    and nothing for me to disclose about a choice I never made.

    **Read my width carefully, and I will say so every time I print.** The
    variance I am built from is the *null* variance: the variance of the
    coefficient under H0, that is, on the assumption the true value is zero.
    That is the right quantity for a hypothesis test, which asks how far from
    zero an estimate is in units of its spread under zero. It is a mismatched
    quantity for a confidence interval, which should use the spread of the
    estimator at the value actually observed. An estimate far from zero has a
    different variance than it would have at zero, and I do not know by how
    much.

    So I am centred on the point estimate but scaled as though the point
    estimate were zero. I am kept because that is exactly what several
    published formulas and several other packages report, which makes me the
    right object for reproducing a published number or comparing against
    another implementation. I am the wrong object for reporting the precision
    of an estimate, and a bootstrap interval is what that calls for.

    **How far off is it?** Measured against percentile bootstrap intervals on
    this package's own datasets, my width runs from about 56% too narrow
    (Fleiss' kappa on Fleiss 1971 and on SR1 prelim) to 29% too wide (Cohen's
    kappa on Gwet cont3x3). The direction depends on the coefficient, so
    there is no correction factor to apply and no safe direction to assume.

    **I warn on construction, and again every time I print.** Unlike
    BootstrapCIResult, which stays quiet about its resample count because one
    pairwise analysis builds C(n,2) of me and fifteen identical warnings
    teach people to filter warnings, there is exactly one of me per analysis,
    and I only exist because somebody asked for a non-default method. The
    config warns too, earlier, about the parameter -- that one says where the
    knob is, this one says how wrong the number is. Different facts are
    useful at the two moments.

    Attributes
    ----------
    variance_source : str
        Which variance produced my half-width. 'null' is the variance under
        H0: coefficient = 0.

    See Also
    --------
    CIResult : my base, which holds the bounds and the level.
    BootstrapCIResult : the resampling alternative, which estimates the
        spread at the observed value rather than at zero.
    """

    variance_source: str

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction, before I freeze."""

        if self.variance_source not in _VALID_VARIANCE_SOURCES:
            raise ValueError(
                f"{_ERROR_PREFIX}: variance_source must be one of "
                f"{list(_VALID_VARIANCE_SOURCES)}, got {self.variance_source!r}"
            )

        # Parent validates the bounds and the level, and coerces them
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Last, so that an invalid interval raises rather than warning first.
        if self.variance_source == "null":
            warnings.warn(
                f"{_WARNING_PREFIX}: this interval is centred on your "
                f"estimate but scaled as though that estimate were zero -- "
                f"its width comes from the variance under H0: kappa = 0. "
                f"Measured against percentile bootstrap intervals on this "
                f"package's datasets, that width runs from about 56% too "
                f"narrow (Fleiss) to 29% too wide (Cohen), and the direction "
                f"depends on the coefficient. Treat it as a reproduction of "
                f"the published formula, not as this estimate's uncertainty.",
                UserWarning,
                stacklevel=5,
            )

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Return a compact string representation of the interval.

        The caveat travels with the number rather than living only in the
        documentation, because the number is what gets copied into a table:

            [0.401, 0.639] (95% CI, null variance -- width assumes kappa = 0,
                            reproduction only)

        I stay on one line. Whoever holds me prints me inside their own
        single-line description, so the full account belongs in the warning
        I raise at construction, not here. What survives the compression is
        the fact that decides whether a reader should use the number at all.
        """
        bounds = (
            f"[{self.lower:.3f}, {self.upper:.3f}]"
            if self.is_defined else "[undefined]"
        )
        caveat = (
            "null variance -- width assumes kappa = 0, reproduction only"
            if self.variance_source == "null" else self.variance_source
        )
        return f"{bounds} ({self.confidence_level:.0%} CI, {caveat})"

    def summary(self) -> dict:
        """
        Return my bounds, level, and variance source in plain types.

        Returns
        -------
        dict
            Keys: 'lower', 'upper' (float), 'alpha', 'confidence_level'
            (float), 'is_defined' (bool), 'variance_source' (str).
        """
        return {
            "lower":            self.lower,
            "upper":            self.upper,
            "alpha":            self.alpha,
            "confidence_level": self.confidence_level,
            "is_defined":       self.is_defined,
            "variance_source":  self.variance_source,
        }
