"""Result for a confidence interval computed from a variance formula."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import ClassVar

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
    formula applied or it did not. I carry no counts and cannot fall back:
    my method is always parametric_null, and when the formula gives no
    bounds, no method ran.

    **Read my width carefully, and I will say so every time I print.** The
    variance I am built from is the *null* variance: the variance of the
    coefficient under H0, that is, on the assumption the true value is zero.
    That is the right quantity for a hypothesis test, which asks how far from
    zero an estimate is in units of its spread under zero. It is a mismatched
    quantity for a confidence interval, which should use the spread of the
    estimator at the value actually observed. An estimate far from zero has a
    different variance than it would have at zero, and I do not know by how
    much.

    The package keeps me as an option to reproduce published analyses that
    build an interval from the null variance, for coefficients such as
    Fleiss' kappa. I am the wrong object for reporting the precision
    of an estimate, and a bootstrap interval is what that calls for.

    **How far off is it?** Against percentile bootstrap intervals on this
    package's datasets, my width runs from about 60% too narrow (Fleiss'
    kappa on SR1 prelim) to more than twice as wide (quadratic-weighted
    Cohen's kappa on Gwet cont3x3). The direction depends on the data, not
    only on the coefficient: Fleiss' kappa is too narrow on most datasets and
    too wide on Gwet cont3x3. So there is no correction factor to apply and
    no safe direction to assume.

    I warn once when I am built, and my printout carries the caveat. The
    config warns earlier, when parametric_null is chosen. There is one of me
    per analysis, so the warning does not repeat.

    Attributes
    ----------
    variance_source : str
        Which variance produced my half-width. 'null' is the variance under
        H0: coefficient = 0.
    variance_reference : str
        Where the variance formula comes from, e.g. "Fleiss, Nee & Landis
        (1979), Eq. 12". I print it.

    See Also
    --------
    CIResult : my base, which holds the bounds and the level.
    BootstrapCIResult : the resampling alternative, which estimates the
        spread at the observed value rather than at zero.
    """

    # The method I can report
    _METHODS: ClassVar[frozenset[str]] = frozenset({"parametric_null"})

    # Inherited from and populated at CIResult parent class
    lower: float
    upper: float
    alpha: float
    method_intended: str
    method_used: str | None

    # Populated at this level
    variance_source: str
    variance_reference: str

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction, before I freeze."""

        # Parent validates and populates lower, upper, and alpha attributes
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc


        if self.variance_source not in _VALID_VARIANCE_SOURCES:
            raise ValueError(
                f"{_ERROR_PREFIX}: variance_source must be one of "
                f"{list(_VALID_VARIANCE_SOURCES)}, got {self.variance_source!r}"
            )
        if not isinstance(self.variance_reference, str) or not self.variance_reference.strip():
            raise ValueError(
                f"{_ERROR_PREFIX}: variance_reference must name the formula's "
                f"source, e.g. 'Fleiss, Nee & Landis (1979), Eq. 12'; got "
                f"{self.variance_reference!r}"
            )

        # Last, so that an invalid interval raises rather than warning first.
        if self.variance_source == "null":
            warnings.warn(
                f"{_WARNING_PREFIX}: this interval is centred on your "
                f"estimate but scaled as though that estimate were zero -- "
                f"its width comes from the variance under H0: kappa = 0. "
                f"Measured against percentile bootstrap intervals on this "
                f"package's datasets, that width runs from about 60% too "
                f"narrow to more than twice as wide, and the direction "
                f"depends on the data. Treat it as a reproduction of "
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

            [0.401, 0.639] (95% CI, null variance of Fleiss, Nee & Landis
                            (1979), Eq. 12 -- width assumes kappa = 0,
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
            f"null variance of {self.variance_reference} -- width assumes "
            f"kappa = 0, reproduction only"
            if self.variance_source == "null"
            else f"{self.variance_source} variance of {self.variance_reference}"
        )
        return f"{bounds} ({self.confidence_level:.0%} CI, {caveat})"

    def summary(self) -> dict:
        """
        Return my bounds, level, and variance source in plain types.

        Returns
        -------
        dict
            Keys: 'lower', 'upper' (float), 'alpha', 'confidence_level'
            (float), 'is_defined' (bool), 'method_intended' (str),
            'method_used' (str or None), 'variance_source',
            'variance_reference' (str).
        """
        return {
            "lower":            self.lower,
            "upper":            self.upper,
            "alpha":            self.alpha,
            "confidence_level": self.confidence_level,
            "is_defined":       self.is_defined,
            "method_intended":  self.method_intended,
            "method_used":      self.method_used,
            "variance_source":  self.variance_source,
            "variance_reference": self.variance_reference,
        }
