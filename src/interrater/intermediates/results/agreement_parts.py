"""The two numbers a chance-corrected coefficient is built from."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

_ERROR_PREFIX = "[AgreementParts] Error"
_TOL = 1e-9
_KINDS = ("agreement", "disagreement")


@dataclass(frozen=True)
class AgreementParts:
    """
    I hold the observed and chance terms a coefficient is built from, and
    rebuild the coefficient from them.

    One number hides how it was reached: kappa = 0.40 can be Po = 0.95 with
    Pe = 0.92 (raters agree on nearly everything, but one category dominates)
    or Po = 0.70 with Pe = 0.50. I keep both terms, so a reader can tell.

    Two kinds:

    - "agreement" (most coefficients): observed agreement Po and chance
      agreement Pe; the coefficient is (Po - Pe) / (1 - Pe).
    - "disagreement" (Krippendorff's alpha): observed disagreement Do and
      expected disagreement De; the coefficient is 1 - Do / De.

    The result holding me checks that its point estimate equals
    ``coefficient()``, so the terms cannot disagree with the number they
    explain.

    Attributes
    ----------
    observed : float
        Po or Do. NaN when it cannot be computed (no subject rated twice).
    chance : float
        Pe or De. NaN under the same condition.
    kind : str
        "agreement" or "disagreement".
    weighted : bool
        Whether a weight matrix gave partial credit to near misses. Weighted
        terms are a different quantity from unweighted ones, so the label
        travels with them.
    """

    observed: float
    chance: float
    kind: str
    weighted: bool

    def __post_init__(self) -> None:
        if self.kind not in _KINDS:
            raise ValueError(
                f"{_ERROR_PREFIX}: kind must be 'agreement' or "
                f"'disagreement', got {self.kind!r}"
            )
        if not isinstance(self.weighted, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: weighted must be bool, "
                f"got {type(self.weighted).__name__}"
            )
        for name in ("observed", "chance"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must be a number, "
                    f"got {type(value).__name__}"
                )
            if math.isnan(value):
                continue
            if self.kind == "agreement" and not (-_TOL <= value <= 1 + _TOL):
                raise ValueError(
                    f"{_ERROR_PREFIX}: {self._labels()[name]} must be in "
                    f"[0, 1], got {value}"
                )
            if self.kind == "disagreement" and value < -_TOL:
                raise ValueError(
                    f"{_ERROR_PREFIX}: {self._labels()[name]} must be >= 0, "
                    f"got {value}"
                )

    # Public API ------------------------------------------------------------- #

    def coefficient(self) -> float:
        """
        The coefficient these terms make.

        NaN when they make none: Pe = 1 (every agreement expected by chance)
        or De = 0 (no disagreement expected), or either term NaN. "= 1" and
        "= 0" use ``np.isclose``, the test the coefficient functions use, so
        a Pe a hair below 1 is NaN here exactly when it is NaN there.
        """
        o, c = float(self.observed), float(self.chance)
        if math.isnan(o) or math.isnan(c):
            return float("nan")
        if self.kind == "agreement":
            return float("nan") if np.isclose(c, 1.0) else (o - c) / (1.0 - c)
        return float("nan") if np.isclose(c, 0.0) else 1.0 - o / c

    def reproduces(self, point: float, tol: float = 1e-9) -> bool:
        """
        Whether ``point`` is the coefficient these terms make (NaN with NaN).
        """
        rebuilt = self.coefficient()
        if math.isnan(rebuilt) or math.isnan(point):
            return math.isnan(rebuilt) and math.isnan(point)
        return abs(rebuilt - point) <= tol * max(1.0, abs(point))

    def __str__(self) -> str:
        labels = self._labels()
        prefix = "weighted " if self.weighted else ""
        return (
            f"{prefix}{labels['observed']} = {self.observed:.3f}, "
            f"{prefix}{labels['chance']} = {self.chance:.3f}"
        )

    def summary(self) -> dict:
        """Return my record in plain types."""
        return {
            "kind": self.kind,
            "observed": float(self.observed),
            "chance": float(self.chance),
            "weighted": self.weighted,
        }

    def _labels(self) -> dict:
        if self.kind == "agreement":
            return {"observed": "Po", "chance": "Pe"}
        return {"observed": "Do", "chance": "De"}
