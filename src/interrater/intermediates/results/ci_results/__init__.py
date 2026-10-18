"""Confidence intervals: the bounds, and an account of how they were produced.

Every class here descends from ``CIResult``, which holds the two bounds and
the level they were taken at. The subclasses differ in what can go wrong on
the way to them, and each describes itself completely -- so whatever holds an
interval can print it without knowing which kind it is.

``BootstrapCIResult`` comes from resampling, and so has failures to report and
a method that may not be the one that was asked for. ``ParametricCIResult``
comes from a variance formula, has nothing to fall back to, and carries a
standing caveat about the variance it was built from.
"""

from .ci_result import CIResult
from .bootstrap_ci_result import BootstrapCIResult
from .parametric_ci_result import ParametricCIResult

__all__ = [
    "CIResult",
    "BootstrapCIResult",
    "ParametricCIResult",
]
