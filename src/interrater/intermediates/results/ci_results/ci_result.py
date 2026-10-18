"""Base class for all confidence interval results."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

_ERROR_PREFIX = "[CIResult] Error"


@dataclass(frozen=True)
class CIResult(ABC):
    """
    I am the base for every confidence interval this package produces.

    I hold the three things every interval has whatever produced it: two
    bounds and the level they were taken at. How they were arrived at --
    resampling, a variance formula, something later -- belongs to my
    subclasses, and so does the account of what could go wrong on the way.

    I carry alpha because an interval without its level is ambiguous, and
    because the thing holding me should not have to consult a config to
    print me.

    I am a frozen dataclass rather than a bare ABC like ``AnalysisResult``,
    because unlike that one I have state of my own to share: my subclasses
    differ in what they add, not in what they all are.

    **I may have no bounds.** When whatever was meant to compute them could
    not run, both are NaN and :attr:`is_defined` is False. That is an outcome
    rather than an error: the point estimate an interval accompanies is
    usually perfectly well defined even when the interval is not, so refusing
    to exist would destroy the useful half of the answer.

    I know nothing about the point estimate and want to know nothing. A
    strongly bias-corrected interval can sit entirely to one side of it,
    which is correct behaviour rather than a contradiction, so whether the
    two are consistent is a question for whoever holds both of us.

    Attributes
    ----------
    lower : float
        Lower bound. NaN when no interval could be computed.
    upper : float
        Upper bound. NaN when no interval could be computed.
    alpha : float
        Significance level the bounds were taken at, e.g. 0.05 for a 95%
        interval.
    """

    lower: float
    upper: float
    alpha: float

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """
        Validate my attributes after construction, before I freeze.

        Structural only, and only what is true of every interval. Each check
        describes something no correct computation could have produced, so
        reaching one means a bug upstream rather than awkward data.
        """

        # Bounds must be floats
        for name, value in (("lower", self.lower), ("upper", self.upper)):
            if not isinstance(value, (int, float, np.floating)) or isinstance(
                value, bool
            ):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must be float, "
                    f"got {type(value).__name__}"
                )

        # alpha must be a float strictly inside (0, 1)
        if not isinstance(self.alpha, (int, float, np.floating)) or isinstance(
            self.alpha, bool
        ):
            raise TypeError(
                f"{_ERROR_PREFIX}: alpha must be float, "
                f"got {type(self.alpha).__name__}"
            )
        if not (0 < self.alpha < 1):
            raise ValueError(
                f"{_ERROR_PREFIX}: alpha must be in (0, 1), got {self.alpha}"
            )

        lower_defined = bool(np.isfinite(self.lower))
        upper_defined = bool(np.isfinite(self.upper))

        # Half an interval is not a weaker interval, it is a broken one
        if lower_defined != upper_defined:
            raise ValueError(
                f"{_ERROR_PREFIX}: bounds must both be defined or both be "
                f"undefined, got ({self.lower}, {self.upper})"
            )
        if lower_defined and self.lower > self.upper:
            raise ValueError(
                f"{_ERROR_PREFIX}: lower bound {self.lower} > upper bound "
                f"{self.upper}"
            )

        # These values arrive from numpy arithmetic, so they arrive as numpy
        # scalars. Coerce once, here, rather than at each place that reads me:
        # a caller reaching for my bounds should get a float, not a float that
        # prints as np.float64(0.4) and carries numpy's semantics with it.
        object.__setattr__(self, "lower", float(self.lower))
        object.__setattr__(self, "upper", float(self.upper))
        object.__setattr__(self, "alpha", float(self.alpha))

    # Public API ------------------------------------------------------------- #

    @property
    def is_defined(self) -> bool:
        """Whether there is an interval at all."""
        return bool(np.isfinite(self.lower) and np.isfinite(self.upper))

    @property
    def confidence_level(self) -> float:
        """1 - alpha, for display. 0.05 reads as 0.95."""
        return 1.0 - self.alpha

    def as_tuple(self) -> tuple[float, float]:
        """My bounds alone, for callers that want only the interval."""
        return (self.lower, self.upper)

    @abstractmethod
    def summary(self) -> dict:
        """
        Return everything I represent, in plain Python types.

        Subclasses include my bounds and level alongside whatever they add
        about how the interval was produced.

        Returns
        -------
        dict
            A plain-type representation of my contents. Keys depend on my
            type.
        """
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: summary() must be implemented by subclass"
        )

    @abstractmethod
    def __str__(self) -> str:
        """
        Return a compact string representation of the interval.

        Every subclass describes itself completely, including anything that
        qualifies how the bounds should be read, so that whoever holds me can
        print me without inspecting my fields or knowing my type.
        """
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: __str__() must be implemented by subclass"
        )
