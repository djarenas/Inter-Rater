"""Result for a bootstrap interval where every draw failed."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ...not_computed import NotComputed
from .bootstrap_ci_result import BootstrapCIResult

_ERROR_PREFIX = "[EmptyBootstrapCIResult] Error"


@dataclass(frozen=True)
class EmptyBootstrapCIResult(BootstrapCIResult):
    """
    I am a bootstrap interval that produced nothing: every draw failed, so
    no method could run and I have no bounds.

    That is an outcome rather than an error. The point estimate I would have
    accompanied may be well defined even when no resample of it is, and I
    still say what was asked for and how many draws failed.

    Attributes
    ----------
    (lower, upper, alpha, method_intended and method_used are inherited from
    CIResult; n_resamples, n_failed and draws from BootstrapCIResult.)
    """

    # Inherited from and populated at CIResult parent class
    lower: float
    upper: float
    alpha: float
    method_intended: str
    method_used: str | None

    # Inherited from and populated at BootstrapCIResult parent class
    n_resamples: int
    n_failed: int
    draws: np.ndarray | NotComputed = field(compare=False, repr=False)

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction, before I freeze."""
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    def _check_method(self) -> None:
        """No method ran, because every draw failed."""
        if self.method_used is not None:
            raise ValueError(
                f"method_used must be None, got {self.method_used!r}"
            )
        if self.n_failed != self.n_resamples:
            raise ValueError(
                f"an empty interval has every draw failed, but "
                f"{self.n_usable} of {self.n_resamples} were usable"
            )
