"""Result for a confidence interval computed from a bootstrap distribution."""

from __future__ import annotations

from abc import abstractmethod
from dataclasses import dataclass, field
from typing import ClassVar

import numpy as np

from ...not_computed import NotComputed
from .ci_result import CIResult

_ERROR_PREFIX = "[BootstrapCIResult] Error"

# Efron & Tibshirani (1986) section 9, p. 72: "B = 1000 is a rough minimum for
# the number of Monte Carlo bootstraps necessary to compute the BC or BC_a
# intervals. Somewhat smaller values, say B = 250, can give a useful percentile
# interval, the difference being that then the constant z_0 need not be
# computed."
#
# These live here rather than beside the CI functions because the interval is
# what reports against them. Nothing enforces them anywhere: too few draws is a
# judgement, it has no defensible cutoff, and an interval computed from 900
# draws is not invalid. It is just worth saying out loud.
MIN_RESAMPLES_PERCENTILE = 250
MIN_RESAMPLES_BCA = 1000

_MIN_RESAMPLES = {
    "percentile": MIN_RESAMPLES_PERCENTILE,
    "bca":        MIN_RESAMPLES_BCA,
}


@dataclass(frozen=True)
class BootstrapCIResult(CIResult):
    """
    I am an interval from a bootstrap, and a clear account of what it was
    computed from.

    I am what a CI method returns and not a panel result: I hold no
    coefficient name and no opinion about the data. Only my BCa child holds a
    point, the one its bias correction was measured against.

    I carry the number of resamples that failed. A resample commonly fails
    when it contains only one category, which drives chance agreement to 1.
    Those draws arrive as NaN and are excluded from my bounds.
    How much that matters depends on whether
    failing correlates with the coefficient's value, which no count can tell
    you, so I report the count and leave the judgement to whoever holds me.

    I record both the method that was asked for and the method that ran. They
    differ when BCa could not apply its own correction and fell back. That is
    a real fallback, not a relabelling: when I say ``percentile`` the bias
    correction was discarded too, and my bounds are plain quantiles.

    I also say when my bounds rest on fewer usable draws than my method is
    normally given -- 1000 for BCa, 250 for percentile, from Efron &
    Tibshirani (1986) section 9. I compare against the usable count rather
    than the requested one, because 1000 draws with 900 failures leaves 100,
    and 100 is what the quantiles were taken over. This is disclosure, not a
    verdict.

    I hold the draws themselves, and check them against my counts: one draw
    per resample, and as many non-finite draws as I counted failed. Each of
    my children checks its own bounds against the draws: Percentile,
    BCA and Empty.

    Attributes
    ----------
    n_resamples : int
        How many draws were handed to the CI method, failures included.
    n_failed : int
        How many of those were unusable -- NaN, which is how an undefined
        coefficient arrives, or infinite.
    draws : np.ndarray | NotComputed
        The draws, shape (n_resamples,), NaN where a resample's coefficient
        was undefined. NotComputed when they were not kept (a permutation
        analysis with keep_draws false), saying why.
    (method_intended and method_used are inherited from CIResult; mine are
    'percentile' or 'bca', and BCa can fall back to percentile.)

    See Also
    --------
    CIResult : my base, which holds the bounds and the level.
    """

    # The methods I can report
    _METHODS: ClassVar[frozenset[str]] = frozenset({"percentile", "bca"})

    # Inherited from and populated at CIResult parent class
    lower: float
    upper: float
    alpha: float
    method_intended: str
    method_used: str | None

    # Populated at this level
    n_resamples:     int
    n_failed:        int
    # Left out of equality and repr: an array has no single truth value,
    # and a thousand numbers do not belong in a one-line repr
    draws:           np.ndarray | NotComputed = field(compare=False, repr=False)

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """
        Validate my attributes after construction.

        Structural only. Each of these describes an interval no correct
        computation could have produced, so reaching one means a bug upstream
        rather than awkward data. Too few usable draws for the method that ran
        is a real concern, but is not one of mine. I hold the counts, so
        whoever holds or prints me can weigh them.
        """
        # Parent validates and populates lower, upper, and alpha attributes
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Counts must be ints, and must agree with each other
        for name, value in (
            ("n_resamples", self.n_resamples), ("n_failed", self.n_failed)
        ):
            if not isinstance(value, (int, np.integer)) or isinstance(value, bool):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must be int, "
                    f"got {type(value).__name__}"
                )
        if self.n_resamples < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_resamples must be >= 1, got {self.n_resamples}"
            )
        if not (0 <= self.n_failed <= self.n_resamples):
            raise ValueError(
                f"{_ERROR_PREFIX}: n_failed must be in [0, {self.n_resamples}], "
                f"got {self.n_failed}"
            )

        # The draws, when kept, are the ones my counts describe
        if not isinstance(self.draws, NotComputed):
            if not isinstance(self.draws, np.ndarray):
                raise TypeError(
                    f"{_ERROR_PREFIX}: draws must be np.ndarray or "
                    f"NotComputed, got {type(self.draws).__name__}"
                )
            if self.draws.ndim != 1:
                raise ValueError(
                    f"{_ERROR_PREFIX}: draws must be 1-dimensional, got "
                    f"shape {self.draws.shape}"
                )
            if self.draws.size != self.n_resamples:
                raise ValueError(
                    f"{_ERROR_PREFIX}: there are {self.draws.size} draws but "
                    f"n_resamples is {self.n_resamples}"
                )
            n_unusable = int(np.count_nonzero(~np.isfinite(self.draws)))
            if n_unusable != self.n_failed:
                raise ValueError(
                    f"{_ERROR_PREFIX}: {n_unusable} draws are non-finite but "
                    f"n_failed is {self.n_failed}"
                )

        # My child checks its method and its bounds
        try:
            self._check_method()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    @abstractmethod
    def _check_method(self) -> None:
        """Check my method, and my bounds against my draws. Each child says."""
        # Just an extra precaution to ensure proper inheritance of the
        # abstract method.
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: _check_method() must be implemented by subclass"
        )

    # Public API ------------------------------------------------------------- #

    @staticmethod
    def usable_draws(draws: np.ndarray) -> np.ndarray:
        """The usable draws: the finite ones, as floats."""
        values = np.asarray(draws, dtype=np.float64)
        return values[np.isfinite(values)]

    def compute_draw_moments(self) -> tuple[float | None, float | None]:
        """
        Mean and bootstrap standard error of my usable draws.

        The standard error divides by B - 1, B being my usable draws (Efron
        & Tibshirani 1986, eq. 2.4). The mean is None with no usable draw; the standard error is
        None with fewer than two. Raises when my draws were not kept.
        """
        if isinstance(self.draws, NotComputed):
            raise ValueError(
                f"{_ERROR_PREFIX}: my draws were not kept "
                f"({self.draws.reason}), so they have no moments"
            )
        usable = self.usable_draws(self.draws)
        mean = float(np.mean(usable)) if usable.size >= 1 else None
        se = float(np.std(usable, ddof=1)) if usable.size >= 2 else None
        return mean, se

    @property
    def n_usable(self) -> int:
        """How many draws my bounds were actually computed from."""
        return self.n_resamples - self.n_failed

    @property
    def recommended_resamples(self) -> int | None:
        """
        How many usable draws the method that ran is normally given.

        None when no method ran -- there is no recommendation for a thing
        that did not happen, and the failure count already says everything.
        """
        if self.method_used is None:
            return None
        return _MIN_RESAMPLES[self.method_used]

    @property
    def below_recommended_resamples(self) -> bool:
        """
        Whether my bounds rest on fewer usable draws than my method expects.

        False when no method ran. Compared against ``n_usable``, since that
        is the number the quantiles were actually taken over.
        """
        minimum = self.recommended_resamples
        return minimum is not None and self.n_usable < minimum

    def __str__(self) -> str:
        """
        Return a compact string representation of the interval.

        Format: [0.481, 0.565] (95% CI, bca, n=1000)

        A fallback is named rather than hidden, and so is a failure count,
        a method that never ran at all, and bounds resting on fewer usable
        draws than the method is normally given:

            [0.481, 0.565] (95% CI, bca, n=1000)
            [0.481, 0.565] (95% CI, percentile [bca requested], n=1000)
            [0.481, 0.565] (95% CI, bca, n=1200, 12 failed)
            [0.339, 0.551] (95% CI, bca, n=50, 50 usable, 1000 recommended)
            [undefined] (95% CI, no method ran [bca requested], n=200, all 200 failed)

        The third line is n=1200 rather than n=1000 deliberately. BCa asked
        for exactly 1000 with any failure at all has fewer than 1000 usable,
        so it would carry the recommendation clause too.
        """
        bounds = (
            f"[{self.lower:.3f}, {self.upper:.3f}]"
            if self.is_defined else "[undefined]"
        )

        meta = [f"{self.confidence_level:.0%} CI"]
        if self.method_used is None:
            meta.append(f"no method ran [{self.method_intended} requested]")
        elif self.degraded:
            meta.append(f"{self.method_used} [{self.method_intended} requested]")
        else:
            meta.append(self.method_used)
        meta.append(f"n={self.n_resamples}")
        if self.n_failed:
            meta.append(
                f"all {self.n_failed} failed" if self.n_usable == 0
                else f"{self.n_failed} failed"
            )
        if self.below_recommended_resamples:
            # Two facts side by side, with no verb between them for a reader
            # to misattach. An earlier wording, "B below 1000 recommended",
            # parsed as a recommendation to stay below 1000.
            meta.append(
                f"{self.n_usable} usable, {self.recommended_resamples} recommended"
            )

        return f"{bounds} ({', '.join(meta)})"

    def summary(self) -> dict:
        """
        Return my bounds, level, counts, methods, and diagnostics in plain types.

        Returns
        -------
        dict
            Keys: 'lower', 'upper' (float), 'alpha', 'confidence_level'
            (float), 'method_intended' (str), 'method_used' (str or None),
            'degraded', 'is_defined', 'below_recommended_resamples' (bool),
            'n_resamples', 'n_failed', 'n_usable' (int),
            'recommended_resamples' (int or None). BCa adds 'acceleration'
            and 'bias_correction' (float).
        """
        result = {
            "lower":            self.lower,
            "upper":            self.upper,
            "alpha":            self.alpha,
            "confidence_level": self.confidence_level,
            "method_intended":  self.method_intended,
            "method_used":      self.method_used,
            "degraded":         self.degraded,
            "is_defined":       self.is_defined,
            "n_resamples":      int(self.n_resamples),
            "n_failed":         int(self.n_failed),
            "n_usable":         self.n_usable,
            "recommended_resamples":      self.recommended_resamples,
            "below_recommended_resamples": self.below_recommended_resamples,
        }
        return result
