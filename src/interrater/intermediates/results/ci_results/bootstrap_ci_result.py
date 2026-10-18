"""Result for a confidence interval computed from a bootstrap distribution."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .ci_result import CIResult

_ERROR_PREFIX = "[BootstrapCIResult] Error"

_VALID_METHODS = ("percentile", "bca")

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
    I am an interval from a bootstrap, and an honest account of what it was
    computed from.

    I am what a CI method returns. I am not an analysis result and do not
    inherit from one: I hold no point estimate, no coefficient name, and no
    opinion about the data. Something that holds me supplies all of that.

    I carry the number of resamples that failed rather than a verdict about
    them. A resample fails when its coefficient is undefined -- a resample
    containing only one category drives chance agreement to 1, and
    (p_o - p_e) / (1 - p_e) has nothing to say. Those draws arrive as NaN and
    are excluded from my bounds. How much that matters depends on whether
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
    verdict: an interval from too few draws is noisy, not wrong, and nothing
    here refuses to produce one. Unlike my other disclosures it does not
    warn, because one analysis builds many of me -- a pairwise run makes one
    per rater pair -- and fifteen identical warnings for one command teaches
    people to silence warnings.

    Attributes
    ----------
    n_resamples : int
        How many draws were handed to the CI method, failures included.
    n_failed : int
        How many of those were unusable -- NaN, which is how an undefined
        coefficient arrives, or infinite.
    method_intended : str
        The method that was asked for: 'percentile' or 'bca'. Set by which
        CI function ran, never read from a config.
    method_used : str | None
        The method that actually produced my bounds, or None when no method
        could be applied at all. None is an absence, not a substitution:
        every resample failed, so there was nothing to take quantiles of.
        Otherwise equal to method_intended, unless BCa fell back.
    acceleration : float | None
        BCa's acceleration constant. None whenever method_used is not 'bca',
        including after a fallback -- an unused diagnostic stored anyway
        invites someone to read it as applied.
    bias_correction : float | None
        BCa's z0. None whenever method_used is not 'bca', for the same reason.

    See Also
    --------
    CIResult : my base, which holds the bounds and the level.
    """

    n_resamples:     int
    n_failed:        int
    method_intended: str
    method_used:     str | None
    acceleration:    float | None = None
    bias_correction: float | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """
        Validate my attributes after construction, before I freeze.

        Structural only. Each of these describes an interval no correct
        computation could have produced, so reaching one means a bug upstream
        rather than awkward data. Too few usable draws for the method that ran
        is a real concern and not one of mine: it takes judgement, it has no
        defensible cutoff, and whoever holds me is the one who can weigh it.
        """

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

        # Something was always asked for
        if self.method_intended not in _VALID_METHODS:
            raise ValueError(
                f"{_ERROR_PREFIX}: method_intended must be one of "
                f"{list(_VALID_METHODS)}, got {self.method_intended!r}"
            )

        # Something ran, or nothing did
        if self.method_used is not None and self.method_used not in _VALID_METHODS:
            raise ValueError(
                f"{_ERROR_PREFIX}: method_used must be None or one of "
                f"{list(_VALID_METHODS)}, got {self.method_used!r}"
            )

        # Percentile has no failure that substitutes another method for it, so
        # it either ran or nothing did. Only BCa can degrade.
        if (
            self.method_intended == "percentile"
            and self.method_used not in ("percentile", None)
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: percentile cannot fall back, but method_used "
                f"is {self.method_used!r}"
            )

        # Three facts that move together: a method runs only when some
        # resample was usable, and bounds exist only when a method ran.
        # Quantiles of nothing do not exist, and quantiles of a non-empty set
        # of finite numbers are finite.
        ran = self.method_used is not None
        bounds_defined = bool(np.isfinite(self.lower) and np.isfinite(self.upper))

        if ran != bounds_defined:
            raise ValueError(
                f"{_ERROR_PREFIX}: bounds are defined if and only if a method "
                f"ran, but method_used is {self.method_used!r} and bounds are "
                f"({self.lower}, {self.upper})"
            )
        if ran != (self.n_usable > 0):
            raise ValueError(
                f"{_ERROR_PREFIX}: a method runs if and only if some resample "
                f"was usable, but method_used is {self.method_used!r} and "
                f"{self.n_usable} of {self.n_resamples} resamples were usable"
            )

        # Diagnostics belong to BCa and only to BCa
        diagnostics = (
            ("acceleration", self.acceleration),
            ("bias_correction", self.bias_correction),
        )
        if self.method_used == "bca":
            for name, value in diagnostics:
                if value is None:
                    raise ValueError(
                        f"{_ERROR_PREFIX}: {name} is required when method_used "
                        f"is 'bca'"
                    )
                if not isinstance(value, (int, float, np.floating)) or isinstance(
                    value, bool
                ):
                    raise TypeError(
                        f"{_ERROR_PREFIX}: {name} must be float, "
                        f"got {type(value).__name__}"
                    )
                if not np.isfinite(value):
                    raise ValueError(
                        f"{_ERROR_PREFIX}: {name} must be finite, got {value}"
                    )
            object.__setattr__(self, "acceleration", float(self.acceleration))
            object.__setattr__(self, "bias_correction", float(self.bias_correction))
        else:
            for name, value in diagnostics:
                if value is not None:
                    raise ValueError(
                        f"{_ERROR_PREFIX}: {name} must be None when method_used "
                        f"is {self.method_used!r}, got {value}"
                    )

        # Parent validates the bounds and the level, and coerces them
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    # Public API ------------------------------------------------------------- #

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

    @property
    def degraded(self) -> bool:
        """
        Whether a different method ran than the one asked for.

        False when nothing ran at all. That is an absence, not a
        substitution, and ``method_used is None`` already says it.
        """
        return (
            self.method_used is not None
            and self.method_intended != self.method_used
        )

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
            'recommended_resamples' (int or None), and 'acceleration',
            'bias_correction' (float) only when method_used is 'bca'.
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
        if self.method_used == "bca":
            result["acceleration"]    = self.acceleration
            result["bias_correction"] = self.bias_correction
        return result
