"""Result for any coefficient that produces a single point estimate."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from ...data_profiles.agreement_profile import AgreementProfile
from ...data_profiles.rater_profile import RaterProfile
from ...not_computed import NotComputed
from .._bootstrap_draws import check_draws_match_interval, draw_moments
from ..agreement_parts import AgreementParts
from ..analysis_result import AnalysisResult
from ..ci_results.ci_result import CIResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile

_ERROR_PREFIX = "[PointEstimateResult] Error"
_TOL = 1e-9

@dataclass(frozen=True)
class PointEstimateResult(AnalysisResult):
    """
    I hold a point estimate, optionally a confidence interval, and optionally
    a bootstrap distribution.

    I am the result for any coefficient that produces a single number: Fleiss'
    kappa, Cohen's kappa, Conger's kappa, Gwet's AC1, Brennan-Prediger,
    Krippendorff's alpha, etc...    
    
    I am immutable (frozen dataclass) because a result is a snapshot of what the data 
    said at one moment and should never change.

    I do not compute my own confidence interval. That is the analyzer's job.
    If confidence_interval is None, no interval was computed (perhaps the user
    asked only for a point estimate, or the data was too sparse to resample).
    When one was computed it arrives as a BootstrapCIResult, which carries not
    only the bounds but which method was asked for, which one ran, and how many
    resamples failed.

    I do not require my point estimate to fall inside my interval. A strongly
    bias-corrected BCa interval can sit entirely to one side of it -- that is
    correct behaviour, not a contradiction, and refusing to exist would be the
    wrong response to it. I warn instead, loudly and three ways, because it is
    also the signature of a degenerate bootstrap and a reader should not have
    to notice it for themselves.

    If bootstrap_distribution is provided, it contains all coefficient values
    computed during resampling. This is useful for visualization, diagnostics,
    and convergence checking. It and the interval are two views of one
    bootstrap, so I refuse them if they disagree: draws come exactly with a
    BootstrapCIResult, one per resample it counted, with as many non-finite
    draws as it counted failed.

    I carry a balance_profile so every result is self-describing: the number
    of subjects, the number of raters, whether the data was balanced, and the
    quality of rater-level information are all available without needing to
    go back to the input data. When balance_profile is None, this structural
    metadata was not available at construction time.

    I also carry what the coefficient is made of, so a reader can tell two
    equal coefficients apart: the observed and chance terms (``parts``), the
    largest value the coefficient could reach with these marginals
    (``kappa_max``), and the agreement profile of the data (prevalence,
    rater marginals, specific agreement). These three are required. Where one
    does not apply to a coefficient it holds NotComputed with the reason,
    never None, so "does not apply" cannot be mistaken for "forgotten".

    I check them against my point estimate: the parts must rebuild it, and it
    must not exceed kappa_max.

    Attributes
    ----------
    point : float
        The coefficient value.
    confidence_interval : CIResult | None
        The interval and an account of how it was produced: bounds, alpha,
        the method asked for against the one that ran, and the number of
        resamples that failed. None if no interval was computed.
    bootstrap_distribution : np.ndarray | NotComputed | None
        Array of shape (n_resamples,) containing coefficient values from each
        bootstrap resample, NaN where a resample's coefficient was undefined.
        Present exactly when confidence_interval is a BootstrapCIResult; None
        otherwise (no interval, or a parametric one). NotComputed when the
        bootstrap ran but its draws were not kept (a permutation analysis
        with keep_draws false), saying so.
    balance_profile : BalanceProfile | None
        Structural metadata about the subjects: n_subjects and how many
        ratings each got. None if not available at construction time.
    rater_profile : RaterProfile | NotComputed
        The raters: names, how much each rated, and how much each pair
        overlapped. NotComputed when the result was built from counts,
        which carry no rater identity.
    parts : AgreementParts | NotComputed
        The observed and chance terms (Po and Pe, or Do and De) the point
        estimate is built from. NotComputed when the coefficient is not one
        such pair, with the reason (Light's kappa is a mean of pair kappas).
    kappa_max : float | NotComputed
        The largest value the coefficient could reach given the marginals it
        corrects for (Cohen 1960). NotComputed when there is no such bound
        for this coefficient, with the reason.
    agreement_profile : AgreementProfile
        Category prevalence, per-rater marginals, specific agreement, and the
        prevalence and bias indices where they apply.
    """

    # What a coefficient is called: its config name, printed name and
    # symbol. Empty here; each coefficient's result declares its own, and
    # every name shown anywhere comes from those (see COEFFICIENT_RESULTS).
    COEFFICIENT: ClassVar[str] = ""
    NAME:        ClassVar[str] = ""
    SYMBOL:      ClassVar[str] = ""

    point: float
    confidence_interval: CIResult | None
    bootstrap_distribution: np.ndarray | NotComputed | None
    balance_profile: BalanceProfile | None
    rater_profile: RaterProfile | NotComputed
    parts: AgreementParts | NotComputed
    kappa_max: float | NotComputed
    agreement_profile: AgreementProfile

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction, before I freeze."""

        # point must be a float
        if not isinstance(self.point, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: point must be float, got {type(self.point).__name__}"
            )

        # If confidence_interval is provided, it must be a BootstrapCIResult.
        # Its own bounds, counts and methods were validated when it was built,
        # so there is nothing here to re-check -- only the relationship between
        # it and my point estimate, which it cannot see.
        if self.confidence_interval is not None:
            if not isinstance(self.confidence_interval, CIResult):
                raise TypeError(
                    f"{_ERROR_PREFIX}: confidence_interval must be a CIResult "
                    f"or None, got {type(self.confidence_interval).__name__}"
                )

            # The point estimate falling outside the interval is permitted --
            # strong bias correction does that legitimately -- but it is also
            # what a degenerate bootstrap looks like, so it is never silent.
            # I say so myself rather than warning: point_outside_interval,
            # summary()["point_outside_interval"], and POINT OUTSIDE INTERVAL
            # in my printed form (and in a permutation table's row). A
            # warning per result would flood a run with many pairs.

        # If bootstrap_distribution is provided, validate it
        if self.bootstrap_distribution is not None and not isinstance(
            self.bootstrap_distribution, NotComputed
        ):
            if not isinstance(self.bootstrap_distribution, np.ndarray):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution must be np.ndarray, NotComputed or None, "
                    f"got {type(self.bootstrap_distribution).__name__}"
                )

            if self.bootstrap_distribution.ndim != 1:
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution must be 1-dimensional, "
                    f"got shape {self.bootstrap_distribution.shape}"
                )

            if len(self.bootstrap_distribution) == 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: bootstrap_distribution must have at least 1 value"
                )

        # My draws and my interval must describe the same bootstrap
        try:
            check_draws_match_interval(
                self.bootstrap_distribution,
                self.confidence_interval,
                draws_name="bootstrap_distribution",
                interval_name="confidence_interval",
            )
        except ValueError as exc:
            raise ValueError(f"{_ERROR_PREFIX}: {exc}") from exc

        # If balance_profile is provided, validate its type
        if self.balance_profile is not None:
            from ...data_profiles.balance_profile import BalanceProfile

            if not isinstance(self.balance_profile, BalanceProfile):
                raise TypeError(
                    f"{_ERROR_PREFIX}: balance_profile must be BalanceProfile or None, "
                    f"got {type(self.balance_profile).__name__}"
                )

        if not isinstance(self.rater_profile, (RaterProfile, NotComputed)):
            raise TypeError(
                f"{_ERROR_PREFIX}: rater_profile must be RaterProfile or "
                f"NotComputed, got {type(self.rater_profile).__name__}"
            )
        # Both profiles describe the same data
        if isinstance(self.rater_profile, RaterProfile) and self.balance_profile is not None:
            if self.rater_profile.n_subjects != self.balance_profile.n_subjects:
                raise ValueError(
                    f"{_ERROR_PREFIX}: rater_profile has "
                    f"{self.rater_profile.n_subjects} subjects but "
                    f"balance_profile has {self.balance_profile.n_subjects}"
                )

        self._check_parts()
        self._check_kappa_max()

        if not isinstance(self.agreement_profile, AgreementProfile):
            raise TypeError(
                f"{_ERROR_PREFIX}: agreement_profile must be AgreementProfile, "
                f"got {type(self.agreement_profile).__name__}"
            )

        self._check_subjects()

    def _check_subjects(self) -> None:
        """Every subject in my data is one the coefficient used."""
        # Ratings and build_from_nij_matrix refuse subjects rated fewer than
        # twice, so none can be left out. If one is, something built me from
        # data that skipped that rule.
        if self.balance_profile is None:
            return
        total = self.balance_profile.n_subjects
        used  = self.agreement_profile.n_subjects
        if used != total:
            raise ValueError(
                f"{_ERROR_PREFIX}: the data has {total:,} subjects, but the "
                f"agreement profile uses {used:,}. Every subject must be rated "
                f"at least twice before a result is built."
            )

    def _check_parts(self) -> None:
        """My parts, when I have them, must rebuild my point estimate."""
        if isinstance(self.parts, NotComputed):
            return
        if not isinstance(self.parts, AgreementParts):
            raise TypeError(
                f"{_ERROR_PREFIX}: parts must be AgreementParts or NotComputed, "
                f"got {type(self.parts).__name__}"
            )
        if not self.parts.reproduces(float(self.point)):
            raise ValueError(
                f"{_ERROR_PREFIX}: parts ({self.parts}) make "
                f"{self.parts.coefficient():.6f}, but the point estimate is "
                f"{self.point:.6f}"
            )

    def _check_kappa_max(self) -> None:
        """kappa_max is in [0, 1], and my point estimate is at most kappa_max."""
        if isinstance(self.kappa_max, NotComputed):
            return
        if isinstance(self.kappa_max, bool) or not isinstance(
            self.kappa_max, (int, float)
        ):
            raise TypeError(
                f"{_ERROR_PREFIX}: kappa_max must be a number or NotComputed, "
                f"got {type(self.kappa_max).__name__}"
            )
        if math.isnan(self.kappa_max):
            return
        # At least 0: ratings independent with these marginals are always
        # possible, and they make the coefficient 0.
        if not (-_TOL <= self.kappa_max <= 1 + _TOL):
            raise ValueError(
                f"{_ERROR_PREFIX}: kappa_max must be in [0, 1], "
                f"got {self.kappa_max}"
            )
        if not math.isnan(self.point) and self.point > self.kappa_max + _TOL:
            raise ValueError(
                f"{_ERROR_PREFIX}: the point estimate {self.point:.6f} exceeds "
                f"kappa_max {self.kappa_max:.6f}, the most these marginals "
                f"allow"
            )

    # Public API ------------------------------------------------------------- #

    @property
    def point_outside_interval(self) -> bool:
        """
        Whether my point estimate falls outside my own interval.

        False when there is no interval, and when the interval exists but has
        no bounds because nothing ran. Legitimate under strong bias
        correction; also the signature of a degenerate bootstrap.
        """
        ci = self.confidence_interval
        if ci is None or not ci.is_defined:
            return False
        return not (ci.lower <= self.point <= ci.upper)

    def __str__(self) -> str:
        """
        Return my printed form: up to three lines.

            Fleiss' κ = 0.493 [0.382, 0.609] (95% CI, bca, n=1000)
              Po = 0.853, Pe = 0.709; κ_max not computed (pooled shares)
              199 subjects; 4 raters; unweighted; exact formula

        The first line is the coefficient's name, my point estimate and my
        interval. The interval describes itself -- I do not inspect it, so
        whatever qualifies it travels with it unchanged, whether that is a
        fallback, a count of failed resamples, or a caveat about the
        variance it was built from. I add only what it cannot know, which is
        where my point estimate sits relative to it ("POINT OUTSIDE
        INTERVAL").

        The second line is what the coefficient is made of, and the third
        the data it was computed on. Each coefficient supplies its own; a
        line it has nothing for is left out. The agreement profile's tables
        print separately: ``print(result.agreement_profile)``.
        """
        name = self._name()
        first = f"{name} {self._headline()}" if name else self._headline()
        lines = [first]
        for line in (self._made_of_line(), self._data_line()):
            if line:
                lines.append(f"  {line}")
        return "\n".join(lines)

    # Printed-form pieces, supplied by each coefficient ---------------------- #

    def _headline(self) -> str:
        """= point [interval] (interval's own description), and where the point sits."""
        parts = [f"= {self.point:.3f}"]

        if self.confidence_interval is not None:
            parts.append(str(self.confidence_interval))
            if self.point_outside_interval:
                parts.append("POINT OUTSIDE INTERVAL")

        return " ".join(parts)

    def _name(self) -> str:
        """The coefficient's printed name: NAME, declared by each coefficient."""
        return self.NAME

    def _kappa_max_note(self) -> str:
        """A few words on kappa_max, to fit on one line."""
        if isinstance(self.kappa_max, NotComputed):
            return f"κ_max not computed ({self.kappa_max.reason})"
        if math.isnan(self.kappa_max):
            return "κ_max undefined"
        return f"κ_max = {self.kappa_max:.3f}"

    def _made_of_line(self) -> str | None:
        """Po and Pe (or Do and De), and kappa_max; None without parts."""
        if isinstance(self.parts, NotComputed):
            return None
        return f"{self.parts}; {self._kappa_max_note()}"

    def _data_line(self) -> str | None:
        """The data the coefficient was computed on; none at this level."""
        return None

    # Shared phrases for the data line -------------------------------------- #

    def _subjects_phrase(self) -> str:
        """'1,617 subjects'. All of them are used: _check_subjects sees to it."""
        return f"{self.balance_profile.n_subjects:,} subjects"

    def _raters_phrase(self) -> str:
        """'4 raters', or 'at least 3 raters' when only counts were given."""
        if isinstance(self.rater_profile, RaterProfile):
            return f"{self.rater_profile.n_raters} raters"
        return f"at least {self.balance_profile.n_raters_lower_bound} raters"

    @staticmethod
    def _weights_phrase(weights_type: str) -> str:
        """'unweighted', or 'linear weights', 'quadratic weights', 'custom weights'."""
        return "unweighted" if weights_type == "identity" else f"{weights_type} weights"

    def summary(self) -> dict:
        """
        Return my point estimate, interval, bootstrap diagnostics, and
        balance metadata in plain types.

        Returns
        -------
        dict
            Keys: 'point' (float); if confidence_interval is not None, 'ci'
            (the interval's own summary dict) plus 'ci_lower', 'ci_upper' and
            'point_outside_interval' hoisted to the top level, since those are
            what a reader scans for; if bootstrap_distribution is not None,
            'bootstrap_draws_mean' and 'bootstrap_se' (the mean and standard
            error of the usable draws, None when too few are usable -- the
            draw counts are in 'ci'; ``{"not_computed": reason}`` when the
            draws were not kept); 'balance' if balance_profile is
            not None; and 'raters', 'parts', 'kappa_max' and
            'agreement_profile', each a plain value or
            ``{"not_computed": reason}``.
        """
        result = {"point": self.point}
        if self.confidence_interval is not None:
            result["ci"] = self.confidence_interval.summary()
            result["ci_lower"] = self.confidence_interval.lower
            result["ci_upper"] = self.confidence_interval.upper
            result["point_outside_interval"] = self.point_outside_interval
        if isinstance(self.bootstrap_distribution, NotComputed):
            result["bootstrap_draws_mean"] = self.bootstrap_distribution.summary()
            result["bootstrap_se"] = self.bootstrap_distribution.summary()
        elif self.bootstrap_distribution is not None:
            mean, se = draw_moments(self.bootstrap_distribution)
            result["bootstrap_draws_mean"] = mean
            result["bootstrap_se"] = se
        if self.balance_profile is not None:
            result["balance"] = self.balance_profile.summary()
        result["raters"] = self.rater_profile.summary()
        result["parts"] = self.parts.summary()
        result["kappa_max"] = (
            self.kappa_max.summary() if isinstance(self.kappa_max, NotComputed)
            else float(self.kappa_max)
        )
        result["agreement_profile"] = self.agreement_profile.summary()
        return result
