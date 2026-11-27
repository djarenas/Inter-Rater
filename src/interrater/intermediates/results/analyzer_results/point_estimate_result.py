"""Result for any coefficient that produces a single point estimate."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

from ...data_profiles.agreement_profile import AgreementProfile
from ...not_computed import NotComputed
from ..agreement_parts import AgreementParts
from ..panel_result import PanelResult
from ..ci_results.bootstrap_ci_result import BootstrapCIResult
from ..ci_results.ci_result import CIResult

if TYPE_CHECKING:
    from ...data_profiles.balance_profile import BalanceProfile
    from ...data_profiles.rater_profile import RaterProfile

_ERROR_PREFIX = "[PointEstimateResult] Error"
_TOL = 1e-9


@dataclass(frozen=True, kw_only=True)
class PointEstimateResult(PanelResult):
    """
    I hold a point estimate, its confidence interval, and the bootstrap
    distribution when the interval is a bootstrap one.

    I am the result for any coefficient that produces a single number: Fleiss'
    kappa, Cohen's kappa, Conger's kappa, Gwet's AC1, Brennan-Prediger,
    Krippendorff's alpha, etc...

    I am immutable (frozen dataclass) because a result is a snapshot of what
    the data said at one moment and should never change.

    I do not compute my own confidence interval. That is the analyzer's job,
    and it always computes one, so I require it. It is a CIResult (a
    BootstrapCIResult or a ParametricCIResult), which carries what is needed
    to understand the interval. When every resample failed it is still a
    BootstrapCIResult, with no bounds.

    I do not require my point estimate to fall inside my interval. A strongly
    bias-corrected BCa interval can sit entirely to one side of it -- that is
    correct behaviour, not a contradiction, and refusing to exist would be the
    wrong response to it. I say so instead, three ways
    (point_outside_interval, summary()["point_outside_interval"], and POINT
    OUTSIDE INTERVAL in my printed form), because it is also the signature of
    a degenerate bootstrap.

    A bootstrap interval holds its draws, every coefficient value computed
    during resampling, useful for visualization, diagnostics and convergence
    checking. It checks them itself: its counts, and its bounds recomputed
    from them.

    My parent, PanelResult, holds and checks what I share with a
    permutation result: the balance and rater profiles (so every result is
    self-describing without going back to the input data), the interval
    method and the weighting.

    I also carry the observed and chance terms (``parts``), the
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
    confidence_interval : CIResult
        The interval and an account of how it was produced: its bounds and
        alpha, and, for a bootstrap interval, the method asked for against
        the one that ran and the number of resamples that failed.
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
    (balance_profile, rater_profile, ci_method, weights_type and
    weight_matrix are inherited from PanelResult.)
    """

    # What a coefficient is called: its config name, printed name and
    # symbol. Empty here; each coefficient's result declares its own, and
    # every name shown anywhere comes from those (see COEFFICIENT_RESULTS).
    COEFFICIENT: ClassVar[str] = ""
    NAME:        ClassVar[str] = ""
    SYMBOL:      ClassVar[str] = ""

    # Whether the coefficient has a weighted form. Each coefficient's result
    # with one sets True; PanelResult refuses weights otherwise.
    _ACCEPTS_WEIGHTS: ClassVar[bool] = False

    # Inherited from parent class PanelResult
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    ci_method: str
    weights_type: str
    weight_matrix: list[list[float]] | None

    # Populated at this level
    point: float
    confidence_interval: CIResult
    parts: AgreementParts | NotComputed
    kappa_max: float | NotComputed
    agreement_profile: AgreementProfile

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes after construction, before I freeze."""

        # The profiles, interval method and weighting, checked by my parent
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # point must be a float
        if not isinstance(self.point, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: point must be float, got {type(self.point).__name__}"
            )

        # confidence_interval must be a CIResult.
        # Its own bounds, counts and methods were validated when it was built,
        # so there is nothing here to re-check -- only the relationship between
        # it and my point estimate, which it cannot see.
        if not isinstance(self.confidence_interval, CIResult):
            raise TypeError(
                f"{_ERROR_PREFIX}: confidence_interval must be a CIResult, "
                f"got {type(self.confidence_interval).__name__}"
            )

        # My interval was asked for with my interval method
        if self.confidence_interval.method_intended != self.ci_method:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method is '{self.ci_method}' but "
                f"confidence_interval was asked for with "
                f"'{self.confidence_interval.method_intended}'"
            )

        # The point estimate falling outside the interval is permitted --
        # strong bias correction does that legitimately -- but it is also
        # what a degenerate bootstrap looks like, so it is never silent.
        # I say so myself rather than warning: point_outside_interval,
        # summary()["point_outside_interval"], and POINT OUTSIDE INTERVAL
        # in my printed form (and in a permutation table's row). A
        # warning per result would flood a run with many pairs.

        # A BCa interval measured its bias correction against my point
        from ..ci_results.bca_bootstrap_ci_result import BCABootstrapCIResult
        if isinstance(self.confidence_interval, BCABootstrapCIResult) and not (
            self.confidence_interval.point == self.point
            or (math.isnan(self.confidence_interval.point) and math.isnan(self.point))
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: the BCa interval was computed at point "
                f"{self.confidence_interval.point}, but my point is "
                f"{self.point}"
            )

        self._check_parts()
        self._check_kappa_max()

        if not isinstance(self.agreement_profile, AgreementProfile):
            raise TypeError(
                f"{_ERROR_PREFIX}: agreement_profile must be AgreementProfile, "
                f"got {type(self.agreement_profile).__name__}"
            )

        self._check_subjects()

    def _accepts_weights(self) -> bool:
        """Whether my coefficient has a weighted form: my _ACCEPTS_WEIGHTS."""
        return self._ACCEPTS_WEIGHTS

    def _check_subjects(self) -> None:
        """Every subject in my data is one the coefficient used."""
        # Ratings refuses subjects rated fewer than twice, so none can be
        # left out. If one is, something built me from
        # data that skipped that rule.
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
                f"{self.parts.compute_coefficient():.6f}, but the point "
                f"estimate is {self.point:.6f}"
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

        False when the interval has no bounds because nothing ran.
        Legitimate under strong bias
        correction; also the signature of a degenerate bootstrap.
        """
        ci = self.confidence_interval
        if not ci.is_defined:
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
    #
    # _subjects_phrase and _raters_phrase come from PanelResult. All my
    # subjects are used: _check_subjects sees to it.

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
            Keys: 'point' (float); 'ci' (the interval's own summary dict) plus 'ci_lower', 'ci_upper' and
            'point_outside_interval' hoisted to the top level, since those are
            what a reader scans for; for a bootstrap interval,
            'bootstrap_draws_mean' and 'bootstrap_se' (the mean and standard
            error of the usable draws, None when too few are usable -- the
            draw counts are in 'ci'; ``{"not_computed": reason}`` when the
            draws were not kept); 'ci_method', 'weights_type',
            'weight_matrix', 'balance' and 'raters' (from PanelResult); and
            'parts', 'kappa_max' and 'agreement_profile', each a plain value
            or ``{"not_computed": reason}``.
        """
        result = {"point": self.point}
        result["ci"] = self.confidence_interval.summary()
        result["ci_lower"] = self.confidence_interval.lower
        result["ci_upper"] = self.confidence_interval.upper
        result["point_outside_interval"] = self.point_outside_interval
        ci = self.confidence_interval
        if isinstance(ci, BootstrapCIResult):
            if isinstance(ci.draws, NotComputed):
                result["bootstrap_draws_mean"] = ci.draws.summary()
                result["bootstrap_se"] = ci.draws.summary()
            else:
                mean, se = ci.compute_draw_moments()
                result["bootstrap_draws_mean"] = mean
                result["bootstrap_se"] = se
        result.update(self._panel_summary())
        result["parts"] = self.parts.summary()
        result["kappa_max"] = (
            self.kappa_max.summary() if isinstance(self.kappa_max, NotComputed)
            else float(self.kappa_max)
        )
        result["agreement_profile"] = self.agreement_profile.summary()
        return result
