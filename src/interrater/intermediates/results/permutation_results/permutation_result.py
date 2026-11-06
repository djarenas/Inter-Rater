"""Base result for permutation analysis (pairwise or subset)."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..analysis_result import AnalysisResult
from ..analyzer_results import result_class_for

if TYPE_CHECKING:
    import pandas as pd

    from ...data_profiles.balance_profile import BalanceProfile
    from ...data_profiles.rater_profile import RaterProfile
    from ..ci_results.bootstrap_ci_result import BootstrapCIResult

_ERROR_PREFIX = "[PermutationResult] Error"

# Shown wherever an "independent" mean interval is shown, so that the text
# and the plots say the same thing about it.
INDEPENDENT_MEAN_CAVEAT = (
    "independent bootstrap: ignores the correlation between combinations "
    "that share raters, so it tends to be too narrow"
)

# A coefficient's printed name, symbol and result class all come from its
# result class (analyzer_results.result_class_for), so nothing here spells them.

# The weightings a result can record; "identity" is unweighted
_WEIGHTS_TYPES = ("identity", "linear", "quadratic", "custom")


def has_weighted_form(coefficient: str) -> bool:
    """Whether the coefficient's own result records a weighting."""
    return "weights_type" in result_class_for(coefficient).__dataclass_fields__


def raters_phrase(rater_profile: RaterProfile) -> str:
    """'4 raters'. A permutation result always knows its raters exactly."""
    return f"{rater_profile.n_raters} raters"


def interval_notes(ci: BootstrapCIResult) -> list[str]:
    """
    What a per-pair or per-subset interval needs said beside its bounds.

    The method that ran, when it is not the one asked for, and failed draws.
    The level, method and draw count are the same for every row, so they
    are said once, in the column heading and the header lines.
    """
    notes = []
    if ci.method_used is None:
        notes.append(f"no method ran [{ci.method_intended} requested]")
    elif ci.degraded:
        notes.append(f"{ci.method_used} [{ci.method_intended} requested]")
    if ci.n_failed:
        notes.append(
            f"all {ci.n_failed} failed" if ci.n_usable == 0
            else f"{ci.n_failed} failed"
        )
    return notes


@dataclass(frozen=True)
class PermutationResult(AnalysisResult):
    """
    I am the base result for permutation analysis across rater subsets.

    I am the permutation-branch counterpart of PointEstimateResult. Where
    PointEstimateResult holds one number for one coefficient, my children
    hold C(n_raters, k) numbers -- one per subset or pair -- plus summary
    statistics and bootstrap intervals. The two branches share
    AnalysisResult as their root but are otherwise independent.

    I hold only what both children genuinely share:
    - Which coefficient was computed
    - How confidence intervals were computed
    - The structural metadata of the dataset

    All statistical substance -- observed values, means, bootstrap
    distributions, per-pair or per-subset CIs -- lives on the children,
    where the vocabulary and semantics belong. Pairwise talks about
    "pairs," subset talks about "subsets." The parent does not impose
    either vocabulary.

    Subclasses
    ----------
    PairwisePermutationResult
        Coefficient computed on every C(n_raters, 2) pair of raters.
    SubsetPermutationResult
        Coefficient computed on every C(n_raters, k) subset of raters.

    Attributes
    ----------
    coefficient_name : str
        Which coefficient was computed: "cohen_kappa", "scotts_pi",
        "conger_kappa", etc.
    ci_method : str
        How confidence intervals were computed: "percentile" or "bca".
        Permutation results have no closed-form null variance (subsets
        sharing a rater are correlated), so "parametric_null" is not
        available.
    balance_profile : BalanceProfile
        Structural metadata about the full dataset (before subsetting):
        n_subjects and how many ratings each got.
    rater_profile : RaterProfile
        The panel's raters: names in position order as the semantics name
        them (injected raters included), how much each rated, and how much
        each pair overlapped. Always present: permutation needs rater
        identity. Pairs and subsets are keyed by position, e.g. (0, 2);
        everything shown to people -- the printed table, to_frame,
        summary, plots -- names them through rater_names.
    weights_type : str
        The weighting every pair or subset was computed with: "identity"
        (unweighted), "linear", "quadratic" or "custom". The mean is over
        those coefficients, so it is weighted the same way. Always
        "identity" for a coefficient with no weighted form.
    weight_matrix : list[list[float]] | None
        The weight matrix that was used, as nested lists; None when
        unweighted. Every pair or subset that records weights must have
        this same matrix.
    """

    # Attribute declaration

    coefficient_name: str
    ci_method: str
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    weights_type: str
    weight_matrix: list[list[float]] | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate permutation result attributes."""

        # Validate coefficient_name
        if not isinstance(self.coefficient_name, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: coefficient_name must be str, "
                f"got {type(self.coefficient_name).__name__}"
            )
        if not self.coefficient_name:
            raise ValueError(
                f"{_ERROR_PREFIX}: coefficient_name must be non-empty"
            )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or 'bca', "
                f"got '{self.ci_method}'. Permutation results have no "
                f"closed-form null variance (subsets are correlated), so "
                f"'parametric_null' is not available."
            )

        # Validate balance_profile
        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None"
            )

        from ...data_profiles.balance_profile import BalanceProfile

        if not isinstance(self.balance_profile, BalanceProfile):
            raise TypeError(
                f"{_ERROR_PREFIX}: balance_profile must be BalanceProfile, "
                f"got {type(self.balance_profile).__name__}"
            )

        try:
            self._check_rater_profile()
            self._check_weights()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

    def _check_weights(self) -> None:
        """The weighting is stated, and a coefficient with no weighted form has none."""
        if self.weights_type not in _WEIGHTS_TYPES:
            raise ValueError(
                f"weights_type must be one of {list(_WEIGHTS_TYPES)}, "
                f"got {self.weights_type!r}"
            )
        if (self.weights_type == "identity") != (self.weight_matrix is None):
            raise ValueError(
                f"weight_matrix must be None exactly when weights_type is "
                f"'identity' (unweighted); got weights_type "
                f"{self.weights_type!r} with "
                f"{'no matrix' if self.weight_matrix is None else 'a matrix'}"
            )
        if self.weight_matrix is not None:
            matrix = np.asarray(self.weight_matrix, dtype=float)
            if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
                raise ValueError(
                    f"weight_matrix must be square, got shape {matrix.shape}"
                )
        if self.weights_type != "identity" and not has_weighted_form(self.coefficient_name):
            raise ValueError(
                f"{self.coefficient_name} has no weighted form, so "
                f"weights_type must be 'identity', got {self.weights_type!r}"
            )

    def _check_rater_profile(self) -> None:
        """A RaterProfile, which checks its own names and counts."""
        from ...data_profiles.rater_profile import RaterProfile

        if not isinstance(self.rater_profile, RaterProfile):
            raise TypeError(
                f"rater_profile must be RaterProfile, got "
                f"{type(self.rater_profile).__name__}"
            )
        # Both profiles describe the same panel
        n_subjects = getattr(self.balance_profile, "n_subjects", None)
        if n_subjects is not None and self.rater_profile.n_subjects != n_subjects:
            raise ValueError(
                f"rater_profile has {self.rater_profile.n_subjects} subjects "
                f"but balance_profile has {n_subjects}"
            )

    @property
    def rater_names(self) -> tuple[str, ...]:
        """The raters in position order: read from rater_profile, not a copy."""
        return self.rater_profile.rater_names

    def _named(self, raters: tuple[int, ...]) -> str:
        """A unit's raters by name, e.g. "ann, cat" for (0, 2)."""
        return ", ".join(self.rater_names[i] for i in raters)

    def _unit_summaries(self, results: dict) -> list[dict]:
        """
        Each pair or subset for summary(): its raters by name, its
        positions, and its own summary (or {"not_computed": reason}).
        """
        return [
            {
                "raters": [self.rater_names[i] for i in raters],
                "positions": list(raters),
                "result": unit.summary(),
            }
            for raters, unit in results.items()
        ]

    # Public API ------------------------------------------------------------- #

    def to_frame(self) -> pd.DataFrame:
        """
        Every pair or subset as a row of a pandas DataFrame, to save or analyse.

        ``result.to_frame().to_csv("pairs.csv")`` writes them all. The
        printed list is made from the same rows, so the two cannot disagree.

        Returns
        -------
        pd.DataFrame
            Columns: raters (by name, "ann, cat"), positions (the key in
            pair_results or subset_results, "0, 2"), estimate, ci_lower,
            ci_upper, point_outside_interval, method_used, n_failed,
            subjects, not_computed.
            Undefined or empty values are NaN, in text columns too, so the
            frame reads back from CSV unchanged. A pair or subset that has
            no result leaves the other columns NaN and says why in
            not_computed; one that has a result has NaN there.
        """
        import pandas as pd

        nan = float("nan")
        return pd.DataFrame(
            [
                {
                    "raters": self._named(raters),
                    "positions": ", ".join(map(str, raters)),
                    "estimate": float(value),
                    "ci_lower": float(ci.lower) if ci is not None else nan,
                    "ci_upper": float(ci.upper) if ci is not None else nan,
                    "point_outside_interval": outside if outside is not None else nan,
                    # None from an interval where no method ran
                    "method_used": (
                        ci.method_used
                        if ci is not None and ci.method_used is not None
                        else nan
                    ),
                    "n_failed": int(ci.n_failed) if ci is not None else nan,
                    "subjects": int(subjects) if subjects is not None else nan,
                    "not_computed": reason if reason is not None else nan,
                }
                for raters, value, ci, outside, subjects, reason in self._units()
            ],
            columns=[
                "raters", "positions", "estimate", "ci_lower", "ci_upper",
                "point_outside_interval", "method_used", "n_failed",
                "subjects", "not_computed",
            ],
        )

    def summary(self) -> dict:
        """
        Return shared permutation metadata in plain Python types.

        Subclasses extend this with their own statistical content.

        Returns
        -------
        dict
            Keys: coefficient_name, ci_method, weights_type,
            weight_matrix, balance, raters.
        """
        return {
            "coefficient_name": self.coefficient_name,
            "ci_method": self.ci_method,
            "weights_type": self.weights_type,
            "weight_matrix": self.weight_matrix,
            "balance": self.balance_profile.summary(),
            "raters": self.rater_profile.summary(),
        }

    # Printed form, shared by pairs and subsets ------------------------------ #

    def _unit_results(self) -> dict:
        """My pairs or subsets: each one's result, or NotComputed."""
        raise NotImplementedError

    def _units(self) -> list[tuple]:
        """
        (raters, value, interval, outside, subjects, reason) for each pair
        or subset, in order. ``outside`` is the unit's own
        point_outside_interval. A unit with a result has no reason; one
        without a result has only its raters and the reason (its subject
        count is in it).
        """
        from ...not_computed import NotComputed

        rows = []
        for raters, unit in self._unit_results().items():
            if isinstance(unit, NotComputed):
                rows.append((raters, float("nan"), None, None, None, unit.reason))
            else:
                rows.append((
                    raters, float(unit.point), unit.confidence_interval,
                    unit.point_outside_interval,
                    unit.balance_profile.n_subjects, None,
                ))
        return rows

    def _check_unit_results(
        self, results: dict, k: int, n_units: int, mean: float, word: str,
    ) -> None:
        """
        My pairs or subsets are the results they claim to be, and the mean
        is theirs.

        Each unit is the coefficient's own result class or NotComputed, keyed
        by a sorted tuple of k rater indices, built with my interval method.
        The mean must be the mean of the units' points that are defined: a
        mean that disagrees with its own units was built from something else.
        """
        from ...not_computed import NotComputed
        from ..analyzer_results.point_estimate_result import PointEstimateResult
        from ..ci_results.bootstrap_ci_result import BootstrapCIResult

        if not isinstance(results, dict):
            raise TypeError(
                f"{word}_results must be a dict, got {type(results).__name__}"
            )
        if len(results) != n_units:
            raise ValueError(
                f"{word}_results has {len(results)} entries but there are "
                f"{n_units} {word}s"
            )
        # Keys are positions into rater_names, so those must be sound first
        self._check_rater_profile()
        expected = result_class_for(self.coefficient_name)
        points = []
        kept = set()  # whether each computed unit kept its draws
        for key, unit in results.items():
            if (
                not isinstance(key, tuple) or len(key) != k
                or not all(isinstance(i, int) for i in key)
                or list(key) != sorted(set(key))
            ):
                raise ValueError(
                    f"{word}_results keys must be sorted tuples of {k} "
                    f"distinct rater indices, got {key!r}"
                )
            if not all(0 <= i < len(self.rater_names) for i in key):
                raise ValueError(
                    f"{word} {key} names a rater position beyond the "
                    f"{len(self.rater_names)} in rater_names"
                )
            if isinstance(unit, NotComputed):
                continue
            if not isinstance(unit, PointEstimateResult):
                raise TypeError(
                    f"{word} {key} must be a result or NotComputed, got "
                    f"{type(unit).__name__}"
                )
            if type(unit) is not expected:
                raise TypeError(
                    f"{word} {key} is a {type(unit).__name__}, but "
                    f"{self.coefficient_name} {word}s are {expected.__name__}"
                )
            if not isinstance(unit.confidence_interval, BootstrapCIResult):
                raise TypeError(
                    f"{word} {key} has no bootstrap interval; the per-{word} "
                    f"bootstrap always runs"
                )
            if unit.ci_method != self.ci_method:
                raise ValueError(
                    f"{word} {key} was built with ci_method "
                    f"'{unit.ci_method}', not '{self.ci_method}'"
                )
            # The mean is weighted as its units are, so they must agree
            if hasattr(unit, "weights_type") and (
                unit.weights_type != self.weights_type
                or unit.weight_matrix != self.weight_matrix
            ):
                raise ValueError(
                    f"{word} {key} was computed with weights_type "
                    f"'{unit.weights_type}' and its own matrix, but this "
                    f"result says '{self.weights_type}'; the mean is "
                    f"weighted as its {word}s are"
                )
            points.append(float(unit.point))
            kept.add(not isinstance(unit.bootstrap_distribution, NotComputed))

        # keep_draws applies to the whole analysis, so a mix is a bug
        if len(kept) > 1:
            raise ValueError(
                f"some {word}s kept their draws and others did not; "
                f"keep_draws applies to every {word}"
            )

        defined = [p for p in points if not math.isnan(p)]
        theirs = sum(defined) / len(defined) if defined else float("nan")
        if math.isnan(theirs) or math.isnan(mean):
            agree = math.isnan(theirs) and math.isnan(mean)
        else:
            agree = math.isclose(theirs, mean, rel_tol=1e-12, abs_tol=1e-12)
        if not agree:
            raise ValueError(
                f"the mean is {mean!r}, but the {word}s' defined points "
                f"average {theirs!r}"
            )

    def _printed(
        self,
        title: str,
        mean: float,
        mean_ci: BootstrapCIResult | None,
        mean_method: str | None,
        unit: str,
        units: str,
        pool: str,
    ) -> str:
        """
        The printed form of a pairwise or subset result.

            Pairwise Cohen's κ, 6 pairs: mean = 0.496 [0.377, 0.605] (...), correlated bootstrap
              ! independent bootstrap: ...           (independent mean only)
              199 subjects; 4 raters; unweighted; each pair on the subjects both raters rated
                raters      κ     95% CI            subjects
                0, 1    0.539    [0.383, 0.682]       199
                ...

        The mean's interval prints itself, so its fallbacks and failures
        travel with it. Up to MAX_ROWS_PRINTED rows are listed; above that,
        a line points to to_frame(), which holds them all.
        """
        from ....contracts import MAX_ROWS_PRINTED

        if mean_ci is None:
            first = f"{title}: mean = {mean:.3f} (no interval for the mean)"
        else:
            first = f"{title}: mean = {mean:.3f} {mean_ci}, {mean_method} bootstrap"
        lines = [first]
        if mean_method == "independent":
            lines.append(f"  ! {INDEPENDENT_MEAN_CAVEAT}")
        bp = self.balance_profile
        # Weighted as the coefficient's own result says it, e.g. 'quadratic weights'
        weights = ""
        if has_weighted_form(self.coefficient_name):
            from ..analyzer_results.point_estimate_result import PointEstimateResult

            weights = f"{PointEstimateResult._weights_phrase(self.weights_type)}; "
        lines.append(
            f"  {bp.n_subjects:,} subjects; {raters_phrase(self.rater_profile)}; "
            f"{weights}each {unit} on {pool}"
        )

        rows = self._units()
        if len(rows) > MAX_ROWS_PRINTED:
            lines.append(
                f"  {len(rows):,} {units}, too many to list here: "
                f"result.to_frame() has them all, e.g. "
                f'result.to_frame().to_csv("{units}.csv")'
            )
            return "\n".join(lines)

        symbol = result_class_for(self.coefficient_name).SYMBOL
        intervals = [row[2] for row in rows if row[2] is not None]
        level = f"{intervals[0].confidence_level:.0%} CI" if intervals else "CI"
        cells = []
        for raters, value, ci, outside, subjects, reason in rows:
            if ci is None:
                cells.append((
                    self._named(raters), "not computed", "", "", reason,
                ))
                continue
            estimate = "undefined" if math.isnan(value) else f"{value:.3f}"
            bounds = (
                f"[{ci.lower:.3f}, {ci.upper:.3f}]" if ci.is_defined
                else "[undefined]"
            )
            # The unit's own flag: said on its row, since no warning is raised
            notes = "; ".join(
                (["POINT OUTSIDE INTERVAL"] if outside else []) + interval_notes(ci)
            )
            cells.append((
                self._named(raters), estimate, bounds,
                f"{subjects:,}", notes,
            ))
        w_raters = max(len("raters"), *(len(c[0]) for c in cells))
        w_value = max(len(symbol), *(len(c[1]) for c in cells))
        w_ci = max(len(level), *(len(c[2]) for c in cells))
        w_subj = max(len("subjects"), *(len(c[3]) for c in cells))
        lines.append(
            f"    {'raters':<{w_raters}}  {symbol:>{w_value}}  "
            f"{level:<{w_ci}}  {'subjects':>{w_subj}}"
        )
        for raters, estimate, bounds, subjects, notes in cells:
            line = (
                f"    {raters:<{w_raters}}  {estimate:>{w_value}}  "
                f"{bounds:<{w_ci}}  {subjects:>{w_subj}}"
            )
            lines.append(f"{line}  {notes}" if notes else line)

        short = sum(1 for ci in intervals if ci.below_recommended_resamples)
        if short:
            recommended = next(
                ci.recommended_resamples for ci in intervals
                if ci.below_recommended_resamples
            )
            lines.append(
                f"  {short} of {len(intervals)} intervals rest on fewer usable "
                f"draws than the {recommended} recommended"
            )
        return "\n".join(lines)
