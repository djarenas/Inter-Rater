"""Base result for permutation analysis (pairwise or subset)."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

from ..analysis_result import AnalysisResult

if TYPE_CHECKING:
    import pandas as pd

    from ...data_profiles.balance_profile import BalanceProfile
    from ..ci_results.bootstrap_ci_result import BootstrapCIResult

_ERROR_PREFIX = "[PermutationResult] Error"

# Shown wherever an "independent" mean interval is shown, so that the text
# and the plots say the same thing about it.
INDEPENDENT_MEAN_CAVEAT = (
    "independent bootstrap: ignores the correlation between combinations "
    "that share raters, so it tends to be too narrow"
)

# Printed name and column symbol for each coefficient
_LABELS = {
    "cohen_kappa": ("Cohen's κ", "κ"),
    "conger_kappa": ("Conger's κ", "κ"),
    "fleiss_kappa": ("Fleiss' κ", "κ"),
    "scotts_pi": ("Scott's π", "π"),
    "gwet_ac1": ("Gwet's AC1", "AC1"),
    "brennan_prediger": ("Brennan-Prediger", "B-P"),
    "krippendorff_alpha": ("Krippendorff's α", "α"),
}


def coefficient_label(name: str) -> tuple[str, str]:
    """The printed name and the column symbol, e.g. ("Cohen's κ", "κ")."""
    return _LABELS.get(name, (name, name))


def raters_phrase(balance_profile: BalanceProfile) -> str:
    """'4 raters', or 'at least 3 raters' when only counts were given."""
    if balance_profile.n_raters is not None and balance_profile.n_raters_exact:
        return f"{balance_profile.n_raters} raters"
    return f"at least {balance_profile.n_raters_lower_bound} raters"


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
        n_subjects, n_raters, balance structure.
    """

    # Attribute declaration

    coefficient_name: str
    ci_method: str
    balance_profile: BalanceProfile

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

    # Public API ------------------------------------------------------------- #

    def to_frame(self) -> pd.DataFrame:
        """
        Every pair or subset as a row of a pandas DataFrame, to save or analyse.

        ``result.to_frame().to_csv("pairs.csv")`` writes them all. The
        printed list is made from the same rows, so the two cannot disagree.

        Returns
        -------
        pd.DataFrame
            Columns: raters ("0, 1" -- positions in the ratings), estimate,
            ci_lower, ci_upper, method_used, n_failed, subjects. Undefined
            values are NaN.
        """
        import pandas as pd

        return pd.DataFrame(
            [
                {
                    "raters": ", ".join(map(str, raters)),
                    "estimate": float(value),
                    "ci_lower": float(ci.lower),
                    "ci_upper": float(ci.upper),
                    "method_used": ci.method_used,
                    "n_failed": int(ci.n_failed),
                    "subjects": int(subjects),
                }
                for raters, value, ci, subjects in self._units()
            ],
            columns=[
                "raters", "estimate", "ci_lower", "ci_upper",
                "method_used", "n_failed", "subjects",
            ],
        )

    def summary(self) -> dict:
        """
        Return shared permutation metadata in plain Python types.

        Subclasses extend this with their own statistical content.

        Returns
        -------
        dict
            Keys: coefficient_name, ci_method, balance.
        """
        return {
            "coefficient_name": self.coefficient_name,
            "ci_method": self.ci_method,
            "balance": self.balance_profile.summary(),
        }

    # Printed form, shared by pairs and subsets ------------------------------ #

    def _units(self) -> list[tuple]:
        """(raters, value, interval, subjects) for each pair or subset, in order."""
        raise NotImplementedError

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
              199 subjects; 4 raters; each pair on the subjects both raters rated
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
        lines.append(
            f"  {bp.n_subjects:,} subjects; {raters_phrase(bp)}; "
            f"each {unit} on {pool}"
        )

        rows = self._units()
        if len(rows) > MAX_ROWS_PRINTED:
            lines.append(
                f"  {len(rows):,} {units}, too many to list here: "
                f"result.to_frame() has them all, e.g. "
                f'result.to_frame().to_csv("{units}.csv")'
            )
            return "\n".join(lines)

        _, symbol = coefficient_label(self.coefficient_name)
        level = f"{rows[0][2].confidence_level:.0%} CI" if rows else "CI"
        cells = []
        for raters, value, ci, subjects in rows:
            estimate = "undefined" if math.isnan(value) else f"{value:.3f}"
            bounds = (
                f"[{ci.lower:.3f}, {ci.upper:.3f}]" if ci.is_defined
                else "[undefined]"
            )
            notes = "; ".join(interval_notes(ci))
            cells.append((
                ", ".join(map(str, raters)), estimate, bounds,
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

        short = sum(1 for _, _, ci, _ in rows if ci.below_recommended_resamples)
        if short:
            recommended = next(
                ci.recommended_resamples for _, _, ci, _ in rows
                if ci.below_recommended_resamples
            )
            lines.append(
                f"  {short} of {len(rows)} intervals rest on fewer usable "
                f"draws than the {recommended} recommended"
            )
        return "\n".join(lines)
