"""What happened to the strings, before anything became a number."""

from __future__ import annotations

from ._wording import counted

from dataclasses import dataclass, field

_ERROR_PREFIX = "[MatchingReport] Error"


@dataclass(frozen=True)
class MatchingReport:
    """
    I count what every cell turned into while it was still a string.

    I am the first half of a cleaning report, covering everything that happens
    before a value becomes a category position: cells that arrived null, cells
    whose text meant nothing, cells that were repaired into a category, and
    cells that were none of those and had to be dealt with by policy.

    I count per rater as well as in total, because the totals hide the thing
    most worth seeing. Four hundred missing cells spread evenly across four
    raters is by design; four hundred concentrated in one rater is
    somebody who stopped halfway, and only the breakdown tells them apart.

    For subjects I keep a distribution rather than a row per subject. With
    sixteen hundred abstracts a per-subject listing is unreadable, and the
    shape is the diagnostic anyway: a long tail of subjects each missing one
    rating means something different from a handful missing all of them.

    Attributes
    ----------
    n_cells:
        Rater cells examined, excluding any subject id column.
    n_null:
        Cells that arrived null. Always missing, whatever the config says --
        an empty spreadsheet cell is not a string and could not be listed.
    missing_matched:
        Per entry in ``missing_values``, how many cells it caught. Entries that
        caught nothing are still listed: a declared marker that never fired is
        worth seeing.
    repairs:
        Per entry in ``value_map``, how many cells it fixed. Empty when no
        repairs were configured, and then omitted from the rendering entirely
        -- somebody who declared no repairs is not asking about them.
    unknown_values:
        Strings that were neither missing nor a category, with how often each
        appeared. Non-empty only when ``unknown_policy`` is ``empty``; under
        ``error`` the clean stopped and there is no report to read.
    per_rater_null, per_rater_unknown:
        The null and unknown counts, keyed by rater name.
    per_rater_missing:
        Cells each rater was left without a rating, whatever the reason:
        null, a missing marker, or an unknown treated as empty.
    subject_missing_distribution:
        How many subjects had how many cells with no rating, before any
        dropping: ``{0: 1204, 1: 380, 2: 33}``.

    "No rating" means one thing everywhere in me: a cell that ended up empty,
    whatever the reason. ``n_missing`` is the sum of the three reasons, and
    the per-rater counts and the subject distribution add up to it. I check
    that when I am built: if they disagree the cleaner has a bug, and a
    report whose own lines contradict each other is not one to cite.
    """

    n_cells:         int
    n_null:          int
    missing_matched: dict[str, int] = field(default_factory=dict)
    repairs:         dict[str, int] = field(default_factory=dict)
    unknown_values:  dict[str, int] = field(default_factory=dict)

    per_rater_null:    dict[str, int] = field(default_factory=dict)
    per_rater_missing: dict[str, int] = field(default_factory=dict)
    per_rater_unknown: dict[str, int] = field(default_factory=dict)

    subject_missing_distribution: dict[int, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Refuse to exist if my own counts do not add up."""
        checks = [
            ("per_rater_missing", sum(self.per_rater_missing.values()), self.n_missing),
            ("per_rater_null", sum(self.per_rater_null.values()), self.n_null),
            ("per_rater_unknown", sum(self.per_rater_unknown.values()), self.n_unknown),
            ("subject_missing_distribution",
             sum(k * n for k, n in self.subject_missing_distribution.items()),
             self.n_missing),
        ]
        for name, got, want in checks:
            # An empty breakdown was not recorded, which is not a disagreement.
            if getattr(self, name) and got != want:
                raise ValueError(
                    f"{_ERROR_PREFIX}: {name} adds up to {got:,}, but the "
                    f"total it breaks down is {want:,}. The counts were built "
                    f"inconsistently; this is a bug, not a property of the data."
                )
        if self.n_missing > self.n_cells:
            raise ValueError(
                f"{_ERROR_PREFIX}: {self.n_missing:,} cells without a rating "
                f"out of {self.n_cells:,} cells."
            )

    @property
    def n_matched_missing(self) -> int:
        """Cells caught by a missing marker."""
        return sum(self.missing_matched.values())

    @property
    def n_missing(self) -> int:
        """Cells with no rating: null, a missing marker, or an unknown treated as empty."""
        return self.n_null + self.n_matched_missing + self.n_unknown

    @property
    def n_repaired(self) -> int:
        return sum(self.repairs.values())

    @property
    def n_unknown(self) -> int:
        return sum(self.unknown_values.values())

    def summary(self) -> dict:
        """Plain types only, so this can be written beside the config."""
        return {
            "n_cells":         self.n_cells,
            "n_null":          self.n_null,
            "n_matched_missing": self.n_matched_missing,
            "n_missing":       self.n_missing,
            "n_repaired":      self.n_repaired,
            "n_unknown":       self.n_unknown,
            "missing_matched": dict(self.missing_matched),
            "repairs":         dict(self.repairs),
            "unknown_values":  dict(self.unknown_values),
            "per_rater_null":    dict(self.per_rater_null),
            "per_rater_missing": dict(self.per_rater_missing),
            "per_rater_unknown": dict(self.per_rater_unknown),
            "subject_missing_distribution": {
                int(k): int(v) for k, v in self.subject_missing_distribution.items()
            },
        }

    def render(self, verbosity: str = "summary") -> str:
        pad   = "    "
        fill  = self.n_missing / self.n_cells if self.n_cells else 0.0
        # The parts add up to the headline: null, marker, unknown.
        parts = [
            f"{self.n_null:,} arrived null",
            f"{self.n_matched_missing:,} matched a missing marker",
        ]
        if self.n_unknown:
            parts.append(f"{self.n_unknown:,} unknown treated as empty")
        # Nothing missing is said in one word; the parts would all be zero.
        no_rating = (
            f"{self.n_missing:,} ({fill:.1%}): " + ", ".join(parts)
            if self.n_missing else "0"
        )
        lines = [
            "Matching",
            f"{pad}cells      : {self.n_cells:,}",
            f"{pad}no rating  : {no_rating}",
        ]

        if self.missing_matched:
            shown = ", ".join(
                f"{k!r} {v:,}" for k, v in self.missing_matched.items()
            )
            lines.append(f"{pad}matched    : {shown}")

        # Nothing configured means nothing to report -- not "0 repairs".
        if self.repairs:
            shown = ", ".join(f"{k!r} {v:,}" for k, v in self.repairs.items())
            lines.append(f"{pad}repaired   : {self.n_repaired:,} ({shown})")

        if self.unknown_values:
            top = sorted(
                self.unknown_values.items(), key=lambda kv: -kv[1]
            )[:3]
            shown = ", ".join(f"{k!r} {v:,}" for k, v in top)
            more  = len(self.unknown_values) - len(top)
            lines.append(
                f"{pad}unknown    : {self.n_unknown:,} treated as empty "
                f"({shown}{f', +{more} more' if more > 0 else ''})"
            )

        if verbosity == "full":
            # Every rater has one cell per row, so rated = rows - no rating.
            rows = self.n_cells // len(self.per_rater_missing) if self.per_rater_missing else 0
            for rater, n in self.per_rater_missing.items():
                extra = self.per_rater_unknown.get(rater, 0)
                lines.append(
                    f"{pad}  {rater}: {rows - n:,} rated, {n:,} with no rating"
                    + (f", {extra:,} of them unknown" if extra else "")
                )
            if self.subject_missing_distribution:
                shown = ", ".join(
                    f"{counted(n, 'subject')} with {k}"
                    for k, n in sorted(self.subject_missing_distribution.items())
                )
                lines.append(f"{pad}  no rating per subject (before drops): {shown}")

        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.render()
