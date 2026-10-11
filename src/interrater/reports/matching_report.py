"""What happened to the strings, before anything became a number."""

from __future__ import annotations

from dataclasses import dataclass, field


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
    per_rater_null, per_rater_missing, per_rater_unknown:
        The same counts, keyed by rater name.
    subject_missing_distribution:
        How many subjects lost how many ratings: ``{0: 1204, 1: 380, 2: 33}``.
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

    @property
    def n_missing(self) -> int:
        """Cells with no rating: nulls plus every string that meant nothing."""
        return self.n_null + sum(self.missing_matched.values())

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
        lines = [
            "Matching",
            f"{pad}cells      : {self.n_cells:,}",
            f"{pad}no rating  : {self.n_missing:,} ({fill:.1%}), "
            f"of which {self.n_null:,} arrived null",
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
            for rater, n in self.per_rater_missing.items():
                extra = self.per_rater_unknown.get(rater, 0)
                lines.append(
                    f"{pad}  {rater}: {n:,} without a rating"
                    + (f", {extra:,} unknown" if extra else "")
                )
            if self.subject_missing_distribution:
                shown = ", ".join(
                    f"{n} subject(s) lost {k}"
                    for k, n in sorted(self.subject_missing_distribution.items())
                )
                lines.append(f"{pad}  by subject: {shown}")

        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.render()
