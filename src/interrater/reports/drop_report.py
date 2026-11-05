"""What was removed for having too little data, and in what order."""

from __future__ import annotations

from ._wording import counted

from dataclasses import dataclass, field

_ERROR_PREFIX = "[DropReport] Error"


@dataclass(frozen=True)
class DropPass:
    """One sweep of the drop loop."""

    number:     int
    n_raters:   int
    n_subjects: int

    @property
    def dropped_nothing(self) -> bool:
        return self.n_raters == 0 and self.n_subjects == 0


@dataclass(frozen=True)
class DropReport:
    """
    I record what was removed for thinness, and how the removing went.

    Dropping is not symmetric: taking out a rater can push a subject below its
    threshold and taking out a subject can push a rater below its own. So
    sweeps repeat until one drops nothing, and I record what that cost.

    I keep every pass separately, because the sequence is the
    diagnostic. A dataset going 1,617 to 1,580 to 1,204 to 340 has cascaded, and
    that is a very different event from a single clean cut to 340 -- yet both
    end at the same number and a total alone cannot tell them apart.

    Attributes
    ----------
    min_ratings_per_rater, min_ratings_per_subject:
        The rule that was applied: the fewest ratings a rater or a subject
        could keep and stay. Recorded so the report says what it enforced,
        not only what it removed -- "nothing dropped" means little without
        the rule nothing fell below.
    passes:
        The sweeps that ran, the last one dropping nothing. Required: the
        drop loop always runs at least one.
    dropped_raters:
        Rater name to how many ratings that rater had when removed.
    dropped_subject_ids:
        Every subject removed, for the full rendering. Summarised as a
        distribution otherwise.
    dropped_subject_distribution:
        How many subjects were removed holding how many ratings:
        ``{0: 12, 1: 380}``.
    """

    min_ratings_per_rater:   int
    min_ratings_per_subject: int
    passes:    tuple[DropPass, ...]

    dropped_raters:      dict[str, int] = field(default_factory=dict)
    dropped_subject_ids: tuple[str | int, ...] = ()
    dropped_subject_distribution: dict[int, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """
        Refuse to exist unless I describe a drop loop that could have run.

        The cleaner sweeps until a pass drops nothing, so a real report has
        at least one pass, numbered from 1, each dropping something but the
        last, and pass totals that add up to what was dropped.
        """
        def fail(message: str) -> None:
            raise ValueError(
                f"{_ERROR_PREFIX}: {message}. The report was built "
                f"inconsistently; this is a bug, not a property of the data."
            )

        if not self.passes:
            fail("no passes, but the drop loop always runs at least one")
        numbers = [p.number for p in self.passes]
        if numbers != list(range(1, len(numbers) + 1)):
            fail(f"passes are numbered {numbers}, not 1 to {len(numbers)}")
        last = self.passes[-1]
        if not last.dropped_nothing:
            fail(
                f"the last pass dropped {counted(last.n_raters, 'rater')} and "
                f"{counted(last.n_subjects, 'subject')}, but sweeps repeat "
                f"until one drops nothing"
            )
        for p in self.passes[:-1]:
            if p.dropped_nothing:
                fail(
                    f"pass {p.number} dropped nothing, so it would have been "
                    f"the last, but {counted(len(self.passes), 'pass', 'passes')} ran"
                )

        totals = [
            ("the passes' raters", sum(p.n_raters for p in self.passes),
             "dropped_raters", self.n_raters_dropped),
            ("the passes' subjects", sum(p.n_subjects for p in self.passes),
             "dropped_subject_ids", self.n_subjects_dropped),
        ]
        if self.dropped_subject_distribution:
            totals.append((
                "dropped_subject_distribution",
                sum(self.dropped_subject_distribution.values()),
                "dropped_subject_ids", self.n_subjects_dropped,
            ))
        for name, got, against, want in totals:
            if got != want:
                fail(f"{name} total {got:,}, but {against} holds {want:,}")

    @property
    def n_passes(self) -> int:
        return len(self.passes)

    @property
    def n_raters_dropped(self) -> int:
        return len(self.dropped_raters)

    @property
    def n_subjects_dropped(self) -> int:
        return len(self.dropped_subject_ids)

    @property
    def dropped_nothing(self) -> bool:
        return not self.dropped_raters and not self.dropped_subject_ids

    def summary(self) -> dict:
        return {
            "min_ratings_per_rater":   self.min_ratings_per_rater,
            "min_ratings_per_subject": self.min_ratings_per_subject,
            "n_passes":  self.n_passes,
            "passes": [
                {"pass": p.number, "raters": p.n_raters, "subjects": p.n_subjects}
                for p in self.passes
            ],
            "n_raters_dropped":   self.n_raters_dropped,
            "n_subjects_dropped": self.n_subjects_dropped,
            "dropped_raters":      dict(self.dropped_raters),
            "dropped_subject_ids": list(self.dropped_subject_ids),
            "dropped_subject_distribution": {
                int(k): int(v)
                for k, v in self.dropped_subject_distribution.items()
            },
        }

    def render(self, verbosity: str = "summary") -> str:
        pad  = "    "
        rule = (
            f"at least {self.min_ratings_per_rater} ratings per rater and "
            f"{self.min_ratings_per_subject} per subject"
        )

        if self.dropped_nothing:
            return "\n".join([
                "Drops",
                f"{pad}rule       : {rule}",
                f"{pad}no raters or subjects dropped",
            ])

        lines = [
            "Drops",
            f"{pad}rule       : {rule} "
            f"(settled after {counted(self.n_passes, 'pass', 'passes')})",
            f"{pad}raters     : {self.n_raters_dropped} removed",
        ]
        if self.dropped_raters:
            shown = ", ".join(
                f"{name} ({counted(n, 'rating')})"
                for name, n in list(self.dropped_raters.items())[:3]
            )
            more = self.n_raters_dropped - min(3, self.n_raters_dropped)
            lines.append(
                f"{pad}             {shown}"
                + (f", +{more} more" if more > 0 else "")
            )
        lines.append(f"{pad}subjects   : {self.n_subjects_dropped:,} removed")

        if self.dropped_subject_distribution:
            shown = ", ".join(
                f"{counted(n, 'subject')} holding {counted(k, 'rating')}"
                for k, n in sorted(self.dropped_subject_distribution.items())
            )
            lines.append(f"{pad}             {shown}")

        if verbosity == "full":
            if len(self.passes) > 1:
                for p in self.passes:
                    removed = (
                        "nothing removed" if p.dropped_nothing else
                        f"removed {counted(p.n_raters, 'rater')} and "
                        f"{counted(p.n_subjects, 'subject')}"
                    )
                    lines.append(f"{pad}  pass {p.number}: {removed}")
            if self.dropped_subject_ids:
                ids = list(self.dropped_subject_ids)
                shown = ", ".join(map(str, ids[:10]))
                more  = len(ids) - min(10, len(ids))
                lines.append(
                    f"{pad}  ids: {shown}"
                    + (f", +{more:,} more" if more > 0 else "")
                )

        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.render()
