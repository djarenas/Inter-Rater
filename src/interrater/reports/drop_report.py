"""What was removed for having too little data, and in what order."""

from __future__ import annotations

from dataclasses import dataclass, field


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
    threshold and taking out a subject can push a rater below its own. So the
    order was a declared choice, and I record which one was made along with
    what it cost.

    Under ``iterate`` I keep every pass separately, because the sequence is the
    diagnostic. A dataset going 1,617 to 1,580 to 1,204 to 340 has cascaded, and
    that is a very different event from a single clean cut to 340 -- yet both
    end at the same number and a total alone cannot tell them apart.

    I also record how the loop ended. Converging means the data settled and the
    result is what the policy asked for. Hitting the pass limit means the
    cascade was still running when it was cut off, and what survived is
    arbitrary rather than principled.

    Attributes
    ----------
    method:
        Which order was used: raters first, subjects first, or iterate.
    passes:
        The sweeps that ran. One entry for the non-iterating methods.
    converged:
        Whether the loop stopped because a pass dropped nothing. False means it
        hit the limit with work still to do, and the result should be treated
        with suspicion.
    dropped_raters:
        Rater name to how many ratings that rater had when removed.
    dropped_subject_ids:
        Every subject removed, for the full rendering. Summarised as a
        distribution otherwise.
    dropped_subject_distribution:
        How many subjects were removed holding how many ratings:
        ``{0: 12, 1: 380}``.
    """

    method:    str
    converged: bool = True
    passes:    tuple[DropPass, ...] = ()

    dropped_raters:      dict[str, int] = field(default_factory=dict)
    dropped_subject_ids: tuple[str | int, ...] = ()
    dropped_subject_distribution: dict[int, int] = field(default_factory=dict)

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
            "method":    self.method,
            "converged": self.converged,
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
        pad = "    "

        if self.dropped_nothing:
            return f"Drops\n{pad}nothing dropped ({self.method})"

        ending = (
            f"converged after {self.n_passes} pass(es)"
            if self.converged
            else f"STOPPED AT THE PASS LIMIT after {self.n_passes}; "
                 f"the cascade had not settled"
        )
        lines = [
            "Drops",
            f"{pad}method     : {self.method} ({ending})",
            f"{pad}raters     : {self.n_raters_dropped} removed",
        ]
        if self.dropped_raters:
            shown = ", ".join(
                f"{name} ({n:,} rating(s))"
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
                f"{n:,} held {k}"
                for k, n in sorted(self.dropped_subject_distribution.items())
            )
            lines.append(f"{pad}             {shown}")

        if verbosity == "full":
            if len(self.passes) > 1:
                for p in self.passes:
                    lines.append(
                        f"{pad}  pass {p.number}: -{p.n_raters} rater(s), "
                        f"-{p.n_subjects:,} subject(s)"
                    )
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
