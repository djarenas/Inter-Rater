"""What a cleaner did, in a form you can read or file away."""

from __future__ import annotations

from ._wording import counted

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from ..semantics.categories import DistanceCategories
from ..semantics.ratings_semantics import RatingsSemantics
from .category_report import CategoryReport
from .drop_report import DropReport
from .matching_report import MatchingReport

if TYPE_CHECKING:  # pragma: no cover
    # Annotation only. A cleaner imports its reports, so importing the config
    # back at runtime would close a cycle between the two packages -- and a
    # report is meant to be readable without the machinery that made it.
    from ..cleaners.cleaner_config import CleanerConfig

_VERBOSITIES = ("summary", "full")


@dataclass(frozen=True)
class CleaningReport:
    """
    I am the receipt for a clean. The soup does not remember the carrots; I do.

    A ``Ratings`` describes its own contents and nothing about where they came
    from, which is what keeps it comparable to any other ratings with the same
    values. Everything that was normalised, repaired, treated as absent or
    dropped along the way lives in me instead, beside the config that decided
    it -- so a methods section can be read off rather than reconstructed.

    I always hold everything. Verbosity is a property of how I am printed,
    never of what I computed: two runs of the same data must produce the same
    record, or a reader cannot tell whether a section is missing because
    nothing happened or because somebody turned it off. Counting is cheap next
    to the work of matching every cell, so there is nothing to buy back by
    computing less.

    I hold two declarations, and the difference between them is the thing I
    exist to make impossible to miss. A cleaner that drops a rater has produced
    data the original declaration no longer describes: five raters were
    declared, four survived, and the Ratings must carry a declaration matching
    what it holds. That is not the cleaner overruling anybody -- it is the
    consequence of a threshold the user set -- but it does mean
    ``ratings.semantics`` is not the object they handed in, and nothing
    downstream will say so again.

    I am frozen. I describe a run that already finished.

    Attributes
    ----------
    config:
        The rules that were applied. Without it my numbers say what happened
        but not why, and the two together are the whole account.
    declared_semantics:
        The rating task as written, before any row was read.
    effective_semantics:
        The rating task the surviving data actually describes. Identical to the
        declared one unless a rater or a category was dropped.
    shape_in:
        Rows and columns of the frame handed to the cleaner.
    shape_out:
        Subjects and raters of the Ratings that came out.
    matching, drops, categories:
        The three stages. See ``MatchingReport``, ``DropReport``,
        ``CategoryReport``.
    """

    # Attributes ---------------------------------------------------------  #

    config:              "CleanerConfig"
    declared_semantics:  RatingsSemantics
    effective_semantics: RatingsSemantics
    shape_in:            tuple[int, int]
    shape_out:           tuple[int, int]
    matching:            MatchingReport
    drops:               DropReport
    categories:          CategoryReport

    # Public API ---------------------------------------------------------  #

    @property
    def n_subjects_lost(self) -> int:
        return self.shape_in[0] - self.shape_out[0]

    @property
    def n_raters_lost(self) -> int:
        return self.drops.n_raters_dropped

    @property
    def survived(self) -> float:
        """Fraction of input rows that made it into the Ratings."""
        return self.shape_out[0] / self.shape_in[0] if self.shape_in[0] else 0.0

    @property
    def semantics_changed(self) -> bool:
        """
        Whether the surviving data describes a different rating task.

        True when a rater or a category did not survive. Subjects never count:
        they are not declared, so removing them narrows the data without
        narrowing the task.
        """
        return self.declared_semantics != self.effective_semantics

    @property
    def k_changed(self) -> bool:
        """
        Whether the number of categories moved.

        Worth its own question, because it is the change with consequences for
        the arithmetic rather than only for the bookkeeping.
        """
        return (
            self.declared_semantics.n_categories
            != self.effective_semantics.n_categories
        )

    def summary(self) -> dict:
        """
        Everything I hold, in plain types.

        Complete at every verbosity: this is the machine-readable record, and
        it would be useless if a display setting could change it. It carries
        the cleaner config that produced it and the package version, so the
        record alone says how the data was cleaned and by which release.

        ``effective_semantics`` comes out in the same shape a semantics YAML
        takes, so somebody who accepts the change can lift it straight into
        their declaration rather than editing by hand.
        """
        try:
            from .. import __version__ as version
        except ImportError:
            version = "unknown"
        return {
            "_package_version": version,
            "cleaner_config":   self.config.to_dict(),
            "shape_in":  list(self.shape_in),
            "shape_out": list(self.shape_out),
            "n_subjects_lost": self.n_subjects_lost,
            "n_raters_lost":   self.n_raters_lost,
            "survived":        round(self.survived, 6),
            "semantics_changed":   self.semantics_changed,
            "declared_semantics":  self.declared_semantics.to_dict(),
            "effective_semantics": self.effective_semantics.to_dict(),
            "matching":   self.matching.summary(),
            "drops":      self.drops.summary(),
            "categories": self.categories.summary(),
        }

    def to_yaml(self, path: str | Path) -> None:
        """Write my summary beside the config that produced it."""
        import yaml

        with open(path, "w", encoding="utf-8") as handle:
            yaml.safe_dump(self.summary(), handle, sort_keys=False)

    # Helper methods -----------------------------------------------------  #

    def _render_semantics(self) -> list[str]:
        """
        The declaration the data ended up describing, when it is not the one
        that was written.

        Silent when nothing changed, on the same rule as every other section:
        somebody whose declaration survived intact is not asking about it.

        The reasons live in the other sections -- which rater, which category,
        and why -- so this reads across into them rather than storing a second
        copy that could disagree.
        """
        if not self.semantics_changed:
            return []

        pad      = "    "
        declared = self.declared_semantics
        effective = self.effective_semantics
        lines = [
            "Semantics",
            f"{pad}declared   : {counted(declared.n_raters, 'rater')} "
            f"{list(declared.raters)}, K={declared.n_categories}",
            f"{pad}effective  : {counted(effective.n_raters, 'rater')} "
            f"{list(effective.raters)}, K={effective.n_categories}",
        ]

        if self.drops.dropped_raters:
            shown = ", ".join(
                f"{name} ({counted(n, 'rating')})"
                for name, n in self.drops.dropped_raters.items()
            )
            lines.append(
                f"{pad}raters     : {shown} fell below "
                f"min_ratings_per_rater and "
                f"{'was' if len(self.drops.dropped_raters) == 1 else 'were'} dropped"
            )

        if self.k_changed:
            lines.append(
                f"{pad}CATEGORIES : {list(self.categories.unused)} was never "
                f"chosen and has been dropped. K is now "
                f"{effective.n_categories}, not the {declared.n_categories} "
                f"you declared."
            )
            lines.append(
                f"{pad}             Gwet's AC1 and Brennan-Prediger put the "
                f"number of categories in their chance term, so values "
                f"computed here will not match the same data with the "
                f"category kept. Kappa is unaffected."
            )
            if isinstance(effective.categories, DistanceCategories):
                lines.append(
                    f"{pad}             The scale carries distances, so every "
                    f"weight has moved as well: categories that were two "
                    f"apart may now be neighbours."
                )

        lines += [
            "",
            f"{pad}Your Ratings carries the effective declaration, not the one "
            f"you wrote.",
            f"{pad}The effective declaration is in summary() in the shape a "
            f"semantics",
            f"{pad}YAML takes, if you want to adopt it.",
        ]
        return lines

    # Display ------------------------------------------------------------  #

    def render(self, verbosity: str = "summary") -> str:
        if verbosity not in _VERBOSITIES:
            raise ValueError(
                f"verbosity must be one of {list(_VERBOSITIES)}, got "
                f"{verbosity!r}."
            )
        rows_in, cols_in   = self.shape_in
        rows_out, cols_out = self.shape_out
        head = [
            "CleaningReport",
            f"    in         : {counted(rows_in, 'row')} x {counted(cols_in, 'column')}",
            f"    out        : {counted(rows_out, 'subject')} x {counted(cols_out, 'rater')} "
            f"({self.survived:.1%} of rows kept)",
            "",
        ]
        body = [
            self.matching.render(verbosity),
            "",
            self.drops.render(verbosity),
            "",
            self.categories.render(verbosity),
        ]
        semantics = self._render_semantics()
        if semantics:
            body += [""] + semantics
        return "\n".join(head + body)

    def __repr__(self) -> str:
        return self.render()
