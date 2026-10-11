"""
Deliberately messy ratings, for exercising a cleaner.

This is not simulated data in the statistical sense. Nothing here has a known
kappa or a rater accuracy you would want to recover. It is clean data that has
then been *damaged* in the specific ways real spreadsheets are damaged, so that
a cleaner can be checked against an answer known in advance.

Every corruption is recorded. A test can therefore assert that the cleaner
found exactly seven unknown codes, not merely that it found some.

The damage is exhaustive rather than realistic: no single export would contain
all of it. What decides whether the file is fatal is the *config*, not the
file. Cleaned permissively -- ``unknown_policy='empty'``,
``unused_categories='keep'`` -- it produces a real Ratings and a report showing
everything that happened. Cleaned with the strict defaults it raises. Two uses,
one fixture.

A config permissive enough to clean this file lives beside it, in
``scruffy_cleaner_config.yaml``. With the strict defaults instead, the same
file raises: first on the unknown codes, and then, were those allowed, on the
category nobody used.

Deliberately absent are the two failures no config can survive: a column name
that does not match the semantics, and a ``value_map`` pointing at a category
that was never declared. Both belong in their own small fixtures, because a
file carrying them can never be cleaned at all.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

# --------------------------------------------------------------- Vocabulary #

SUBJECT_COL = "subject_id"
RATERS      = ("ann", "bob", "cara", "dev", "eve")
CATEGORIES  = ("include", "unsure", "exclude")

# ------------------------------------------------------------- Damage kinds #

CASE       = "case"
"""Right category, wrong capitalisation. Normalisation recovers these."""

WHITESPACE = "whitespace"
"""Right category, stray spaces. Normalisation recovers these."""

REPAIRABLE = "repairable"
"""A shorthand the rating task did not declare, fixable by value_map."""

UNKNOWN    = "unknown"
"""A response nobody declared and no repair covers. A real decision."""

MISSING    = "missing"
"""One of the several ways a spreadsheet says nothing."""

NULL       = "null"
"""A truly empty cell. Not a string, so it can never be listed in a config."""

CASE_VARIANTS = {
    "include": ("Include", "INCLUDE"),
    "unsure":  ("Unsure", "UNSURE"),
    "exclude": ("Exclude", "EXCLUDE"),
}

WHITESPACE_VARIANTS = {
    "include": (" include", "include "),
    "unsure":  (" unsure", "unsure  "),
    "exclude": ("exclude ", "  exclude"),
}

REPAIRABLE_CODES = {"inc": "include", "exc": "exclude", "unk": "unsure"}
"""Shorthand a rater used. The value_map that fixes these is in the docstring
of ``suggested_config`` below."""

UNKNOWN_CODES = ("maybe", "mabye", "?", "TBD")
"""'maybe' was never offered; 'mabye' is a typo of nothing."""

MISSING_STRINGS = ("", "   ", "NA", "N/A", "NOT SHOWN", "-")
"""Blank, whitespace-only, two spellings of not-available, an explicit
not-shown marker, and a dash."""


# ----------------------------------------------------------------- Manifest #


@dataclass(frozen=True)
class Damage:
    """One corrupted cell, and what it used to be."""

    row:      int
    rater:    str
    original: str
    written:  object
    kind:     str


@dataclass
class ScruffyRatings:
    """
    Messy ratings alongside the truth about how they were messed up.

    Attributes
    ----------
    frame:
        What a cleaner is handed. A subject id column and one column per rater.
    clean_frame:
        The same table before damage. Not given to the cleaner; used to check
        that normalisation recovers the original category.
    damage:
        Every corruption applied, in the order applied.
    quit_rater:
        The rater who stops partway through, if one was made to.
    quit_after:
        The row at which they stopped.
    thin_subjects:
        Row indices deliberately left with too few ratings to survive a
        minimum.
    unused_category:
        A category removed from the data entirely, if one was.
    duplicated_ids:
        Subject ids made to appear twice, if any were.
    """

    frame:           pd.DataFrame
    clean_frame:     pd.DataFrame
    damage:          list[Damage]
    quit_rater:      str | None = None
    quit_after:      int | None = None
    thin_subjects:   tuple[int, ...] = ()
    unused_category: str | None = None
    duplicated_ids:  tuple[str, ...] = ()

    # ---- Counting, so tests can assert exactly ---- #

    def count(self, kind: str) -> int:
        return sum(1 for d in self.damage if d.kind == kind)

    def cells(self, kind: str) -> list[Damage]:
        return [d for d in self.damage if d.kind == kind]

    @property
    def n_recoverable(self) -> int:
        """Damage that normalisation alone should undo."""
        return self.count(CASE) + self.count(WHITESPACE)

    @property
    def n_absent(self) -> int:
        """Cells a cleaner should end up treating as no rating."""
        return self.count(MISSING) + self.count(NULL)

    def summary(self) -> str:
        rows, cols = self.frame.shape
        lines = [
            "ScruffyRatings(",
            f"    shape       = {rows} subject(s) x {cols - 1} rater(s)",
            f"    categories  = {list(CATEGORIES)}",
            f"    case        = {self.count(CASE)}",
            f"    whitespace  = {self.count(WHITESPACE)}",
            f"    repairable  = {self.count(REPAIRABLE)}",
            f"    unknown     = {self.count(UNKNOWN)}",
            f"    missing     = {self.count(MISSING)}",
            f"    null        = {self.count(NULL)}",
            f"    total dmg   = {len(self.damage)} / {rows * (cols - 1)} cells",
        ]
        if self.quit_rater:
            lines.append(
                f"    quit        = {self.quit_rater} after row {self.quit_after}"
            )
        if self.thin_subjects:
            lines.append(f"    thin rows   = {list(self.thin_subjects)}")
        if self.unused_category:
            lines.append(f"    unused      = {self.unused_category!r}")
        if self.duplicated_ids:
            lines.append(f"    duplicated  = {list(self.duplicated_ids)}")
        lines.append(")")
        return "\n".join(lines)


# ---------------------------------------------------------------- Generator #


def make_scruffy_ratings(
    n_subjects:   int = 40,
    n_case:       int = 8,
    n_whitespace: int = 8,
    n_repairable: int = 6,
    n_unknown:    int = 7,
    n_missing:    int = 10,
    n_null:       int = 4,
    quit_rater:   str | None = "eve",
    quit_after:   int = 12,
    n_thin:       int = 3,
    drop_category: str | None = "unsure",
    n_duplicate:  int = 2,
    seed:         int = 0,
) -> ScruffyRatings:
    """
    Build a messy wide ratings frame with a manifest of its own defects.

    Prevalence is roughly uniform, which makes the numbers easy to reason about
    while debugging. Real screening data is heavily skewed; that belongs in a
    separate fixture, because skew changes what the coefficients do rather than
    what the cleaner does.

    Note that 'unsure' is an ordinary category, not an absence. It was offered
    and a rater chose it, so it counts like any other response. Absence is what
    the missing and null damage produce.

    Parameters
    ----------
    quit_rater, quit_after:
        Make one rater stop rating partway through, which is what a per-rater
        breakdown exists to reveal: four hundred absent cells spread evenly is
        by design, four hundred in one column is somebody who left.
    n_thin:
        Subjects left with a single rating, so a minimum has something to drop.
    drop_category:
        A declared category removed from the data entirely, so the unused check
        has something to find. Set None to leave every category used.
    n_duplicate:
        Subject ids made to appear twice, which is usually an export bug.

    Examples
    --------
    >>> data = make_scruffy_ratings(seed=0)
    >>> data.count("unknown")
    7
    """
    rng = np.random.default_rng(seed)

    usable = [c for c in CATEGORIES if c != drop_category]
    clean  = pd.DataFrame(
        rng.choice(usable, size=(n_subjects, len(RATERS))), columns=list(RATERS)
    )
    clean.insert(0, SUBJECT_COL, [f"s{i:03d}" for i in range(n_subjects)])

    frame  = clean.copy()
    damage: list[Damage] = []

    def record(row: int, rater: str, written: object, kind: str) -> None:
        damage.append(
            Damage(row, rater, str(clean.at[row, rater]), written, kind)
        )
        frame.at[row, rater] = written

    # ---- Cell-level damage, on cells nothing else will touch ---- #

    protected = set(range(n_thin)) | set(range(quit_after, n_subjects)) \
        if quit_rater else set(range(n_thin))

    free = [
        (r, c)
        for r in range(n_subjects)
        for c in RATERS
        if r >= n_thin and not (quit_rater and c == quit_rater and r >= quit_after)
    ]
    rng.shuffle(free)

    plan = (
        [(CASE, CASE_VARIANTS)] * n_case
        + [(WHITESPACE, WHITESPACE_VARIANTS)] * n_whitespace
        + [(REPAIRABLE, None)] * n_repairable
        + [(UNKNOWN, None)] * n_unknown
        + [(MISSING, None)] * n_missing
        + [(NULL, None)] * n_null
    )
    if len(plan) > len(free):
        raise ValueError(
            f"asked for {len(plan)} damaged cells but only {len(free)} are "
            f"free once the quitting rater and thin subjects are reserved."
        )

    reverse_repair = {v: k for k, v in REPAIRABLE_CODES.items()}

    for (row, rater), (kind, variants) in zip(free, plan):
        original = str(clean.at[row, rater])
        if kind in (CASE, WHITESPACE):
            written = str(rng.choice(variants[original]))
        elif kind == REPAIRABLE:
            written = reverse_repair[original]
        elif kind == UNKNOWN:
            written = str(rng.choice(UNKNOWN_CODES))
        elif kind == MISSING:
            written = str(rng.choice(MISSING_STRINGS))
        else:
            written = np.nan
        record(row, rater, written, kind)

    # ---- A rater who stops ---- #

    if quit_rater:
        for row in range(quit_after, n_subjects):
            record(row, quit_rater, "", MISSING)

    # ---- Subjects too thin to keep ---- #

    for row in range(n_thin):
        for rater in RATERS[1:]:
            record(row, rater, "", MISSING)

    # ---- Ids that repeat ---- #

    duplicated: list[str] = []
    for i in range(n_duplicate):
        target = n_subjects - 1 - i
        frame.at[target, SUBJECT_COL] = frame.at[i + n_thin, SUBJECT_COL]
        duplicated.append(str(frame.at[target, SUBJECT_COL]))

    return ScruffyRatings(
        frame           = frame,
        clean_frame     = clean,
        damage          = damage,
        quit_rater      = quit_rater,
        quit_after      = quit_after if quit_rater else None,
        thin_subjects   = tuple(range(n_thin)),
        unused_category = drop_category,
        duplicated_ids  = tuple(duplicated),
    )
