"""Small wording helpers shared by the reports."""

from __future__ import annotations


def counted(n: int, singular: str, plural: str | None = None) -> str:
    """``counted(1, "rater")`` is '1 rater'; ``counted(1617, "row")`` is '1,617 rows'."""
    word = singular if n == 1 else (plural or singular + "s")
    return f"{n:,} {word}"
