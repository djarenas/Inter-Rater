"""
String work for the cleaner.

Free functions over plain strings. Nothing here knows what a rating task is
or holds state. The config calls into it to check its own consistency before
any data exists; the cleaner calls the same functions on every cell. One rule,
one implementation -- a config validated under one interpretation and matched
under another would be worse than no validation at all.
"""

from __future__ import annotations

from collections.abc import Sequence


def fold(value: str, *, case_sensitive: bool, strip_whitespace: bool) -> str:
    """
    Normalise a string the way it will be matched.

    Under the defaults ``" NA "`` becomes ``"na"``, so a whitespace variant
    need not appear in every list and a value cannot match a category
    case-insensitively while missing the missing-values list.

    Named for ``str.casefold``, which is the comparison-grade lowercase: it
    handles cases ``lower()`` does not, such as German ``ß`` folding to ``ss``.
    """
    if strip_whitespace:
        value = value.strip()
    return value if case_sensitive else value.casefold()


def fold_collisions(
    values: Sequence[str], *, case_sensitive: bool, strip_whitespace: bool
) -> list[tuple[str, ...]]:
    """
    The groups of different strings that ``fold`` makes the same.

    ``["A", "a", "B"]`` gives ``[("A", "a")]`` under case-insensitive
    matching and ``[]`` under case-sensitive. Repeats of the very same string
    are not a collision. Groups keep the order the strings were given in, and
    come back in the order their first member was given.

    Two declared names in one group would be matched as one: every answer
    for either would be coded as the same category.
    """
    groups: dict[str, list[str]] = {}
    for v in values:
        members = groups.setdefault(
            fold(v, case_sensitive=case_sensitive, strip_whitespace=strip_whitespace),
            [],
        )
        if v not in members:
            members.append(v)
    return [tuple(m) for m in groups.values() if len(m) > 1]


def describe_folding(*, case_sensitive: bool, strip_whitespace: bool) -> str:
    """One sentence saying what ``fold`` ignores, for error messages."""
    ignored = [
        what for what, on in (
            ("case", not case_sensitive),
            ("spaces at either end", strip_whitespace),
        ) if on
    ]
    if not ignored:
        return "Matching is exact."
    return f"Matching ignores {' and '.join(ignored)}."


def fold_all(
    values: Sequence[str], *, case_sensitive: bool, strip_whitespace: bool
) -> list[str]:
    """``fold`` over an iterable, in order."""
    return [
        fold(v, case_sensitive=case_sensitive, strip_whitespace=strip_whitespace)
        for v in values
    ]


# ------------------------------------------------------- The boundary #


def build_lookup(
    category_names: tuple[str, ...],
    missing_values: tuple[str, ...],
    value_map:      dict[str, str],
    empty:          int,
    *,
    case_sensitive:   bool,
    strip_whitespace: bool,
) -> tuple[dict[str, int], dict[str, str]]:
    """
    One folded string to one code. This is where strings become integers.

    Returns the lookup and, alongside it, a note of where each entry came
    from -- ``'category'``, ``'missing'`` or ``'repair'`` -- so a report can
    attribute a cell to the rule that caught it rather than only counting
    outcomes.

    Repairs are added last and win. Somebody who writes a repair for a string
    that is also a category name has said something contradictory, but the more
    specific instruction is the one they wrote deliberately.
    """
    def folded(v: str) -> str:
        return fold(v, case_sensitive=case_sensitive, strip_whitespace=strip_whitespace)

    lookup: dict[str, int] = {}
    origin: dict[str, str] = {}

    for index, name in enumerate(category_names):
        lookup[folded(name)] = index
        origin[folded(name)] = "category"

    for marker in missing_values:
        lookup[folded(marker)] = empty
        origin[folded(marker)] = "missing"

    for written, meant in value_map.items():
        lookup[folded(written)] = category_names.index(meant)
        origin[folded(written)] = "repair"

    return lookup, origin
