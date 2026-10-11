"""What a ratings array means, declared before the data is seen."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from ..contracts import MIN_RATERS, MIN_SUBJECTS
from .categories import Categories, build_categories

_ERROR_PREFIX = "[RatingsSemantics] Error"


class RatingsSemantics:
    """
    I say what a ratings array means. A numpy array does not know what it is.

    My biggest contribution is that I can be built from a YAML file: one
    declaration that constructs me, is versioned with the analysis, and doubles
    as the documentation of what the rating task was.

    I name the columns, I own the categories that give the numbers meaning, and
    I record what the data was expected to look like. I hold no ratings and do
    no counting.

    I am written before the data is read, not derived from it. That is the
    point of me: if a spreadsheet arrives with five rater columns when four
    raters were hired, the fifth is a mistake to be reported, not a rater to
    be invented. Declaring the rating task up front is what turns a surprise into an
    error message.

    I am reusable. The same four raters working a second dataset need the
    same me, which is why I know nothing about any particular dataset -- no
    subject ids, no counts, nothing a cleaner discovered. My ``n_subjects`` is
    an expectation you may state and I may check, never a field anything fills
    in later.

    Attributes
    ----------
    raters:
        Column labels, in column order, never sorted. Given a count instead of
        names I generate ``rater0, rater1, ...``, so a label is file-relative:
        ``rater0`` means "the first rater column of whatever was fed in", not a
        person who carries across rating tasks. Numbering starts at zero so a label
        matches its column index.
    categories:
        What the integers in the array mean, and how far apart they are. See
        ``Categories``.
    n_subjects:
        How many rows the data should have, or None for no expectation. Stated,
        it is checked; absent, nothing is checked. Never populated after the
        fact -- that would tie me to one dataset and end my reuse.
    unique_subjects:
        Whether repeated subject ids are an error. True by default, because a
        repeated id can be frequently due to an export bug. Set False for a rating task that
        deliberately re-shows subjects to measure intra-rater consistency.
        Bootstrap draws repeat rows too, but they bypass validation entirely
        and never consult me.

    Examples
    --------
    >>> from interrater.semantics import UnorderedCategories
    >>> sem = RatingsSemantics(
    ...     raters     = ["ann", "bob", "cara", "dev"],
    ...     categories = UnorderedCategories(("include", "unsure", "exclude")),
    ... )
    >>> sem.n_raters
    4
    >>> RatingsSemantics(raters=3, categories=sem.categories).raters
    ('rater0', 'rater1', 'rater2')
    """

    # Attributes -------------------------------------------------------------  #

    raters:          tuple[str, ...]
    # Always a tuple, even when a count was passed: __init__ generates the
    # names immediately, so nothing downstream ever sees an int.
    categories:      Categories
    n_subjects:      int | None
    unique_subjects: bool

    # Constructor ------------------------------------------------------------  #

    def __init__(
        self,
        raters:          Sequence[str] | int,
        categories:      Categories,
        n_subjects:      int | None = None,
        unique_subjects: bool = True,
    ) -> None:
        self._check_constructor_inputs(raters, categories, n_subjects, unique_subjects)

        self.raters          = self._as_rater_names(raters)
        self.categories      = categories
        self.n_subjects      = n_subjects
        self.unique_subjects = unique_subjects
        self._validate()

    # Class methods ----------------------------------------------------------  #

    @classmethod
    def from_dict(cls, spec: dict) -> RatingsSemantics:
        """
        Build from a plain mapping, as parsed from YAML.

        ``raters`` is either a list of names or an integer count.
        ``categories`` is a mapping dispatched by ``kind``; see
        ``build_categories``.
        """
        missing = {"raters", "categories"} - set(spec)
        if missing:
            raise KeyError(
                f"{_ERROR_PREFIX}: config is missing {sorted(missing)}. A "
                f"rating task must declare who rated and what they could say."
            )

        # build_categories names itself but not the config it was reached
        # through, and somebody debugging a YAML file needs to know which
        # block was wrong.
        try:
            categories = build_categories(spec["categories"])
        except (KeyError, TypeError, ValueError) as error:
            raise type(error)(f"{_ERROR_PREFIX}: {error}") from error

        return cls(
            raters          = spec["raters"],
            categories      = categories,
            n_subjects      = spec.get("n_subjects"),
            unique_subjects = spec.get("unique_subjects", True),
        )

    @classmethod
    def build_from_yaml(cls, path: str | Path) -> RatingsSemantics:
        """
        Build from a YAML file.

        The file is the rating task declaration, versioned alongside the analysis:

        .. code-block:: yaml

            raters: [ann, bob, cara, dev]
            categories:
              kind: ordered
              names: [mild, moderate, severe]
              weighting: quadratic
            n_subjects: 1617
            unique_subjects: true
        """
        import yaml

        with open(path, "r", encoding="utf-8") as handle:
            spec = yaml.safe_load(handle)
        if not isinstance(spec, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: {path} must contain a YAML mapping, got "
                f"{type(spec).__name__}."
            )
        return cls.from_dict(spec)

    # Public API -------------------------------------------------------------  #

    @property
    def n_raters(self) -> int:
        return len(self.raters)

    @property
    def n_categories(self) -> int:
        """K, owned by my categories. Never inferred from any array."""
        return self.categories.n_categories

    def to_dict(self) -> dict:
        """
        The mapping ``from_dict`` would rebuild me from.

        The inverse of construction. A cleaner that had to drop a rater cannot
        hand back a corrected declaration unless I can write myself out, and a
        user should never have to reconstruct one by hand from a report.

        Every field is emitted, including ``n_subjects: null``. The result is
        meant to be read and edited, and a field that exists should be visible
        rather than inferred from its absence.

        Generated rater names are emitted as names. Somebody who declared a
        count of three gets ``[rater0, rater1, rater2]`` back, which is what
        those names always were.
        """
        return {
            "raters":          list(self.raters),
            "categories":      self.categories.to_dict(),
            "n_subjects":      self.n_subjects,
            "unique_subjects": self.unique_subjects,
        }

    def index_of_rater(self, rater: str) -> int:
        """Which column holds ``rater``."""
        try:
            return self.raters.index(rater)
        except ValueError:
            raise KeyError(
                f"{_ERROR_PREFIX}: has no rater {rater!r}. Known: "
                f"{list(self.raters)}."
            ) from None

    # Comparison and display -------------------------------------------------  #

    def __eq__(self, other: object) -> bool:
        """
        Two of me are equal when we declare the same rating task.

        ``other`` is annotated ``object`` and returns ``NotImplemented`` rather
        than raising, because that is the protocol: Python then tries the
        reflected comparison and falls back to False. Raising here would break
        ``sem in [a, b]``, dict comparison, and any assertion that compares me
        with something else -- all of which are supposed to answer False, not
        explode. A test helper is the place for a loud, diffing comparison.
        """
        if not isinstance(other, RatingsSemantics):
            return NotImplemented
        return (
            self.raters == other.raters
            and self.categories == other.categories
            and self.n_subjects == other.n_subjects
            and self.unique_subjects == other.unique_subjects
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        expected = "any" if self.n_subjects is None else f"{self.n_subjects:,}"
        return (
            f"RatingsSemantics(\n"
            f"    raters          = {list(self.raters)}\n"
            f"    categories      = {type(self.categories).__name__}"
            f"{list(self.categories.names)}\n"
            f"    n_subjects      = {expected}\n"
            f"    unique_subjects = {self.unique_subjects}\n"
            f")"
        )

    # Helper methods ---------------------------------------------------------  #

    @staticmethod
    def _check_constructor_inputs(
        raters:          Sequence[str] | int,
        categories:      Categories,
        n_subjects:      int | None,
        unique_subjects: bool,
    ) -> None:
        """
        Check that what was handed in is the right kind of thing.

        Separate from ``_validate``, which asks whether what I hold makes
        sense: whether there are enough raters, whether any name repeats. This
        asks the earlier question -- whether the arguments are even the types
        they claim to be -- and runs first, because a bad type turns those
        later checks into confusing failures rather than clear ones.

        Three traps make this less obvious than it looks.

        A ``str`` is a ``Sequence``, so ``raters="ann"`` would pass a plain
        Sequence check and quietly become three raters called 'a', 'n' and 'n'.
        It is caught by name.

        A ``bool`` is an ``int``, so ``raters=True`` would pass an int check
        and mean one rater. Both integer fields exclude bools explicitly.

        ``isinstance`` takes a tuple of types, never a list.
        """
        name = "raters"

        if isinstance(raters, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: raters must be a sequence of names or a "
                f"count, got the string {raters!r}. A string is a sequence of "
                f"characters, so this would have become "
                f"{len(raters)} rater(s) named "
                f"{[c for c in raters[:3]]}. Wrap it in a list if you meant "
                f"one name -- though one rater is not an agreement task."
            )

        if isinstance(raters, bool) or not isinstance(raters, (int, Sequence)):
            raise TypeError(
                f"{_ERROR_PREFIX}: raters must be a sequence of names or an "
                f"integer count, got {type(raters).__name__}."
            )

        if not isinstance(raters, int):
            bad = [r for r in raters if not isinstance(r, str)]
            if bad:
                raise TypeError(
                    f"{_ERROR_PREFIX}: raters must hold strings, got "
                    f"{[type(r).__name__ for r in bad[:3]]} "
                    f"(e.g. {bad[:3]}). Rater names label columns; a number "
                    f"here is usually a header row that was read as data."
                )

        if not isinstance(categories, Categories):
            raise TypeError(
                f"{_ERROR_PREFIX}: categories must be a Categories, got "
                f"{type(categories).__name__}. Use UnorderedCategories for a "
                f"nominal scheme, OrderedCategories for a ranked one, or "
                f"WeightedCategories to supply distances directly."
            )

        if n_subjects is not None and (
            isinstance(n_subjects, bool) or not isinstance(n_subjects, int)
        ):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_subjects must be an integer or None, got "
                f"{type(n_subjects).__name__}. None means no expectation."
            )

        if not isinstance(unique_subjects, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: unique_subjects must be a boolean, got "
                f"{type(unique_subjects).__name__}."
            )

    @staticmethod
    def _as_rater_names(raters: Sequence[str] | int) -> tuple[str, ...]:
        """A count becomes generated names; names are taken as given."""
        if isinstance(raters, int):
            if raters < MIN_RATERS:
                raise ValueError(
                    f"{_ERROR_PREFIX}: needs at least {MIN_RATERS} raters, "
                    f"got {raters}. Agreement is undefined for a single rater."
                )
            return tuple(f"rater{i}" for i in range(raters))
        return tuple(str(name) for name in raters)

    def _validate(self) -> None:
        """
        Whether what I hold makes sense.

        Types were checked before anything was stored; this asks the later
        question of whether the values work together.
        """
        if len(self.raters) < MIN_RATERS:
            raise ValueError(
                f"{_ERROR_PREFIX}: needs at least {MIN_RATERS} raters, got "
                f"{len(self.raters)}. Agreement is undefined for a single "
                f"rater."
            )

        seen: set[str] = set()
        dupes = [r for r in self.raters if r in seen or seen.add(r)]
        if dupes:
            raise ValueError(
                f"{_ERROR_PREFIX}: raters contains {len(dupes)} duplicate "
                f"name(s) "
                f"(e.g. {dupes[:3]}). Two columns naming the same rater would "
                f"be counted as independent by every coefficient."
            )

        blank = [i for i, r in enumerate(self.raters) if not r.strip()]
        if blank:
            raise ValueError(
                f"{_ERROR_PREFIX}: raters has {len(blank)} blank name(s) at "
                f"position(s) "
                f"{blank[:3]}. Give every column a name, or pass a count and "
                f"let me generate rater0, rater1, ..."
            )

        if self.n_subjects is not None and self.n_subjects < MIN_SUBJECTS:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_subjects must be at least {MIN_SUBJECTS} "
                f"if "
                f"stated, got {self.n_subjects}. Pass None for no expectation."
            )

