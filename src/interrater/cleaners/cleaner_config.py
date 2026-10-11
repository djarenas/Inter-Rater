"""Every cleaning decision, declared in one place and versioned with the rating task."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from ._string_utils import fold
from .frame_layout import FrameLayout
from .insufficient_data_policy import InsufficientDataPolicy

_ERROR_PREFIX = "[CleanerConfig] Error"
_DEFAULT_MISSING_VALUES: tuple[str, ...] = ("",) # The only spelling of nothing safe to assume.
_UNKNOWN_POLICIES = ("error", "empty")
_UNUSED_POLICIES  = ("error", "keep", "drop")


class CleanerConfig:
    """
    I am every decision a cleaner makes, written down before it runs.

    I can serve as part of the methods section of a paper in machine-readable form. Which strings
    meant nothing, which were repaired and into what, what happened to a
    response nobody declared, how thin the data was allowed to get -- all of it
    stated in one file that is versioned with the analysis rather than
    reconstructed from memory later.

    My constructor serves to validate that a YAML (or dictionary) has the required 
    information to build a configuration object for the cleaner. 

    Since as a configuration object I know nothing about the data, it is left for the cleaner
    to check me against semantics to see if some choices are sensible. For example, whether ``value_map`` 
    sends a string to a category that actually exists is a question about a rating task.

    Attributes
    ----------
    frame_layout:
        How to read the frame. See ``FrameLayout``. Mainly deals with whether the dataframe
        has headers to check/keep and/or a column with subjectIDs.
    case_sensitive:
        Whether string matching distinguishes case. Governs everything matched:
        category names, ``missing_values``, and ``value_map`` keys alike, so
        there is no way for a value to match one and miss another.
    strip_whitespace:
        Whether to trim before matching. On, so ``" NA "`` catches without
        needing a whitespace variant in every list.
    missing_values:
        Strings that mean no rating. A separate list rather than entries in
        ``value_map`` because absence is the decision most worth seeing at a
        glance. Nulls are always missing regardless of what is listed here --
        an empty spreadsheet cell arrives as a null, not as a string, so it
        could not be listed anyway.
    value_map:
        Repairs: a string as it appears in the file, mapped to the category it
        was meant to be.
    unknown_policy:
        What to do with a value that survives normalisation and is still not a
        category. ``error`` by default, because a response nobody declared is
        usually a sign the data is not what its owner thinks it is.
    unused_categories:
        What to do when a declared category was never chosen. ``error`` by
        default, for the same reason. ``keep`` is the safe alternative: K is a
        fact about what raters could have said, and dropping a category would
        move any coefficient whose chance term counts them.
    insufficient_data_drops:
        How thin the data may get. See ``InsufficientDataPolicy``.
    """

    # Attributes -------------------------------------------------------------  #

    frame_layout:            FrameLayout
    case_sensitive:          bool
    strip_whitespace:        bool
    missing_values:          tuple[str, ...]
    value_map:               dict[str, str]
    unknown_policy:          str
    unused_categories:       str
    insufficient_data_drops: InsufficientDataPolicy

    # Constructor ------------------------------------------------------------  #

    def __init__(
        self,
        frame_layout:            FrameLayout | None = None,
        case_sensitive:          bool = False,
        strip_whitespace:        bool = True,
        missing_values:          Sequence[str] | None = None,
        value_map:               dict[str, str] | None = None,
        unknown_policy:          str = "error",
        unused_categories:       str = "error",
        insufficient_data_drops: InsufficientDataPolicy | None = None,
    ) -> None:

        self.frame_layout = frame_layout if frame_layout is not None else FrameLayout()
        self.case_sensitive    = case_sensitive
        self.strip_whitespace  = strip_whitespace
        self.missing_values    = tuple(
            _DEFAULT_MISSING_VALUES if missing_values is None
            else missing_values
        )
        self.value_map = dict(value_map) if value_map is not None else {}
        self.unknown_policy    = unknown_policy
        self.unused_categories = unused_categories
        self.insufficient_data_drops = (
            insufficient_data_drops if insufficient_data_drops is not None
            else InsufficientDataPolicy()
        )

        self._validate()

    # Validation -------------------------------------------------------------  #

    def _validate(self) -> None:
        self._validate_delegates()
        self._validate_flags()
        self._validate_missing_values()
        self._validate_value_map()
        self._validate_no_double_claim()
        self._validate_policies()

    def _validate_delegates(self) -> None:
        if not isinstance(self.frame_layout, FrameLayout):
            raise TypeError(
                f"{_ERROR_PREFIX}.frame_layout must be a FrameLayout, got "
                f"{type(self.frame_layout).__name__}."
            )
        if not isinstance(self.insufficient_data_drops, InsufficientDataPolicy):
            raise TypeError(
                f"{_ERROR_PREFIX}.insufficient_data_drops must be an "
                f"InsufficientDataPolicy, got "
                f"{type(self.insufficient_data_drops).__name__}."
            )

    def _validate_flags(self) -> None:
        for field in ("case_sensitive", "strip_whitespace"):
            if not isinstance(getattr(self, field), bool):
                raise TypeError(
                    f"{_ERROR_PREFIX}.{field} must be a boolean, got "
                    f"{type(getattr(self, field)).__name__}."
                )

    def _validate_missing_values(self) -> None:
        bad = [v for v in self.missing_values if not isinstance(v, str)]
        if bad:
            raise TypeError(
                f"{_ERROR_PREFIX}.missing_values must hold strings, got "
                f"{[type(v).__name__ for v in bad[:3]]}. A null cell is always "
                f"missing and cannot be listed here."
            )

    def _validate_value_map(self) -> None:
        if not isinstance(self.value_map, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}.value_map must be a mapping of file string "
                f"to category name, got {type(self.value_map).__name__}."
            )
        for key, value in self.value_map.items():
            if not isinstance(key, str) or not isinstance(value, str):
                raise TypeError(
                    f"{_ERROR_PREFIX}.value_map must map string to string, but "
                    f"{key!r} -> {value!r} is "
                    f"{type(key).__name__} -> {type(value).__name__}. To mark "
                    f"a string as absent, list it in missing_values instead."
                )
            if not key:
                raise ValueError(
                    f"{_ERROR_PREFIX}.value_map has an empty key. An empty "
                    f"cell is handled by missing_values, not here."
                )

    def _validate_no_double_claim(self) -> None:
        """
        A string cannot both mean nothing and mean something.

        Folded before comparing, using the same rule the cleaner will apply to
        data, so a config that validates here cannot behave differently there.
        """
        missing = {self._fold(v) for v in self.missing_values}
        claimed = sorted(
            key for key in self.value_map if self._fold(key) in missing
        )
        if claimed:
            raise ValueError(
                f"{_ERROR_PREFIX} claims {len(claimed)} string(s) twice "
                f"(e.g. {claimed[:3]}): they appear in missing_values and are "
                f"also mapped by value_map. Decide whether each means no "
                f"rating or a category, and list it once."
                + ("" if self.case_sensitive else
                   " Matching is case-insensitive, so entries differing only "
                   "in case collide.")
            )

    def _validate_policies(self) -> None:
        if self.unknown_policy not in _UNKNOWN_POLICIES:
            raise ValueError(
                f"{_ERROR_PREFIX}.unknown_policy must be one of "
                f"{list(_UNKNOWN_POLICIES)}, got {self.unknown_policy!r}. "
                f"'error' stops on a response nobody declared; 'empty' treats "
                f"it as no rating."
            )
        if self.unused_categories not in _UNUSED_POLICIES:
            raise ValueError(
                f"{_ERROR_PREFIX}.unused_categories must be one of "
                f"{list(_UNUSED_POLICIES)}, got "
                f"{self.unused_categories!r}."
            )

    def _fold(self, value: str) -> str:
        """My matching rule, applied to one string."""
        return fold(
            value,
            case_sensitive   = self.case_sensitive,
            strip_whitespace = self.strip_whitespace,
        )

    # Classmethods -----------------------------------------------------------  #

    @classmethod
    def from_dict(cls, spec: dict) -> CleanerConfig:
        spec   = dict(spec)
        layout = spec.pop("frame_layout", None)
        drops  = spec.pop("insufficient_data_drops", None)
        known  = {
            "case_sensitive", "strip_whitespace", "missing_values",
            "value_map", "unknown_policy", "unused_categories",
        }
        unknown = sorted(set(spec) - known)
        if unknown:
            raise ValueError(
                f"{_ERROR_PREFIX} has unknown key(s) {unknown}. Known: "
                f"{sorted(known)}, frame_layout, insufficient_data_drops."
            )

        try:
            frame_layout = FrameLayout.from_dict(layout) if layout else None
        except (KeyError, TypeError, ValueError) as error:
            raise type(error)(f"{_ERROR_PREFIX}: {error}") from error

        try:
            insufficient = (
                InsufficientDataPolicy.from_dict(drops) if drops else None
            )
        except (KeyError, TypeError, ValueError) as error:
            raise type(error)(f"{_ERROR_PREFIX}: {error}") from error

        return cls(
            frame_layout            = frame_layout,
            insufficient_data_drops = insufficient,
            **spec,
        )

    @classmethod
    def build_from_yaml(cls, path: str | Path) -> CleanerConfig:
        """
        Build from a YAML file:

        .. code-block:: yaml
        EXAMPLE:
            frame_layout:
              use_column_names: true
              subject_id_col: abstract_id
            case_sensitive: false
            missing_values: ["", "NA", "NOT SHOWN"]
            value_map:
              inc: include, 1
              exc: exclude, -1
            unknown_policy: error
            unused_categories: error
            insufficient_data_drops:
              min_ratings_per_rater: 2
              min_ratings_per_subject: 2
              method: raters_first
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

    # Comparison and display -------------------------------------------------  #

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CleanerConfig):
            return NotImplemented
        return all(
            getattr(self, field) == getattr(other, field)
            for field in (
                "frame_layout", "case_sensitive", "strip_whitespace",
                "missing_values", "value_map", "unknown_policy",
                "unused_categories", "insufficient_data_drops",
            )
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        drops = self.insufficient_data_drops
        return (
            f"CleanerConfig(\n"
            f"    frame_layout      = {self.frame_layout}\n"
            f"    case_sensitive    = {self.case_sensitive}\n"
            f"    strip_whitespace  = {self.strip_whitespace}\n"
            f"    missing_values    = {list(self.missing_values)}\n"
            f"    value_map         = {len(self.value_map)} repair(s)\n"
            f"    unknown_policy    = {self.unknown_policy!r}\n"
            f"    unused_categories = {self.unused_categories!r}\n"
            f"    insufficient_data_drops = per_rater>={drops.min_ratings_per_rater}, "
            f"per_subject>={drops.min_ratings_per_subject}, "
            f"{drops.method!r}\n"
            f")"
        )