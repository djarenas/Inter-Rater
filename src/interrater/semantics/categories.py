"""What a rater was allowed to say, and how far apart those answers are."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np

_MIN_CATEGORIES = 2


class Categories(ABC):
    """
    I am the set of responses a rater was allowed to give.

    I am the rating task's design, not a summary of the data. I list what was offered,
    never what was chosen. If I name five grades and nobody earned a D, I still
    have five categories, because the raters could have said D and did not --
    and any coefficient whose chance term divides by the number of categories
    would move if I quietly shrank. This may not be relevant to all agreement statistics,
    but it is relevant to some.

    I own K. Nothing else in the package may infer it: the largest code present
    in a dataset is a lower bound on K and never K itself, and that gap is
    exactly where a bootstrap draw goes wrong when a rare category happens to
    miss a resample.

    I hold no ratings and do no counting. My position in a list is the integer
    a ratings array stores, so my order is my encoding.

    Attributes
    ----------
    names:
        The allowed responses, in encoding order. Position i is code i.
        Always text: a number such as ``1`` or ``-1`` (what YAML gives for
        unquoted ``1`` and ``-1``) is stored as ``'1'`` and ``'-1'``, because
        ratings are matched as strings. A boolean is refused rather than
        converted: YAML reads unquoted ``yes``/``no`` as True/False, and
        ``'True'`` is not what anyone meant to declare.
    labels:
        Optional display names. Cosmetic; never used in arithmetic.
    """

    # Attributes -------------------------------------------------------------  #

    names:  tuple[str, ...]
    labels: dict[str, str]

    # Constructor ------------------------------------------------------------  #

    def __init__(
        self,
        names:  object,
        labels: dict[str, str] | None = None,
    ) -> None:
        self.names  = tuple(_as_name(n, type(self).__name__) for n in names)
        self.labels = {
            _as_name(k, type(self).__name__): v
            for k, v in dict(labels or {}).items()
        }
        self._validate()

    def _validate(self) -> None:
        """Every check I make. Subclasses extend, never replace."""
        name = type(self).__name__

        if len(self.names) < _MIN_CATEGORIES:
            raise ValueError(
                f"{name} needs at least {_MIN_CATEGORIES} categories, got "
                f"{len(self.names)}. A single category cannot produce "
                f"agreement or disagreement."
            )

        seen: set[str] = set()
        dupes = [n for n in self.names if n in seen or seen.add(n)]
        if dupes:
            raise ValueError(
                f"{name}.names contains {len(dupes)} duplicate entr(ies) "
                f"(e.g. {dupes[:3]}). Position is the integer encoding, so a "
                f"duplicate would make two codes mean the same thing."
            )

        unknown = sorted(set(self.labels) - set(self.names))
        if unknown:
            raise ValueError(
                f"{name}.labels names {len(unknown)} categor(ies) that are not "
                f"in this set (e.g. {unknown[:3]}). Add them to names, or drop "
                f"them from labels."
            )

    # Public API -------------------------------------------------------------  #

    @property
    def n_categories(self) -> int:
        """K. What the raters could have said, not what they did say."""
        return len(self.names)

    def index_of(self, category: str) -> int:
        """The integer a ratings array stores for ``category``."""
        try:
            return self.names.index(category)
        except ValueError:
            raise KeyError(
                f"{type(self).__name__} has no category {category!r}. Known: "
                f"{list(self.names)}."
            ) from None

    def label_for(self, category: str) -> str:
        """Display name, falling back to the category itself."""
        return self.labels.get(category, category)

    def to_dict(self) -> dict:
        """
        The mapping ``build_categories`` would rebuild me from.

        The inverse of construction, so a rating task's declaration can be
        written back out -- which is what lets a cleaner offer corrected
        semantics for the user's records rather than asking them to edit by
        hand. Subclasses add whatever they own.
        """
        spec: dict = {"kind": self._KIND, "names": list(self.names)}
        if self.labels:
            spec["labels"] = dict(self.labels)
        return spec

    def _eq_extras(self, other: "Categories") -> bool:
        """
        Whatever a subclass owns beyond type and names.

        Names alone are not identity. Two ordered scales over the same names
        with different weighting produce different distances between every
        pair, and comparing equal would make a round trip through to_dict look
        correct while silently losing the weighting.
        """
        return True

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Categories):
            return NotImplemented
        if type(self) is not type(other) or self.names != other.names:
            return False
        return self._eq_extras(other)

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(\n"
            f"    n_categories = {self.n_categories}\n"
            f"    names        = {list(self.names)}\n"
            f")"
        )


class UnorderedCategories(Categories):
    """
    I am categories with no ranking and no notion of distance.

    Nominal. 'include', 'unsure' and 'exclude' are three different answers, and
    nothing in me says 'unsure' sits between the other two -- even when a human
    reading the words would assume it does. If you want that assumption in the
    arithmetic, say so with OrderedCategories.

    Every disagreement I describe is equally wrong, which is what unweighted
    kappa already computes. A weighted coefficient handed me would reduce to
    the unweighted one, so analysers that need distance should refuse me rather
    than silently return the same number under a different name.
    """

    _KIND = "unordered"

    names:  tuple[str, ...]
    labels: dict[str, str]


class DistanceCategories(Categories, ABC):
    """
    I am categories that know how wrong each disagreement is.

    My children disagree about where that knowledge comes from -- derived from
    rank, or supplied outright -- but both can answer the same question, so a
    weighted coefficient can take either and need not care which.

    I exist because that question has exactly one form: a K x K matrix of
    agreement weights, 1.0 on the diagonal, falling towards 0.0 as two
    categories grow further apart.
    """

    @abstractmethod
    def weight_matrix(self) -> np.ndarray:
        """
        (K, K) float. Agreement weight between every pair of categories.

        1.0 means indistinguishable, 0.0 means maximally different. The
        diagonal is 1.0 and the matrix is symmetric: how wrong a disagreement
        is does not depend on which rater said which.
        """


class OrderedCategories(DistanceCategories):
    """
    I am categories on a ranked scale, and I derive distance from position.

    'mild' < 'moderate' < 'severe'. My order is not merely an encoding, it is a
    claim: adjacent categories are more alike than distant ones, so a rater
    who says 'mild' where another said 'moderate' has disagreed less than one
    who said 'severe'.

    I compute weights from rank distance, either linearly or quadratically.
    Quadratic is the common default in the medical literature and forgives
    near-misses more readily; linear penalises distance evenly. Neither is more
    correct, so I make you say which -- the choice changes the coefficient and
    ought to appear in a methods section.

    Attributes
    ----------
    weighting:
        'linear' or 'quadratic'. How rank distance becomes disagreement.
    """

    _KIND = "ordered"
    _ERROR_PREFIX = "[OrderedCategories] Error"

    names:     tuple[str, ...]
    labels:    dict[str, str]
    weighting: str

    _ALLOWED = ("linear", "quadratic")

    def __init__(
        self,
        names:     object,
        weighting: str = "quadratic",
        labels:    dict[str, str] | None = None,
    ) -> None:
        self.weighting = weighting
        super().__init__(names, labels)

    def _validate(self) -> None:
        try:
            super()._validate()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{self._ERROR_PREFIX}: {exc}") from exc

        if self.weighting not in self._ALLOWED:
            raise ValueError(
                f"{type(self).__name__}.weighting must be one of "
                f"{list(self._ALLOWED)}, got {self.weighting!r}. Linear "
                f"penalises rank distance evenly; quadratic forgives "
                f"near-misses. The choice changes the coefficient."
            )

    def to_dict(self) -> dict:
        spec = super().to_dict()
        spec["weighting"] = self.weighting
        return spec

    def _eq_extras(self, other: "OrderedCategories") -> bool:
        return self.weighting == other.weighting

    def weight_matrix(self) -> np.ndarray:
        k        = self.n_categories
        rank     = np.arange(k)
        distance = np.abs(rank[:, None] - rank[None, :]).astype(np.float64)
        if self.weighting == "quadratic":
            distance = distance ** 2
        return 1.0 - distance / distance.max()


class WeightedCategories(DistanceCategories):
    """
    I am categories with distances supplied outright, not inferred.

    I do not require a ranking. 'cat' and 'dog' may be closer to each other
    than either is to 'car' without any of the three being greater or less than
    another, and no ordered scale can express that.

    Use me when the similarity you mean cannot be recovered from an order:
    grouped diagnoses, a coding frame where two labels are near-synonyms, or a
    scale whose steps are genuinely uneven. Otherwise OrderedCategories is
    fewer numbers to get wrong.

    Attributes
    ----------
    weights:
        (K, K) float. Symmetric, 1.0 on the diagonal, entries in [0, 1].
    """

    _KIND = "weighted"
    _ERROR_PREFIX = "[WeightedCategories] Error"

    names:   tuple[str, ...]
    labels:  dict[str, str]
    weights: np.ndarray

    def __init__(
        self,
        names:   object,
        weights: np.ndarray,
        labels:  dict[str, str] | None = None,
    ) -> None:
        self.weights = np.asarray(weights, dtype=np.float64)
        super().__init__(names, labels)

    def _validate(self) -> None:
        try:
            super()._validate()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{self._ERROR_PREFIX}: {exc}") from exc
        
        name = type(self).__name__
        k    = self.n_categories

        if self.weights.shape != (k, k):
            raise ValueError(
                f"{name}.weights must be ({k}, {k}) to match its {k} "
                f"categories, got {self.weights.shape}."
            )
        if not np.allclose(self.weights, self.weights.T):
            raise ValueError(
                f"{name}.weights must be symmetric. How wrong a disagreement "
                f"is cannot depend on which rater said which."
            )
        if not np.allclose(np.diag(self.weights), 1.0):
            raise ValueError(
                f"{name}.weights must have 1.0 on the diagonal: a category "
                f"agrees perfectly with itself. Got "
                f"{np.diag(self.weights).round(3).tolist()}."
            )
        if self.weights.min() < 0.0 or self.weights.max() > 1.0:
            raise ValueError(
                f"{name}.weights must lie in [0, 1], where 1.0 is "
                f"indistinguishable and 0.0 maximally different. Got range "
                f"[{self.weights.min():.3f}, {self.weights.max():.3f}]. If you "
                f"have distances rather than similarities, subtract them from "
                f"1 after scaling to [0, 1]."
            )

    def to_dict(self) -> dict:
        """
        Includes the whole matrix, as nested lists.

        Verbose -- K rows of K floats -- but there is no shorthand to emit. A
        matrix that could be derived from a rule would have been an
        OrderedCategories.
        """
        spec = super().to_dict()
        spec["weights"] = [[float(x) for x in row] for row in self.weights]
        return spec

    def _eq_extras(self, other: "WeightedCategories") -> bool:
        return bool(np.array_equal(self.weights, other.weights))

    def weight_matrix(self) -> np.ndarray:
        return self.weights


# Names -----------------------------------------------------------------------  #

def _as_name(value: object, owner: str) -> str:
    """
    A category name as text.

    Numbers become their text, so YAML's unquoted ``1`` and ``-1`` mean the
    ``'1'`` and ``'-1'`` a ratings file holds. Booleans are refused: YAML turns
    unquoted ``yes``, ``no``, ``on`` and ``off`` into True and False, and the
    original word cannot be recovered from that. Anything else is refused too.

    I raise ValueError, not TypeError, because ``build_categories`` reads a
    TypeError from a constructor as a misspelt key.
    """
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(
            f"{owner}: category name {value!r} is a boolean. YAML reads "
            f"unquoted yes/no/on/off as True/False; quote the names, e.g. "
            f"['yes', 'no']."
        )
    if isinstance(value, str):
        return value
    if isinstance(value, (int, float, np.integer, np.floating)):
        return str(value)
    raise ValueError(
        f"{owner}: category name {value!r} must be text or a number, got "
        f"{type(value).__name__}."
    )


# Factory ---------------------------------------------------------------------  #

_KINDS = {
    cls._KIND: cls
    for cls in (UnorderedCategories, OrderedCategories, WeightedCategories)
}
"""Built from each class's own token, so ``kind`` has one definition."""


def build_categories(
    spec: dict, base_dir: str | Path | None = None
) -> Categories:
    """
    Build the right Categories from a plain mapping, as parsed from YAML.

    ``kind`` selects the class; the remaining keys are its arguments. Kept here
    rather than with the config that calls it, because this is the module that
    knows what subclasses exist.

    A ``weighted`` scheme may give its matrix inline as ``weights`` or in a
    separate YAML file as ``weights_file``, never both. A relative
    ``weights_file`` is read from ``base_dir``, which ``RatingsSemantics``
    sets to the folder of the semantics file; without one, from the current
    working directory.

    Examples
    --------
    >>> build_categories({"kind": "unordered", "names": ["yes", "no"]})
    UnorderedCategories(
        n_categories = 2
        names        = ['yes', 'no']
    )
    """
    spec = dict(spec)
    kind = spec.pop("kind", None)
    if kind not in _KINDS:
        raise ValueError(
            f"categories.kind must be one of {sorted(_KINDS)}, got {kind!r}. "
            f"Use 'unordered' for a nominal scheme, 'ordered' for a ranked "
            f"scale, or 'weighted' to supply a similarity matrix directly."
        )
    if "names" not in spec:
        raise KeyError(
            f"categories of kind {kind!r} needs 'names': the responses a rater "
            f"was allowed to give, in encoding order."
        )
    if "weights_file" in spec:
        spec["weights"] = _load_weights_file(
            spec.pop("weights_file"), spec, kind, base_dir
        )
    try:
        return _KINDS[kind](**spec)
    except TypeError as error:
        raise TypeError(
            f"categories of kind {kind!r} got unexpected keys "
            f"{sorted(set(spec) - {'names'})}. {error}"
        ) from None


def _load_weights_file(
    file: str | Path, spec: dict, kind: str, base_dir: str | Path | None
) -> list:
    """
    Read a weighted scheme's matrix from its own YAML file.

    The file holds exactly ``names`` and ``matrix``. The names are there so
    that the matrix cannot be read against the wrong categories: a bare matrix
    depends on row order, and reordering the categories in the semantics would
    otherwise move every weight onto the wrong pair without a word.

    I return the matrix as parsed. ``WeightedCategories`` then validates it
    exactly as it would an inline one.
    """
    import yaml

    if kind != "weighted":
        raise ValueError(
            f"categories.weights_file is only allowed with kind 'weighted', "
            f"got kind {kind!r}. An ordered scale derives its weights from "
            f"'weighting'; an unordered one has none."
        )
    if "weights" in spec:
        raise ValueError(
            "categories gives both 'weights' and 'weights_file'. Give the "
            "matrix in one place only, so there is one answer to what it is."
        )

    path = Path(file)
    if not path.is_absolute():
        path = Path(base_dir if base_dir is not None else ".") / path
    if path.suffix.lower() not in {".yaml", ".yml"}:
        raise ValueError(
            f"categories.weights_file must be a YAML file (.yaml or .yml), "
            f"got {path.name!r}."
        )

    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
    except FileNotFoundError as error:
        raise FileNotFoundError(
            f"categories.weights_file not found: {path}. A relative path is "
            f"read from the folder of the semantics file."
        ) from error

    if not isinstance(data, dict) or set(data) != {"names", "matrix"}:
        found = sorted(data) if isinstance(data, dict) else type(data).__name__
        raise ValueError(
            f"categories.weights_file {path} must hold exactly 'names' and "
            f"'matrix', got {found}."
        )

    file_names     = [str(n) for n in data["names"]]
    declared_names = [str(n) for n in spec["names"]]
    if file_names != declared_names:
        raise ValueError(
            f"categories.weights_file {path} is written for categories "
            f"{file_names}, but the semantics declares {declared_names}. "
            f"They must match exactly, in the same order, because row and "
            f"column i of the matrix belong to category i."
        )
    return data["matrix"]
