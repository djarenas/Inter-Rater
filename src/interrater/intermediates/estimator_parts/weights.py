"""Weights for ordinal data: how much each disagreement costs."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_WEIGHTS_ERROR = "[Weights] Error"

_KINDS = ("identity", "linear", "quadratic", "custom")


@dataclass(frozen=True, eq=False)
class Weights:
    """
    I say how much credit a disagreement still earns.

    A numpy array of size K by K. Entry (i, j) is the credit given when one
    rater says category i and another says category j: 1.0 for the same
    category, 0.0 for maximally different, something between for a
    disagreement that is only partly wrong.

    Unweighted agreement is not a separate idea -- it is me with an identity
    matrix, where every disagreement costs everything. That is why nothing in
    this package has a weighted and an unweighted code path: there is one path,
    and identity is the default position on it.

    I only mean anything when the categories are ordered or otherwise have a
    notion of distance. 'cat', 'dog' and 'car' can be given weights saying the
    first two are closer, but a nominal scale with no such claim should use
    identity and say nothing it cannot defend.

    I remember which kind I am. A report that says a coefficient was weighted
    should be able to say *how*, and 'quadratic' is a great deal more useful to
    a reader than a printed matrix.

    Attributes
    ----------
    matrix:
        A numpy array of size K by K, holding floats in [0, 1], symmetric, with
        1.0 on the diagonal.
    kind:
        How the matrix was made: identity, linear, quadratic, or custom.

    Examples
    --------
    >>> w = Weights.quadratic(3)
    >>> w.matrix[0].round(3).tolist()
    [1.0, 0.75, 0.0]
    >>> Weights.identity(3).is_identity
    True
    """

    # Attributes ---------------------------------------------------------  #

    matrix: np.ndarray
    kind:   str = "custom"

    # Constructor --------------------------------------------------------  #

    def __post_init__(self) -> None:
        object.__setattr__(self, "matrix", np.asarray(self.matrix, dtype=np.float64))
        self._validate()

    def _validate(self) -> None:
        matrix = self.matrix

        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ValueError(
                f"{_WEIGHTS_ERROR}: matrix must be square -- one row and one "
                f"column per category -- got shape {matrix.shape}."
            )

        if matrix.shape[0] < 2:
            raise ValueError(
                f"{_WEIGHTS_ERROR}: needs at least 2 categories, got "
                f"{matrix.shape[0]}."
            )

        if not np.allclose(matrix, matrix.T):
            raise ValueError(
                f"{_WEIGHTS_ERROR}: matrix must be symmetric. How much a "
                f"disagreement costs cannot depend on which rater is named "
                f"first."
            )

        if not np.allclose(np.diag(matrix), 1.0):
            raise ValueError(
                f"{_WEIGHTS_ERROR}: matrix must have 1.0 on the diagonal -- a "
                f"category agrees perfectly with itself. Got "
                f"{np.diag(matrix).round(3).tolist()}."
            )

        if matrix.min() < 0.0 or matrix.max() > 1.0:
            raise ValueError(
                f"{_WEIGHTS_ERROR}: matrix must lie in [0, 1], where 1.0 is "
                f"full credit and 0.0 is none. Got range "
                f"[{matrix.min():.3f}, {matrix.max():.3f}]. If you have "
                f"distances rather than credits, subtract them from 1 after "
                f"scaling."
            )

        if self.kind not in _KINDS:
            raise ValueError(
                f"{_WEIGHTS_ERROR}: kind must be one of {list(_KINDS)}, got "
                f"{self.kind!r}."
            )

    # Class methods ------------------------------------------------------  #

    @classmethod
    def identity(cls, n_categories: int) -> Weights:
        """
        Every disagreement costs everything. The unweighted case.
        """
        return cls(np.eye(n_categories), kind="identity")

    @classmethod
    def linear(cls, n_categories: int) -> Weights:
        """
        Credit falls in proportion to how far apart two categories sit.

        For an ordered scale only: position is taken to be rank, so categories
        one step apart are treated as more alike than categories three steps
        apart. Penalises distance evenly, where quadratic forgives near misses.
        Neither is more correct and the choice changes the coefficient, which
        is why it has to be made rather than defaulted.
        """
        rank     = np.arange(n_categories)
        distance = np.abs(rank[:, None] - rank[None, :]).astype(np.float64)
        return cls(1.0 - distance / distance.max(), kind="linear")

    @classmethod
    def quadratic(cls, n_categories: int) -> Weights:
        """
        Credit falls with the square of the distance between categories.

        The usual choice in the medical literature. Forgives adjacent
        disagreements much more readily than linear weighting does -- on a
        four-point scale, neighbours keep 0.89 of their credit against 0.67.
        """
        rank     = np.arange(n_categories)
        distance = np.abs(rank[:, None] - rank[None, :]).astype(np.float64) ** 2
        return cls(1.0 - distance / distance.max(), kind="quadratic")

    @classmethod
    def from_type(cls, weights_type: str, n_categories: int) -> Weights | None:
        """
        Factory method to build Weights from a type string.

        Parameters
        ----------
        weights_type : str
            The type of weights: "identity", "linear", "quadratic", or "custom".
        n_categories : int
            Number of categories (not used for "custom").

        Returns
        -------
        Weights | None
            The constructed Weights object, or None if weights_type is "identity"
            (since identity weights are the default and don't need an object).

        Raises
        ------
        ValueError
            If weights_type is "custom" (use from_file or pass matrix directly).
            If weights_type is unknown.
        """
        if weights_type == "identity":
            return None  # Identity is the default; no Weights needed
        elif weights_type == "linear":
            return cls.linear(n_categories)
        elif weights_type == "quadratic":
            return cls.quadratic(n_categories)
        elif weights_type == "custom":
            raise ValueError(
                f"{_WEIGHTS_ERROR}: from_type does not support 'custom'. "
                "Use from_file() or pass a matrix directly to Weights()."
            )
        else:
            raise ValueError(
                f"{_WEIGHTS_ERROR}: unknown weights_type '{weights_type}'. "
                f"Expected one of: 'identity', 'linear', 'quadratic', 'custom'."
            )

    @classmethod
    def from_file(cls, path: str) -> Weights:
        """
        Load custom weights from a file.

        Supports YAML, JSON, CSV, and plain text formats. Auto-detects format
        from file extension. Loaded weights default to kind="custom".

        Parameters
        ----------
        path : str
            Path to weights file. Extension must be one of:
            .yaml, .yml (YAML), .json (JSON), .csv (CSV), .txt (plain text)

        Returns
        -------
        Weights
            Loaded and validated weights matrix with kind="custom".

        Raises
        ------
        FileNotFoundError
            If file does not exist.
        ValueError
            If file extension not recognized, or matrix is invalid.
        IOError
            If file cannot be read.

        Examples
        --------
        Load from YAML:
        >>> w = Weights.from_file("weights.yaml")

        Load from CSV:
        >>> w = Weights.from_file("weights.csv")
        """
        import pathlib

        # Get file extension
        path_obj = pathlib.Path(path)
        if not path_obj.exists():
            raise FileNotFoundError(f"{_WEIGHTS_ERROR}: file not found: {path}")

        ext = path_obj.suffix.lower()

        try:
            if ext in {".yaml", ".yml"}:
                matrix = cls._load_yaml(path)
            elif ext == ".json":
                matrix = cls._load_json(path)
            elif ext == ".csv":
                matrix = cls._load_csv(path)
            elif ext == ".txt":
                matrix = cls._load_text(path)
            else:
                raise ValueError(
                    f"{_WEIGHTS_ERROR}: unknown file extension '{ext}'. "
                    f"Expected one of: .yaml, .yml, .json, .csv, .txt"
                )
        except (FileNotFoundError, ValueError):
            raise
        except Exception as e:
            raise ValueError(
                f"{_WEIGHTS_ERROR}: failed to load weights from {path}: {e}"
            ) from e

        # Build and return Weights (validation happens in __post_init__)
        return cls(matrix, kind="custom")

    @staticmethod
    def _load_yaml(path: str) -> np.ndarray:
        """Load matrix from YAML file."""
        import yaml

        with open(path) as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError(
                f"YAML must contain a dict with 'matrix' key, got {type(data).__name__}"
            )

        if "matrix" not in data:
            raise ValueError("YAML must contain 'matrix' key")

        matrix = np.array(data["matrix"], dtype=np.float64)
        return matrix

    @staticmethod
    def _load_json(path: str) -> np.ndarray:
        """Load matrix from JSON file."""
        import json

        with open(path) as f:
            data = json.load(f)

        if not isinstance(data, dict):
            raise ValueError(
                f"JSON must contain an object with 'matrix' key, got {type(data).__name__}"
            )

        if "matrix" not in data:
            raise ValueError("JSON must contain 'matrix' key")

        matrix = np.array(data["matrix"], dtype=np.float64)
        return matrix

    @staticmethod
    def _load_csv(path: str) -> np.ndarray:
        """Load matrix from CSV file (space, tab, or comma separated)."""
        matrix = np.loadtxt(path, delimiter=None, dtype=np.float64)
        if matrix.ndim != 2:
            raise ValueError(
                f"CSV must be 2-dimensional (K×K matrix), got shape {matrix.shape}"
            )
        return matrix

    @staticmethod
    def _load_text(path: str) -> np.ndarray:
        """Load matrix from plain text file (whitespace separated)."""
        matrix = np.loadtxt(path, dtype=np.float64)
        if matrix.ndim != 2:
            raise ValueError(
                f"Text file must be 2-dimensional (K×K matrix), got shape {matrix.shape}"
            )
        return matrix

    # Public API ---------------------------------------------------------  #

    @property
    def n_categories(self) -> int:
        return int(self.matrix.shape[0])

    @property
    def is_identity(self) -> bool:
        """Whether I would leave a coefficient unweighted."""
        return bool(np.allclose(self.matrix, np.eye(self.n_categories)))

    # Comparison and display ---------------------------------------------  #

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Weights):
            return NotImplemented
        return np.array_equal(self.matrix, other.matrix)

    __hash__ = None  # type: ignore[assignment]
