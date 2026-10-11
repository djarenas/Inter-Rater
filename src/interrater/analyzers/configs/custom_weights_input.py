"""Configuration for custom weight matrices."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ...intermediates.estimator_parts import Weights

_ERROR_PREFIX = "[CustomWeightsInput] Error"


@dataclass(frozen=True)
class CustomWeightsInput:
    """
    I describe how to load or construct custom weights.

    I am NOT a Weights object; I am configuration that tells the analyzer
    how to build one. Provide either a matrix directly (inline) or a file_path
    to load from. Exactly one must be set.

    Attributes
    ----------
    kind : str, default "custom"
        The type identifier. Always "custom" for this class.
    matrix : np.ndarray | None, default None
        The K×K weight matrix, inline. Values in [0, 1], symmetric, 1.0 diagonal.
        If provided, file_path must be None.
    file_path : str | None, default None
        Path to a file containing the weight matrix (YAML, CSV, JSON, or plain text).
        If provided, matrix must be None. File will be loaded via Weights.from_file().
    """

    kind: str = "custom"
    matrix: np.ndarray | None = None
    file_path: str | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate that exactly one of matrix or file_path is set."""

        name = "kind"
        if not isinstance(getattr(self, name), str):
            raise TypeError(
                f"{_ERROR_PREFIX}: {name} must be a string, "
                f"got {type(getattr(self, name)).__name__}"
            )

        # Check that kind is "custom"
        if self.kind != "custom":
            raise ValueError(
                f"{_ERROR_PREFIX}: kind must be 'custom', got '{self.kind}'"
            )

        # Exactly one of matrix or file_path must be set
        has_matrix = self.matrix is not None
        has_path = self.file_path is not None

        if not (has_matrix ^ has_path):  # XOR: exactly one true
            raise ValueError(
                f"{_ERROR_PREFIX}: exactly one of matrix or file_path must be set. "
                f"Got matrix={has_matrix}, file_path={has_path}"
            )

        # If matrix provided, validate it's a numpy array
        if has_matrix and not isinstance(self.matrix, np.ndarray):
            raise TypeError(
                f"{_ERROR_PREFIX}: matrix must be np.ndarray, "
                f"got {type(self.matrix).__name__}"
            )

        # If file_path provided, validate it's a string
        if has_path and not isinstance(self.file_path, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: file_path must be str, "
                f"got {type(self.file_path).__name__}"
            )

    # Public API ------------------------------------------------------------- #

    def build_weights(self) -> Weights:
        """
        Build and return a Weights object from this input.

        Returns
        -------
        Weights
            The constructed weights matrix with validation applied.

        Raises
        ------
        ValueError
            If matrix or file cannot be loaded/validated.
        """
        from ...intermediates.estimator_parts import Weights

        if self.matrix is not None:
            # Build from inline matrix
            return Weights(self.matrix, kind="custom")
        else:
            # Load from file
            if self.file_path is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: file_path is None but matrix also None. "
                    f"This should not happen after __post_init__ validation."
                )
            return Weights.from_file(self.file_path)

    def summary(self) -> dict:
        """Return configuration for audit trail."""
        result = {
            "kind": self.kind,
            "file_path": self.file_path,
        }
        # Don't include matrix in summary (large, not needed in audit trail)
        if self.matrix is not None:
            result["matrix_shape"] = tuple(self.matrix.shape)
        return result
