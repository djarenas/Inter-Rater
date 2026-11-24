"""Configuration for subset permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass

from .permutation_config import PermutationConfig
from ...contracts import SUBSET_PERMUTATION_COEFFICIENTS

_ERROR_PREFIX = "[SubsetPermutationConfig] Error"

# Refused by name: listed among the subset coefficients until 2.1.0
_SCOTTS_PI_NOT_FOR_SUBSETS = (
    "Scott's pi is defined for exactly two raters, so it has no subsets: "
    "with more raters it is Fleiss' kappa under another name. Use "
    "coefficient_to_compute: fleiss_kappa for subsets, or pairwise "
    "permutation for Scott's pi on every pair."
)


@dataclass(frozen=True)
class SubsetPermutationConfig(PermutationConfig):
    """
    I am the configuration for subset permutation analysis.

    I configure a coefficient computed for all C(n_raters, k) combinations of
    k raters.
    This is useful for understanding how agreement scales with panel size or
    for analyzing specific subgroups of raters.

    subset_k controls which subset sizes to evaluate:
    - Providing a single integer (e.g. 3) tells the analyzer to run every
      subset of that size.
    - Providing a list of integers (e.g. [2, 3]) signals to run multiple
      subset sizes.
    - Providing the string "all" tells the analyzer to run every k from 2
      through n_raters (which exact subsets are run is decided by the
      analyzer at run time, when it has the data in front of it).

    I enforce that:
    - coefficient_to_compute is in SUBSET_PERMUTATION_COEFFICIENTS
    - config_of_coefficient is the correct config type for that coefficient
    - subset_k is provided and structurally valid (>= 2 for every k)

    The analysis class is responsible for:
    - Normalizing k: checking whether it is an int or [int], with [int] as
      the final form
    - Expanding "all" to list(range(2, n_raters + 1))
    - Validating each k <= n_raters

    Attributes
    ----------
    coefficient_to_compute : str
        Must be one of: "fleiss_kappa", "gwet_ac1", "brennan_prediger",
        "krippendorff_alpha". Not "scotts_pi": it is two-rater only.
    config_of_coefficient : AnalyzerConfig
        The nested config matching the coefficient
        (FleissKappaConfig, GwetAC1Config, etc.)
    subset_k : int | list[int] | str | None, default None
        Must be provided (not None). A single int, a list of ints
        (each >= 2), or "all". The analysis class expands "all" and
        validates each k against n_raters.
    """

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate subset permutation configuration."""

        # Scott's pi is refused by name, before the parent's generic message
        if self.coefficient_to_compute == "scotts_pi":
            raise ValueError(f"{_ERROR_PREFIX}: {_SCOTTS_PI_NOT_FOR_SUBSETS}")

        # Parent validation
        try:
            super().__post_init__()
        except (TypeError, ValueError) as e:
            raise type(e)(
                f"{_ERROR_PREFIX}: parent validation failed: {e}"
            ) from e

        # Require subset_k for subset permutation
        if self.subset_k is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: subset_k is required for subset "
                f"permutation. Specify the number of raters per subset "
                f"(e.g., subset_k=3, subset_k=[2, 3], or subset_k='all')"
            )

        # Validate subset_k by type
        self._validate_subset_k()

        # Validate that coefficient is subset-capable
        if self.coefficient_to_compute not in SUBSET_PERMUTATION_COEFFICIENTS:
            raise ValueError(
                f"{_ERROR_PREFIX}: '{self.coefficient_to_compute}' does not "
                f"support subset permutation. Allowed coefficients: "
                f"{sorted(SUBSET_PERMUTATION_COEFFICIENTS)}"
            )

        # Validate that config_of_coefficient matches the coefficient
        expected = self._resolve_config_class(self.coefficient_to_compute)
        if not isinstance(self.config_of_coefficient, expected):
            raise TypeError(
                f"{_ERROR_PREFIX}: coefficient_to_compute="
                f"'{self.coefficient_to_compute}' requires "
                f"config_of_coefficient to be {expected.__name__}, "
                f"got {type(self.config_of_coefficient).__name__}"
            )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_dict(cls, data_dict: dict) -> SubsetPermutationConfig:
        """
        Build a subset permutation config from a plain dictionary.

        Does everything the parent's ``_prepare_dict`` does (validate
        inputs, normalize YAML aliases, hydrate bootstrap_mean_config),
        then hydrates ``config_of_coefficient`` into the correct
        AnalyzerConfig subclass based on ``coefficient_to_compute``.

        If ``coefficient_to_compute`` is unknown, hydration is skipped and
        ``__post_init__`` rejects it with its own clear error message.

        Parameters
        ----------
        data_dict : dict
            Configuration as plain Python types. Mutated in place.

        Returns
        -------
        SubsetPermutationConfig
            Instantiated config with all nested objects hydrated
            and validation applied.

        Raises
        ------
        TypeError
            If data_dict is not a dict.
        ValueError
            If data_dict is None or validation fails.
        """
        cls._prepare_dict(data_dict)

        # Hydrate config_of_coefficient if it is still a raw dict
        coefficient = data_dict.get("coefficient_to_compute")
        config_class = cls._resolve_config_class(coefficient)

        if config_class is not None:
            coc = data_dict.get("config_of_coefficient")
            if isinstance(coc, dict):
                data_dict["config_of_coefficient"] = (
                    config_class.build_from_dict(coc)
                )

        return cls(**data_dict)

    # Private helpers -------------------------------------------------------- #

    def _validate_subset_k(self) -> None:
        """
        Validate subset_k type and values.

        Accepts int (>= 2), list[int] (non-empty, all >= 2, no duplicates),
        or the string "all". Anything else raises TypeError or ValueError.

        Upper-bound validation (k <= n_raters) and expansion of "all" happen
        at analysis time when the data is available.
        """
        k = self.subset_k

        if isinstance(k, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: subset_k must be int, list[int], or "
                f"'all', got {type(k).__name__}"
            )

        if isinstance(k, int):
            if k < 2:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subset_k must be >= 2, got {k}"
                )

        elif isinstance(k, list):
            if len(k) == 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: subset_k list must not be empty"
                )
            for i, val in enumerate(k):
                if isinstance(val, bool) or not isinstance(val, int):
                    raise TypeError(
                        f"{_ERROR_PREFIX}: subset_k list elements must be "
                        f"int, got {type(val).__name__} at index {i}"
                    )
                if val < 2:
                    raise ValueError(
                        f"{_ERROR_PREFIX}: subset_k values must be >= 2, "
                        f"got {val} at index {i}"
                    )
            if len(k) != len(set(k)):
                raise ValueError(
                    f"{_ERROR_PREFIX}: subset_k list must not contain "
                    f"duplicates, got {k}"
                )

        elif isinstance(k, str):
            if k != "all":
                raise ValueError(
                    f"{_ERROR_PREFIX}: subset_k string must be 'all', "
                    f"got '{k}'"
                )

        else:
            raise TypeError(
                f"{_ERROR_PREFIX}: subset_k must be int, list[int], or "
                f"'all', got {type(k).__name__}"
            )
