"""Configuration for pairwise permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ...contracts import PAIRWISE_PERMUTATION_COEFFICIENTS
from .permutation_config import PermutationConfig

if TYPE_CHECKING:
    from .analyzer_config import AnalyzerConfig
    from .permutation_bootstrap_mean_config import (
        PermutationBootstrapMeanConfig,
    )

_ERROR_PREFIX = "[PairwisePermutationConfig] Error"


@dataclass(frozen=True)
class PairwisePermutationConfig(PermutationConfig):
    """
    I am the configuration for pairwise permutation analysis.

    I configure a coefficient computed for each of the C(n_raters, 2) pairs
    of raters. This shows which pairs agree less than others, which can point
    to raters who are not interchangeable. A rater whose pairs disagree more
    is not necessarily wrong: they may have more training or more context.
    Low agreement suggests the raters are not interchangeable, not which of
    them is right.

    I enforce that:
    - coefficient_to_compute is in PAIRWISE_PERMUTATION_COEFFICIENTS
    - config_of_coefficient is the correct config type for that coefficient
    - subset_k must be None (pairwise is fixed at k=2)
    """
    # Attribute declaration

    # Inherited from PermutationConfig
    coefficient_to_compute: str
    # Checked against coefficient_to_compute
    config_of_coefficient: AnalyzerConfig
    # Must be None: pairwise is k = 2 (a leftover from a subset config is
    # refused)
    subset_k: int | None
    bootstrap_mean_config: PermutationBootstrapMeanConfig | None
    balanced_data_only: bool
    keep_draws: bool
    max_evaluations: int

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate pairwise permutation configuration."""

        # Call parent validation first
        try:
            super().__post_init__()
        except (TypeError, ValueError) as e:
            raise type(e)(
                f"{_ERROR_PREFIX}: parent validation failed: {e}"
            ) from e

        # Reject subset_k for pairwise (pairwise is always k=2, implicit)
        if self.subset_k is not None:
            raise ValueError(
                f"{_ERROR_PREFIX}: subset_k must be None for pairwise "
                f"permutation (pairwise is always k=2), "
                f"got subset_k={self.subset_k}"
            )

        # Validate that coefficient is pairwise-capable. Every coefficient
        # is today, so this only guards a future one that cannot run on pairs.
        coefficient = self.coefficient_to_compute
        if coefficient not in PAIRWISE_PERMUTATION_COEFFICIENTS:
            raise ValueError(
                f"{_ERROR_PREFIX}: '{self.coefficient_to_compute}' does "
                f"not support pairwise permutation. Allowed coefficients: "
                f"{sorted(PAIRWISE_PERMUTATION_COEFFICIENTS)}"
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
    def build_from_dict(cls, data_dict: dict) -> PairwisePermutationConfig:
        """
        Build a pairwise permutation config from a plain dictionary.

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
        PairwisePermutationConfig
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
