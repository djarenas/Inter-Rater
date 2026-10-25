"""Configuration for subset permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass

from .permutation_config import PermutationConfig
from ...contracts import SUBSET_PERMUTATION_COEFFICIENTS

_ERROR_PREFIX = "[SubsetPermutationConfig] Error"


@dataclass(frozen=True)
class SubsetPermutationConfig(PermutationConfig):
    """
    I am the configuration for subset permutation analysis.

    I compute a coefficient for all C(n_raters, k) combinations of k raters.
    This is useful for understanding how agreement scales with panel size or
    for analyzing specific subgroups of raters.

    subset_k controls which subset sizes to evaluate:
    - A single int (e.g. 3) runs one subset size.
    - A list of ints (e.g. [2, 3]) runs multiple subset sizes.
    - The string "all" runs every k from 2 through n_raters (expanded
      at analysis time when the data is available).

    I enforce that:
    - coefficient_to_compute is in SUBSET_PERMUTATION_COEFFICIENTS
    - config_of_coefficient is the correct config type for that coefficient
    - subset_k is provided and structurally valid (>= 2 for every k)
    - if bootstrap_mean_config uses "independent", per-subset bootstrap
      must be enabled (n_resamples >= 2)

    The analysis class is responsible for:
    - Normalizing int to [int]
    - Expanding "all" to list(range(2, n_raters + 1))
    - Validating each k <= n_raters

    Attributes
    ----------
    coefficient_to_compute : str
        Must be one of: "fleiss_kappa", "scotts_pi", "gwet_ac1",
        "brennan_prediger", "krippendorff_alpha"
    config_of_coefficient : AnalyzerConfig
        The nested config matching the coefficient
        (FleissKappaConfig, ScottsConfig, GwetAC1Config, etc.)
    subset_k : int | list[int] | str | None, default None
        Must be provided (not None). A single int, a list of ints
        (each >= 2), or "all". The analysis class expands "all" and
        validates each k against n_raters.
    """

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate subset permutation configuration."""

        # Call parent validation first
        try:
            super().__post_init__()
        except (TypeError, ValueError) as e:
            raise type(e)(
                f"{_ERROR_PREFIX}: parent validation failed: {e}"
            ) from e

        # Require subset_k for subset permutation
        if self.subset_k is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: subset_k is required for subset permutation. "
                f"Specify the number of raters per subset (e.g., subset_k=3, "
                f"subset_k=[2, 3], or subset_k='all')"
            )

        # Validate subset_k by type
        self._validate_subset_k()

        # Validate that coefficient is subset-capable
        if self.coefficient_to_compute not in SUBSET_PERMUTATION_COEFFICIENTS:
            raise ValueError(
                f"{_ERROR_PREFIX}: '{self.coefficient_to_compute}' does not support "
                f"subset permutation. Allowed coefficients: "
                f"{sorted(SUBSET_PERMUTATION_COEFFICIENTS)}"
            )

        # Validate that config_of_coefficient matches the coefficient
        # Map coefficient names to expected config class names
        expected_config_map = {
            "fleiss_kappa": "FleissKappaConfig",
            "scotts_pi": "ScottsConfig",
            "gwet_ac1": "GwetAC1Config",
            "brennan_prediger": "BrennanPredigerConfig",
            "krippendorff_alpha": "KrippendorffAlphaConfig",
        }

        expected_config_class = expected_config_map.get(
            self.coefficient_to_compute
        )
        if expected_config_class is None:
            # This shouldn't happen if contracts are correct, but guard it
            raise ValueError(
                f"{_ERROR_PREFIX}: Unknown coefficient '{self.coefficient_to_compute}'"
            )

        actual_config_class = type(self.config_of_coefficient).__name__

        if actual_config_class != expected_config_class:
            raise TypeError(
                f"{_ERROR_PREFIX}: coefficient_to_compute='{self.coefficient_to_compute}' "
                f"requires config_of_coefficient to be {expected_config_class}, "
                f"got {actual_config_class}"
            )

        # Cross-validate: "independent" mean bootstrap requires per-subset
        # bootstrap enabled in config_of_coefficient.
        # The "independent" method assembles the mean from independently-
        # drawn per-subset bootstrap values, so each subset must be
        # bootstrapped. The "correlated" method does its own subject
        # resampling and does not depend on per-subset bootstrap.
        if (
            self.bootstrap_mean_config is not None
            and self.bootstrap_mean_config.method == "independent"
            and self.config_of_coefficient.n_resamples < 2
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: bootstrap_mean_config method "
                f"'independent' requires per-subset bootstrap enabled "
                f"in config_of_coefficient (n_resamples >= 2), "
                f"but n_resamples={self.config_of_coefficient.n_resamples}"
            )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_dict(cls, data_dict: dict) -> SubsetPermutationConfig:
        """
        Build a subset permutation config from a plain dictionary.

        I do everything the parent's ``_prepare_dict`` does (validate
        inputs, normalize YAML aliases, hydrate bootstrap_mean_config),
        then I hydrate ``config_of_coefficient`` into the correct
        AnalyzerConfig subclass based on ``coefficient_to_compute``.

        If ``coefficient_to_compute`` is not in the subset map, I skip
        hydration and let ``__post_init__`` reject it with its own clear
        error message.

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

        # Map subset coefficients to their config classes
        from .fleiss_kappa_config import FleissKappaConfig
        from .scotts_config import ScottsConfig
        from .gwet_ac1_config import GwetAC1Config
        from .brennan_prediger_config import BrennanPredigerConfig
        from .krippendorff_alpha_config import KrippendorffAlphaConfig

        config_class_map = {
            "fleiss_kappa": FleissKappaConfig,
            "scotts_pi": ScottsConfig,
            "gwet_ac1": GwetAC1Config,
            "brennan_prediger": BrennanPredigerConfig,
            "krippendorff_alpha": KrippendorffAlphaConfig,
        }

        # Hydrate config_of_coefficient if it is still a raw dict
        coefficient = data_dict.get("coefficient_to_compute")
        config_class = config_class_map.get(coefficient)

        if config_class is not None:
            coc = data_dict.get("config_of_coefficient")
            if isinstance(coc, dict):
                data_dict["config_of_coefficient"] = (
                    config_class.build_from_dict(coc)
                )

        return cls(**data_dict)

    def summary(self) -> dict:
        """Return configuration (version and nested config from parent)."""
        return super().summary()

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
                if isinstance(val, bool):
                    raise TypeError(
                        f"{_ERROR_PREFIX}: subset_k list elements must be "
                        f"int, got {type(val).__name__} at index {i}"
                    )
                if not isinstance(val, int):
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
