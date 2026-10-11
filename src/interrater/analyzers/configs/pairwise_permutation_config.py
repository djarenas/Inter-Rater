"""Configuration for pairwise permutation analysis."""

from __future__ import annotations

from dataclasses import dataclass

from ...contracts import PAIRWISE_PERMUTATION_COEFFICIENTS
from .permutation_config import PermutationConfig
from .analyzer_config import AnalyzerConfig
from .permutation_bootstrap_mean_config import PermutationBootstrapMeanConfig

_ERROR_PREFIX = "[PairwisePermutationConfig] Error"


@dataclass(frozen=True)
class PairwisePermutationConfig(PermutationConfig):
    """
    I am the configuration for pairwise permutation analysis.

    I compute a coefficient for all C(n_raters, 2) pairs of raters. This is
    useful for detecting which raters may be problematic or for understanding
    how each pair of raters agrees independently.

    I enforce that:
    - coefficient_to_compute is in PAIRWISE_PERMUTATION_COEFFICIENTS
    - config_of_coefficient is the correct config type for that coefficient
    - subset_k must be None (pairwise is fixed at k=2)
    - if bootstrap_mean_config uses the "independent" method,
      config_of_coefficient must have per-pair bootstrap enabled
      (n_resamples >= 2)
    """
    # Attribute declaration

    # Inherited from PermutationConfig
    coefficient_to_compute: str
    config_of_coefficient: AnalyzerConfig # In this object, checked that you want an analyzer that is consistent with pairwise permutations
    subset_k: int | None  # Checks that this value is None. Subsets should not happen in Pairwise permutations, but a user may put by mistake from a leftover config
    bootstrap_mean_config: PermutationBootstrapMeanConfig | None

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
                f"permutation (pairwise is always k=2). "
                f"Got subset_k={self.subset_k}"
            )

        # Validate that coefficient is pairwise-capable
        if self.coefficient_to_compute not in PAIRWISE_PERMUTATION_COEFFICIENTS:
            raise ValueError(
                f"{_ERROR_PREFIX}: '{self.coefficient_to_compute}' does "
                f"not support pairwise permutation. Allowed coefficients: "
                f"{sorted(PAIRWISE_PERMUTATION_COEFFICIENTS)}"
            )

        # Validate that config_of_coefficient matches the coefficient
        expected_config_map = {
            "cohen_kappa": "CohenKappaConfig",
            "scotts_pi": "ScottsConfig",
            "conger_kappa": "CongerKappaConfig",
            "gwet_ac1": "GwetAC1Config",
            "fleiss_kappa": "FleissKappaConfig",
            "brennan_prediger": "BrennanPredigerConfig",
            "krippendorff_alpha": "KrippendorffAlphaConfig",
        }

        expected_config_class = expected_config_map[self.coefficient_to_compute]
        actual_config_class = type(self.config_of_coefficient).__name__

        if actual_config_class != expected_config_class:
            raise TypeError(
                f"{_ERROR_PREFIX}: coefficient_to_compute="
                f"'{self.coefficient_to_compute}' requires "
                f"config_of_coefficient to be {expected_config_class}, "
                f"got {actual_config_class}"
            )

        # Cross-validate: "independent" mean bootstrap requires per-pair
        # bootstrap enabled in config_of_coefficient.
        # The "independent" method assembles the mean from independently-
        # drawn per-pair bootstrap values, so each pair must be
        # bootstrapped. The "correlated" method does its own subject
        # resampling and does not depend on per-pair bootstrap.
        if (
            self.bootstrap_mean_config is not None
            and self.bootstrap_mean_config.method == "independent"
            and self.config_of_coefficient.n_resamples < 2
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: bootstrap_mean_config method "
                f"'independent' requires per-pair bootstrap enabled "
                f"in config_of_coefficient (n_resamples >= 2), "
                f"but n_resamples={self.config_of_coefficient.n_resamples}"
            )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_dict(cls, data_dict: dict) -> PairwisePermutationConfig:
        """
        Build a pairwise permutation config from a plain dictionary.

        I do everything the parent's ``_prepare_dict`` does (validate
        inputs, normalize YAML aliases, hydrate bootstrap_mean_config),
        then I hydrate ``config_of_coefficient`` into the correct
        AnalyzerConfig subclass based on ``coefficient_to_compute``.

        If ``coefficient_to_compute`` is not in the pairwise map, I skip
        hydration and let ``__post_init__`` reject it with its own clear
        error message.

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

        # Map pairwise coefficients to their config classes
        from .cohen_kappa_config import CohenKappaConfig
        from .conger_kappa_config import CongerKappaConfig
        from .scotts_config import ScottsConfig
        from .gwet_ac1_config import GwetAC1Config
        from .fleiss_kappa_config import FleissKappaConfig
        from .brennan_prediger_config import BrennanPredigerConfig
        from .krippendorff_alpha_config import KrippendorffAlphaConfig

        config_class_map = {
            "cohen_kappa": CohenKappaConfig,
            "scotts_pi": ScottsConfig,
            "conger_kappa": CongerKappaConfig,
            "gwet_ac1": GwetAC1Config,
            "fleiss_kappa": FleissKappaConfig,
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
        """Return configuration (inherited from parent)."""
        return super().summary()
