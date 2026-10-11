"""Configuration for permutation analysis across rater combinations."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import TYPE_CHECKING

from .permutation_bootstrap_mean_config import PermutationBootstrapMeanConfig

if TYPE_CHECKING:
    from .analyzer_config import AnalyzerConfig

_ERROR_PREFIX = "[PermutationConfig] Error"


@dataclass(frozen=True)
class PermutationConfig:
    """
    I am the configuration for permutation analysis.

    Permutation analysis computes a coefficient across multiple rater
    combinations (pairs, subsets). I hold three things:

    - The name of the coefficient to compute (as a string)
    - The nested AnalyzerConfig that specifies all analytical decisions
      for computing that coefficient (including per-combination bootstrap)
    - An optional PermutationBootstrapMeanConfig that controls whether
      and how the mean across combinations is bootstrapped

    I validate that the coefficient name is one of the known, allowed
    coefficients. Subclasses (PairwisePermutationConfig,
    SubsetPermutationConfig) further validate that the coefficient is
    appropriate for their specific permutation mode, that
    config_of_coefficient is the right type, and cross-validate
    bootstrap_mean_config against config_of_coefficient when the
    "independent" method requires per-combination bootstrap.

    Attributes
    ----------
    coefficient_to_compute : str
        The name of the coefficient: "fleiss_kappa", "cohen_kappa",
        "scotts_pi", "conger_kappa", etc. Must be one of the known
        coefficients in contracts.
    config_of_coefficient : AnalyzerConfig
        The nested config (FleissKappaConfig, CohenKappaConfig, etc.) that
        specifies analytical decisions: n_resamples, alpha, random_state,
        ci_method, weights_type, etc. This controls per-combination
        bootstrap (resampling subjects within a combination).
    subset_k : int | list[int] | str | None, default None
        Subset size(s) for subset permutation. Accepts a single int,
        a list of ints, or "all". Only used if this config is used by
        SubsetPermutationConfig (which validates fully).
        PairwisePermutationConfig will reject this in __post_init__.
    bootstrap_mean_config : PermutationBootstrapMeanConfig | None, default None
        Configuration for bootstrapping the mean coefficient across
        combinations. None means no mean-level bootstrap. When present,
        controls the procedure ("correlated" or "independent"),
        n_bootstraps, ci_method, alpha, and random_state for the mean.
    """

    # Attribute declaration

    # Required (no defaults)
    coefficient_to_compute: str
    config_of_coefficient: AnalyzerConfig

    # With defaults
    subset_k: int | list | str | None = None
    bootstrap_mean_config: PermutationBootstrapMeanConfig | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate permutation configuration."""
        from ... import contracts as _c

        from .analyzer_config import AnalyzerConfig

        # Get all known permutation coefficients
        all_permutation_coefficients = (
            _c.PAIRWISE_PERMUTATION_COEFFICIENTS
            | _c.SUBSET_PERMUTATION_COEFFICIENTS
        )

        name = "coefficient_to_compute"
        if not isinstance(getattr(self, name), str):
            raise TypeError(
                f"{_ERROR_PREFIX}: {name} must be a string, "
                f"got {type(getattr(self, name)).__name__}"
            )

        # Validate coefficient_to_compute is a known coefficient
        if self.coefficient_to_compute not in all_permutation_coefficients:
            raise ValueError(
                f"{_ERROR_PREFIX}: coefficient_to_compute must be one of "
                f"{sorted(all_permutation_coefficients)}, "
                f"got '{self.coefficient_to_compute}'"
            )

        # Validate that config_of_coefficient is an AnalyzerConfig
        if not isinstance(self.config_of_coefficient, AnalyzerConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config_of_coefficient must be an "
                f"AnalyzerConfig (FleissKappaConfig, CohenKappaConfig, "
                f"etc.), got {type(self.config_of_coefficient).__name__}"
            )

        # Validate subset_k type if provided.
        # Accept int, list, or str ("all"). Detailed validation
        # (values >= 2, no dupes, "all" only) is in SubsetPermutationConfig.
        if self.subset_k is not None:
            if isinstance(self.subset_k, bool):
                raise TypeError(
                    f"{_ERROR_PREFIX}: subset_k must be int, list[int], "
                    f"'all', or None, "
                    f"got {type(self.subset_k).__name__}"
                )          
            if not isinstance(self.subset_k, (int, list, str)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: subset_k must be int, list[int], "
                    f"'all', or None, "
                    f"got {type(self.subset_k).__name__}"
                )

        # Validate bootstrap_mean_config type if provided
        if self.bootstrap_mean_config is not None:
            if not isinstance(
                self.bootstrap_mean_config, PermutationBootstrapMeanConfig
            ):
                raise TypeError(
                    f"{_ERROR_PREFIX}: bootstrap_mean_config must be "
                    f"PermutationBootstrapMeanConfig or None, "
                    f"got {type(self.bootstrap_mean_config).__name__}"
                )


    # Public API ------------------------------------------------------------- #

    @classmethod
    def _prepare_dict(cls, data_dict: dict) -> dict:
        """
        Validate inputs and hydrate shared nested configs.

        I do the work that is common to all PermutationConfig children:
        validate that data_dict is a dict, normalize YAML key aliases,
        and hydrate bootstrap_mean_config from a sub-dict into a
        PermutationBootstrapMeanConfig object.

        I do NOT hydrate config_of_coefficient -- that requires knowing
        which AnalyzerConfig subclass to use, which is the child's job.

        Parameters
        ----------
        data_dict : dict
            Raw configuration dict. Mutated in place.

        Returns
        -------
        dict
            The same dict, with bootstrap_mean_config hydrated.

        Raises
        ------
        TypeError
            If data_dict is not a dict.
        ValueError
            If data_dict is None.
        """
        if data_dict is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: _prepare_dict: data_dict must be "
                f"a dictionary. Got None type"
            )
        if not isinstance(data_dict, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: _prepare_dict: data_dict must be "
                f"a dictionary. Got {type(data_dict).__name__}"
            )

        from ...yaml_utils import normalize_yaml_key

        # Normalize alias → canonical key
        normalize_yaml_key(data_dict, "bootstrap_mean_config")

        # Hydrate: convert sub-dict to typed config object
        if "bootstrap_mean_config" in data_dict:
            bmc = data_dict["bootstrap_mean_config"]
            if isinstance(bmc, dict):
                data_dict["bootstrap_mean_config"] = (
                    PermutationBootstrapMeanConfig(**bmc)
                )

        return data_dict

    @classmethod
    def build_from_dict(cls, data_dict: dict) -> PermutationConfig:
        """
        Build a config from a plain dictionary.

        I normalize YAML key aliases and hydrate nested objects.
        Hydration means converting a plain sub-dict into its typed
        config object -- e.g. a ``bootstrap_mean`` dict becomes a
        ``PermutationBootstrapMeanConfig``.

        Children override this to also hydrate config_of_coefficient
        into the correct AnalyzerConfig subclass, calling
        ``cls._prepare_dict(data_dict)`` for the shared work.

        Parameters
        ----------
        data_dict : dict
            Configuration as plain Python types. Mutated in place.

        Returns
        -------
        PermutationConfig
            Instantiated config object with validation applied.

        Raises
        ------
        TypeError
            If data_dict is not a dict, or if config_of_coefficient
            is still a raw dict (not hydrated by a child class).
        ValueError
            If data_dict is None or validation fails.
        """
        cls._prepare_dict(data_dict)
        return cls(**data_dict)

    @classmethod
    def build_from_yaml(cls, path: str) -> PermutationConfig:
        """
        Load configuration from a YAML file.

        I read the file, validate that all keys are recognized, and
        delegate to ``build_from_dict`` for alias normalization and
        hydration. Children that override ``build_from_dict`` get
        correct dispatch via ``cls`` without overriding this method.

        Parameters
        ----------
        path : str
            Path to YAML file.

        Returns
        -------
        PermutationConfig
            Instantiated config object with validation applied.

        Raises
        ------
        FileNotFoundError
            If the YAML file does not exist.
        TypeError
            If the YAML contains unknown keys.
        ValueError
            If configuration values fail validation.
        """
        import yaml

        from ...yaml_utils import YAML_KEY_ALIASES

        with open(path) as f:
            data_dict = yaml.safe_load(f)

        # Check for unknown YAML keys (allow aliases too)
        valid_fields = {f.name for f in fields(cls)}
        for aliases in YAML_KEY_ALIASES.values():
            valid_fields.update(aliases)
        unknown = set(data_dict) - valid_fields
        if unknown:
            raise TypeError(
                f"{_ERROR_PREFIX}: build_from_yaml: Unknown YAML "
                f"key(s): {unknown}. Valid keys are: "
                f"{sorted(valid_fields)}. "
                f"Check your YAML configuration file."
            )

        return cls.build_from_dict(data_dict)

    def summary(self) -> dict:
        """
        Return configuration including nested coefficient config.

        Returns
        -------
        dict
            Keys: coefficient_to_compute, subset_k,
            config_of_coefficient (nested summary),
            bootstrap_mean_config (nested summary or None).
        """
        result = {
            "coefficient_to_compute": self.coefficient_to_compute,
            "subset_k": self.subset_k,
            "config_of_coefficient": self.config_of_coefficient.summary(),
            "bootstrap_mean_config": (
                self.bootstrap_mean_config.summary()
                if self.bootstrap_mean_config is not None
                else None
            ),
        }
        return result
