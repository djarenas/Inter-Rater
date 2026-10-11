"""Base configuration for all rater injectors."""

from __future__ import annotations

from dataclasses import dataclass, fields

_ERROR_PREFIX = "[InjectorConfig] Error"


@dataclass(frozen=True)
class InjectorConfig:
    """
    I am the base configuration for rater injection strategies.

    I hold parameters that every injector needs. Leaf configs inherit me
    and add strategy-specific fields (prevalences, source rater index,
    noise rate, etc.).

    I am frozen -- once built, nothing changes. The injector reads me
    but never writes me.

    Every analysis includes package version in its audit trail for
    reproducibility.

    Attributes
    ----------
    n_new_raters : int
        How many synthetic raters to inject. Each draws independently.
    random_state : int or None
        Seed for reproducibility. None means unseeded (non-reproducible).
        Every injector uses this to initialise its RNG.
    """

    n_new_raters: int = 1
    random_state: int | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate base fields."""
        # n_new_raters: positive int, no bools
        if isinstance(self.n_new_raters, bool) or not isinstance(
            self.n_new_raters, int
        ):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_new_raters must be int, "
                f"got {type(self.n_new_raters).__name__}"
            )
        if self.n_new_raters < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_new_raters must be >= 1, "
                f"got {self.n_new_raters}"
            )

        # random_state: int or None, no bools
        if self.random_state is not None:
            if isinstance(self.random_state, bool) or not isinstance(
                self.random_state, int
            ):
                raise TypeError(
                    f"{_ERROR_PREFIX}: random_state must be int or None, "
                    f"got {type(self.random_state).__name__}"
                )
            if self.random_state < 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: random_state must be non-negative, "
                    f"got {self.random_state}"
                )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_dict(cls, data_dict: dict) -> InjectorConfig:
        """
        Build a config from a plain dictionary.

        I validate the dict, and delegate to ``cls(**data_dict)`` so
        subclass constructors handle their own fields. Child classes
        override me to hydrate nested sub-dicts, calling
        ``super().build_from_dict(data_dict)`` after their hydration.

        Parameters
        ----------
        data_dict : dict
            Configuration as plain Python types.

        Returns
        -------
        InjectorConfig
            Instantiated config object with validation applied.

        Raises
        ------
        TypeError
            If data_dict is not a dict.
        ValueError
            If data_dict is None or configuration values fail
            validation.

        Examples
        --------
        >>> config = DistributionInjectorConfig.build_from_dict(
        ...     {"prevalences": [0.5, 0.3, 0.2], "n_new_raters": 2}
        ... )
        """
        if data_dict is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: build_from_dict: data_dict must be "
                f"a dictionary. Got None type"
            )
        if not isinstance(data_dict, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: build_from_dict: data_dict must be "
                f"a dictionary. Got {type(data_dict).__name__}"
            )

        return cls(**data_dict)

    @classmethod
    def build_from_yaml(cls, path: str) -> InjectorConfig:
        """
        Load configuration from a YAML file.

        I read the file, validate that all keys are recognized, and
        delegate to ``build_from_dict`` for hydration.

        Parameters
        ----------
        path : str
            Path to YAML file containing injector configuration.

        Returns
        -------
        InjectorConfig
            Instantiated config object with validation applied.

        Raises
        ------
        FileNotFoundError
            If the YAML file does not exist.
        TypeError
            If the YAML contains unknown keys.
        ValueError
            If configuration values fail validation.

        Examples
        --------
        >>> config = DistributionInjectorConfig.build_from_yaml(
        ...     "distribution_config.yaml"
        ... )
        """
        import yaml

        with open(path) as f:
            data_dict = yaml.safe_load(f)

        # Check for unknown YAML keys
        valid_fields = {f.name for f in fields(cls)}
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
        Return configuration as a dictionary for documentation and audit.

        Child classes override me, calling ``super().summary()`` to get
        the base dict and extending it with their own fields.

        Returns
        -------
        dict
            Plain-type representation of all configuration fields,
            including package version for reproducibility.
        """
        result = {
            "n_new_raters": self.n_new_raters,
            "random_state": self.random_state,
        }

        # Add package version to audit trail for reproducibility
        try:
            from ... import __version__

            result["_package_version"] = __version__
        except ImportError:
            result["_package_version"] = "unknown"

        return result
