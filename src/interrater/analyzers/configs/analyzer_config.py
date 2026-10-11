"""Base configuration for all inter-rater agreement analyzers."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[AnalyzerConfig] Error"

@dataclass(frozen=True)
class AnalyzerConfig:
    """
    I am the base configuration for all inter-rater agreement analyzers.

    Every analyzer needs these decisions: how many bootstrap resamples, what
    confidence level, which CI method, and what random seed for
    reproducibility. I enforce these choices upfront so they're documented in
    YAML and auditable in every analysis.

    I also hold a placeholder for convergence diagnostics -- the ability to
    verify that your n_resamples is large enough that CI bounds have
    stabilized. This is computationally expensive but statistically rigorous,
    and it is what separates responsible analysis from wishful thinking.

    Every analysis includes package version in its audit trail for
    reproducibility.

    Attributes
    ----------
    n_resamples : int, default 1000
        Number of bootstrap draws to take. Higher is more stable but slower.
    alpha : float, default 0.05
        Significance level. For alpha=0.05, compute 95% CI via [0.025, 0.975].
    random_state : int | None, default 42
        Random seed for reproducibility. None for unseeded randomness.
    ci_method : str, default "percentile"
        Method for computing confidence interval: "percentile" or "bca"
        (bias-corrected and accelerated).
    convergence_test : ConvergenceTest | None, default None
        Optional convergence diagnostics. When set, verifies that n_resamples
        is large enough for stable CI bounds. None means skip (for speed in
        exploratory work).
    hypothesis_test_config : HypothesisTestConfig | None, default None
        Optional configuration for hypothesis testing against H₀: coefficient = 0.
        When set, defines the method (permutation or parametric_null), permutation
        scheme, and number of permutations. None means no hypothesis test.
    injection : tuple[InjectorConfig, ...] | None, default None
        Injection pipeline -- a sequence of injection steps applied to
        the Ratings before analysis. Each step adds synthetic raters.
        Steps are applied in order (chaining). None means no injection.
        Accepts a single InjectorConfig, a list, or a tuple -- all are
        normalized to a tuple in __post_init__.
        In YAML, this is a list of dicts, each with a ``type`` key that
        dispatches to the appropriate injector config class.
    """

    n_resamples: int = 1000
    alpha: float = 0.05
    random_state: int | None = 42
    ci_method: str = "percentile"
    convergence_test: ConvergenceTest | None = None
    hypothesis_test_config: HypothesisTestConfig | None = None
    injection: tuple[InjectorConfig] | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate configuration after construction."""

        # Check the bool python bug. After around Python2.3, bools were made to inherit
        # int for backwards compatibility 
        for name in ("n_resamples", "alpha", "random_state"):
            if isinstance(getattr(self, name), bool):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} cannot be bool, "
                    f"got {type(getattr(self, name)).__name__}"
                )

        # Validate n_resamples
        if not isinstance(self.n_resamples, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_resamples must be int, "
                f"got {type(self.n_resamples).__name__}"
            )
        if self.n_resamples < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_resamples must be >= 1, got {self.n_resamples}"
            )

        # Validate alpha
        if not isinstance(self.alpha, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: alpha must be float, "
                f"got {type(self.alpha).__name__}"
            )
        if not (0 < self.alpha < 1):
            raise ValueError(
                f"{_ERROR_PREFIX}: alpha must be in (0, 1), "
                f"got {self.alpha}"
            )

        # Validate random_state
        if self.random_state is not None:
            if not isinstance(self.random_state, int):
                raise TypeError(
                    f"{_ERROR_PREFIX}: random_state must be int or None, "
                    f"got {type(self.random_state).__name__}"
                )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be a string, "
                f"got {type(self.ci_method).__name__}"
            )

        if self.ci_method not in {"percentile", "bca", "parametric_null"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile', 'bca', or 'parametric_null', "
                f"got '{self.ci_method}'"
            )

        # Validate convergence_test type
        if self.convergence_test is not None:
            from ..convergence_test import ConvergenceTest

            if not isinstance(self.convergence_test, ConvergenceTest):
                raise TypeError(
                    f"{_ERROR_PREFIX}: convergence_test must be ConvergenceTest or None, "
                    f"got {type(self.convergence_test).__name__}"
                )

        # Validate hypothesis_test_config type
        if self.hypothesis_test_config is not None:
            from .hypothesis_test_config import HypothesisTestConfig

            if not isinstance(self.hypothesis_test_config, HypothesisTestConfig):
                raise TypeError(
                    f"{_ERROR_PREFIX}: hypothesis_test_config must be HypothesisTestConfig or None, "
                    f"got {type(self.hypothesis_test_config).__name__}"
                )

        # Normalize and validate injection
        if self.injection is not None:
            from ...injectors.configs import InjectorConfig

            # Normalize: single config → tuple, list → tuple
            if isinstance(self.injection, InjectorConfig):
                object.__setattr__(
                    self, "injection", (self.injection,)
                )
            elif isinstance(self.injection, list):
                object.__setattr__(
                    self, "injection", tuple(self.injection)
                )

            # Validate type is tuple after normalization
            if not isinstance(self.injection, tuple):
                raise TypeError(
                    f"{_ERROR_PREFIX}: injection must be a tuple, list, "
                    f"or single InjectorConfig (or None), "
                    f"got {type(self.injection).__name__}"
                )

            # Validate non-empty
            if len(self.injection) == 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: injection cannot be empty. "
                    f"Provide at least one injection step, or None "
                    f"to skip injection."
                )

            # Validate each element
            for i, step in enumerate(self.injection):
                if not isinstance(step, InjectorConfig):
                    raise TypeError(
                        f"{_ERROR_PREFIX}: injection[{i}] must be "
                        f"an InjectorConfig, "
                        f"got {type(step).__name__}"
                    )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_dict(cls, data_dict: dict) -> AnalyzerConfig:
        """
        Build a config from a plain dictionary.

        I normalize YAML key aliases and hydrate nested objects.
        Hydration means converting a plain sub-dict into its typed
        config object -- e.g. a ``hypothesis_test`` dict becomes a
        ``HypothesisTestConfig``.

        Parameters
        ----------
        data_dict : dict
            Configuration as plain Python types. Mutated in place
            (alias keys are renamed, nested sub-dicts are replaced
            with hydrated config objects).

        Returns
        -------
        AnalyzerConfig
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
        >>> config = CohenKappaConfig.build_from_dict(
        ...     {"n_resamples": 2000, "alpha": 0.05}
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

        from ...yaml_utils import normalize_yaml_key
        from .hypothesis_test_config import HypothesisTestConfig

        # Normalize alias → canonical key
        normalize_yaml_key(data_dict, "hypothesis_test_config")

        # Hydrate: convert sub-dict to typed config object
        if "hypothesis_test_config" in data_dict:
            ht_data = data_dict["hypothesis_test_config"]
            if isinstance(ht_data, dict):
                data_dict["hypothesis_test_config"] = (
                    HypothesisTestConfig(**ht_data)
                )

        # Hydrate injection steps: dict or list of dicts → tuple of InjectorConfigs
        if "injection" in data_dict:
            inj = data_dict["injection"]
            # Normalize single dict → list (YAML shorthand for one step)
            if isinstance(inj, dict):
                inj = [inj]
            if isinstance(inj, list) and len(inj) > 0 and isinstance(inj[0], dict):
                from ...injectors.registry import hydrate_injection_step

                data_dict["injection"] = tuple(
                    hydrate_injection_step(step) for step in inj
                )

        return cls(**data_dict)

    @classmethod
    def build_from_yaml(cls, path: str) -> AnalyzerConfig:
        """
        Load configuration from a YAML file.

        I read the file, validate that all keys are recognized, and
        delegate to ``build_from_dict`` for alias normalization and
        hydration.

        Parameters
        ----------
        path : str
            Path to YAML file containing analyzer configuration.

        Returns
        -------
        AnalyzerConfig
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
        >>> config = FleissKappaConfig.build_from_yaml("fleiss_config.yaml")
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
        Return configuration as a dictionary for documentation and audit.

        Returns
        -------
        dict
            Plain-type representation of all configuration fields, including
            package version for reproducibility.
        """
        result = {
            "n_resamples": self.n_resamples,
            "alpha": self.alpha,
            "random_state": self.random_state,
            "ci_method": self.ci_method,
        }

        if self.convergence_test is not None:
            result["convergence_test"] = self.convergence_test.summary()

        if self.hypothesis_test_config is not None:
            result["hypothesis_test"] = self.hypothesis_test_config.summary()

        if self.injection is not None:
            result["injection"] = [
                step.summary() for step in self.injection
            ]

        # Add package version to audit trail for reproducibility
        try:
            from ... import __version__

            result["_package_version"] = __version__
        except ImportError:
            result["_package_version"] = "unknown"

        return result
