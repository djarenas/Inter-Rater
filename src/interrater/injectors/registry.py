"""Registry mapping YAML type strings to injector config and class pairs."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .configs import InjectorConfig
    from .rater_injector import RaterInjector

_ERROR_PREFIX = "[InjectorRegistry] Error"


def _build_registry() -> dict[str, tuple[type[InjectorConfig], type[RaterInjector]]]:
    """
    I build the registry lazily to avoid circular imports.

    Returns
    -------
    dict[str, tuple[type[InjectorConfig], type[RaterInjector]]]
        Maps YAML type strings to (ConfigClass, InjectorClass) pairs.
    """
    from .configs import (
        DistributionInjectorConfig,
        ExtremalInjectorConfig,
        PooledInjectorConfig,
        RaterBootstrapInjectorConfig,
        WithinSubjectInjectorConfig,
    )
    from .distribution_injector import DistributionInjector
    from .extremal_injector import ExtremalInjector
    from .pooled_injector import PooledInjector
    from .rater_bootstrap_injector import RaterBootstrapInjector
    from .within_subject_injector import WithinSubjectInjector

    return {
        "distribution": (DistributionInjectorConfig, DistributionInjector),
        "within_subject": (WithinSubjectInjectorConfig, WithinSubjectInjector),
        "pooled": (PooledInjectorConfig, PooledInjector),
        "rater_bootstrap": (RaterBootstrapInjectorConfig, RaterBootstrapInjector),
        "extremal": (ExtremalInjectorConfig, ExtremalInjector),
    }


def get_registry() -> dict[str, tuple[type[InjectorConfig], type[RaterInjector]]]:
    """
    I return the injector registry, building it on first call.

    The registry maps YAML type strings (e.g. ``"distribution"``,
    ``"rater_bootstrap"``) to ``(ConfigClass, InjectorClass)`` pairs.

    - Config side uses the first element for YAML hydration.
    - Runner side uses both to instantiate the injector from its config.

    Returns
    -------
    dict[str, tuple[type[InjectorConfig], type[RaterInjector]]]
        The registry mapping.
    """
    global _REGISTRY
    if _REGISTRY is None:
        _REGISTRY = _build_registry()
    return _REGISTRY


def hydrate_injection_step(step_dict: dict) -> InjectorConfig:
    """
    I convert a single YAML injection step dict into a typed config.

    I pop the ``type`` key, look it up in the registry, and delegate
    to the config class's ``build_from_dict``.

    Parameters
    ----------
    step_dict : dict
        A single injection step from YAML, e.g.
        ``{"type": "rater_bootstrap", "source_raters": [0, 1], ...}``.
        Mutated in place (``type`` key is popped).

    Returns
    -------
    InjectorConfig
        The typed config object (e.g. RaterBootstrapInjectorConfig).

    Raises
    ------
    ValueError
        If ``type`` key is missing or not a recognized injector type.
    TypeError
        If step_dict is not a dict.
    """
    if not isinstance(step_dict, dict):
        raise TypeError(
            f"{_ERROR_PREFIX}: each injection step must be a dict, "
            f"got {type(step_dict).__name__}"
        )

    if "type" not in step_dict:
        raise ValueError(
            f"{_ERROR_PREFIX}: each injection step must have a "
            f"'type' key. Got keys: {sorted(step_dict.keys())}"
        )

    registry = get_registry()
    step_type = step_dict.pop("type")

    if step_type not in registry:
        raise ValueError(
            f"{_ERROR_PREFIX}: unknown injection type '{step_type}'. "
            f"Valid types: {sorted(registry.keys())}"
        )

    config_cls, _ = registry[step_type]
    return config_cls.build_from_dict(step_dict)


def get_injector_class(config: InjectorConfig) -> type[RaterInjector]:
    """
    I return the injector class for a given config instance.

    I build a reverse mapping from ConfigClass → InjectorClass on
    first call, then look up the config's type. This is used at
    runtime by ``Analyzer.apply_injection()`` to go from a typed
    config object to the injector that knows how to execute it.

    Parameters
    ----------
    config : InjectorConfig
        A typed config instance (e.g. ExtremalInjectorConfig).

    Returns
    -------
    type[RaterInjector]
        The injector class (e.g. ExtremalInjector).

    Raises
    ------
    TypeError
        If config is not an InjectorConfig.
    ValueError
        If no injector is registered for this config type.
    """
    from .configs import InjectorConfig as InjectorConfigBase

    if not isinstance(config, InjectorConfigBase):
        raise TypeError(
            f"{_ERROR_PREFIX}: get_injector_class: config must be "
            f"an InjectorConfig, got {type(config).__name__}"
        )

    global _REVERSE_REGISTRY
    if _REVERSE_REGISTRY is None:
        registry = get_registry()
        _REVERSE_REGISTRY = {
            config_cls: injector_cls
            for config_cls, injector_cls in registry.values()
        }

    config_type = type(config)
    if config_type not in _REVERSE_REGISTRY:
        raise ValueError(
            f"{_ERROR_PREFIX}: get_injector_class: no injector "
            f"registered for config type {config_type.__name__}. "
            f"Registered config types: "
            f"{sorted(c.__name__ for c in _REVERSE_REGISTRY)}"
        )

    return _REVERSE_REGISTRY[config_type]


# Module-level caches
_REGISTRY: dict | None = None
_REVERSE_REGISTRY: dict | None = None
