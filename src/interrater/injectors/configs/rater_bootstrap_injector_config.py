"""Configuration for rater bootstrap injection."""

from __future__ import annotations

from dataclasses import dataclass

from .injector_config import InjectorConfig

_ERROR_PREFIX = "[RaterBootstrapInjectorConfig] Error"


@dataclass(frozen=True)
class RaterBootstrapInjectorConfig(InjectorConfig):
    """
    I configure a rater cloned from one or more existing raters' columns.

    The injected rater is built by pooling the selected source raters'
    non-EMPTY ratings and resampling with replacement. The new rater has
    the same marginal distribution as the source subset but no
    subject-specific information -- their habits are cloned, their
    judgments are scrambled.

    With a single source rater, this clones one person's biases. With
    all raters as sources, this collapses into the PooledInjector
    (without smoothing or exclusion knobs).

    "What if we added someone with rater k's biases but who can't
    actually distinguish between subjects?"

    Attributes
    ----------
    source_raters : tuple of int
        Column indices of the real raters to clone from. Their
        non-EMPTY ratings are pooled, then resampled with replacement.
        Bounds checking happens at generate time in the injector,
        not here.
    n_new_raters : int
        Inherited from InjectorConfig.
    random_state : int or None
        Inherited from InjectorConfig.
    """

    source_raters: tuple[int, ...] = (0,)

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate rater-bootstrap-specific fields."""
        try:
            super().__post_init__()
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # source_raters: non-empty tuple of non-negative ints, no duplicates
        if not isinstance(self.source_raters, tuple):
            raise TypeError(
                f"{_ERROR_PREFIX}: source_raters must be a tuple, "
                f"got {type(self.source_raters).__name__}"
            )
        if len(self.source_raters) == 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: source_raters cannot be empty. "
                f"Provide at least one rater index to clone from."
            )
        for i, val in enumerate(self.source_raters):
            if isinstance(val, bool) or not isinstance(val, int):
                raise TypeError(
                    f"{_ERROR_PREFIX}: source_raters[{i}] must be int, "
                    f"got {type(val).__name__}"
                )
            if val < 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: source_raters[{i}] must be "
                    f"non-negative, got {val}"
                )
        if len(self.source_raters) != len(set(self.source_raters)):
            raise ValueError(
                f"{_ERROR_PREFIX}: source_raters contains duplicates"
            )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_dict(cls, data_dict: dict) -> RaterBootstrapInjectorConfig:
        """
        Build a config from a plain dictionary.

        I convert source_raters from a list (as YAML provides) to a
        tuple, then delegate to the base class.

        Parameters
        ----------
        data_dict : dict
            Configuration as plain Python types.

        Returns
        -------
        RaterBootstrapInjectorConfig
            Instantiated config object with validation applied.
        """
        # Hydrate: list -> tuple for source_raters
        if "source_raters" in data_dict:
            val = data_dict["source_raters"]
            if isinstance(val, list):
                data_dict["source_raters"] = tuple(val)

        return super().build_from_dict(data_dict)

    def summary(self) -> dict:
        """
        Return configuration as a dictionary for documentation and audit.

        Returns
        -------
        dict
            Base fields plus source_raters.
        """
        result = super().summary()
        result["source_raters"] = self.source_raters
        return result
