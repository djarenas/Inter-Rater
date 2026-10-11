"""Configuration for pooled rater injection."""

from __future__ import annotations

from dataclasses import dataclass

from .injector_config import InjectorConfig

_ERROR_PREFIX = "[PooledInjectorConfig] Error"


@dataclass(frozen=True)
class PooledInjectorConfig(InjectorConfig):
    """
    I configure a rater who samples from the entire pool of all ratings.

    At every subject, each injected rater randomly selects one rating
    from anywhere in the data -- any rater, any subject. This preserves
    the empirical global prevalences but destroys all subject-level
    structure.

    With smoothing > 0, I blend the empirical pool distribution with a
    uniform distribution over categories. At smoothing=0.0 you get the
    pure empirical marginals; at smoothing=1.0 the rater draws uniformly
    regardless of the data.

    With exclude_raters, I leave specific raters out of the pool before
    sampling. Useful for leave-one-out sensitivity analysis.

    "What if we added someone who uses the right categories at the right
    rates but doesn't actually look at the subjects?"

    Attributes
    ----------
    smoothing : float
        Blend between empirical pool distribution (0.0) and uniform
        over categories (1.0). Default 0.0 means pure empirical.
    exclude_raters : tuple of int
        Rater indices to exclude from the pool. Empty tuple means
        use all raters. Bounds checking happens at generate time
        in the injector, not here.
    n_new_raters : int
        Inherited from InjectorConfig.
    random_state : int or None
        Inherited from InjectorConfig.
    """

    smoothing: float = 0.0
    exclude_raters: tuple[int, ...] = ()

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate pooled-specific fields."""
        try:
            super().__post_init__()
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # smoothing: float in [0.0, 1.0], no bools
        if isinstance(self.smoothing, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: smoothing must be a number, "
                f"got bool"
            )
        if not isinstance(self.smoothing, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: smoothing must be a number, "
                f"got {type(self.smoothing).__name__}"
            )
        if not (0.0 <= self.smoothing <= 1.0):
            raise ValueError(
                f"{_ERROR_PREFIX}: smoothing must be in [0.0, 1.0], "
                f"got {self.smoothing}"
            )

        # exclude_raters: tuple of non-negative ints, no duplicates
        if not isinstance(self.exclude_raters, tuple):
            raise TypeError(
                f"{_ERROR_PREFIX}: exclude_raters must be a tuple, "
                f"got {type(self.exclude_raters).__name__}"
            )
        for i, val in enumerate(self.exclude_raters):
            if isinstance(val, bool) or not isinstance(val, int):
                raise TypeError(
                    f"{_ERROR_PREFIX}: exclude_raters[{i}] must be int, "
                    f"got {type(val).__name__}"
                )
            if val < 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: exclude_raters[{i}] must be "
                    f"non-negative, got {val}"
                )
        if len(self.exclude_raters) != len(set(self.exclude_raters)):
            raise ValueError(
                f"{_ERROR_PREFIX}: exclude_raters contains duplicates"
            )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_dict(cls, data_dict: dict) -> PooledInjectorConfig:
        """
        Build a config from a plain dictionary.

        I convert exclude_raters from a list (as YAML provides) to a
        tuple, then delegate to the base class.

        Parameters
        ----------
        data_dict : dict
            Configuration as plain Python types.

        Returns
        -------
        PooledInjectorConfig
            Instantiated config object with validation applied.
        """
        # Hydrate: list -> tuple for exclude_raters
        if "exclude_raters" in data_dict:
            val = data_dict["exclude_raters"]
            if isinstance(val, list):
                data_dict["exclude_raters"] = tuple(val)

        return super().build_from_dict(data_dict)

    def summary(self) -> dict:
        """
        Return configuration as a dictionary for documentation and audit.

        Returns
        -------
        dict
            Base fields plus smoothing and exclude_raters.
        """
        result = super().summary()
        result["smoothing"] = self.smoothing
        result["exclude_raters"] = self.exclude_raters
        return result
