"""Configuration for extremal rater injection."""

from __future__ import annotations

from dataclasses import dataclass

from .injector_config import InjectorConfig

_ERROR_PREFIX = "[ExtremalInjectorConfig] Error"

_VALID_STRATEGIES = ("majority", "minority", "median", "contrarian")


@dataclass(frozen=True)
class ExtremalInjectorConfig(InjectorConfig):
    """
    I configure a deterministic rater who picks an extremal category.

    At each subject, the injected rater looks at how the real raters
    voted and picks a category based on prevalence rank:

    - ``"majority"``: the most common category (mode). Simulates a
      conformist who always goes with the crowd.
    - ``"minority"``: the least common category that was actually
      used. Simulates a dissenter who picks the unpopular option
      someone did choose.
    - ``"median"``: the category with middle-ranked prevalence among
      the categories that appeared at this subject. Works for nominal
      (unordered) data -- it ranks by count, not by category index.
    - ``"contrarian"``: a category that no real rater used at this
      subject. The true contrarian -- picks where nobody went. When
      all categories were used (no unused category available), falls
      back to minority behavior.

    Ties are broken randomly using ``random_state``.

    "What if we added someone who always agrees with the majority?"
    "What if we added a true contrarian who picks what nobody picked?"

    Attributes
    ----------
    strategy : str
        One of ``"majority"``, ``"minority"``, ``"median"``, or
        ``"contrarian"``. Default is ``"majority"``.
    n_new_raters : int
        Inherited from InjectorConfig.
    random_state : int or None
        Inherited from InjectorConfig. Used for tie-breaking.
    """

    strategy: str = "majority"

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate extremal-specific fields."""
        try:
            super().__post_init__()
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # strategy: must be one of the valid strings
        if not isinstance(self.strategy, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: strategy must be a string, "
                f"got {type(self.strategy).__name__}"
            )
        if self.strategy not in _VALID_STRATEGIES:
            raise ValueError(
                f"{_ERROR_PREFIX}: strategy must be one of "
                f"{_VALID_STRATEGIES}, got {self.strategy!r}"
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """
        Return configuration as a dictionary for documentation and audit.

        Returns
        -------
        dict
            Base fields plus strategy.
        """
        result = super().summary()
        result["strategy"] = self.strategy
        return result
