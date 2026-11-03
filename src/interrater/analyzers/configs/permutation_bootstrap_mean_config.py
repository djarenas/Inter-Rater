"""Configuration for bootstrapping the mean across rater combinations."""

from __future__ import annotations

from dataclasses import dataclass, field

_ERROR_PREFIX = "[PermutationBootstrapMeanConfig] Error"

_METHODS = ("correlated", "independent")

# Removed in 2.1.0, and where each now comes from in config_of_coefficient
_REMOVED = {
    "n_bootstraps": "n_resamples",
    "ci_method": "ci_method",
    "alpha": "alpha",
    "random_state": "random_state",
}


@dataclass(frozen=True)
class PermutationBootstrapMeanConfig:
    """
    I choose how the mean coefficient across rater combinations (pairs or
    subsets) is bootstrapped.

    I am bundled as one object so that ``None`` on the parent
    PermutationConfig means "no mean-level bootstrap." When I am present,
    the permutation analyzer bootstraps the mean and builds its interval.

    I hold only the method. The draw count, interval method, level and
    seed come from the parent's ``config_of_coefficient``, the same ones
    the per-combination intervals use.

    Two bootstrap procedures are supported:

    - "correlated": resample subjects once per draw, compute all
      combinations on the same resampled data, then average. Preserves
      the correlation structure between combinations that share a rater.
    - "independent": draw one value from each combination's own
      bootstrap draws, then average. It ignores the correlation between
      combinations that share raters, so its interval tends to be too
      narrow. It needs at least 2 draws per combination and
      ci_method='percentile'; the parent config checks both.

    Attributes
    ----------
    method : str
        Bootstrap procedure: "correlated" or "independent".
    n_bootstraps, ci_method, alpha, random_state : None
        Removed in 2.1.0. Kept as fields only so that a call or a config
        that still sets one is told where it moved, instead of failing on
        an unknown argument. Any value is refused.
    """

    method: str

    # Removed in 2.1.0; any value is refused with where it moved
    n_bootstraps: None = field(default=None, repr=False, compare=False)
    ci_method: None = field(default=None, repr=False, compare=False)
    alpha: None = field(default=None, repr=False, compare=False)
    random_state: None = field(default=None, repr=False, compare=False)

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate the method."""
        if not isinstance(self.method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: method must be str, "
                f"got {type(self.method).__name__}"
            )
        if self.method not in _METHODS:
            raise ValueError(
                f"{_ERROR_PREFIX}: method must be 'correlated' or "
                f"'independent', got '{self.method}'"
            )

        removed = [name for name in _REMOVED if getattr(self, name) is not None]
        if removed:
            moved = ", ".join(f"{name} -> {_REMOVED[name]}" for name in removed)
            raise TypeError(
                f"{_ERROR_PREFIX}: {removed} were removed in 2.1.0. The "
                f"mean's draw count, interval method, level and seed are the "
                f"ones in config_of_coefficient ({moved}), the same ones the "
                f"per-combination intervals use. Delete these from the "
                f"bootstrap mean config; the only key is 'method'."
            )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def from_dict(cls, data: dict) -> PermutationBootstrapMeanConfig:
        """
        Build from a plain dict, as read from YAML.

        Raises
        ------
        TypeError
            If data is not a dict, or holds a key removed in 2.1.0 (the
            message says where it moved).
        ValueError
            If it holds any other key besides ``method``.
        """
        if not isinstance(data, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: expected a dict, got {type(data).__name__}"
            )
        # Removed keys go to the constructor, which says where they moved
        unknown = sorted(set(data) - {"method"} - set(_REMOVED))
        if unknown:
            raise ValueError(
                f"{_ERROR_PREFIX}: unknown key(s) {unknown}; the only key is "
                f"'method'. The draw count, interval method, level and seed "
                f"come from config_of_coefficient."
            )
        return cls(**data)

    def summary(self) -> dict:
        """
        Return configuration in plain Python types.

        Returns
        -------
        dict
            Keys: method.
        """
        return {"method": self.method}
