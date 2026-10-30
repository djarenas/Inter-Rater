"""Configuration for bootstrapping the mean across rater combinations."""

from __future__ import annotations

from dataclasses import dataclass

_ERROR_PREFIX = "[PermutationBootstrapMeanConfig] Error"

_METHODS = ("correlated", "independent")


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
    """

    method: str

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

    # Public API ------------------------------------------------------------- #

    @classmethod
    def from_dict(cls, data: dict) -> PermutationBootstrapMeanConfig:
        """
        Build from a plain dict, as read from YAML.

        Raises
        ------
        TypeError
            If data is not a dict.
        ValueError
            If it holds any key other than ``method``.
        """
        if not isinstance(data, dict):
            raise TypeError(
                f"{_ERROR_PREFIX}: expected a dict, got {type(data).__name__}"
            )
        unknown = sorted(set(data) - {"method"})
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
