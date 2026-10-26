"""Configuration for bootstrapping the mean across rater combinations."""

from __future__ import annotations

from dataclasses import dataclass

_ERROR_PREFIX = "[PermutationBootstrapMeanConfig] Error"


@dataclass(frozen=True)
class PermutationBootstrapMeanConfig:
    """
    I am the configuration for bootstrapping the mean coefficient across
    rater combinations (pairs or subsets).

    I am bundled as one object so that ``None`` on the parent
    PermutationConfig means "no mean-level bootstrap." When I am present,
    the permutation analyzer will bootstrap the mean and produce CIs.

    Two bootstrap procedures are supported:

    - "correlated": resample subjects once per draw, compute all
      combinations on the same resampled data, then average. Preserves
      the correlation structure between combinations that share a rater.
    - "independent": each combination is bootstrapped separately (via
      the nested config_of_coefficient's bootstrap settings, which
      always run), and the mean is assembled from those independent
      draws. Needs at least 2 draws per combination, so the parent
      config refuses it when config_of_coefficient.n_resamples < 2.

    For the "independent" method, reproducible per-combination seeds
    are derived from ``random_state`` using ``np.random.SeedSequence``.

    Attributes
    ----------
    method : str
        Bootstrap procedure: "correlated" or "independent".
    n_bootstraps : int
        Number of bootstrap draws.
    ci_method : str
        How CIs are computed: "percentile" or "bca". Default
        "percentile". Permutation results have no closed-form null
        variance, so "parametric_null" is not available. "bca" needs
        method "correlated" (see the check in __post_init__).
    alpha : float
        Significance level for CIs. Default 0.05 (95% CI).
    random_state : int | None
        Seed for reproducibility. None for unseeded. Default None.
    """

    # Attribute declaration

    # Required (no defaults)
    method: str
    n_bootstraps: int

    # With defaults
    ci_method: str = "percentile"
    alpha: float = 0.05
    random_state: int | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate permutation bootstrap mean configuration."""

        # Validate method
        if not isinstance(self.method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: method must be str, "
                f"got {type(self.method).__name__}"
            )
        if self.method not in {"correlated", "independent"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: method must be 'correlated' or "
                f"'independent', got '{self.method}'"
            )

        # Guard against the bool backwards compatibility with int
        for name in ("alpha", "random_state", "n_bootstraps"):
            if isinstance(getattr(self, name), bool):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must not be bool, "
                    f"got {type(getattr(self, name)).__name__}"
                )

        # Validate n_bootstraps
        if not isinstance(self.n_bootstraps, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_bootstraps must be int, "
                f"got {type(self.n_bootstraps).__name__}"
            )
        if self.n_bootstraps < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_bootstraps must be >= 1, "
                f"got {self.n_bootstraps}"
            )

        # Validate ci_method
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method not in {"percentile", "bca"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be 'percentile' or "
                f"'bca', got '{self.ci_method}'. Permutation results "
                f"have no closed-form null variance, so "
                f"'parametric_null' is not available."
            )

        # BCa's acceleration comes from a jackknife that must leave out the
        # same units the bootstrap resamples. The independent mean resamples
        # each combination separately and mixes the draws, so no single
        # leave-one-out matches it.
        if self.method == "independent" and self.ci_method == "bca":
            raise ValueError(
                f"{_ERROR_PREFIX}: method 'independent' cannot use "
                f"ci_method='bca'. BCa's acceleration needs a jackknife that "
                f"leaves out the units the bootstrap resamples; the "
                f"independent mean resamples each combination separately "
                f"and mixes the draws, so no jackknife matches it. Use "
                f"ci_method='percentile', or method='correlated'."
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

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """
        Return configuration in plain Python types.

        Returns
        -------
        dict
            Keys: method, n_bootstraps, ci_method, alpha, random_state.
        """
        return {
            "method": self.method,
            "n_bootstraps": self.n_bootstraps,
            "ci_method": self.ci_method,
            "alpha": self.alpha,
            "random_state": self.random_state,
        }
