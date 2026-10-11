"""Configuration for hypothesis testing."""

from __future__ import annotations

from dataclasses import dataclass

_ERROR_PREFIX = "[HypothesisTestConfig] Error"

_VALID_METHODS = {"permutation", "parametric_null"}
_VALID_SCHEMES = {"column_shuffle", "within_subject", "pool_and_redistribute"}


@dataclass(frozen=True)
class HypothesisTestConfig:
    """
    I am the configuration for hypothesis testing against H₀: coefficient = 0.

    I define the method, permutation scheme, and parameters for hypothesis
    testing. I am stored as an optional attribute on AnalyzerConfig and read
    by each leaf analyzer's hypothesis_test() method.

    Two methods are supported:
    - "permutation": build null distribution by shuffling ratings. Requires
      Ratings data (not NijMatrix). The scheme parameter controls how
      ratings are shuffled.
    - "parametric_null": use an analytical null variance formula. Coefficient-
      specific (e.g., Fleiss 1979 Eq. 12 for Fleiss' kappa). Works with
      both Ratings and NijMatrix.

    Three permutation schemes are supported (only used when method="permutation"):
    - "column_shuffle": shuffle each rater's column independently. Preserves
      each rater's marginal distribution. Tests: "Do raters agree more than
      expected given their individual biases?"
    - "within_subject": shuffle category assignments within each subject's row.
      Preserves the number of raters per subject. Tests: "Is the pattern
      within subjects non-random?"
    - "pool_and_redistribute": pool all ratings and redistribute randomly.
      Preserves global category proportions only. Tests: "Is there ANY
      structure in the data?"

    Attributes
    ----------
    method : str
        The hypothesis testing method ("permutation" or "parametric_null").
    scheme : str
        The permutation scheme (only used when method="permutation").
    n_permutations : int
        Number of permutations to generate (only used when method="permutation").
    alpha : float
        The significance level for the test.
    random_state : int
        Random seed for reproducibility of permutation tests.
    """

    method: str = "permutation"
    scheme: str = "column_shuffle"
    n_permutations: int = 5000
    alpha: float = 0.05
    random_state: int = 42

    def __post_init__(self) -> None:
        """Validate hypothesis test configuration."""
        for name in ("method", "scheme"):
            if not isinstance(getattr(self, name), str):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must be a string, "
                    f"got {type(getattr(self, name)).__name__}"
                )

        if self.method not in _VALID_METHODS:
            raise ValueError(
                f"{_ERROR_PREFIX}: method must be one of {_VALID_METHODS}, "
                f"got '{self.method}'."
            )

        # Validate scheme (only meaningful for permutation, but always validated)
        if self.scheme not in _VALID_SCHEMES:
            raise ValueError(
                f"{_ERROR_PREFIX}: scheme must be one of {_VALID_SCHEMES}, "
                f"got '{self.scheme}'."
            )

        # Check the bool python bug. After around Python2.3, bools were made to inherit
        # int for backwards compatibility 
        for name in ("n_permutations", "random_state", "alpha"):
            if isinstance(getattr(self, name), bool):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} cannot be bool, "
                    f"got {type(getattr(self, name)).__name__}"
                )

        # Validate number of permutations
        if not isinstance(self.n_permutations, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: n_permutations must be int, "
                f"got {type(self.n_permutations).__name__}"
            )
        if self.n_permutations < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_permutations must be >= 1, "
                f"got {self.n_permutations}"
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
        if not isinstance(self.random_state, int):
            raise TypeError(
                f"{_ERROR_PREFIX}: random_state must be an integer, "
                f"got {type(self.random_state).__name__}."
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration for audit trail."""
        return {
            "method": self.method,
            "scheme": self.scheme,
            "n_permutations": self.n_permutations,
            "alpha": self.alpha,
            "random_state": self.random_state,
        }
