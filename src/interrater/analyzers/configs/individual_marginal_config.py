"""Configuration for individual marginal coefficients (Cohen, Conger, Light)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .analyzer_config import AnalyzerConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from .convergence_test_config import ConvergenceTestConfig
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[IndividualMarginalConfig] Error"


@dataclass(frozen=True)
class IndividualMarginalConfig(AnalyzerConfig):
    """
    I am the base configuration for coefficients that keep rater identity.

    Individual marginal coefficients (Cohen's kappa, Conger's kappa, Light's
    kappa) maintain which rater said what, enabling per-rater marginal
    probabilities.

    I do not allow rater_jackknife (leave-one-rater-out diagnostics): it is
    not supported in this version yet, for any coefficient.
    SubsetPermutationAnalysis with subset_k: n_raters - 1 gives
    leave-one-rater-out today.

    Attributes
    ----------
    rater_jackknife : bool, default False
        Must be False: leave-one-rater-out diagnostics are not supported in
        this version yet. SubsetPermutationAnalysis with subset_k:
        n_raters - 1 gives the same numbers meanwhile.
    (all other fields inherited from AnalyzerConfig)
    """

    # Attribute declaration

    # Inherited from parent class AnalyzerConfig
    n_resamples: int
    alpha: float
    random_state: int
    ci_method: str
    convergence_test_config: ConvergenceTestConfig | None
    hypothesis_test_config: HypothesisTestConfig | None
    injection: tuple[InjectorConfig, ...] | None
    # Populated at this level (although not supported yet)
    rater_jackknife: bool = False

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate individual marginal configuration."""

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate rater_jackknife
        if not isinstance(self.rater_jackknife, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: rater_jackknife must be bool, "
                f"got {type(self.rater_jackknife).__name__}"
            )

        # Rater jackknife not currently supported for individual marginals
        if self.rater_jackknife:
            raise ValueError(
                f"{_ERROR_PREFIX}: rater_jackknife (leave-one-rater-out) is not "
                f"supported in this version yet. Set rater_jackknife: false. "
                f"SubsetPermutationAnalysis with subset_k: n_raters - 1 gives "
                f"the coefficient without each rater in turn."
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration including rater_jackknife."""
        result = super().summary()
        result.update({
            "rater_jackknife": self.rater_jackknife,
        })
        return result
