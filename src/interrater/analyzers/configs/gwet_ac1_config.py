"""Configuration for Gwet's AC1 analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .no_weight_group_config import NoWeightGroupConfig

if TYPE_CHECKING:
    from .hypothesis_test_config import HypothesisTestConfig
    from ..convergence_test import ConvergenceTest
    from ...injectors import InjectorConfig

_ERROR_PREFIX = "[GwetAC1Config] Error"


@dataclass(frozen=True)
class GwetAC1Config(NoWeightGroupConfig):
    """
    I am the configuration for Gwet's AC1 analysis.

    Gwet's AC1 is an unweighted, chance-corrected agreement coefficient for
    nominal data with any number of raters (n >= 2). It addresses the "kappa
    paradox" -- the counter-intuitive collapse of Cohen/Fleiss kappa under
    extreme prevalence -- by replacing kappa's marginal-based chance term with
    a uniform-distribution-based term that reflects the information content
    of the category set.

    AC1 uses a chance-agreement probability:

        Pe(gamma) = Tw / (K - 1)

    where Tw = 1 - sum(pj^2) measures category heterogeneity. When prevalence
    is skewed (one category dominates), Tw is small, Pe(gamma) is small, and
    AC1 stays high -- matching the intuition that raters who mostly agree
    really do agree, even when categories are imbalanced.

    For weighted (ordinal) analysis, use AC2 (GwetAC2Config) instead.

    I inherit all fields from NoWeightGroupConfig and add no new ones.
    The unweighted constraint is enforced by the analyzer, which refuses
    ratings whose semantics declares ordered or weighted categories.

    References
    ----------
    - Gwet, K. L. (2008). Computing inter-rater reliability and its variance
      in the presence of high agreement. British Journal of Mathematical and
      Statistical Psychology, 61(1), 29--48.
    - Gwet, K. L. (2014). Handbook of Inter-Rater Reliability (4th ed.).
      Advanced Analytics, LLC.
    - Wongpakaran, N., et al. (2013). A comparison of Cohen's Kappa and
      Gwet's AC1 when calculating inter-rater reliability coefficients:
      a study conducted with personality disorder samples.
      BMC Medical Research Methodology, 13, 61.

    Attributes
    ----------
    balanced_data_only : bool, default False
        AC1 works on both balanced and unbalanced data. The formula
        naturally generalizes to varying raters per subject.
    (all other fields inherited from GroupCoefficientConfig and AnalyzerConfig)
    """

    # Attribute declaration

    # Inherited from parent class AnalyzerConfig
    n_resamples: int 
    alpha: float 
    random_state: int | None
    ci_method: str
    convergence_test: ConvergenceTest | None
    hypothesis_test_config: HypothesisTestConfig | None
    injection: tuple[InjectorConfig] | None
    # Inherited from parent class GroupCoefficientConfig
    balanced_data_only: bool
    rater_jackknife: bool
    # Inherited from parent class NoWeightGroupConfig (validated there to stay at defaults)

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Gwet AC1 configuration."""

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # AC1 has no parametric null variance formula, so reject it as ci_method
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method='parametric_null' is not available "
                f"for Gwet's AC1. The parametric null variance formula is specific "
                f"to Fleiss' kappa (Fleiss 1979, Eq. 12). "
                f"Use 'percentile' or 'bca' instead."
            )

        # Also reject parametric_null as hypothesis test method
        if (
            self.hypothesis_test_config is not None
            and self.hypothesis_test_config.method == "parametric_null"
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis test method='parametric_null' is not "
                f"available for Gwet's AC1. There is no closed-form null variance "
                f"for AC1. Use method='permutation' instead."
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration (inherits parent's summary)."""
        return super().summary()
