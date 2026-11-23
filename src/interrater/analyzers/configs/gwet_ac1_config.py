"""Configuration for Gwet's AC1 analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .group_coefficient_config import GroupCoefficientConfig

if TYPE_CHECKING:
    from .hypothesis_test_config import HypothesisTestConfig
    from ..convergence_test import ConvergenceTest
    from ...injectors import InjectorConfig

_ERROR_PREFIX = "[GwetAC1Config] Error"


@dataclass(frozen=True)
class GwetAC1Config(GroupCoefficientConfig):
    """
    I am the configuration for Gwet's AC1 analysis.

    Gwet's AC1 is a chance-corrected agreement coefficient for categorical
    data with two or more raters. It addresses the low values of Cohen's and
    Fleiss' kappa under skewed prevalence. Its chance term (Gwet 2008,
    Equation 38) is:

        Pe = sum_j p_j (1 - p_j) / (K - 1)

    When prevalence is skewed (one category dominates), Pe is small, so AC1
    stays close to Po.

    With ordered or weighted categories in the semantics, the analyzer
    computes AC2, the weighted AC1. Weights are not configured here.

    I inherit all fields from GroupCoefficientConfig and add no new ones.

    References
    ----------
    - Gwet, K. L. (2008). Computing inter-rater reliability and its variance
      in the presence of high agreement. British Journal of Mathematical and
      Statistical Psychology, 61(1), 29--48.
    - Gwet, K. L. (2021). Handbook of Inter-Rater Reliability (5th ed.),
      Vol. 1: Analysis of Categorical Ratings. AgreeStat Analytics.
      Chapter 4, AC2.

    Attributes
    ----------
    balanced_data_only : bool, default False
        AC1 works on both balanced and unbalanced data. The formula
        naturally generalizes to varying raters per subject.
    (all other fields inherited from GroupCoefficientConfig and
    AnalyzerConfig)
    """

    # Attribute declaration

    # Inherited from parent class AnalyzerConfig
    n_resamples: int
    alpha: float
    random_state: int
    ci_method: str
    convergence_test: ConvergenceTest | None
    hypothesis_test_config: HypothesisTestConfig | None
    injection: tuple[InjectorConfig, ...] | None
    # Inherited from parent class GroupCoefficientConfig
    balanced_data_only: bool
    rater_jackknife: bool

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Gwet AC1 configuration."""

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # No parametric null is implemented for AC1, so refuse it
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method='parametric_null' is not "
                f"available for Gwet's AC1. The package does not implement a "
                f"parametric null for Gwet's AC1. Use 'percentile' or 'bca' "
                f"instead."
            )

        # Also reject parametric_null as hypothesis test method
        if (
            self.hypothesis_test_config is not None
            and self.hypothesis_test_config.method == "parametric_null"
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis test method='parametric_null' is not "
                f"available for Gwet's AC1. The package does not implement a "
                f"parametric null for it. Use method='permutation' instead."
            )
