"""Configuration for Brennan-Prediger analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .group_coefficient_config import GroupCoefficientConfig

if TYPE_CHECKING:
    from ...injectors import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[BrennanPredigerConfig] Error"


@dataclass(frozen=True)
class BrennanPredigerConfig(GroupCoefficientConfig):
    """
    I am the configuration for Brennan-Prediger analysis.

    Brennan and Prediger's coefficient defines chance agreement as a uniform
    guess among K categories, making it independent of the observed category
    proportions.

    The weighted form replaces 1/K with sum(W) / K^2, where W is the weight
    matrix. The chance term remains independent of what the raters actually
    did (the hallmark of Brennan-Prediger) and the weighted coefficient
    is the unweighted one with identity W as a special case. The formula is
    the one in Gwet's irrCAC 1.4 (``bp.coeff.raw``).

    Weights are declared in the semantics of the ratings (categories kind
    'ordered' or 'weighted'). The analyzer reads the weight matrix from there
    at construction time.

    There is no parametric null variance formula for Brennan-Prediger, so
    ``ci_method='parametric_null'`` and hypothesis test
    ``method='parametric_null'`` are both rejected. Use bootstrap CIs
    (percentile or BCa) and permutation tests instead.

    References
    ----------
    - Brennan, R. L., & Prediger, D. J. (1981). Coefficient kappa: Some
      uses, misuses, and alternatives. Educational and Psychological
      Measurement, 41(3), 687--699.
    - Gwet, K. L. (2014). Handbook of Inter-Rater Reliability (4th ed.).
      Advanced Analytics, LLC.
    - Gwet, K. L. (2026). irrCAC: Computing the Extent of Agreement among
      Raters with Chance-Corrected Agreement Coefficient (CAC). R package
      version 1.4.

    Attributes
    ----------
    balanced_data_only : bool, default False
        Brennan-Prediger works on both balanced and unbalanced data.
        The formula naturally generalizes to varying raters per subject.
    (all other fields inherited from GroupCoefficientConfig and parents)
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
        """Validate Brennan-Prediger configuration."""

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # B-P has no parametric null variance formula, so reject it as ci_method
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method='parametric_null' is not available "
                f"for Brennan-Prediger. There is no closed-form null variance for "
                f"this coefficient. "
                f"Use 'percentile' or 'bca' instead."
            )

        # Also reject parametric_null as hypothesis test method
        if (
            self.hypothesis_test_config is not None
            and self.hypothesis_test_config.method == "parametric_null"
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis test method='parametric_null' is not "
                f"available for Brennan-Prediger. There is no closed-form null "
                f"variance for this coefficient. Use method='permutation' instead."
            )
