"""Configuration for Scott's pi analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .no_weight_group_config import NoWeightGroupConfig

if TYPE_CHECKING:
    from ...injectors import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[ScottsConfig] Error"


@dataclass(frozen=True)
class ScottsConfig(NoWeightGroupConfig):
    """
    I am the configuration for Scott's pi analysis.

    Scott's pi (Scott 1955) measures agreement between exactly two raters on
    nominal (unweighted) categorical data, corrected for chance using pooled
    marginal proportions.

    **Why this class exists:**

    Scott's pi is historically important and appears in many systematic reviews.
    It is defined only for two raters. Fleiss (1971) generalized it to many
    raters by keeping the same observed agreement and chance formulas. At two
    raters, Fleiss' kappa *is* Scott's pi -- they are algebraically identical.

    This class exists to:
    - Honor Scott's priority and the many papers citing him
    - Provide a clear, two-rater-only entry point (no confusion about whether
      to use Scott or Fleiss for many-rater data)
    - Document the nominal-data-only constraint explicitly

    For many raters, use FleissKappaConfig instead (it reduces to Scott's pi
    at n=2 and is the established generalization).

    I am unweighted. The analyzer refuses ratings whose semantics declares
    ordered or weighted categories, with a ValueError.

    Because Scott's pi is algebraically identical to Fleiss' kappa at two
    raters, the parametric null variance formula (Fleiss 1979, Eq. 12) is
    available here -- ci_method='parametric_null' is allowed.

    References
    ----------
    Scott, W. A. (1955). Reliability of content analysis: The case of
        nominal scale coding. Public Opinion Quarterly, 19(3), 321--325.
    Fleiss, J. L. (1971). Measuring nominal scale agreement among many
        raters. Psychological Bulletin, 76(5), 378--382.
    Gwet, K. L. (2014). Handbook of Inter-Rater Reliability (4th ed.).
        Advanced Analytics, LLC.

    Attributes
    ----------
    (all fields inherited from NoWeightGroupConfig and parents)
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
    # Inherited from parent class NoWeightGroupConfig (validated to stay at defaults)

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Scott's pi configuration."""

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Scott's pi equals Fleiss' kappa at two raters, so it inherits the
        # same closed-form null variance and the same caveat about using it
        # for an interval. Called last, and from here rather than from the
        # base class -- see AnalyzerConfig._warn_parametric_null.
        if self.ci_method == "parametric_null":
            self._warn_parametric_null()

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration (inherits parent's summary)."""
        return super().summary()
