"""Configuration for Krippendorff's alpha analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .weightable_group_config import WeightableGroupConfig

if TYPE_CHECKING:
    from ...injectors import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[KrippendorffAlphaConfig] Error"


@dataclass(frozen=True)
class KrippendorffAlphaConfig(WeightableGroupConfig):
    """
    I am the configuration for Krippendorff's alpha analysis.

    Krippendorff's alpha measures disagreement rather than agreement:

        alpha = 1 - Do / De

    where Do is the observed disagreement pooled across all pairable rating
    pairs, and De is the expected disagreement under random permutation of
    the same ratings. This is not the same expression as the kappa family's
    ``(Po - Pe) / (1 - Pe)`` rearranged -- the two use different units
    (pairs pooled vs. subjects averaged) and different chance models
    (finite-sample permutation vs. population rates).

    Alpha handles missing data natively. There is no assumption that every
    subject was rated the same number of times or by the same raters. A
    subject rated once contributes no pairs and therefore contributes
    nothing -- it is not dropped by a rule, it simply has nothing to offer.
    That is why alpha is the coefficient usually reached for when the data
    is ragged.

    The weighted form uses a cost function d_ab = 1 - w_ab entering both
    the observed and expected disagreement terms identically. Identity
    weights give the nominal (unweighted) case. The weight matrix encodes
    how much credit near-misses receive, exactly as in the Fleiss and B-P
    weighted forms.

    There is no parametric null variance formula for Krippendorff's alpha.
    Krippendorff recommends bootstrapping for inference (Hayes &
    Krippendorff 2007). ``ci_method='parametric_null'`` and hypothesis
    test ``method='parametric_null'`` are both rejected. Use bootstrap CIs
    (percentile or BCa) and permutation tests instead.

    References
    ----------
    - Krippendorff, K. (2004). Content Analysis: An Introduction to Its
      Methodology (2nd ed.), chapter 11. Sage.
    - Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
      Departmental Papers, Annenberg School for Communication, University
      of Pennsylvania.
    - Hayes, A. F., & Krippendorff, K. (2007). Answering the call for a
      standard reliability measure for coding data. Communication Methods
      and Measures, 1(1), 77--89.

    Attributes
    ----------
    balanced_data_only : bool, default False
        Krippendorff's alpha works on both balanced and unbalanced data.
        Unbalanced data is its natural habitat -- subjects with varying
        numbers of raters are handled by pair-pooling rather than
        subject-averaging.
    weights_type : str, default "identity"
        Supports "identity" (unweighted/nominal), "linear", "quadratic",
        "custom". The cost of confusing categories a and b is 1 - w_ab.
    (all other fields inherited from WeightableGroupConfig and parents)
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
    # Inherited from parent class WeightableGroupConfig
    weights_type: str = "identity"
    custom_weights_path: str | None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Krippendorff's alpha configuration."""
        # Call parent validation (WeightableGroupConfig validates weights)
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Alpha has no parametric null variance formula, so reject it as ci_method
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method='parametric_null' is not available "
                f"for Krippendorff's alpha. There is no closed-form null "
                f"variance formula for this coefficient. Krippendorff "
                f"recommends bootstrapping (Hayes & Krippendorff 2007). "
                f"Use 'percentile' or 'bca' instead."
            )

        # Also reject parametric_null as hypothesis test method
        if (
            self.hypothesis_test_config is not None
            and self.hypothesis_test_config.method == "parametric_null"
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis test method='parametric_null' is not "
                f"available for Krippendorff's alpha. There is no closed-form "
                f"null variance for this coefficient. Use method='permutation' "
                f"instead."
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration (inherits parent's summary)."""
        return super().summary()
