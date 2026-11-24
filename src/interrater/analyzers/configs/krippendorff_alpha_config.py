"""Configuration for Krippendorff's alpha analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .group_coefficient_config import GroupCoefficientConfig

if TYPE_CHECKING:
    from ...injectors import InjectorConfig
    from .convergence_test_config import ConvergenceTestConfig
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[KrippendorffAlphaConfig] Error"


@dataclass(frozen=True)
class KrippendorffAlphaConfig(GroupCoefficientConfig):
    """
    I am the configuration for Krippendorff's alpha analysis.

    Krippendorff's alpha measures disagreement rather than agreement
    (Krippendorff 2011):

        alpha = 1 - Do / De

    where Do is the observed disagreement pooled across all pairable rating
    pairs, and De is the expected disagreement under random permutation of
    the same ratings. Written with Po = 1 - Do and Pe = 1 - De, this is the
    kappa family's ``(Po - Pe) / (1 - Pe)``, but the two compute Po and Pe
    differently.

    Alpha handles missing data natively. There is no assumption that every
    subject was rated the same number of times or by the same raters. A
    subject rated once contributes no pairs and therefore contributes
    nothing.

    The weighted form uses a cost function d_ab = 1 - w_ab entering both
    the observed and expected disagreement terms identically. Identity
    weights give the nominal (unweighted) case. The weight matrix encodes
    how much credit near-misses receive, exactly as in the Fleiss and B-P
    weighted forms.

    The package does not implement a parametric null for Krippendorff's alpha.
    ``ci_method='parametric_null'`` and hypothesis test
    ``method='parametric_null'`` are both rejected. Use bootstrap CIs
    (percentile or BCa) and permutation tests instead.

    References
    ----------
    - Krippendorff, K. (2011). Computing Krippendorff's alpha-reliability.
      Departmental Papers, Annenberg School for Communication, University
      of Pennsylvania.

    Attributes
    ----------
    balanced_data_only : bool, default False
        Krippendorff's alpha works on both balanced and unbalanced data.
    (all other fields inherited from GroupCoefficientConfig and parents)
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
    # Inherited from parent class GroupCoefficientConfig
    balanced_data_only: bool
    rater_jackknife: bool

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Krippendorff's alpha configuration."""
        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # No parametric null is implemented for alpha, so refuse it
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method='parametric_null' is not available "
                f"for Krippendorff's alpha. The package does not implement a "
                f"parametric null for it. Use 'percentile' or 'bca' instead."
            )

        # Also reject parametric_null as hypothesis test method
        if (
            self.hypothesis_test_config is not None
            and self.hypothesis_test_config.method == "parametric_null"
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis test method='parametric_null' is not "
                f"available for Krippendorff's alpha. The package implements "
                f"no parametric null for it. Use method='permutation' instead."
            )
