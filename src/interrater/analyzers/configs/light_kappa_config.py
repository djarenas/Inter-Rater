"""Configuration for Light's kappa analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .no_weight_individual_config import NoWeightIndividualConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig


_ERROR_PREFIX = "[LightKappaConfig] Error"


@dataclass(frozen=True)
class LightKappaConfig(NoWeightIndividualConfig):
    """
    I am the configuration for Light's kappa analysis.

    Light's kappa (Light 1971) is the unweighted mean of all C(n, 2) pairwise
    Cohen's kappas. It summarises how well raters agree on average across every
    possible pair, without pooling marginals as Fleiss' kappa does.

    Light's kappa keeps rater identity through Cohen's kappa: each pair uses
    its own marginals. The mean then compresses the C(n, 2) values into a
    single summary. For per-pair detail, spread, and per-rater diagnostics,
    use the pairwise permutation analyzer -- Light's kappa is the headline
    number, not the diagnostic.

    At exactly two raters, Light's kappa reduces to Cohen's unweighted kappa
    (there is only one pair).

    I do not support weighting. Light (1971) defined the coefficient as the
    mean of unweighted pairwise kappas. The analyzer refuses ratings whose
    semantics declares ordered or weighted categories.

    There is no closed-form null variance for Light's kappa. The C(n, 2)
    pairwise kappas are correlated (pairs sharing a rater are not independent),
    so the variance of the mean cannot be derived from the individual null
    variances. Use bootstrap ("percentile" or "bca") for confidence intervals.

    References
    ----------
    Light, R. J. (1971). Measures of response agreement for qualitative
        data: Some generalizations and alternatives. Psychological
        Bulletin, 76(5), 365--377.

    Attributes
    ----------
    balanced_data_only : bool, default True
        Must be True. The field exists so that every coefficient config
        accepts the same keys, and a config moved here from Fleiss' kappa
        is told why False does not apply instead of failing on an unknown
        key. Light (1971) defines the coefficient for complete data: the
        analyzer requires every rater to rate every subject, which is
        stricter than balanced.
    (all other fields inherited from NoWeightIndividualConfig and parents)
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
    # Inherited from parent class IndividualMarginalConfig
    rater_jackknife: bool
    # Inherited from parent class NoWeightIndividualConfig (validated to stay at defaults)
    # Populated at this level
    # Must be True. Exists so that every coefficient config accepts the same
    # keys, and a config moved here is told why False does not apply.
    balanced_data_only: bool = True

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Light's kappa configuration."""

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        if not isinstance(self.balanced_data_only, bool):
            raise TypeError(f"{_ERROR_PREFIX}: balanced_data_only must be bool")

        if not self.balanced_data_only:
            raise ValueError(
                f"{_ERROR_PREFIX}: balanced_data_only: false is not available "
                f"for Light's kappa. It has no approximation for missing "
                f"ratings: every rater must rate every subject (Light 1971). "
                f"For data with missing ratings use Fleiss' kappa with "
                f"balanced_data_only: false, Gwet's AC1, Brennan-Prediger or "
                f"Krippendorff's alpha, or Cohen's kappa pair by pair in a "
                f"pairwise permutation analysis."
            )

        # Light has no parametric null variance formula, so reject it as ci_method
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method='parametric_null' is not available "
                f"for Light's kappa. The C(n, 2) pairwise kappas are correlated "
                f"(pairs sharing a rater are not independent), so no closed-form "
                f"null variance exists. Use 'percentile' or 'bca' instead."
            )

        # Also reject parametric_null as hypothesis test method
        if (
            self.hypothesis_test_config is not None
            and self.hypothesis_test_config.method == "parametric_null"
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis test method='parametric_null' is not "
                f"available for Light's kappa. The pairwise kappas are correlated, "
                f"so no closed-form null variance exists. "
                f"Use method='permutation' instead."
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration (inherits parent's summary)."""
        return super().summary()
