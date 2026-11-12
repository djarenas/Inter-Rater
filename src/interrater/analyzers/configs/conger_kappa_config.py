"""Configuration for Conger's kappa analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .no_weight_individual_config import NoWeightIndividualConfig

if TYPE_CHECKING:
    from ...injectors.configs import InjectorConfig
    from ..convergence_test import ConvergenceTest
    from .hypothesis_test_config import HypothesisTestConfig

_ERROR_PREFIX = "[CongerKappaConfig] Error"

@dataclass(frozen=True)
class CongerKappaConfig(NoWeightIndividualConfig):
    """
    I am the configuration for Conger's kappa analysis.

    Conger's kappa (Conger 1980) generalises Cohen's kappa to any number of
    raters (n >= 2) while preserving individual marginals. Where Fleiss' kappa
    pools marginal distributions across raters, Conger's kappa keeps each
    rater's own category distribution in the chance term -- the same principle
    as Cohen's, extended beyond a single pair.

    This makes Conger's kappa the multi-rater member of the individual-marginal
    family. At exactly two raters, Conger's kappa reduces to Cohen's unweighted
    kappa.

    I do not support weighting. Conger (1980) defined the coefficient for
    nominal data only; there is no standard weighted generalisation. The
    analyzer refuses ratings whose semantics declares ordered or weighted
    categories.

    There is no known closed-form null variance for Conger's kappa, so
    ci_method="parametric_null" is not available. Use bootstrap
    ("percentile" or "bca") for confidence intervals.


    References
    ----------
    Conger, A. J. (1980). Integration and generalization of kappas for
        multiple raters. Psychological Bulletin, 88(2), 322--328.

    Attributes
    ----------
    balanced_data_only : bool, default True
        Must be True. The field exists so that every coefficient config
        accepts the same keys, and a config moved here from Fleiss' kappa
        is told why False does not apply instead of failing on an unknown
        key. Conger's kappa has no approximation for missing ratings: the
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
        """Validate Conger's kappa configuration."""

        # Call parent validation
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        if not isinstance (self.balanced_data_only, bool):
            raise TypeError(f"{_ERROR_PREFIX}: balanced_data_only must be bool")

        if not self.balanced_data_only:
            raise ValueError(
                f"{_ERROR_PREFIX}: balanced_data_only: false is not available "
                f"for Conger's kappa. It has no approximation for missing "
                f"ratings: every rater must rate every subject (Conger 1980). "
                f"For data with missing ratings use Fleiss' kappa with "
                f"balanced_data_only: false, Gwet's AC1, Brennan-Prediger or "
                f"Krippendorff's alpha, or Cohen's kappa pair by pair in a "
                f"pairwise permutation analysis."
            )


        # Conger has no parametric null variance formula, so reject it as ci_method
        if self.ci_method == "parametric_null":
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method='parametric_null' is not available "
                f"for Conger's kappa. There is no closed-form null variance for "
                f"Conger's kappa. Use 'percentile' or 'bca' instead."
            )

        # Also reject parametric_null as hypothesis test method
        if (
            self.hypothesis_test_config is not None
            and self.hypothesis_test_config.method == "parametric_null"
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: hypothesis test method='parametric_null' is not "
                f"available for Conger's kappa. There is no closed-form null "
                f"variance. Use method='permutation' instead."
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration, with balanced_data_only as the group configs report it."""
        result = super().summary()
        result["balanced_data_only"] = self.balanced_data_only
        return result
