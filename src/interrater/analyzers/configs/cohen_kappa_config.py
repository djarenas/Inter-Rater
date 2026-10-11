"""Configuration for Cohen's kappa analysis."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

from .weightable_individual_config import WeightableIndividualConfig

if TYPE_CHECKING:
    from .hypothesis_test_config import HypothesisTestConfig
    from ..convergence_test import ConvergenceTest
    from ...injectors import InjectorConfig

_ERROR_PREFIX = "[CohenKappaConfig] Error"
_WARNING_PREFIX = "[CohenKappaConfig] Warning"

@dataclass(frozen=True)
class CohenKappaConfig(WeightableIndividualConfig):
    """
    I am the configuration for Cohen's kappa analysis.

    Cohen's kappa (Cohen 1960) measures agreement between exactly two raters
    on categorical data, correcting for chance using each rater's individual
    marginal distribution. Unlike Scott's pi or Fleiss' kappa, which pool
    marginals across raters, Cohen's kappa preserves rater identity -- it
    compares what rater A actually said with what rater B actually said.

    I support weighting for ordinal data (Cohen 1968). With weights_type
    set to "linear" or "quadratic", disagreements between adjacent categories
    are penalised less than disagreements between distant ones. The unweighted
    form (weights_type="identity") treats all disagreements equally.

    Cohen's kappa has a known asymptotic null variance (Fleiss, Cohen &
    Everitt 1969), so ci_method="parametric_null" is available here.

    References
    ----------
    Cohen, J. (1960). A coefficient of agreement for nominal scales.
        Educational and Psychological Measurement, 20(1), 37--46.
    Cohen, J. (1968). Weighted kappa: Nominal scale agreement provision
        for scaled disagreement or partial credit. Psychological Bulletin,
        70(4), 213--220.
    Fleiss, J. L., Cohen, J., & Everitt, B. S. (1969). Large sample
        standard errors of kappa and weighted kappa. Psychological
        Bulletin, 72(5), 323--327.

    Attributes
    ----------
    balanced_data_only: bool, default True
        Must always be True. The variable exists in the class so that if someone feeds a YAML with a true value
        the config takes it and does not complain.
    (all other fields inherited from WeightableIndividualConfig and parents)
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
    # Inherited from parent class WeightableIndividualConfig
    weights_type: str 
    custom_weights_path: str | None
    # Populated at this level: Changing default 
    balanced_data_only: bool = True

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate Cohen's kappa configuration."""

        # Call validation by parent classes
        try:
            super().__post_init__()
        except (ValueError, TypeError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        if not isinstance (self.balanced_data_only, bool):
            raise TypeError(f"{_ERROR_PREFIX}: balanced_data_only must be bool")

        if not self.balanced_data_only:
            # Bypass the frozen dataclass and warn
            object.__setattr__(self, "balanced_data_only", True)
            warnings.warn(
                f"{_WARNING_PREFIX}: Cohen's kappa can only be calculated "
                f"with balanced data. Ignoring balanced_data_only=False.",
                UserWarning,
                stacklevel=2,
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """Return configuration (inherits parent's summary)."""
        return super().summary()
