"""Configuration objects for inter-rater agreement analyzers.

Every analyzer is configured via a frozen dataclass that declares all
analytical decisions upfront: bootstrap resamples, confidence level, CI method,
random seed, and optional convergence diagnostics. Configurations load from
YAML and are fully auditable.

Hierarchy:
  AnalyzerConfig (base)
    ├── GroupCoefficientConfig
    │   ├── FleissKappaConfig
    │   ├── ScottsConfig
    │   ├── GwetAC1Config
    │   ├── BrennanPredigerConfig
    │   └── KrippendorffAlphaConfig
    │
    ├── IndividualMarginalConfig
    │   ├── CohenKappaConfig
    │   ├── CongerKappaConfig
    │   └── LightKappaConfig
    │
    └── PermutationConfig (abstract; holds nested AnalyzerConfig)
        ├── PairwisePermutationConfig
        └── SubsetPermutationConfig

Supporting types:
  HypothesisTestConfig (hypothesis testing against H₀: coefficient = 0)
  PermutationBootstrapMeanConfig (mean-level bootstrap for permutation)
"""

from .analyzer_config import AnalyzerConfig
from .group_coefficient_config import GroupCoefficientConfig
from .fleiss_kappa_config import FleissKappaConfig
from .scotts_config import ScottsConfig
from .individual_marginal_config import IndividualMarginalConfig
from .gwet_ac1_config import GwetAC1Config
from .brennan_prediger_config import BrennanPredigerConfig
from .krippendorff_alpha_config import KrippendorffAlphaConfig
from .cohen_kappa_config import CohenKappaConfig
from .conger_kappa_config import CongerKappaConfig
from .light_kappa_config import LightKappaConfig
from .hypothesis_test_config import HypothesisTestConfig
from .permutation_config import PermutationConfig
from .permutation_bootstrap_mean_config import PermutationBootstrapMeanConfig
from .pairwise_permutation_config import PairwisePermutationConfig
from .subset_permutation_config import SubsetPermutationConfig

__all__ = [
    "AnalyzerConfig",
    "GroupCoefficientConfig",
    "FleissKappaConfig",
    "ScottsConfig",
    "GwetAC1Config",
    "BrennanPredigerConfig",
    "KrippendorffAlphaConfig",
    "IndividualMarginalConfig",
    "CohenKappaConfig",
    "CongerKappaConfig",
    "LightKappaConfig",
    "HypothesisTestConfig",
    "PermutationConfig",
    "PermutationBootstrapMeanConfig",
    "PairwisePermutationConfig",
    "SubsetPermutationConfig",
]
