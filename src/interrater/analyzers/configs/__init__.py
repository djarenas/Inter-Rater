"""Configuration objects for inter-rater agreement analyzers.

Every analyzer is configured via a frozen dataclass that declares all
analytical decisions upfront: bootstrap resamples, confidence level, CI method,
random seed, and optional convergence diagnostics. Configurations load from
YAML and are fully auditable.

Hierarchy:
  AnalyzerConfig (base)
    ├── GroupCoefficientConfig
    │   ├── NoWeightGroupConfig (AC1, Scott)
    │   │   ├── GwetAC1Config
    │   │   └── ScottsConfig
    │   │
    │   └── WeightableGroupConfig (Fleiss, BP, Krippendorff)
    │       ├── FleissKappaConfig
    │       ├── BrennanPredigerConfig
    │       └── KrippendorffAlphaConfig
    │
    ├── IndividualMarginalConfig
    │   ├── NoWeightIndividualConfig (Conger, Light)
    │   │   ├── CongerKappaConfig
    │   │   └── LightKappaConfig
    │   │
    │   └── WeightableIndividualConfig (Cohen)
    │       └── CohenKappaConfig
    │
    └── PermutationConfig (holds nested AnalyzerConfig)
        ├── PairwisePermutationConfig
        └── SubsetPermutationConfig

Supporting types:
  HypothesisTestConfig (hypothesis testing against H₀: coefficient = 0)
  PermutationBootstrapMeanConfig (mean-level bootstrap for permutation)
"""

from .analyzer_config import AnalyzerConfig
from .group_coefficient_config import GroupCoefficientConfig
from .no_weight_group_config import NoWeightGroupConfig
from .weightable_group_config import WeightableGroupConfig
from .fleiss_kappa_config import FleissKappaConfig
from .scotts_config import ScottsConfig
from .individual_marginal_config import IndividualMarginalConfig
from .no_weight_individual_config import NoWeightIndividualConfig
from .gwet_ac1_config import GwetAC1Config
from .brennan_prediger_config import BrennanPredigerConfig
from .krippendorff_alpha_config import KrippendorffAlphaConfig
from .weightable_individual_config import WeightableIndividualConfig
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
    "WeightableGroupConfig",
    "NoWeightGroupConfig",
    "FleissKappaConfig",
    "ScottsConfig",
    "GwetAC1Config",
    "BrennanPredigerConfig",
    "KrippendorffAlphaConfig",
    "IndividualMarginalConfig",
    "NoWeightIndividualConfig",
    "WeightableIndividualConfig",
    "CohenKappaConfig",
    "CongerKappaConfig",
    "LightKappaConfig",
    "HypothesisTestConfig",
    "PermutationConfig",
    "PermutationBootstrapMeanConfig",
    "PairwisePermutationConfig",
    "SubsetPermutationConfig",
]
