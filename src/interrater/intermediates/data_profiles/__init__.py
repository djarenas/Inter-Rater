"""Data structure profiles: balance, overlap, and other data metadata.

Profiles describe the structure and properties of ratings data. They exist
because visualization and interpretation depend on understanding data shape,
not just coefficient values.
"""

from .agreement_profile import AgreementProfile
from .balance_profile import BalanceProfile
from .prevalence_bias_index import PrevalenceBiasIndex
from .rater_shared_matrix import RaterSharedMatrix

__all__ = [
    "AgreementProfile",
    "BalanceProfile",
    "PrevalenceBiasIndex",
    "RaterSharedMatrix",
]
