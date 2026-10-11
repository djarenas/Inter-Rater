"""Data structure profiles: balance, overlap, and other data metadata.

Profiles describe the structure and properties of ratings data. They exist
because visualization and interpretation depend on understanding data shape,
not just coefficient values.
"""

from .balance_profile import BalanceProfile
from .rater_shared_matrix import RaterSharedMatrix

__all__ = ["BalanceProfile", "RaterSharedMatrix"]
