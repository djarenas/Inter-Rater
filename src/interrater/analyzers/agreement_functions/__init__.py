"""
Coefficient arithmetic. Called by analysers, not by users directly.

Each function takes an ``EstimatorInput`` and returns a float. None of them
validate: whether a coefficient applies to the data is the analyser's question,
asked once before any of this runs.
"""

from .group_agreement import (
    agreement_corrected_by_chance,
    brennan_prediger,
    category_proportions,
    chance_brennan_prediger,
    chance_fleiss,
    chance_gwet_ac1,
    fleiss_kappa,
    gwet_ac1,
    observed_agreement,
    percent_agreement,
    scotts_pi,
)

from .individual_marginals import (
    chance_cohen,
    chance_conger,
    cohen_kappa,
    conger_kappa,
    rater_marginals,
    shared_rows,
)
from .krippendorff import (
    expected_disagreement,
    krippendorff_alpha,
    observed_disagreement,
    pairable_ratings,
)
from .rater_combinations import (
    leave_one_rater_out,
    mean_over_same_subsets,
    mean_over_subsets,
    over_pairs,
    over_subsets,
    rater_subsets,
    subjects_per_leave_one_out,
    subjects_per_subset,
)

__all__ = [
    "observed_agreement",
    "category_proportions",
    "agreement_corrected_by_chance",
    "chance_fleiss",
    "chance_gwet_ac1",
    "chance_brennan_prediger",
    "percent_agreement",
    "fleiss_kappa",
    "scotts_pi",
    "gwet_ac1",
    "brennan_prediger",
    "rater_marginals",
    "shared_rows",
    "chance_cohen",
    "chance_conger",
    "cohen_kappa",
    "conger_kappa",
    "rater_subsets",
    "subjects_per_subset",
    "over_subsets",
    "over_pairs",
    "leave_one_rater_out",
    "subjects_per_leave_one_out",
    "mean_over_subsets",
    "mean_over_same_subsets",
    "pairable_ratings",
    "observed_disagreement",
    "expected_disagreement",
    "krippendorff_alpha",
]
