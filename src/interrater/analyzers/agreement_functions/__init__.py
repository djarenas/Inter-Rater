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
    cohen_kappa_max,
    max_observed_agreement,
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
    "cohen_kappa_max",
    "max_observed_agreement",
    "chance_conger",
    "cohen_kappa",
    "conger_kappa",
    "pairable_ratings",
    "observed_disagreement",
    "expected_disagreement",
    "krippendorff_alpha",
]
