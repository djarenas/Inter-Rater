"""Base analyzer for individual marginal coefficients."""

from __future__ import annotations

from abc import ABC
from typing import TYPE_CHECKING

import numpy as np

from ..contracts import EMPTY
from .analyzer import Analyzer

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.data_profiles import BalanceProfile, RaterProfile
    from ..intermediates.estimator_parts.estimator_input import EstimatorInput
    from .configs import IndividualMarginalConfig

_ERROR_PREFIX = "[IndividualMarginalAnalyzer] Error"

# Who defined each coefficient for complete data, for the refusal message
_COMPLETE_DATA_SOURCES = {
    "Conger's kappa": "Conger 1980",
    "Light's kappa": "Light 1971",
}


def check_complete(category_index: np.ndarray, coefficient: str) -> None:
    """
    Refuse ratings where some rater did not rate some subject.

    Conger's and Light's kappas are defined for complete data: every rater
    rates every subject. One place for the rule and its message, used by
    their analyzers and by the permutation analyses.

    Parameters
    ----------
    category_index : np.ndarray
        The ratings, EMPTY where a rater did not rate a subject.
    coefficient : str
        "Conger's kappa" or "Light's kappa", as the message names it.

    Raises
    ------
    ValueError
        If any rating is missing, saying how many and what to use instead.
    """
    missing = category_index == EMPTY
    if not missing.any():
        return
    raise ValueError(
        f"{coefficient} needs every rater to rate every subject "
        f"({_COMPLETE_DATA_SOURCES[coefficient]}); these data have "
        f"{int(missing.sum()):,} missing ratings in "
        f"{int(missing.any(axis=1).sum()):,} subjects. For data with missing "
        f"ratings use Fleiss' kappa (exact when every subject has the same "
        f"number of ratings, an approximation otherwise with "
        f"balanced_data_only: false), Gwet's AC1, Brennan-Prediger or "
        f"Krippendorff's alpha, or Cohen's kappa pair by pair in a pairwise "
        f"permutation analysis."
    )


class IndividualMarginalAnalyzer(Analyzer, ABC):
    """
    I am the abstract base for individual marginal coefficients (Cohen, Conger, Light).

    Individual marginal coefficients keep rater identity: they read which rater
    said what, and use each rater's own category distribution in the chance
    term. I am the sibling of GroupAnalyzer, whose inheritors do not need rater identity.

    My inheritors' chance terms use each rater's own category rates, so their
    coefficient functions read the ratings array, where each column is one
    rater. Data available only as counts per subject fit a group coefficient
    instead (Fleiss, Scott, Gwet's AC1, Brennan-Prediger, Krippendorff).

    Attributes
    ----------
    ratings : Ratings
        The ratings data. Individual marginal coefficients read which rater
        gave which rating.
    config : IndividualMarginalConfig
        Must be an instance of IndividualMarginalConfig or a subclass
        (CohenKappaConfig, CongerKappaConfig, LightKappaConfig).
    balance_profile : BalanceProfile
        Describes the data structure (from Analyzer).
    rater_profile : RaterProfile
        Computed from ratings: who the raters are, how much each
        rated, how much each pair overlapped (from Analyzer).
    weights : Weights | None
        From the semantics (from Analyzer).
    estimator_input : EstimatorInput
        Carrying category_index, n_categories, and weights (from Analyzer).
    """

    # Attribute declaration

    # Populated at parent class Analyzer
    ratings: Ratings
    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    estimator_input: EstimatorInput
    # Checked at this level
    config: IndividualMarginalConfig

    # Constructor ------------------------------------------------------------ #

    def __init__(
        self, ratings: Ratings, config: IndividualMarginalConfig
    ) -> None:
        """
        Run the parent's checks, then check the config type and the weights.

        Parameters
        ----------
        ratings : Ratings
            The ratings data to analyze.
        config : IndividualMarginalConfig
            The analysis configuration for an individual marginal coefficient.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object, or config object is not from an inheritor
            of IndividualMarginalConfig.
        ValueError
            If the semantics declares weights and my coefficient has no
            weighted form.
        """
        # Call parent __init__ (validates ratings and config, applies the
        # injection, builds the profiles, the weights and the EstimatorInput)
        try:
            super().__init__(ratings, config)
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # Validate config type
        self._validate_individual_config()

        # Weights arrive with the semantics; refuse them if unusable here
        self._refuse_weights_if_unsupported()

    # Private helpers -------------------------------------------------------- #

    def _validate_individual_config(self) -> None:
        """
        Validate that self.config is an IndividualMarginalConfig or subclass.

        Raises
        ------
        TypeError
            If config is not IndividualMarginalConfig.
        """
        from .configs import IndividualMarginalConfig

        if not isinstance(self.config, IndividualMarginalConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be IndividualMarginalConfig or "
                f"subclass (CohenKappaConfig, CongerKappaConfig, "
                f"LightKappaConfig). "
                f"Got {type(self.config).__name__}"
            )
