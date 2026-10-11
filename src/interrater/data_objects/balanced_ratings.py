"""Ratings where every subject was rated the same number of times."""

from __future__ import annotations

import numpy as np

from .. import array_utils as _a
from . import _ratings_validation as _v
from .ratings import Ratings
from ..semantics.ratings_semantics import RatingsSemantics

_ERROR_PREFIX = "[BalancedRatings] Error"

class BalancedRatings(Ratings):
    """
    I am Ratings where every subject was rated by the same *number* of raters.

    Not necessarily the same raters. Subject 1 may have gone to Ann and Bob
    while subject 2 went to Cara and Dev; what matters is that both were seen
    twice. This is the ordinary shape of a rating task where each subject is
    assigned to a fixed-size subset of the rater pool.

    What I buy you is exact Fleiss' kappa. The formula needs a constant number
    of raters per subject and never needs the same raters, so this guarantee --
    not the stronger one of a complete rectangle -- is precisely what it asks
    for. Given only Ratings, an approximation using an averaged rater count is
    the best available, and that is an approximation nobody should have to
    discover by reading the source.

    I take the same three things my parent does and check one more. Nothing
    narrows me for you: hand me an array and I will tell you whether it is
    what I claim.
    """

    # Attributes -------------------------------------------------------------  #

    # Inherited from parent class (Ratings)
    category_index: np.ndarray
    subject_ids:    np.ndarray
    semantics:      RatingsSemantics

    def _validate(self) -> None:
        try:
            super()._validate()
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e

        _v.validate_balanced(self.category_index, type(self).__name__)

    # Public API ---------------------------------------------------------  #

    @property
    def raters_per_subject(self) -> int:
        """
        The constant every subject shares -- the n exact Fleiss' kappa needs.

        Computed, not held, for the same reason my parent holds nothing it was
        not given.
        """
        return _a.raters_per_subject(self.category_index)

    @classmethod
    def build_from_ratings(cls, ratings: Ratings) -> BalancedRatings:
        """
        Convert a Ratings object to BalancedRatings.

        Takes a Ratings object and returns a BalancedRatings object with the
        same data. The conversion will validate that the data is actually
        balanced (constant raters per subject).

        Parameters
        ----------
        ratings : Ratings
            A Ratings object to convert.

        Returns
        -------
        BalancedRatings
            A BalancedRatings object wrapping the same data.

        Raises
        ------
        ValueError
            If the ratings are not actually balanced (varying raters per subject).
        """
        return cls(
            category_index=ratings.category_index,
            semantics=ratings.semantics,
            subject_ids=ratings.subject_ids,
        )
