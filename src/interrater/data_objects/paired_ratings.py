"""Ratings from exactly two raters, with no gaps between them."""

from __future__ import annotations

from . import _ratings_validation as _v
from .balanced_ratings import BalancedRatings

_ERROR_PREFIX = "[PairedRatings] Error"

class PairedRatings(BalancedRatings):
    """
    I am two raters and the subjects they both rated.

    Cohen's kappa can consume me and nothing narrower is needed. A coefficient
    defined for two raters can simply be handed an object that cannot exist
    with three, so the precondition lives in the signature and no analyser
    needs to check it.

    Unlike my parent class, I have no empties; for two raters it means both rated every subject I
    hold. That is why extracting me from a wider dataset means first reducing
    to the subjects the pair share: the overlap is where a paired comparison
    is defined, and taking it is part of what makes me rather than something
    done to me afterwards. It also makes me crossed, and therefore balanced,
    which is why I sit where I do.

    A pair sharing few subjects produces a kappa that is valid and worth
    little. I do not police that -- somebody may have no better data and still
    need the number.

    I know my raters by name, through my semantics, but not which columns of a
    wider dataset I came from. Whoever extracted me knows that, and it is
    their bookkeeping rather than my memory.
    """

    def _validate(self) -> None:
        # My own two checks run before the inherited chain, because for a pair
        # they are the direct statement of what is wrong. Three raters should
        # hear "exactly 2 raters", and a gap should hear "no gaps" -- not a
        # complaint about balance, which for two raters is only a symptom of
        # the gap.
        _v.validate_exactly_two_raters(self.category_index, type(self).__name__)
        _v.validate_crossed(self.category_index, type(self).__name__)
        try:
            super()._validate()
        except (TypeError, ValueError) as e:
            raise type(e)(f"{_ERROR_PREFIX}: {str(e)}") from e
