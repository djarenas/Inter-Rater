"""How thin the data may get before rows and columns are dropped."""

from __future__ import annotations

from ..contracts import MIN_RATINGS_PER_RATER, MIN_RATINGS_PER_SUBJECT


class InsufficientDataPolicy:
    """
    I say how much data is too little, and what to do about it.

    There is extra detail to handle the fact that dropping is not symmetric.
    Removing a rater can push a subject below its
    threshold, and removing a subject can push a rater below its own, so one
    sweep is not enough. Sweeps repeat until one drops nothing. What is left is
    the largest data meeting both minimums, and it is the same whichever axis
    is swept first, so the order is not a choice. The cascade can run a dataset
    down to very little -- the report shows each pass so that is visible as a
    sequence rather than a single surprising total.

    My minimums may be raised but never lowered: the floors in ``contracts``
    are what makes a row or column able to participate at all, and below them a
    subject has nobody to disagree with. Raising a minimum is a rating task saying it
    wants more evidence than the arithmetic strictly needs.

    Attributes
    ----------
    min_ratings_per_rater:
        How many subjects a rater must have rated to be kept.
    min_ratings_per_subject:
        How many raters must have rated a subject for it to be kept.
    """

    # Attributes -------------------------------------------------------------  #

    min_ratings_per_rater:   int
    min_ratings_per_subject: int

    # Constructor ------------------------------------------------------------  #

    def __init__(
        self,
        min_ratings_per_rater:   int = MIN_RATINGS_PER_RATER,
        min_ratings_per_subject: int = MIN_RATINGS_PER_SUBJECT,
        method:                  object = None,
    ) -> None:
        # Removed in 2.1.0; kept as a name so it gets this error, not
        # "unexpected keyword".
        if method is not None:
            raise TypeError(
                f"{type(self).__name__}.method was removed in 2.1.0 (got "
                f"{method!r}). Dropping always repeats until nothing more falls "
                f"below the minimums: the largest data meeting both, whatever "
                f"the order. Delete the method line."
            )
        self.min_ratings_per_rater   = min_ratings_per_rater
        self.min_ratings_per_subject = min_ratings_per_subject
        self._validate()

    def _validate(self) -> None:
        name = type(self).__name__

        for field, value, floor in (
            ("min_ratings_per_rater", self.min_ratings_per_rater,
             MIN_RATINGS_PER_RATER),
            ("min_ratings_per_subject", self.min_ratings_per_subject,
             MIN_RATINGS_PER_SUBJECT),
        ):
            if not isinstance(value, int) or isinstance(value, bool):
                raise TypeError(
                    f"{name}.{field} must be an integer, got "
                    f"{type(value).__name__}."
                )
            if value < floor:
                raise ValueError(
                    f"{name}.{field} is {value}, below the floor of {floor}. "
                    f"That floor is what makes a rater or subject able to take "
                    f"part in any comparison at all; it can be raised as a "
                    f"deliberate choice, never lowered."
                )

    # Class methods ----------------------------------------------------------  #

    @classmethod
    def from_dict(cls, spec: dict) -> InsufficientDataPolicy:
        # method is passed on so it gets the removed-in-2.1.0 error.
        unknown = sorted(set(spec) - {
            "min_ratings_per_rater", "min_ratings_per_subject", "method"
        })
        if unknown:
            raise ValueError(
                f"insufficient_data_drops has unknown key(s) {unknown}. Known: "
                f"min_ratings_per_rater, min_ratings_per_subject."
            )
        return cls(**spec)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, InsufficientDataPolicy):
            return NotImplemented
        return (
            self.min_ratings_per_rater == other.min_ratings_per_rater
            and self.min_ratings_per_subject == other.min_ratings_per_subject
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return (
            f"InsufficientDataPolicy(per_rater="
            f"{self.min_ratings_per_rater}, per_subject="
            f"{self.min_ratings_per_subject})"
        )
