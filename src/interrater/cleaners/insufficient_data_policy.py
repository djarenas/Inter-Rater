"""How thin the data may get before rows and columns are dropped."""

from __future__ import annotations

from ..contracts import MIN_RATINGS_PER_RATER, MIN_RATINGS_PER_SUBJECT

MAX_DROP_PASSES: int = 100
"""Safety rail on ``method: iterate``.

Each pass only removes rows and columns, so iteration always converges and
normally stops after two or three. Reaching this limit means something
pathological is happening, and a result truncated mid-cascade is arbitrary --
which is why the report says which way the loop ended rather than only how many
passes it took.
"""


class InsufficientDataPolicy:
    """
    I say how much data is too little, and what to do about it.

    There is extra detail to handle the fact that dropping is not symmetric. 
    Removing a rater can push a subject below its
    threshold, and removing a subject can push a rater below its own, so the
    order changes the answer. I make that order a declared choice rather than
    an accident of implementation.

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
    method:
        ``raters_first`` drops raters, then subjects, once. ``subjects_first``
        reverses that. ``iterate`` repeats until a pass drops nothing, which
        can cascade a dataset down to very little -- the report shows each pass
        so that is visible as a sequence rather than a single surprising total.
    """

    # Attributes -------------------------------------------------------------  #

    min_ratings_per_rater:   int
    min_ratings_per_subject: int
    method:                  str

    _METHODS = ("raters_first", "subjects_first", "iterate")

    # Constructor ------------------------------------------------------------  #

    def __init__(
        self,
        min_ratings_per_rater:   int = MIN_RATINGS_PER_RATER,
        min_ratings_per_subject: int = MIN_RATINGS_PER_SUBJECT,
        method:                  str = "raters_first",
    ) -> None:
        self.min_ratings_per_rater   = min_ratings_per_rater
        self.min_ratings_per_subject = min_ratings_per_subject
        self.method                  = method
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

        if self.method not in self._METHODS:
            raise ValueError(
                f"{name}.method must be one of {list(self._METHODS)}, got "
                f"{self.method!r}. Dropping is not symmetric, so the order is "
                f"a choice you have to make."
            )

    # Class methods ----------------------------------------------------------  #

    @classmethod
    def from_dict(cls, spec: dict) -> InsufficientDataPolicy:
        unknown = sorted(set(spec) - {
            "min_ratings_per_rater", "min_ratings_per_subject", "method"
        })
        if unknown:
            raise ValueError(
                f"insufficient_data_drops has unknown key(s) {unknown}. Known: "
                f"min_ratings_per_rater, min_ratings_per_subject, method."
            )
        return cls(**spec)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, InsufficientDataPolicy):
            return NotImplemented
        return (
            self.min_ratings_per_rater == other.min_ratings_per_rater
            and self.min_ratings_per_subject == other.min_ratings_per_subject
            and self.method == other.method
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return (
            f"InsufficientDataPolicy(per_rater="
            f"{self.min_ratings_per_rater}, per_subject="
            f"{self.min_ratings_per_subject}, method={self.method!r})"
        )
