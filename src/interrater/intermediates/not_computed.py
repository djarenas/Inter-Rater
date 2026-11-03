"""A value that was deliberately not computed, and why."""

from __future__ import annotations

from dataclasses import dataclass

_ERROR_PREFIX = "[NotComputed] Error"


@dataclass(frozen=True)
class NotComputed:
    """
    I stand in for a value that does not apply here, and say why.

    ``None`` cannot tell "does not apply" from "was never filled in", so a
    field that is legitimately absent holds me instead, with a reason a
    reader can act on -- for example, that a bias index needs each rater's
    ratings and the result was built from category counts only.

    Three states, kept apart:

    - a number: computed;
    - NaN: computed, but undefined for this data (kappa when Pe = 1);
    - NotComputed: not applicable here, with the reason.

    I refuse to be used as a truth value. ``if result.kappa_max:`` would be
    true for me and false for a computed 0.0, which is wrong both ways; check
    ``isinstance(value, NotComputed)`` instead.

    Attributes
    ----------
    reason : str
        Why the value was not computed. Must not be empty.
    """

    reason: str

    def __post_init__(self) -> None:
        if not isinstance(self.reason, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: reason must be str, "
                f"got {type(self.reason).__name__}"
            )
        if not self.reason.strip():
            raise ValueError(
                f"{_ERROR_PREFIX}: reason must say why the value was not "
                f"computed; got an empty string"
            )

    def __bool__(self) -> bool:
        raise TypeError(
            f"{_ERROR_PREFIX}: a NotComputed value has no truth value "
            f"({self.reason}). Check isinstance(value, NotComputed) instead."
        )

    def __str__(self) -> str:
        return f"not computed: {self.reason}"

    def summary(self) -> dict:
        """
        Return my record in plain types, for writing to disk.

        Returns
        -------
        dict
            ``{"not_computed": reason}``.
        """
        return {"not_computed": self.reason}
