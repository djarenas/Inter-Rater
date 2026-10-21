"""Which declared categories were actually used."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CategoryReport:
    """
    I say which of the declared categories anybody chose.

    A category nobody chose is not automatically wrong. If the scheme offered
    five grades and no student earned a D, that is a fact about the cohort. But
    it is also the shape a mistake takes when a rating task is not what its owner
    thinks it is, which is why the default policy stops rather than shrugs.

    Keeping an unused category is the safe answer, because K is a fact about
    what raters could have said. Dropping one changes any coefficient whose
    chance term counts categories -- Gwet's AC1 and Brennan-Prediger both --
    while leaving kappa untouched, so the damage is quiet and easy to miss.

    I also record the rare categories, which the unused check does not catch
    and which are more dangerous than the empty ones. A category chosen once in
    three hundred cells passes any test for being used, and every statistic
    involving it is noise.

    Attributes
    ----------
    counts:
        Category name to how many ratings it received.
    unused:
        Categories with no ratings at all.
    policy:
        What ``unused_categories`` was set to.
    action:
        What that policy did. Under ``error`` there is no report to read, and
        the cleaner never removes a category, so this is always 'kept'. Kept
        as a field so that report summaries keep the same keys.
    rare_threshold:
        Below this many ratings a category is called rare.
    """

    counts:         dict[str, int]
    policy:         str
    action:         str = "kept"
    rare_threshold: int = 10
    unused:         tuple[str, ...] = ()

    @property
    def rare(self) -> dict[str, int]:
        """Used, but so seldom that anything involving them is noise."""
        return {
            name: n
            for name, n in self.counts.items()
            if 0 < n < self.rare_threshold
        }

    def summary(self) -> dict:
        return {
            "counts":         dict(self.counts),
            "unused":         list(self.unused),
            "rare":           dict(self.rare),
            "rare_threshold": self.rare_threshold,
            "policy":         self.policy,
            "action":         self.action,
        }

    def render(self, verbosity: str = "summary") -> str:
        pad   = "    "
        total = sum(self.counts.values()) or 1
        lines = ["Categories"]

        if verbosity == "full":
            for name, n in self.counts.items():
                lines.append(f"{pad}{name:<14} {n:>8,}  ({n / total:.1%})")
        else:
            shown = ", ".join(
                f"{name} {n / total:.0%}" for name, n in self.counts.items()
            )
            lines.append(f"{pad}used       : {shown}")

        if self.unused:
            lines.append(
                f"{pad}unused     : {list(self.unused)} ({self.action})"
            )
        if self.rare:
            shown = ", ".join(f"{k} ({v})" for k, v in self.rare.items())
            lines.append(
                f"{pad}rare       : {shown} -- fewer than "
                f"{self.rare_threshold} ratings; statistics involving these "
                f"are not estimable"
            )
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.render()
