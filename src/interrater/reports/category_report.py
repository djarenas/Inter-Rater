"""Which declared categories were actually used."""

from __future__ import annotations

from ._wording import counted

from dataclasses import dataclass, field


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
    three hundred cells passes any test for being used, and every estimate
    involving it is imprecise.

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
    labels:
        Category name to its declared label, shown beside the name when it
        says something the name does not.

    The counts are of the cleaned data: what survived matching and dropping,
    which is what every coefficient is computed on.
    """

    counts:         dict[str, int]
    policy:         str
    action:         str = "kept"
    rare_threshold: int = 10
    unused:         tuple[str, ...] = ()
    labels:         dict[str, str] = field(default_factory=dict)

    @property
    def rare(self) -> dict[str, int]:
        """Used, but so seldom that estimates involving them are imprecise."""
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
        rare  = self.rare

        # Rare and unused are marked on the category's own row, so nothing
        # else in this section sits where a category name would.
        def mark(name: str) -> str:
            if name in rare:
                return " (rare)"
            if name in self.unused:
                return f" (unused, {self.action})"
            return ""

        def named(name: str) -> str:
            label = self.labels.get(name, "")
            # A label only differing in case says nothing the name does not.
            if label and label.casefold() != name.casefold():
                return f"{name} ({label})"
            return name

        n_ratings = sum(self.counts.values())
        lines = [f"Categories (cleaned data: {counted(n_ratings, 'rating')})"]
        if verbosity == "full":
            width = max([14] + [len(named(n)) + 1 for n in self.counts])
            for name, n in self.counts.items():
                lines.append(
                    f"{pad}{named(name):<{width}} {n:>8,}  ({n / total:.1%}){mark(name)}"
                )
        else:
            shown = ", ".join(
                f"{name} {n / total:.0%}{mark(name)}"
                for name, n in self.counts.items()
            )
            lines.append(f"{pad}share      : {shown}")

        if rare:
            lines += [
                "",
                f"{pad}(rare = fewer than {self.rare_threshold} ratings; "
                f"estimates involving it are imprecise)",
            ]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.render()
