"""The draws of the independent mean, built while the units are bootstrapped."""

from __future__ import annotations

import numpy as np

_ERROR_PREFIX = "[IndependentMeanDraws] Error"


class IndependentMeanDraws:
    """
    I am the draws of the independent mean (Computation C, method
    "independent").

    Each of my draws averages one value picked at random from each defined
    unit's bootstrap draws (Computation B). Units are added one at a time, in
    unit order, so B can discard a unit's draws once added.

    My picks come from one stream seeded with config_of_coefficient's
    random_state, so the order fixes the result; a unit out of order is
    refused as a bug. The same seed gives the same draws whatever keep_draws
    says.

    A NaN pick makes that draw NaN, so it counts as failed. With no defined
    unit, every draw is NaN.
    """

    def __init__(
        self, units: list, n_resamples: int, random_state: int | None
    ) -> None:
        """
        Start an empty running total.

        Parameters
        ----------
        units : list
            The defined units (pairs or subsets, as rater-position tuples),
            in the order they will be added.
        n_resamples : int
            How many draws to build.
        random_state : int | None
            Seed for the picks: config_of_coefficient's random_state.
        """
        self._units = list(units)
        self._defined = set(self._units)
        self._added = 0
        self._rng = np.random.default_rng(random_state)
        self._total = np.zeros(n_resamples, dtype=np.float64)

    def add(self, unit: tuple, draws: np.ndarray) -> None:
        """Pick from a unit's draws; units that are not defined are skipped."""
        if unit not in self._defined:
            return
        # The picks come off one stream, so the order fixes the result.
        if self._added >= len(self._units) or unit != self._units[self._added]:
            raise ValueError(
                f"{_ERROR_PREFIX}: independent mean got unit {unit} out of "
                f"order; this is a bug, not a property of the data."
            )
        picks = self._rng.integers(0, draws.size, size=self._total.size)
        self._total += draws[picks]
        self._added += 1

    def draws(self) -> np.ndarray:
        """The mean's draws, once every defined unit has been added."""
        if self._added != len(self._units):
            raise ValueError(
                f"{_ERROR_PREFIX}: independent mean has {self._added} of "
                f"{len(self._units)} units; this is a bug, not a property "
                f"of the data."
            )
        if not self._units:
            return np.full(self._total.size, np.nan)
        return self._total / len(self._units)
