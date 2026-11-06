"""Plotter for a BalanceProfile: how many raters each subject got."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt

from ..intermediates.data_profiles import BalanceProfile

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.estimator_parts import NijMatrix
    from ..intermediates.results import PermutationResult, PointEstimateResult

_ERROR_PREFIX = "[BalanceProfilePlotter] Error"


@dataclass(frozen=True)
class BalanceProfilePlotter:
    """
    I plot a BalanceProfile: how many subjects got 1, 2, 3, ... raters.

    One bar for every count from the fewest raters any subject got to the
    most, so a count no subject had shows as a zero bar rather than a gap.
    Each bar has its count on top: next to 1,614 subjects with two raters,
    3 subjects with three would otherwise be invisible. The title says
    whether the data is balanced.

    A BalanceProfile is subject-side only, so a counts matrix has one too:
    I can be built from Ratings, an NijMatrix, or any result. How much each
    rater rated is RaterProfilePlotter.draw_counts.

    Attributes
    ----------
    balance_profile : BalanceProfile
        The subjects to plot.
    """

    # Attribute declaration
    balance_profile: BalanceProfile

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes."""
        if not isinstance(self.balance_profile, BalanceProfile):
            raise TypeError(
                f"{_ERROR_PREFIX}: balance_profile must be a BalanceProfile, "
                f"got {type(self.balance_profile).__name__}"
            )

    @classmethod
    def build_from_ratings(cls, ratings: Ratings) -> BalanceProfilePlotter:
        """
        Build me from a Ratings object.

        Parameters
        ----------
        ratings : Ratings
            The ratings whose subjects to plot.

        Returns
        -------
        BalanceProfilePlotter
        """
        return cls(BalanceProfile.build_from_ratings(ratings))

    @classmethod
    def build_from_nij_matrix(cls, nij: NijMatrix) -> BalanceProfilePlotter:
        """
        Build me from a counts matrix.

        Parameters
        ----------
        nij : NijMatrix
            The counts matrix (n_subjects x n_categories).

        Returns
        -------
        BalanceProfilePlotter
        """
        return cls(BalanceProfile.build_from_nij_matrix(nij))

    @classmethod
    def build_from_result(
        cls, result: PointEstimateResult | PermutationResult,
    ) -> BalanceProfilePlotter:
        """
        Build me from a coefficient's result or a permutation result.

        Parameters
        ----------
        result : PointEstimateResult | PermutationResult
            Any result that carries a balance profile.

        Returns
        -------
        BalanceProfilePlotter

        Raises
        ------
        TypeError
            If result is neither kind of result, or has no balance profile.
        """
        from ..intermediates.results import PermutationResult, PointEstimateResult

        if not isinstance(result, (PointEstimateResult, PermutationResult)):
            raise TypeError(
                f"{_ERROR_PREFIX}: result must be a PointEstimateResult or a "
                f"PermutationResult, got {type(result).__name__}"
            )
        return cls(result.balance_profile)

    # Public API ------------------------------------------------------------- #

    def draw(self, ax: plt.Axes | None = None) -> plt.Figure:
        """
        Draw the number of subjects for each count of raters.

        Parameters
        ----------
        ax : plt.Axes | None
            Axes to draw on. If None, a new Figure and Axes are created.

        Returns
        -------
        plt.Figure
            The figure containing the plot.
        """
        profile = self.balance_profile
        counts = profile.subject_counts
        lowest, highest = min(counts), max(counts)
        raters = list(range(lowest, highest + 1))
        subjects = [counts.get(k, 0) for k in raters]

        if ax is None:
            fig, ax = plt.subplots(figsize=(max(5.0, 0.5 * len(raters) + 3.0), 4.5))
        else:
            fig = ax.figure

        bars = ax.bar(raters, subjects, color="steelblue", zorder=2)
        ax.bar_label(bars, labels=[f"{s:,}" for s in subjects], fontsize=8, padding=2)

        # ------------------------------------------------------------------ #
        # Labels
        # ------------------------------------------------------------------ #

        def _raters(k: int) -> str:
            return f"{k} rater" if k == 1 else f"{k} raters"

        if profile.is_balanced:
            title = f"Balanced: every subject has {_raters(profile.raters_per_subject)}"
        else:
            title = f"Unbalanced: subjects have {lowest} to {_raters(highest)}"

        ax.set_xticks(raters)
        ax.set_ylim(0, max(subjects) * 1.15)
        ax.set_title(f"{title}\n{profile.n_subjects:,} subjects", fontsize=12)
        ax.set_xlabel("Raters per subject", fontsize=10)
        ax.set_ylabel("Subjects", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()

        return fig
