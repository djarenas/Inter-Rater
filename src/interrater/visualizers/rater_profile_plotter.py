"""Plotter for a RaterProfile: who shared subjects with whom, and how much each rated."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle

from ..contracts import MIN_OVERLAP_HARD
from ..intermediates.data_profiles import RaterProfile
from ..intermediates.not_computed import NotComputed

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates.results import PermutationResult, PointEstimateResult

_ERROR_PREFIX = "[RaterProfilePlotter] Error"

# What a hatched cell means: the one rule gives such a pair no coefficient
THIN_PAIR_LABEL = f"fewer than {MIN_OVERLAP_HARD} shared: no coefficient"


@dataclass(frozen=True)
class RaterProfilePlotter:
    """
    I plot a RaterProfile, in two views.

    draw_overlap is a heatmap of how many subjects each pair of raters both
    rated. Every pairwise coefficient rests on one of those cells, and on
    partially crossed data a thin pair's coefficient means little without
    it. The diagonal is each rater's own total; it is white, outside the
    colour scale, so the colours compare pairs only. A pair sharing fewer
    than MIN_OVERLAP_HARD subjects is hatched: the analysis gives it no
    coefficient at all.

    draw_counts is a bar per rater: how many subjects each one rated, in
    the declared order (the heatmap's order), with a dashed line at all
    subjects.

    I need rater identity, so I am built from Ratings or from a result that
    has a RaterProfile. A result built from counts has none, and I say so
    when built rather than when drawn.

    Attributes
    ----------
    rater_profile : RaterProfile
        The raters to plot.
    annotate : bool, default True
        Write each cell's count in the overlap heatmap. Turn off on large
        panels, where the numbers get too small to read.
    """

    # Attribute declaration
    rater_profile: RaterProfile
    annotate: bool = True

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate my attributes."""
        if not isinstance(self.rater_profile, RaterProfile):
            raise TypeError(
                f"{_ERROR_PREFIX}: rater_profile must be a RaterProfile, "
                f"got {type(self.rater_profile).__name__}"
            )
        if not isinstance(self.annotate, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: annotate must be bool, "
                f"got {type(self.annotate).__name__}"
            )

    @classmethod
    def build_from_ratings(cls, ratings: Ratings, annotate: bool = True) -> RaterProfilePlotter:
        """
        Build me from a Ratings object.

        Parameters
        ----------
        ratings : Ratings
            The ratings whose raters to plot.
        annotate : bool, default True
            See the class docstring.

        Returns
        -------
        RaterProfilePlotter
        """
        return cls(RaterProfile.build_from_ratings(ratings), annotate=annotate)

    @classmethod
    def build_from_result(
        cls,
        result: PointEstimateResult | PermutationResult,
        annotate: bool = True,
    ) -> RaterProfilePlotter:
        """
        Build me from a coefficient's result or a permutation result.

        Parameters
        ----------
        result : PointEstimateResult | PermutationResult
            Any result that carries a rater profile.
        annotate : bool, default True
            See the class docstring.

        Returns
        -------
        RaterProfilePlotter

        Raises
        ------
        TypeError
            If result is neither kind of result.
        ValueError
            If the result was built from counts, which carry no rater
            identity.
        """
        from ..intermediates.results import PermutationResult, PointEstimateResult

        if not isinstance(result, (PointEstimateResult, PermutationResult)):
            raise TypeError(
                f"{_ERROR_PREFIX}: result must be a PointEstimateResult or a "
                f"PermutationResult, got {type(result).__name__}"
            )
        if isinstance(result.rater_profile, NotComputed):
            raise ValueError(
                f"{_ERROR_PREFIX}: this result has no rater profile "
                f"({result.rater_profile.reason}). Build the plotter from "
                f"the Ratings instead."
            )
        return cls(result.rater_profile, annotate=annotate)

    # Public API ------------------------------------------------------------- #

    def draw_overlap(self, ax: plt.Axes | None = None) -> plt.Figure:
        """
        Draw the heatmap of subjects shared by each pair of raters.

        Parameters
        ----------
        ax : plt.Axes | None
            Axes to draw on. If None, a new Figure and Axes are created.

        Returns
        -------
        plt.Figure
            The figure containing the plot.
        """
        profile = self.rater_profile
        shared = profile.rater_overlap.shared
        n = profile.n_raters
        off_diagonal = ~np.eye(n, dtype=bool)

        if ax is None:
            size = max(4.5, 0.55 * n + 2.5)
            fig, ax = plt.subplots(figsize=(size + 1.2, size))
        else:
            fig = ax.figure

        # ------------------------------------------------------------------ #
        # Colours: pairs only; the diagonal is white, outside the scale
        # ------------------------------------------------------------------ #

        values = shared.astype(float)
        values[~off_diagonal] = np.nan
        top = int(shared[off_diagonal].max()) if n > 1 else 0
        cmap = plt.get_cmap("Blues").with_extremes(bad="white")
        image = ax.imshow(
            np.ma.masked_invalid(values), cmap=cmap,
            vmin=0, vmax=max(top, 1),
        )
        fig.colorbar(image, ax=ax, label="Subjects both rated")

        # Outline the diagonal cells so they read as cells, not gaps
        for i in range(n):
            ax.add_patch(Rectangle(
                (i - 0.5, i - 0.5), 1, 1,
                fill=False, edgecolor="lightgray", linewidth=0.8,
            ))

        # ------------------------------------------------------------------ #
        # Pairs the analysis gives no coefficient
        # ------------------------------------------------------------------ #

        thin = off_diagonal & (shared < MIN_OVERLAP_HARD)
        for i, j in zip(*np.nonzero(thin)):
            ax.add_patch(Rectangle(
                (j - 0.5, i - 0.5), 1, 1,
                fill=False, hatch="///", edgecolor="gray", linewidth=0,
            ))
        if thin.any():
            ax.legend(
                handles=[Patch(facecolor="white", edgecolor="gray",
                               hatch="///", label=THIN_PAIR_LABEL)],
                loc="lower left", bbox_to_anchor=(0, 1.08),
                fontsize=8, frameon=False,
            )

        # ------------------------------------------------------------------ #
        # Counts in the cells
        # ------------------------------------------------------------------ #

        if self.annotate:
            for i in range(n):
                for j in range(n):
                    dark = off_diagonal[i, j] and shared[i, j] > 0.6 * top
                    ax.text(
                        j, i, f"{int(shared[i, j]):,}",
                        ha="center", va="center", fontsize=8,
                        color="white" if dark else "black",
                    )

        # ------------------------------------------------------------------ #
        # Labels
        # ------------------------------------------------------------------ #

        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(profile.rater_names, rotation=45, ha="right", fontsize=9)
        ax.set_yticklabels(profile.rater_names, fontsize=9)
        ax.set_title("Subjects shared by each pair of raters", fontsize=12)
        ax.set_xlabel(
            f"{profile.n_subjects:,} subjects; diagonal: each rater's own total",
            fontsize=10,
        )
        fig.tight_layout()

        return fig

    def draw_counts(self, ax: plt.Axes | None = None) -> plt.Figure:
        """
        Draw a bar per rater: how many subjects each one rated.

        Parameters
        ----------
        ax : plt.Axes | None
            Axes to draw on. If None, a new Figure and Axes are created.

        Returns
        -------
        plt.Figure
            The figure containing the plot.
        """
        profile = self.rater_profile
        n = profile.n_raters

        if ax is None:
            fig, ax = plt.subplots(figsize=(max(5.0, 0.45 * n + 2.0), 4.5))
        else:
            fig = ax.figure

        # In declared order, the heatmap's order, so the two views line up
        ax.bar(range(n), profile.n_per_rater, color="steelblue", zorder=2)
        ax.axhline(
            profile.n_subjects, color="black", linestyle="--", linewidth=1.2,
            zorder=3, label=f"all {profile.n_subjects:,} subjects",
        )

        ax.set_xticks(range(n))
        ax.set_xticklabels(profile.rater_names, rotation=45, ha="right", fontsize=9)
        ax.set_ylim(0, profile.n_subjects * 1.25)
        ax.set_title("Subjects each rater rated", fontsize=12)
        ax.set_ylabel("Subjects", fontsize=10)
        ax.legend(fontsize=8, loc="upper right")
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()

        return fig
