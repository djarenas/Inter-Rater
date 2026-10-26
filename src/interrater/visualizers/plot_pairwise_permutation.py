"""Visualization for pairwise permutation analysis results."""

from __future__ import annotations

import textwrap
from collections import defaultdict
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from ..intermediates.results.permutation_results.permutation_result import (
    INDEPENDENT_MEAN_CAVEAT,
)

if TYPE_CHECKING:
    from ..intermediates.results import PairwisePermutationResult

_ERROR_PREFIX = "[plot_pairwise_permutation] Error"


def plot_pairwise_permutation(
    result: PairwisePermutationResult,
    rater_names: list[str] | None = None,
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """
    I plot pairwise permutation results as a grouped-by-rater scatter.

    Each rater gets a column on the x-axis. Every pair is plotted under
    both raters it belongs to, with a label showing the partner. Each
    partner rater is assigned a unique color and marker shape so that
    pairs are visually distinguishable (including for color-blind users).
    The overall mean is drawn as a horizontal reference line, and its
    bootstrap CI (if available) as a shaded band. Per-pair bootstrap CIs
    appear as vertical whiskers when present.

    A pair can be bootstrapped and still have no interval: every draw
    failed, so no method ran. Those get no whisker, which on its own would
    look exactly like a pair that was never bootstrapped, so their marker
    is drawn unfilled and the legend says what that means. A plot that
    silently omits an interval is the visual form of the same mistake the
    result object was built to stop.

    Parameters
    ----------
    result : PairwisePermutationResult
        The result object from a pairwise permutation analysis.
    rater_names : list[str] | None
        Human-readable names for the raters. If None, raters are labeled
        by their 0-based index ("Rater 0", "Rater 1", ...). Length must
        match the number of distinct raters found in the pair keys.
    ax : plt.Axes | None
        Axes to draw on. If None, a new Figure and Axes are created.

    Returns
    -------
    plt.Figure
        The figure containing the plot.

    Raises
    ------
    TypeError
        If result is not a PairwisePermutationResult.
    ValueError
        If rater_names length does not match the number of raters.
    """
    # Lazy import to avoid circular dependency at module level
    from ..intermediates.results import PairwisePermutationResult

    if not isinstance(result, PairwisePermutationResult):
        raise TypeError(
            f"{_ERROR_PREFIX}: result must be PairwisePermutationResult, "
            f"got {type(result).__name__}"
        )

    # ------------------------------------------------------------------ #
    # Resolve rater identities
    # ------------------------------------------------------------------ #

    # Collect all unique rater indices from pair keys
    all_rater_indices = sorted(
        {idx for pair in result.pairwise_values for idx in pair}
    )
    n_raters = len(all_rater_indices)

    if rater_names is not None:
        if len(rater_names) != n_raters:
            raise ValueError(
                f"{_ERROR_PREFIX}: rater_names has {len(rater_names)} "
                f"entries but there are {n_raters} raters in the result"
            )
        name_map = dict(zip(all_rater_indices, rater_names))
    else:
        name_map = {idx: f"Rater {idx}" for idx in all_rater_indices}

    # ------------------------------------------------------------------ #
    # Group pairs by rater
    # ------------------------------------------------------------------ #

    # For each rater, collect (partner_index, pair_key, value)
    rater_pairs = defaultdict(list)
    for pair, value in result.pairwise_values.items():
        r_a, r_b = pair
        rater_pairs[r_a].append((r_b, pair, value))
        rater_pairs[r_b].append((r_a, pair, value))

    # Sort each rater's partners by index for consistent ordering
    for rater_idx in rater_pairs:
        rater_pairs[rater_idx].sort(key=lambda x: x[0])

    # ------------------------------------------------------------------ #
    # Create figure
    # ------------------------------------------------------------------ #

    if ax is None:
        fig, ax = plt.subplots(figsize=(max(6, n_raters * 1.5), 5))
    else:
        fig = ax.get_figure()

    # ------------------------------------------------------------------ #
    # Draw mean reference line and bootstrap CI band
    # ------------------------------------------------------------------ #

    mean_val = result.mean_over_pairs
    ax.axhline(
        mean_val, color="black", linewidth=1.0, linestyle="--",
        label=f"Mean = {mean_val:.3f}", zorder=1,
    )

    mean_ci = result.bootstrap_ci_of_mean
    if mean_ci is not None and mean_ci.is_defined:
        label = (
            f"Mean {mean_ci.confidence_level:.0%} CI "
            f"[{mean_ci.lower:.3f}, {mean_ci.upper:.3f}]"
        )
        if result.bootstrap_mean_method == "independent":
            label += "\n" + textwrap.fill(INDEPENDENT_MEAN_CAVEAT, 43)
        ax.axhspan(
            mean_ci.lower, mean_ci.upper, alpha=0.15, color="steelblue",
            label=label,
            zorder=0,
        )

    # ------------------------------------------------------------------ #
    # Assign color + marker per rater (for partner identity)
    # ------------------------------------------------------------------ #

    _COLORS = [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf",
    ]
    _MARKERS = ["o", "s", "D", "^", "v", "h", "p", "*", "X", "P"]

    rater_style = {}
    for i, idx in enumerate(all_rater_indices):
        rater_style[idx] = {
            "color": _COLORS[i % len(_COLORS)],
            "marker": _MARKERS[i % len(_MARKERS)],
        }

    # ------------------------------------------------------------------ #
    # Plot points with jitter within each rater column
    # ------------------------------------------------------------------ #

    # x positions: one integer per rater
    rater_x = {idx: i for i, idx in enumerate(all_rater_indices)}

    # Jitter spread within a column
    jitter_width = 0.15  # total half-width of jitter

    # Set when some pair was bootstrapped but produced no interval
    any_undefined_ci = False

    for rater_idx in all_rater_indices:
        partners = rater_pairs[rater_idx]
        n_pts = len(partners)
        x_center = rater_x[rater_idx]

        # Evenly space points within the jitter band
        if n_pts == 1:
            offsets = [0.0]
        else:
            offsets = np.linspace(-jitter_width, jitter_width, n_pts)

        for j, (partner_idx, pair_key, value) in enumerate(partners):
            x_pos = x_center + offsets[j]

            # Skip NaN values
            if np.isnan(value):
                continue

            # Style from partner rater
            style = rater_style[partner_idx]

            # Per-pair CI whisker. An interval can exist and still have no
            # bounds -- every draw failed, so no method ran. There is
            # nothing to draw for one, so the marker carries the news
            # instead.
            pair_ci = result.bootstrap_cis_per_pair[pair_key]
            has_interval = pair_ci.is_defined
            if not has_interval:
                any_undefined_ci = True
            else:
                ci_lo, ci_hi = pair_ci.lower, pair_ci.upper
                ax.plot(
                    [x_pos, x_pos], [ci_lo, ci_hi],
                    color="gray", linewidth=1.0, zorder=2,
                )
                # Small caps on whiskers
                cap_w = 0.04
                ax.plot(
                    [x_pos - cap_w, x_pos + cap_w], [ci_lo, ci_lo],
                    color="gray", linewidth=1.0, zorder=2,
                )
                ax.plot(
                    [x_pos - cap_w, x_pos + cap_w], [ci_hi, ci_hi],
                    color="gray", linewidth=1.0, zorder=2,
                )

            # Point — color and marker identify the partner rater, and an
            # unfilled face says the interval was attempted and came back
            # undefined.
            ax.plot(
                x_pos, value,
                marker=style["marker"], color=style["color"],
                markerfacecolor=style["color"] if has_interval else "none",
                markersize=7, linestyle="None",
                zorder=3,
            )

            # Small text label showing partner
            partner_label = name_map[partner_idx]
            ax.annotate(
                partner_label,
                (x_pos, value),
                textcoords="offset points",
                xytext=(6, 2),
                fontsize=6,
                color="dimgray",
                zorder=4,
            )

    # ------------------------------------------------------------------ #
    # Axes formatting
    # ------------------------------------------------------------------ #

    ax.set_xticks(range(n_raters))
    ax.set_xticklabels(
        [name_map[idx] for idx in all_rater_indices],
        fontsize=9,
    )
    ax.set_xlabel("Rater", fontsize=10)
    ax.set_ylabel(result.coefficient_name, fontsize=10)

    # Title
    method_str = ""
    if result.bootstrap_mean_method is not None:
        method_str = f" ({result.bootstrap_mean_method} bootstrap)"
    ax.set_title(
        f"Pairwise {result.coefficient_name}{method_str}",
        fontsize=12,
    )

    # Y-axis: default [0, 1] for kappa scale, widen if data falls outside
    all_values = [v for v in result.pairwise_values.values() if not np.isnan(v)]
    y_min = 0.0
    y_max = 1.0
    if all_values:
        data_min = min(all_values)
        data_max = max(all_values)
        # Include CI bounds if present. Undefined intervals are skipped --
        # min() against a NaN returns whichever argument came first, so one
        # of them would set the axis limits at random.
        for pair_ci in result.bootstrap_cis_per_pair.values():
            if pair_ci.is_defined:
                data_min = min(data_min, pair_ci.lower)
                data_max = max(data_max, pair_ci.upper)
        if mean_ci is not None and mean_ci.is_defined:
            data_min = min(data_min, mean_ci.lower)
            data_max = max(data_max, mean_ci.upper)
        if data_min < 0:
            y_min = data_min - 0.05
        if data_max > 1:
            y_max = data_max + 0.05
    ax.set_ylim(y_min, y_max)

    # Proxy entry so an unfilled marker is readable without the docstring
    if any_undefined_ci:
        ax.plot(
            [], [], marker="o", color="dimgray", markerfacecolor="none",
            linestyle="None", markersize=7,
            label="Open marker: no CI (every draw failed)",
        )

    ax.legend(fontsize=8, loc="best")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    return fig
