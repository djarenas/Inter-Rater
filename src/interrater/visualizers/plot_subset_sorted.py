"""Visualization for subset permutation results -- sorted dot plot."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    from ..intermediates.analysis_results import SubsetPermutationResult

_ERROR_PREFIX = "[plot_subset_sorted] Error"


def plot_subset_sorted(
    result: SubsetPermutationResult,
    rater_names: list[str] | None = None,
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """
    I plot subset permutation results as a sorted dot plot.

    Subsets are sorted by coefficient value (ascending, left to right).
    Each subset is a dot at its coefficient value. Per-subset bootstrap
    CIs appear as vertical whiskers when present. The overall mean is
    drawn as a horizontal reference line, and its bootstrap CI (if
    available) as a shaded band.

    X-tick labels show the rater composition of each subset (e.g.
    "A, B, C"). When n_subsets exceeds a readability threshold, labels
    are rotated or dropped entirely.

    Parameters
    ----------
    result : SubsetPermutationResult
        The result object from a subset permutation analysis (single k).
    rater_names : list[str] | None
        Human-readable names for the raters. If None, raters are labeled
        by their 0-based index ("0", "1", ...). Length must match the
        number of distinct raters found in the subset keys.
    ax : plt.Axes | None
        Axes to draw on. If None, a new Figure and Axes are created.

    Returns
    -------
    plt.Figure
        The figure containing the plot.

    Raises
    ------
    TypeError
        If result is not a SubsetPermutationResult.
    ValueError
        If rater_names length does not match the number of raters.
    """
    # Lazy import to avoid circular dependency at module level
    from ..intermediates.analysis_results import SubsetPermutationResult

    if not isinstance(result, SubsetPermutationResult):
        raise TypeError(
            f"{_ERROR_PREFIX}: result must be SubsetPermutationResult, "
            f"got {type(result).__name__}"
        )

    # ------------------------------------------------------------------ #
    # Resolve rater names
    # ------------------------------------------------------------------ #

    all_rater_indices = sorted(
        {idx for subset in result.subset_values for idx in subset}
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
        name_map = {idx: str(idx) for idx in all_rater_indices}

    # ------------------------------------------------------------------ #
    # Sort subsets by coefficient value
    # ------------------------------------------------------------------ #

    sorted_subsets = sorted(
        result.subset_values.items(),
        key=lambda item: (
            item[1] if not np.isnan(item[1]) else float("-inf")
        ),
    )
    subset_keys = [s[0] for s in sorted_subsets]
    subset_vals = np.array([s[1] for s in sorted_subsets])
    n_subsets = len(subset_keys)

    # ------------------------------------------------------------------ #
    # Build subset labels
    # ------------------------------------------------------------------ #

    subset_labels = [
        ", ".join(name_map[idx] for idx in subset)
        for subset in subset_keys
    ]

    # ------------------------------------------------------------------ #
    # Create figure
    # ------------------------------------------------------------------ #

    width = max(6, n_subsets * 0.4)
    width = min(width, 20)  # cap at reasonable size

    if ax is None:
        fig, ax = plt.subplots(figsize=(width, 5))
    else:
        fig = ax.get_figure()

    x_positions = np.arange(n_subsets)

    # ------------------------------------------------------------------ #
    # Draw mean reference line and bootstrap CI band
    # ------------------------------------------------------------------ #

    mean_val = result.mean_over_subsets
    ax.axhline(
        mean_val, color="black", linewidth=1.0, linestyle="--",
        label=f"Mean = {mean_val:.3f}", zorder=1,
    )

    if result.bootstrap_ci_of_mean is not None:
        ci_lower, ci_upper = result.bootstrap_ci_of_mean
        ax.axhspan(
            ci_lower, ci_upper, alpha=0.15, color="steelblue",
            label=f"Mean 95% CI [{ci_lower:.3f}, {ci_upper:.3f}]",
            zorder=0,
        )

    # ------------------------------------------------------------------ #
    # Draw per-subset CI whiskers
    # ------------------------------------------------------------------ #

    if result.bootstrap_cis_per_subset is not None:
        for i, subset in enumerate(subset_keys):
            ci = result.bootstrap_cis_per_subset.get(subset)
            if ci is not None:
                ci_lo, ci_hi = ci
                if np.isnan(ci_lo) or np.isnan(ci_hi):
                    continue
                ax.plot(
                    [i, i], [ci_lo, ci_hi],
                    color="gray", linewidth=1.0, zorder=2,
                )
                # Small caps on whiskers
                cap_w = 0.15
                ax.plot(
                    [i - cap_w, i + cap_w], [ci_lo, ci_lo],
                    color="gray", linewidth=1.0, zorder=2,
                )
                ax.plot(
                    [i - cap_w, i + cap_w], [ci_hi, ci_hi],
                    color="gray", linewidth=1.0, zorder=2,
                )

    # ------------------------------------------------------------------ #
    # Draw dots
    # ------------------------------------------------------------------ #

    # Color: subsets below mean in muted red, above in muted blue
    colors = []
    for val in subset_vals:
        if np.isnan(val):
            colors.append("lightgray")
        elif val < mean_val:
            colors.append("#d66")
        else:
            colors.append("#4682b4")

    ax.scatter(
        x_positions, subset_vals,
        c=colors, s=40, edgecolors="black", linewidths=0.5,
        zorder=3,
    )

    # ------------------------------------------------------------------ #
    # X-axis labels -- adapt to n_subsets
    # ------------------------------------------------------------------ #

    if n_subsets <= 20:
        ax.set_xticks(x_positions)
        ax.set_xticklabels(subset_labels, fontsize=8)
        if n_subsets > 10:
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    elif n_subsets <= 50:
        ax.set_xticks(x_positions)
        ax.set_xticklabels(
            subset_labels, fontsize=6, rotation=90, ha="center",
        )
    else:
        # Too many subsets -- drop labels entirely
        ax.set_xticks([])
        ax.set_xlabel(
            f"Subsets (n={n_subsets}, sorted by value)", fontsize=10,
        )

    if n_subsets <= 50:
        ax.set_xlabel("Subset (raters)", fontsize=10)

    # ------------------------------------------------------------------ #
    # Y-axis and title
    # ------------------------------------------------------------------ #

    ax.set_ylabel(result.coefficient_name, fontsize=10)

    method_str = ""
    if result.bootstrap_mean_method is not None:
        method_str = f" ({result.bootstrap_mean_method} bootstrap)"
    ax.set_title(
        f"Subset k={result.subset_k} {result.coefficient_name}"
        f" — sorted by value{method_str}",
        fontsize=12,
    )

    # Y-axis bounds: default [0, 1], widen if data falls outside
    all_values = [v for v in subset_vals if not np.isnan(v)]
    y_min = 0.0
    y_max = 1.0
    if all_values:
        data_min = min(all_values)
        data_max = max(all_values)
        # Include CI bounds if present
        if result.bootstrap_cis_per_subset is not None:
            for ci_lo, ci_hi in result.bootstrap_cis_per_subset.values():
                if not np.isnan(ci_lo):
                    data_min = min(data_min, ci_lo)
                if not np.isnan(ci_hi):
                    data_max = max(data_max, ci_hi)
        if result.bootstrap_ci_of_mean is not None:
            data_min = min(data_min, result.bootstrap_ci_of_mean[0])
            data_max = max(data_max, result.bootstrap_ci_of_mean[1])
        if data_min < 0:
            y_min = data_min - 0.05
        if data_max > 1:
            y_max = data_max + 0.05
    ax.set_ylim(y_min, y_max)

    ax.legend(fontsize=8, loc="best")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    return fig
