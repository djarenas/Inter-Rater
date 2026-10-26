"""Visualization for multi-k subset results -- strip-and-trend plot."""

from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from ..intermediates.results.permutation_results.permutation_result import (
    INDEPENDENT_MEAN_CAVEAT,
)

if TYPE_CHECKING:
    from ..intermediates.results import MultiKSubsetResult

_ERROR_PREFIX = "[plot_subset_strip] Error"


def plot_subset_strip(
    result: MultiKSubsetResult,
    rater_names: list[str] | None = None,
    ax: plt.Axes | None = None,
) -> plt.Figure:
    """
    I plot multi-k subset results as a strip-and-trend chart.

    Each k value gets a column on the x-axis. Within each column, every
    subset's coefficient is plotted as a jittered dot, showing the spread
    across rater combinations at that panel size. The mean at each k is
    drawn as a larger diamond marker, connected across k values with a
    line to show how agreement scales with panel size. Bootstrap CIs of
    the mean appear as vertical whiskers at each k.

    When n_subsets at a given k is small (<=15), subset dots can be
    individually distinguished. When large, dots become semi-transparent
    to show density without overplotting.

    A k can be bootstrapped and still have no interval on its mean: every
    draw failed, so no method ran. That k gets no whisker, which alone
    would look like a k that was never bootstrapped, so its mean diamond
    is drawn hollow and the legend says so. Large k is where this is
    likeliest -- the more raters a subset demands, the fewer subjects
    reach two of them -- which is exactly the end of the trend line a
    reader is drawing conclusions from.

    Parameters
    ----------
    result : MultiKSubsetResult
        The result object from a multi-k subset permutation analysis.
    rater_names : list[str] | None
        Human-readable names for the raters. If None, raters are labeled
        by their 0-based index. Used for tooltip-style annotations when
        n_subsets is small. Length must match the number of distinct
        raters found in the subset keys.
    ax : plt.Axes | None
        Axes to draw on. If None, a new Figure and Axes are created.

    Returns
    -------
    plt.Figure
        The figure containing the plot.

    Raises
    ------
    TypeError
        If result is not a MultiKSubsetResult.
    ValueError
        If rater_names length does not match the number of raters.
    """
    # Lazy import to avoid circular dependency at module level
    from ..intermediates.results import MultiKSubsetResult

    if not isinstance(result, MultiKSubsetResult):
        raise TypeError(
            f"{_ERROR_PREFIX}: result must be MultiKSubsetResult, "
            f"got {type(result).__name__}"
        )

    # ------------------------------------------------------------------ #
    # Resolve rater names
    # ------------------------------------------------------------------ #

    # Collect all rater indices across all k results
    all_rater_indices = set()
    for k in result.ks:
        per_k = result[k]
        for subset in per_k.subset_values:
            all_rater_indices.update(subset)
    all_rater_indices = sorted(all_rater_indices)
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
    # Create figure
    # ------------------------------------------------------------------ #

    n_ks = len(result.ks)
    width = max(6, n_ks * 2.0)
    width = min(width, 16)

    if ax is None:
        fig, ax = plt.subplots(figsize=(width, 5))
    else:
        fig = ax.get_figure()

    # ------------------------------------------------------------------ #
    # Collect means and CIs for the trend line
    # ------------------------------------------------------------------ #

    x_positions = np.arange(n_ks)
    means = []
    ci_lowers = []
    ci_uppers = []
    has_any_ci = False

    # True where a bootstrap ran for that k and produced no interval, as
    # distinct from no bootstrap at all. Both leave NaN bounds and so draw
    # no whisker; only this tells them apart afterwards.
    mean_ci_undefined = []

    for k in result.ks:
        per_k = result[k]
        means.append(per_k.mean_over_subsets)
        mean_ci = per_k.bootstrap_ci_of_mean
        mean_ci_undefined.append(
            mean_ci is not None and not mean_ci.is_defined
        )
        if mean_ci is not None and mean_ci.is_defined:
            ci_lowers.append(mean_ci.lower)
            ci_uppers.append(mean_ci.upper)
            has_any_ci = True
        else:
            ci_lowers.append(np.nan)
            ci_uppers.append(np.nan)

    means = np.array(means)
    ci_lowers = np.array(ci_lowers)
    ci_uppers = np.array(ci_uppers)

    # ------------------------------------------------------------------ #
    # Draw strip dots per k
    # ------------------------------------------------------------------ #

    jitter_width = 0.2  # half-width of jitter band

    # Track all plotted values for y-axis bounds
    all_plotted_values = []

    for i, k in enumerate(result.ks):
        per_k = result[k]
        values = list(per_k.subset_values.values())
        subsets = list(per_k.subset_values.keys())
        n_subsets = len(values)

        # Filter NaN for plotting
        valid_mask = [not np.isnan(v) for v in values]
        valid_values = [v for v, m in zip(values, valid_mask) if m]
        valid_subsets = [s for s, m in zip(subsets, valid_mask) if m]
        all_plotted_values.extend(valid_values)

        if not valid_values:
            continue

        # Jitter x positions
        rng_jitter = np.random.default_rng(seed=42 + k)
        jitter = rng_jitter.uniform(
            -jitter_width, jitter_width, size=len(valid_values)
        )
        x_jittered = i + jitter

        # Adapt appearance to n_subsets
        if n_subsets <= 15:
            alpha = 0.7
            dot_size = 30
        elif n_subsets <= 50:
            alpha = 0.4
            dot_size = 20
        else:
            alpha = 0.2
            dot_size = 12

        ax.scatter(
            x_jittered, valid_values,
            s=dot_size, color="#888888", alpha=alpha,
            edgecolors="none", zorder=2,
        )

        # Label individual subsets when few enough to read
        if n_subsets <= 8:
            # Sort by value to stagger labels and avoid overlap
            indexed = sorted(
                enumerate(zip(x_jittered, valid_values, valid_subsets)),
                key=lambda t: t[1][1],  # sort by value
            )

            # Compute staggered y-offsets with minimum gap
            min_gap_pts = 12  # minimum gap in points
            label_y_offsets = _stagger_offsets(
                [val for _, (_, val, _) in indexed],
                min_gap_pts,
                ax,
            )

            for rank, (j, (x_j, val, subset)) in enumerate(indexed):
                label = ", ".join(name_map[idx] for idx in subset)
                y_off = label_y_offsets[rank]
                ax.annotate(
                    label,
                    (x_j, val),
                    textcoords="offset points",
                    xytext=(8, y_off),
                    fontsize=6,
                    color="dimgray",
                    zorder=4,
                )

    # ------------------------------------------------------------------ #
    # Draw mean CI whiskers
    # ------------------------------------------------------------------ #

    if has_any_ci:
        for i in range(n_ks):
            if not np.isnan(ci_lowers[i]):
                ax.plot(
                    [i, i], [ci_lowers[i], ci_uppers[i]],
                    color="#4682b4", linewidth=2.0, zorder=3,
                    solid_capstyle="round",
                )
                # Caps
                cap_w = 0.08
                ax.plot(
                    [i - cap_w, i + cap_w],
                    [ci_lowers[i], ci_lowers[i]],
                    color="#4682b4", linewidth=2.0, zorder=3,
                )
                ax.plot(
                    [i - cap_w, i + cap_w],
                    [ci_uppers[i], ci_uppers[i]],
                    color="#4682b4", linewidth=2.0, zorder=3,
                )
                all_plotted_values.extend(
                    [ci_lowers[i], ci_uppers[i]]
                )

    # ------------------------------------------------------------------ #
    # Draw mean trend line and diamonds
    # ------------------------------------------------------------------ #

    # Trend line connecting means
    if n_ks > 1:
        ax.plot(
            x_positions, means,
            color="black", linewidth=1.0, linestyle="-",
            zorder=4, alpha=0.5,
        )

    # Diamond markers for means. A hollow diamond is a mean whose interval
    # was attempted and came back undefined -- the mean itself is still a
    # real number, so the marker stays on the trend line rather than
    # vanishing from it.
    mean_faces = [
        "none" if undefined else "black" for undefined in mean_ci_undefined
    ]
    mean_edges = [
        "black" if undefined else "white" for undefined in mean_ci_undefined
    ]
    mean_widths = [
        1.2 if undefined else 0.8 for undefined in mean_ci_undefined
    ]

    ax.scatter(
        x_positions, means,
        s=80, c=mean_faces, marker="D",
        edgecolors=mean_edges, linewidths=mean_widths,
        zorder=5,
        label=f"Mean (n={n_ks} k values)",
    )

    # Annotate mean values
    for i, m in enumerate(means):
        ax.annotate(
            f"{m:.3f}",
            (i, m),
            textcoords="offset points",
            xytext=(0, 10),
            fontsize=8,
            ha="center",
            fontweight="bold",
            zorder=6,
        )

    # ------------------------------------------------------------------ #
    # X-axis
    # ------------------------------------------------------------------ #

    # Build x-tick labels with subset count
    tick_labels = []
    for k in result.ks:
        per_k = result[k]
        noun = "subset" if per_k.n_subsets == 1 else "subsets"
        tick_labels.append(f"k={k}\n({per_k.n_subsets} {noun})")

    ax.set_xticks(x_positions)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.set_xlabel("Subset size (k)", fontsize=10)

    # ------------------------------------------------------------------ #
    # Y-axis and title
    # ------------------------------------------------------------------ #

    ax.set_ylabel(result.coefficient_name, fontsize=10)

    # Determine bootstrap method (should be same across all k)
    methods = set()
    for k in result.ks:
        per_k = result[k]
        if per_k.bootstrap_mean_method is not None:
            methods.add(per_k.bootstrap_mean_method)
    method_str = ""
    if methods:
        method_str = f" ({', '.join(sorted(methods))} bootstrap)"

    ax.set_title(
        f"{result.coefficient_name} across subset sizes"
        f"{method_str}",
        fontsize=12,
    )

    # The whiskers carry no legend entry of their own; an independent mean
    # gets one, so the caveat sits next to the intervals it qualifies.
    if "independent" in methods and has_any_ci:
        ax.plot(
            [], [], color="#4682b4", linewidth=2.0,
            label="Mean CI: " + textwrap.fill(INDEPENDENT_MEAN_CAVEAT, 43),
        )

    # Y-axis bounds: default [0, 1], widen if data falls outside
    y_min = 0.0
    y_max = 1.0
    if all_plotted_values:
        data_min = min(all_plotted_values)
        data_max = max(all_plotted_values)
        if data_min < 0:
            y_min = data_min - 0.05
        if data_max > 1:
            y_max = data_max + 0.05
    ax.set_ylim(y_min, y_max)

    # Proxy entry so a hollow diamond is readable without the docstring
    if any(mean_ci_undefined):
        ax.scatter(
            [], [], s=80, marker="D", facecolors="none",
            edgecolors="black", linewidths=1.2,
            label="Hollow: no CI (every draw failed)",
        )

    ax.legend(fontsize=8, loc="best")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    return fig


# Private helpers ---------------------------------------------------------- #


def _stagger_offsets(
    sorted_values: list[float],
    min_gap_pts: float,
    ax: plt.Axes,
) -> list[float]:
    """
    I compute vertical label offsets that prevent overlap.

    Given y-values sorted ascending, I walk through them and push each
    label's offset so that consecutive labels are at least min_gap_pts
    apart in display coordinates. When values are well-separated, all
    offsets are 0 (labels sit next to their dot). When values cluster,
    labels fan out vertically.

    Parameters
    ----------
    sorted_values : list[float]
        Y-values in ascending order.
    min_gap_pts : float
        Minimum vertical distance between consecutive labels in points.
    ax : plt.Axes
        Axes, used to convert data coordinates to display coordinates.

    Returns
    -------
    list[float]
        Vertical offset in points for each label.
    """
    n = len(sorted_values)
    if n == 0:
        return []
    if n == 1:
        return [0.0]

    # Convert data y-values to display points for gap measurement
    # We need the transform from data to display coordinates
    fig = ax.get_figure()
    fig.canvas.draw()  # ensure transform is up to date

    display_ys = []
    for val in sorted_values:
        _, dy = ax.transData.transform((0, val))
        display_ys.append(dy)

    # Compute adjusted display positions with minimum gap
    adjusted = [display_ys[0]]
    for i in range(1, n):
        prev = adjusted[-1]
        curr = display_ys[i]
        if curr - prev < min_gap_pts:
            adjusted.append(prev + min_gap_pts)
        else:
            adjusted.append(curr)

    # Convert back to offsets in points relative to actual position
    offsets = []
    for i in range(n):
        offsets.append(adjusted[i] - display_ys[i])

    return offsets
