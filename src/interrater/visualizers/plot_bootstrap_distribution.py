"""Visualization of one coefficient's bootstrap distribution and its interval."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    from ..intermediates.results import PointEstimateResult

_ERROR_PREFIX = "[plot_bootstrap_distribution] Error"

# Display names for the result classes; anything else falls back to the
# class name, so a new coefficient still plots.
_COEFFICIENT_NAMES = {
    "FleissKappaResult": "Fleiss' kappa",
    "ScottsPiResult": "Scott's pi",
    "GwetsAC1Result": "Gwet's AC1",
    "BrennanPredigerResult": "Brennan-Prediger",
    "KrippendorffAlphaResult": "Krippendorff's alpha",
    "CohenKappaResult": "Cohen's kappa",
    "CongerKappaResult": "Conger's kappa",
    "LightKappaResult": "Light's kappa",
}

_METHOD_NAMES = {"bca": "BCa", "percentile": "percentile"}

# Up to this many distinct values, the draws get one bar per value.
_MAX_DISCRETE_VALUES = 60


def plot_bootstrap_distribution(
    result: PointEstimateResult,
    ax: plt.Axes | None = None,
    bins: int | str = "auto",
) -> plt.Figure:
    """
    I plot a coefficient's bootstrap distribution with its interval.

    The histogram is the usable draws; draws whose coefficient was undefined
    are left out and counted in the caption, as the interval counted them.
    A solid line marks the point estimate and the shaded band the interval,
    labelled with the method that actually produced it -- a BCa that fell
    back says 'percentile', because it is one.

    When the interval is BCa I also show what BCa did. The plain percentile
    interval from the same draws is drawn as a dashed outline, so the shift
    BCa made is visible, and a box gives its two corrections in words:

    - z0, the bias correction: whether the draws sit centred on the point
      estimate. It is read from the share of draws below the point; 0 means
      half, positive means more than half, so the interval moves up.
    - a, the acceleration: one sixth of the skewness of the
      leave-one-subject-out values (Efron & Tibshirani 1986). It says
      whether the coefficient's spread changes with the coefficient;
      positive leans the interval toward higher values, negative toward
      lower.

    Parameters
    ----------
    result : PointEstimateResult
        A coefficient's result from ``estimate()``, with a bootstrap
        interval (ci_method 'percentile' or 'bca').
    ax : plt.Axes | None
        Axes to draw on. If None, a new Figure and Axes are created.
    bins : int | str
        Passed to ``Axes.hist``. Default 'auto'.

    Returns
    -------
    plt.Figure
        The figure containing the plot.

    Raises
    ------
    TypeError
        If result is not a PointEstimateResult.
    ValueError
        If the result has no bootstrap draws to plot: a parametric_null
        interval has none, and an interval where every draw failed has
        nothing to show.
    """
    from ..intermediates.results import PointEstimateResult
    from ..intermediates.results.ci_results import BootstrapCIResult

    if not isinstance(result, PointEstimateResult):
        raise TypeError(
            f"{_ERROR_PREFIX}: result must be a coefficient's result from "
            f"estimate(), got {type(result).__name__}"
        )

    ci = result.confidence_interval
    if not isinstance(ci, BootstrapCIResult) or result.bootstrap_distribution is None:
        raise ValueError(
            f"{_ERROR_PREFIX}: this result has no bootstrap distribution. "
            f"A parametric_null interval is computed from a formula, not "
            f"from draws; use ci_method='percentile' or 'bca' to get one."
        )

    draws = np.asarray(result.bootstrap_distribution, dtype=float)
    usable = draws[np.isfinite(draws)]
    if usable.size == 0:
        raise ValueError(
            f"{_ERROR_PREFIX}: every one of the {draws.size} bootstrap draws "
            f"was undefined, so there is no distribution to plot."
        )

    name = _COEFFICIENT_NAMES.get(type(result).__name__, type(result).__name__)

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4.5))
    else:
        fig = ax.figure

    # ------------------------------------------------------------------ #
    # Draws, point, interval
    # ------------------------------------------------------------------ #

    bins, discrete = _bins_for(usable, bins)
    distinct = _distinct(usable)

    ax.hist(usable, bins=bins, color="lightgray", edgecolor="gray",
            linewidth=0.5, zorder=1)

    level = f"{ci.confidence_level:.0%}"
    if ci.method_used is not None:
        ax.axvspan(
            ci.lower, ci.upper, alpha=0.15, color="steelblue", zorder=0,
            label=f"{level} {_METHOD_NAMES[ci.method_used]} interval "
                  f"[{ci.lower:.3f}, {ci.upper:.3f}]",
        )

    ax.axvline(result.point, color="black", linewidth=1.5, zorder=3,
               label=f"point estimate {result.point:.3f}")

    if ci.method_used == "bca":
        _draw_bca_view(ax, usable, ci)
    elif ci.method_intended == "bca":
        # Asked for BCa, got percentile: say so where the corrections would be.
        _text_box(ax, "BCa could not be applied;\nthis is the percentile interval.")

    # ------------------------------------------------------------------ #
    # Labels
    # ------------------------------------------------------------------ #

    caption = f"{ci.n_usable:,} of {ci.n_resamples:,} draws usable"
    if ci.below_recommended_resamples:
        caption += (
            f" -- below the {ci.recommended_resamples:,} recommended "
            f"for {_METHOD_NAMES[ci.method_used]}"
        )
    if discrete:
        at_point = float(np.mean(np.round(usable, 10) == round(result.point, 10)))
        caption += (
            f"\n{distinct.size} distinct values; {at_point:.1%} of draws "
            f"equal the point estimate"
        )

    ax.set_title(f"{name}: bootstrap distribution", fontsize=12)
    ax.set_xlabel(f"{name}\n{caption}", fontsize=10)
    ax.set_ylabel("Draws", fontsize=10)
    # Headroom, so the legend and the BCa box sit above the bars.
    ax.set_ylim(top=ax.get_ylim()[1] * 1.45)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    return fig


# Helpers --------------------------------------------------------------- #


def _distinct(values: np.ndarray) -> np.ndarray:
    """
    The distinct values, treating floats equal to 10 decimals as one.

    The same coefficient reached through different tables can differ in the
    last bit: on 2,000 shuffles of SR1 prelim, 105 'distinct' Fleiss values
    were 51 real ones.
    """
    return np.unique(np.round(values, 10))


def _bins_for(values: np.ndarray, bins: int | str) -> tuple:
    """
    Histogram bins, and whether the values were treated as discrete.

    A coefficient on small data takes few distinct values -- Brennan-
    Prediger depends only on raw agreement -- and automatic bins then leave
    empty gaps between them. With bins='auto' and few distinct values, each
    value gets its own bar instead. Shared with plot_null_distribution.
    """
    distinct = _distinct(values)
    if not (isinstance(bins, str) and bins == "auto"
            and 1 < distinct.size <= _MAX_DISCRETE_VALUES):
        return bins, False
    half_gap = np.min(np.diff(distinct)) / 2
    edges = np.concatenate((
        [distinct[0] - half_gap],
        (distinct[:-1] + distinct[1:]) / 2,
        [distinct[-1] + half_gap],
    ))
    return edges, True


def _draw_bca_view(ax: plt.Axes, usable: np.ndarray, ci) -> None:
    """The percentile interval for comparison, and BCa's corrections in words."""
    alpha = ci.alpha
    p_lower = float(np.percentile(usable, (alpha / 2) * 100))
    p_upper = float(np.percentile(usable, (1 - alpha / 2) * 100))
    for x in (p_lower, p_upper):
        ax.axvline(x, color="steelblue", linestyle="--", linewidth=1.0, zorder=2)
    # One legend entry for the pair of dashed lines.
    ax.plot([], [], color="steelblue", linestyle="--", linewidth=1.0,
            label=f"percentile interval, for comparison "
                  f"[{p_lower:.3f}, {p_upper:.3f}]")

    from scipy.stats import norm

    z0 = ci.bias_correction
    a = ci.acceleration
    share_below = float(norm.cdf(z0))
    if abs(a) < 1e-3:
        lean = "no lean"
    elif a > 0:
        lean = "leans toward higher values"
    else:
        lean = "leans toward lower values"
    _text_box(
        ax,
        f"BCa corrections\n"
        f"z0 = {z0:+.3f}  ({share_below:.0%} of draws below the point)\n"
        f"a  = {a:+.4f}  ({lean})",
    )


def _text_box(ax: plt.Axes, text: str) -> None:
    ax.text(
        0.98, 0.97, text, transform=ax.transAxes, ha="right", va="top",
        fontsize=8, family="monospace",
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray", alpha=0.9),
        zorder=4,
    )
