"""Visualization of a permutation test: the null distribution and the observed value."""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np

from .plot_bootstrap_distribution import _bins_for, _text_box

if TYPE_CHECKING:
    from ..intermediates.results import HypothesisTestResult

_ERROR_PREFIX = "[plot_null_distribution] Error"

# Display names, keyed by HypothesisTestResult.coefficient.
_COEFFICIENT_NAMES = {
    "fleiss_kappa": "Fleiss' kappa",
    "scotts_pi": "Scott's pi",
    "gwet_ac1": "Gwet's AC1",
    "brennan_prediger": "Brennan-Prediger",
    "krippendorff_alpha": "Krippendorff's alpha",
    "cohen_kappa": "Cohen's kappa",
    "conger_kappa": "Conger's kappa",
    "light_kappa": "Light's kappa",
}


def plot_null_distribution(
    result: HypothesisTestResult,
    ax: plt.Axes | None = None,
    bins: int | str = "auto",
) -> plt.Figure:
    """
    I plot a permutation test: what agreement looks like by chance, and yours.

    The histogram is the coefficient on every usable shuffle of the ratings
    -- the null distribution. A solid line marks the observed value. The
    shuffles counted toward the p-value are coloured: those at or above the
    observed value for 'greater', at or below for 'less', and for
    'two-sided' the smaller of the two tails, since the p-value is twice
    that tail.

    A box gives the p-value with the counts it was computed from, so the
    plot shows the arithmetic and not only its answer:
    (count + 1) / (usable + 1), doubled for two-sided and capped at 1
    (Phipson & Smyth 2010).

    Parameters
    ----------
    result : HypothesisTestResult
        From ``hypothesis_test()`` with method='permutation'.
    ax : plt.Axes | None
        Axes to draw on. If None, a new Figure and Axes are created.
    bins : int | str
        Passed to ``Axes.hist``. Default 'auto'; with few distinct values,
        each value gets its own bar.

    Returns
    -------
    plt.Figure
        The figure containing the plot.

    Raises
    ------
    TypeError
        If result is not a HypothesisTestResult.
    ValueError
        If the test was parametric_null, which has a formula rather than
        shuffles, or the result does not record the counts its p-value
        came from.
    """
    from ..intermediates.results import HypothesisTestResult

    if not isinstance(result, HypothesisTestResult):
        raise TypeError(
            f"{_ERROR_PREFIX}: result must be a HypothesisTestResult from "
            f"hypothesis_test(), got {type(result).__name__}"
        )
    if result.method != "permutation" or result.null_distribution is None:
        raise ValueError(
            f"{_ERROR_PREFIX}: only a permutation test has a null "
            f"distribution to plot. A parametric_null test uses a formula; "
            f"its z and p are in the result."
        )
    if result.n_usable is None or result.alternative is None:
        raise ValueError(
            f"{_ERROR_PREFIX}: this result does not record the counts its "
            f"p-value came from (n_usable, n_at_or_above, n_at_or_below, "
            f"alternative), so the plot cannot show them."
        )

    label = _COEFFICIENT_NAMES.get(result.coefficient, result.coefficient)
    draws = np.asarray(result.null_distribution, dtype=float)
    usable = draws[np.isfinite(draws)]

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4.5))
    else:
        fig = ax.figure

    # ------------------------------------------------------------------ #
    # Which tail counted toward p
    # ------------------------------------------------------------------ #

    n = result.n_usable
    greater = (result.n_at_or_above + 1) / (n + 1)
    less = (result.n_at_or_below + 1) / (n + 1)

    # Ties within the same tolerance Analyzer._permutation_p_value counts
    # with, so the coloured bars add up to the count in the box.
    tol = 100 * np.finfo(float).eps * max(abs(result.point), 1.0)
    at_or_above = usable >= result.point - tol
    at_or_below = usable <= result.point + tol
    if result.alternative == "greater" or (
        result.alternative == "two-sided" and greater <= less
    ):
        tail, side = at_or_above, "above"
    else:
        tail, side = at_or_below, "below"

    # ------------------------------------------------------------------ #
    # Histogram, split into the counted tail and the rest
    # ------------------------------------------------------------------ #

    bins, _ = _bins_for(usable, bins)
    edges = np.histogram_bin_edges(usable, bins=bins)
    ax.hist(
        [usable[~tail], usable[tail]], bins=edges, stacked=True,
        color=["lightgray", "firebrick"], edgecolor="gray", linewidth=0.5,
        label=["shuffles", f"shuffles at or {side} the observed value"],
        zorder=1,
    )
    ax.axvline(result.point, color="black", linewidth=1.5, zorder=3,
               label=f"observed {result.point:.3f}")

    # ------------------------------------------------------------------ #
    # The p-value and its arithmetic
    # ------------------------------------------------------------------ #

    count = result.n_at_or_above if side == "above" else result.n_at_or_below
    formula = f"({count:,} + 1) / ({n:,} + 1)"
    if result.alternative == "two-sided":
        formula = f"2 x {formula}"
        if 2 * min(greater, less) > 1:
            formula += ", capped at 1"
    decision = "reject" if result.reject_null else "fail to reject"
    _text_box(
        ax,
        f"{result.alternative} p = {result.p_value:.4f}\n"
        f"  = {formula}\n"
        f"{decision} at alpha = {result.alpha}",
    )

    # ------------------------------------------------------------------ #
    # Labels
    # ------------------------------------------------------------------ #

    caption = f"{n:,} of {result.n_permutations:,} shuffles usable"
    ax.set_title(f"{label}: permutation test", fontsize=12)
    ax.set_xlabel(f"{label}\n{caption}", fontsize=10)
    ax.set_ylabel("Shuffles", fontsize=10)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.45)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()

    return fig
