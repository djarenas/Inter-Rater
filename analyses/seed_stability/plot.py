"""
Figures for the seed-stability analysis.

Reads results.csv and writes three figures, each with two panels -- percentile
on the left, BCa on the right, so every figure reads the same way:

  fig1_sd_lower   SD of the lower bound against B, with tolerance lines
  fig2_sd_upper   the same for the upper bound
  fig3_collapse   SD normalised by interval width, with the fitted c/sqrt(B)

Axes are log-log throughout. On linear axes the whole curve crushes into the
left corner; on log-log the 1/sqrt(B) law is a straight line, which is the
thing a reader is being asked to check. The x ticks are the actual resample
counts rather than powers of ten, because those are the values someone is
choosing between.

Colour identifies the dataset and marker shape repeats that identification,
so the figures survive greyscale printing and colour-vision deficiency. The
palette is the six-colour Okabe-Ito set, which passes the all-pairs checks in
scripts/validate_palette.js; its worst CVD pair sits in the 6-8 band, which is
admissible only with a secondary encoding, hence the markers.

Nothing here recomputes a measurement: every number is read from the CSV, and
the only arithmetic is the standard deviation across seeds and the fitted
constant.

Run from the repository root, after measure.py:

    python analyses/seed_stability/plot.py

Dependencies: numpy, pandas, matplotlib, all three runtime dependencies of
the package, so a working install already has them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")           # no display needed; write files only
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CSV = HERE / "results.csv"
FIGDIR = HERE / "figures"

TOLERANCES = [0.005, 0.01]
METHODS = ["percentile", "bca"]
METHOD_LABEL = {"percentile": "Percentile", "bca": "BCa"}

# Okabe-Ito, six slots. Validated on the all-pairs list; yellow and black are
# left out (yellow fails the contrast floor on a light surface).
COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00"]
MARKERS = ["o", "s", "^", "D", "v", "X"]


def load() -> pd.DataFrame:
    """One row per (dataset, method, B, bound): the SD across seeds."""
    if not CSV.is_file():
        print(f"error: {CSV} not found. Run measure.py first.", file=sys.stderr)
        raise SystemExit(1)

    df = pd.read_csv(CSV, comment="#")

    rows = []
    for (ds, method, B), g in df.groupby(["dataset", "ci_method", "B"]):
        width = float(g.width.mean())
        for bound in ("lower", "upper"):
            rows.append({
                "dataset": ds,
                "method":  method,
                "B":       int(B),
                "bound":   bound,
                "sd":      float(g[bound].std(ddof=1)),
                "width":   width,
                "n_seeds": len(g),
            })
    return pd.DataFrame(rows)


def style(ax, xs: list[int]) -> None:
    """Log-log axes with recessive furniture and literal x ticks."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels([str(x) for x in xs])
    ax.minorticks_off()
    ax.grid(True, which="major", lw=0.4, alpha=0.3, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
        ax.spines[side].set_alpha(0.5)
    ax.tick_params(labelsize=8, length=3, width=0.6)


def fitted_constant(sub: pd.DataFrame) -> float:
    """c in SD = c * width / sqrt(B), fitted with the exponent held at -1/2."""
    return float(np.exp(np.mean(np.log(sub.sd / sub.width) + 0.5 * np.log(sub.B))))


def width_labels(data: pd.DataFrame, xs: list[int]) -> dict[str, str]:
    """
    Legend labels carrying each dataset's interval width.

    The width is what orders the lines in the un-normalised figures, so
    printing it lets a reader see that before the normalised figure proves
    it. One number serves both panels: percentile and BCa widths differ by
    at most 0.3% on these datasets, and the drift from B=250 to B=4000 is
    1-2%, so the value quoted is the most converged one -- percentile at
    the largest B.
    """
    at_max = data[(data.method == "percentile") & (data.B == xs[-1])]
    width = at_max.groupby("dataset").width.first()
    return {ds: f"{ds}  (w {width[ds]:.3f})" for ds in width.index}


def bound_figure(data: pd.DataFrame, bound: str, order: list[str],
                 xs: list[int], outstem: str, n_seeds: int) -> None:
    """SD of one bound against B, one panel per method."""
    sub = data[data.bound == bound]
    labels = width_labels(data, xs)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)

    for ax, method in zip(axes, METHODS):
        m = sub[sub.method == method]
        for i, ds in enumerate(order):
            s = m[m.dataset == ds].sort_values("B")
            ax.plot(s.B, s.sd, color=COLOURS[i], marker=MARKERS[i],
                    markersize=5.5, lw=1.6, label=labels[ds], zorder=3,
                    markeredgecolor="white", markeredgewidth=0.6)

        for tol in TOLERANCES:
            ax.axhline(tol, color="#444444", lw=0.9, ls=(0, (4, 3)),
                       alpha=0.65, zorder=2)
            ax.annotate(f"SD = {tol}", xy=(xs[-1], tol),
                        xytext=(-2, 3), textcoords="offset points",
                        ha="right", va="bottom", fontsize=7.5,
                        color="#444444")

        style(ax, xs)
        ax.set_title(METHOD_LABEL[method], fontsize=10.5, pad=8)
        ax.set_xlabel("resamples $B$", fontsize=9)

    axes[0].set_ylabel(f"SD of the {bound} bound across seeds", fontsize=9)

    # Spanning barely one decade, a log axis prints two labels and nothing
    # between them, so the ticks are named explicitly.
    yt = [0.001, 0.002, 0.005, 0.01, 0.02]
    axes[0].set_yticks(yt)
    axes[0].set_yticklabels([f"{v:g}" for v in yt])

    axes[1].legend(fontsize=8, frameon=False, loc="lower left",
                   bbox_to_anchor=(1.02, 0.0), handlelength=2.4)

    fig.suptitle(
        f"Seed-to-seed variability of the {bound} 95% bound "
        f"(Fleiss' kappa, {n_seeds} seeds per point)",
        fontsize=11, y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 0.84, 0.94))
    save(fig, outstem)


def collapse_figure(data: pd.DataFrame, order: list[str], xs: list[int],
                    outstem: str, n_seeds: int) -> None:
    """SD / width against B, with the fitted law drawn through it."""
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)
    grid = np.logspace(np.log10(xs[0] * 0.92), np.log10(xs[-1] * 1.08), 100)
    line_legends = []

    for ax, method in zip(axes, METHODS):
        m = data[data.method == method]
        c = fitted_constant(m)

        fit_line, = ax.plot(grid, c / np.sqrt(grid), color="#444444", lw=1.1,
                            ls=(0, (5, 3)), zorder=2,
                            label=f"${c:.3f}/\\sqrt{{B}}$")

        for i, ds in enumerate(order):
            s = m[m.dataset == ds]
            for bound, filled in (("lower", True), ("upper", False)):
                t = s[s.bound == bound].sort_values("B")
                ax.plot(t.B, t.sd / t.width, color=COLOURS[i],
                        marker=MARKERS[i], markersize=5.0, lw=0,
                        zorder=3, markeredgecolor=COLOURS[i],
                        markeredgewidth=0.9,
                        markerfacecolor=COLOURS[i] if filled else "white",
                        label=ds if bound == "lower" else None)

        style(ax, xs)
        ax.set_title(f"{METHOD_LABEL[method]}  —  $c$ = {c:.3f}",
                     fontsize=10.5, pad=8)
        ax.set_xlabel("resamples $B$", fontsize=9)
        # Only the fitted line goes in this legend. Left to itself, matplotlib
        # collects every labelled artist and repeats all six datasets here,
        # on top of the data, duplicating the panel legend on the right.
        line_legends.append(
            ax.legend(handles=[fit_line], fontsize=8.5, frameon=False,
                      loc="upper right")
        )

    axes[0].set_ylabel("SD of bound  /  interval width", fontsize=9)

    # A log axis this narrow prints one decade label and nothing else, so the
    # ticks are named explicitly.
    yt = [0.01, 0.015, 0.02, 0.03, 0.05]
    axes[0].set_yticks(yt)
    axes[0].set_yticklabels([f"{v:g}" for v in yt])

    handles = [
        plt.Line2D([], [], color=COLOURS[i], marker=MARKERS[i], lw=0,
                   markersize=5.0, label=ds, markeredgecolor=COLOURS[i])
        for i, ds in enumerate(order)
    ] + [
        plt.Line2D([], [], color="#444444", marker="o", lw=0, markersize=5.0,
                   markerfacecolor="#444444", label="lower bound"),
        plt.Line2D([], [], color="#444444", marker="o", lw=0, markersize=5.0,
                   markerfacecolor="white", label="upper bound"),
    ]
    # Keep the fitted-line legend: a second legend() call on the same axes
    # replaces the first unless it is re-added as an artist.
    axes[1].add_artist(line_legends[1])
    axes[1].legend(handles=handles, fontsize=8, frameon=False,
                   loc="lower left", bbox_to_anchor=(1.02, 0.0),
                   handlelength=1.6)

    fig.suptitle(
        f"Normalised by interval width, the datasets collapse onto one line "
        f"per method ({n_seeds} seeds per point)",
        fontsize=11, y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 0.84, 0.94))
    save(fig, outstem)


def save(fig, stem: str) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 200})):
        path = FIGDIR / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", **kw)
        print(f"  wrote {path.relative_to(HERE.parent.parent)}")
    plt.close(fig)


def main() -> int:
    data = load()
    xs = sorted(int(b) for b in data.B.unique())
    n_seeds = int(data.n_seeds.iloc[0])

    # Fixed dataset order, so colour and marker mean the same thing in every
    # figure. Sorted by how noisy the dataset is, which makes the legend read
    # in the same order as the lines.
    order = (
        data[(data.method == "percentile") & (data.B == xs[-1])]
        .groupby("dataset").sd.max().sort_values().index.tolist()
    )

    print(f"{len(data)} cells from {CSV.name}; "
          f"B = {xs}; {n_seeds} seeds per cell")
    print(f"datasets, least to most noisy: {', '.join(order)}\n")

    bound_figure(data, "lower", order, xs, "fig1_sd_lower", n_seeds)
    bound_figure(data, "upper", order, xs, "fig2_sd_upper", n_seeds)
    collapse_figure(data, order, xs, "fig3_collapse", n_seeds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
