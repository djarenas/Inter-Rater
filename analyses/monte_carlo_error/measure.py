"""
What governs the Monte Carlo error of a bootstrap percentile interval?

The companion analysis (analyses/seed_stability/) measured, on six real
datasets, that the standard deviation of each bound of a 95% percentile
interval behaves as

    SD  ~=  0.69 * width / sqrt(B)

That constant is *not* a property of percentile intervals. It is a
property of percentile intervals taken from an approximately normal
bootstrap distribution, which is what those six kappa datasets happened
to produce. Published as a general rule it would be wrong, in one case
below by a factor of forty.

What does generalise is the standard error of a sample quantile,

    SD of the p-quantile  =  sqrt( p (1 - p) / B )  /  f(q_p)

where f is the density of the bootstrap distribution at that quantile.
Everything follows from it: short tails mean a high density at the
quantile and therefore a small error, while a skewed distribution has a
dense end and a sparse end and so its two bounds are not equally
reproducible.

This script demonstrates that on distributions whose densities are known
in closed form, so the prediction can be checked rather than asserted.
For each shape it reports the simulated SD of both bounds, the SD the
formula predicts, and their ratio. A ratio near 1 everywhere is the
result; the accompanying spread in the normalised constant c is what the
0.69 shorthand would have gotten wrong.

**It deliberately does not import interrater.** The claim is about the
percentile method, not about this implementation, so nothing here depends
on the package being correct -- which is what makes it usable as evidence
for the package.

Run from the repository root:

    python analyses/monte_carlo_error/measure.py

About a minute. Dependencies: numpy, scipy, pandas.
"""

from __future__ import annotations

import csv
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats

# -- Settings --------------------------------------------------------------- #

ALPHA = 0.05            # 95% interval, so bounds at the 2.5th and 97.5th
P = ALPHA / 2
B_VALUES = [1000, 4000]  # two values, to show the 1/sqrt(B) rate also holds
                         # when the distribution is not normal
REPS = 4000              # independent bootstrap distributions per cell
SEED = 20261018

HERE = Path(__file__).resolve().parent
OUT = HERE / "results.csv"

# Frozen scipy distributions, so the density at the quantile is exact
# rather than estimated. These stand in for bootstrap distributions of
# differing shape; none is claimed to arise from a real estimator.
SHAPES = [
    ("normal",                   stats.norm()),
    ("t(8), heavy tails",        stats.t(8)),
    ("uniform, short tails",     stats.uniform()),
    ("beta(1.4, 8), skewed",     stats.beta(1.4, 8)),
    ("exponential",              stats.expon()),
    ("lognormal(0, 0.6)",        stats.lognorm(0.6)),
]


# -- The prediction --------------------------------------------------------- #

def predicted_sd(dist, p: float, B: int) -> float:
    """
    The standard error of the p-quantile of B draws from ``dist``.

    sqrt(p(1-p)/B) / f(q_p) -- the classical asymptotic standard error of
    a sample quantile. Nothing here is specific to bootstrapping; a
    bootstrap distribution is just the distribution being sampled.
    """
    return float(np.sqrt(p * (1 - p) / B) / dist.pdf(dist.ppf(p)))


def normal_constant(p: float) -> float:
    """
    The 0.6815 that the normal case reduces to.

    For a normal with sd sigma, f(q_p) = phi(z_p)/sigma and the interval
    width is 2|z_p|sigma, so sigma cancels and the standard error becomes
    a pure number times width/sqrt(B).
    """
    z = stats.norm.ppf(p)
    return float(np.sqrt(p * (1 - p)) / stats.norm.pdf(z) / (2 * abs(z)))


# -- Main ------------------------------------------------------------------- #

def main() -> int:
    rng = np.random.default_rng(SEED)
    c_normal = normal_constant(P)

    print(f"95% percentile interval, {REPS} replications per cell, "
          f"seed {SEED}")
    print(f"normal-case constant from theory: {c_normal:.4f}\n")

    rows = []
    for B in B_VALUES:
        print(f"{'=' * 94}\nB = {B}\n{'=' * 94}")
        print(f"{'distribution':<24} {'skew':>6} "
              f"{'sim lo':>8} {'pred lo':>8} {'ratio':>6}   "
              f"{'sim up':>8} {'pred up':>8} {'ratio':>6}   "
              f"{'c lo':>6} {'c up':>6}")
        print("-" * 94)

        for name, dist in SHAPES:
            # Truth, from the closed-form density
            pred_lo = predicted_sd(dist, P, B)
            pred_up = predicted_sd(dist, 1 - P, B)
            width_true = float(dist.ppf(1 - P) - dist.ppf(P))
            skew = float(dist.stats("s"))

            # Simulation: REPS independent "bootstrap distributions"
            draws = dist.rvs(size=(REPS, B), random_state=rng)
            los, ups = np.quantile(draws, [P, 1 - P], axis=1)

            sim_lo = float(los.std(ddof=1))
            sim_up = float(ups.std(ddof=1))
            width_sim = float(ups.mean() - los.mean())

            c_lo = sim_lo * np.sqrt(B) / width_sim
            c_up = sim_up * np.sqrt(B) / width_sim

            print(f"{name:<24} {skew:>6.2f} "
                  f"{sim_lo:>8.5f} {pred_lo:>8.5f} {sim_lo/pred_lo:>6.3f}   "
                  f"{sim_up:>8.5f} {pred_up:>8.5f} {sim_up/pred_up:>6.3f}   "
                  f"{c_lo:>6.3f} {c_up:>6.3f}", flush=True)

            rows.append({
                "distribution":   name,
                "skew":           skew,
                "B":              B,
                "reps":           REPS,
                "width_true":     width_true,
                "width_simulated": width_sim,
                "sd_lower_simulated":  sim_lo,
                "sd_lower_predicted":  pred_lo,
                "sd_upper_simulated":  sim_up,
                "sd_upper_predicted":  pred_up,
                "c_lower":        c_lo,
                "c_upper":        c_up,
                "c_normal_theory": c_normal,
            })
        print()

    df = pd.DataFrame(rows)
    ratios = pd.concat([
        df.sd_lower_simulated / df.sd_lower_predicted,
        df.sd_upper_simulated / df.sd_upper_predicted,
    ])
    cs = pd.concat([df.c_lower, df.c_upper])

    print(f"{'=' * 94}")
    print(f"simulated / predicted, over all {len(ratios)} bound-cells:")
    print(f"    mean {ratios.mean():.4f}   sd {ratios.std():.4f}   "
          f"min {ratios.min():.4f}   max {ratios.max():.4f}")
    print(f"the normalised constant c, over the same cells:")
    print(f"    min {cs.min():.3f}   max {cs.max():.3f}   "
          f"ratio max/min {cs.max()/cs.min():.1f}x")
    print()
    print("Read together: the quantile formula predicts every shape, while a")
    print("single constant spans the range above and so describes none of")
    print("them except the normal.")

    with OUT.open("w", newline="") as fh:
        for line in (
            f"# numpy {np.__version__}",
            f"# scipy {scipy.__version__}",
            f"# pandas {pd.__version__}",
            f"# python {platform.python_version()} on {platform.system()}",
            f"# generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
            f"# alpha {ALPHA}, p {P}, B values {B_VALUES}, {REPS} reps, seed {SEED}",
            f"# normal-case constant from theory {c_normal:.6f}",
            f"# predicted SD = sqrt(p(1-p)/B) / f(q_p), f from the closed-form density",
            f"# does not import interrater: this is about the percentile method",
            f"# produced by analyses/monte_carlo_error/measure.py",
        ):
            fh.write(line + "\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\nwrote {len(rows)} rows to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
