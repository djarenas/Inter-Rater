"""
How reproducible is a bootstrap confidence interval?

Everything about a bootstrap interval is fixed by the data except which
resamples get drawn. Fix the data, fix B, fix the method, vary only the
seed, and whatever movement remains is Monte Carlo noise -- the part of the
published interval that is an artefact of the draw rather than a property
of the ratings.

This script measures that movement. For each dataset, each method and each
B it builds the same 95% interval under 100 different seeds and records
every bound, so the standard deviation of the lower limit and of the upper
limit can be computed from the saved values rather than taken on trust.

Why it matters: B = 1000 is near-universal in applied work, and it comes
from one sentence in Efron & Tibshirani (1986) section 9 -- "B = 1000 is a
rough minimum" -- rather than from any analysis of a particular dataset.
That is a reasonable default and a poor justification. The output here lets
a reader pick a tolerance they can defend and read off the B that meets it
on data resembling their own.

**Both CI methods are swept, because they are not interchangeable here.**
The same passage in Efron & Tibshirani gives 250 for percentile and 1000
for BCa, the stated reason being that percentile need not compute z0. BCa
therefore estimates an extra quantity from the same draws, and takes its
quantiles at levels that shift with the seed rather than at fixed 2.5 and
97.5. Numbers measured on one method should not be quoted for the other.

Coefficient: Fleiss' kappa, the group coefficient these datasets support
(3-4 raters each). balanced_data_only=False so the two unbalanced datasets
are included; on those, Fleiss' kappa is the approximated form, and the
`balanced` column records which is which.

Output: results.csv, one row per (dataset, ci_method, B, seed). Aggregation
is left to the plotting script, so a different statistic can be computed
later without re-running anything. `ci_method` is what was asked for and
`method_used` is what ran -- BCa can fall back to percentile, and a
fallback must not be averaged in with the method that was requested.

Run from the repository root:

    python analyses/seed_stability/measure.py

Roughly 45 minutes. Dependencies: the package, numpy, pandas. No matplotlib.
"""

from __future__ import annotations

import csv
import platform
import sys
import time
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import interrater
from interrater import (
    CleanerConfig, FleissKappaAnalyzer, RatingsCleaner, RatingsSemantics,
)
from interrater.analyzers.configs import FleissKappaConfig


# -- What gets measured ----------------------------------------------------- #

CI_METHODS = ["percentile", "bca"]

B_VALUES = [250, 500, 1000, 2000, 4000]

# 100 seeds. Written out rather than generated, because a reproducibility
# study whose own seed list is itself generated has merely moved the
# question somewhere less visible.
SEEDS = [
    1, 2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
    71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
    113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
    173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
    229, 233, 239, 241, 251, 257, 263, 269, 271, 277,
    281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
    349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
    409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
    463, 467, 479, 487, 491, 499, 503, 509, 521, 523,
]

ALPHA = 0.05

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
EXAMPLE = REPO / "example_data"
OUT = HERE / "results.csv"

# (label, review, folder, data file, semantics yaml, cleaner yaml)
DATASETS = [
    ("SR1 prelim", "SR1", "SR1_OSF2018", "SR1_prelim_data.csv",
     "SR1_prelim_semantics.yaml", "SR1_prelim_cleaner_config.yaml"),
    ("SR1 main", "SR1", "SR1_OSF2018", "SR1_main_data.txt",
     "SR1_main_semantics.yaml", "SR1_main_cleaner_config.yaml"),
    ("SR2", "SR2", "SR2_JSRC2019", "SR2_data.txt",
     "SR2_semantics.yaml", "SR2_cleaner_config.yaml"),
    ("SR3 FIS Psych", "SR3", "SR3_JGIM2019", "SR_FIS_Psych_data.txt",
     "SR_FIS_Psych_semantics.yaml", "SR_FIS_Psych_cleaner_config.yaml"),
    ("SR4 FIS CVD", "SR4", "SR4_FIS_CVD", "SR_FIS_CVD_data.txt",
     "SR_FIS_CVD_semantics.yaml", "SR_FIS_CVD_cleaner_config.yaml"),
    ("SR5 Coc Cardio", "SR5", "SR5_CocCardio", "SR_cocaine_cardio_data.txt",
     "SR_cocaine_cardio_semantics.yaml", "SR_cocaine_cardio_cleaner_config.yaml"),
]


# -- Loading ---------------------------------------------------------------- #
#
# A copy of the loader in tests/test_example_data/conftest.py. Duplicated
# deliberately: a script should not import from the test suite, and the
# alternative -- shipping a loader in the package that only these scripts
# use -- would put it in every user's install for no reason.
#
# What the two copies share is not really code. It is one fact about the
# data folder: .csv files ship with a header, .txt files are tab-separated
# and headerless. That convention is written down in example_data/README.md
# under "Folder layout", which is the file to change if a dataset ever
# arrives in a third format -- and then both copies of this function.

def load_ratings(folder: str, data: str, semantics: str, cleaner: str):
    """Load a shipped dataset through the cleaning pipeline."""
    d = EXAMPLE / folder
    sem = RatingsSemantics.build_from_yaml(str(d / semantics))
    cfg = CleanerConfig.build_from_yaml(str(d / cleaner))

    path = d / data
    # .csv files ship with a header; .txt files do not.
    df = pd.read_csv(
        path,
        sep="," if path.suffix == ".csv" else "\t",
        header=0 if path.suffix == ".csv" else None,
        dtype=str,
    )
    return RatingsCleaner(df, sem, config=cfg).ratings


# -- Provenance ------------------------------------------------------------- #

def provenance() -> list[str]:
    """
    Lines recording what produced the numbers below them.

    The numpy version is here for a specific reason: numpy does not promise
    that ``np.random.Generator`` yields the same stream across versions --
    that guarantee covers only the legacy ``RandomState``. So these exact
    bounds are reproducible under the numpy recorded here and need not be
    under another. The standard deviations are a different matter: they
    describe Monte Carlo noise at a given B, and that does not depend on
    which stream produced it.
    """
    return [
        f"# interrater {interrater.__version__}",
        f"# imported from {Path(interrater.__file__).resolve().parent}",
        f"# numpy {np.__version__}",
        f"# pandas {pd.__version__}",
        f"# python {platform.python_version()} on {platform.system()}",
        f"# generated {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"# coefficient fleiss_kappa, alpha {ALPHA}",
        f"# ci methods {CI_METHODS}",
        f"# balanced_data_only False",
        f"# B values {B_VALUES}",
        f"# {len(SEEDS)} seeds, first {SEEDS[0]} last {SEEDS[-1]}",
        f"# produced by analyses/seed_stability/measure.py",
    ]


# -- Main ------------------------------------------------------------------- #

def check_which_interrater() -> None:
    """
    Say which installed copy of the package is being measured.

    A bare ``import interrater`` finds whatever is on sys.path, which is
    not necessarily the source tree beside this script -- an old install,
    or a stray copy elsewhere, wins silently and the results then describe
    code nobody is looking at. That has happened during this project more
    than once, so the resolved path is printed every run and a mismatch is
    called out rather than left to be noticed.
    """
    imported = Path(interrater.__file__).resolve().parent
    expected = (REPO / "src" / "interrater").resolve()

    print(f"interrater {interrater.__version__} from {imported}")
    if imported != expected:
        print(
            f"\n  WARNING: that is not this repository's copy, which is at\n"
            f"           {expected}\n"
            f"  The numbers below will describe the imported copy above.\n"
            f"  To measure this repository instead:\n"
            f"           pip install -e .        (or)\n"
            f"           PYTHONPATH=src python {Path(__file__).relative_to(REPO)}\n",
            file=sys.stderr,
        )


def main() -> int:
    if not EXAMPLE.is_dir():
        print(f"error: {EXAMPLE} not found. Run from the repository root.",
              file=sys.stderr)
        return 1

    check_which_interrater()

    # The null-variance warning is not relevant here -- every interval this
    # script builds is a bootstrap one -- but the unbalanced datasets emit
    # their own notices, and thousands of them would bury the progress
    # output.
    warnings.simplefilter("ignore", UserWarning)

    print("loading datasets ...", flush=True)
    loaded = []
    for label, review, folder, data, sem, cln in DATASETS:
        ratings = load_ratings(folder, data, sem, cln)
        loaded.append((label, review, ratings))
        print(f"  {label:<16} {ratings.category_index.shape[0]:>5} subjects, "
              f"{ratings.category_index.shape[1]} raters", flush=True)

    total = len(loaded) * len(CI_METHODS) * len(B_VALUES) * len(SEEDS)
    print(f"\n{total} interval computations "
          f"({len(loaded)} datasets x {len(CI_METHODS)} methods x "
          f"{len(B_VALUES)} B values x {len(SEEDS)} seeds)\n", flush=True)

    rows = []
    fallbacks = Counter()
    started = time.perf_counter()
    done = 0

    for label, review, ratings in loaded:
        for method in CI_METHODS:
            for B in B_VALUES:
                t0 = time.perf_counter()
                points, lowers = set(), []

                for seed in SEEDS:
                    result = FleissKappaAnalyzer(
                        ratings,
                        FleissKappaConfig(
                            ci_method=method,
                            n_resamples=B,
                            alpha=ALPHA,
                            balanced_data_only=False,
                            random_state=seed,
                        ),
                    ).estimate()

                    ci = result.confidence_interval
                    points.add(round(float(result.point), 12))
                    lowers.append(ci.lower)

                    if ci.method_used != method:
                        fallbacks[(label, method, B, ci.method_used)] += 1

                    rows.append({
                        "dataset":     label,
                        "review":      review,
                        "balanced":    result.balance_profile.is_balanced,
                        "n_subjects":  result.balance_profile.n_subjects,
                        "n_raters":    result.balance_profile.n_raters,
                        "ci_method":   method,
                        "B":           B,
                        "seed":        seed,
                        "point":       result.point,
                        "lower":       ci.lower,
                        "upper":       ci.upper,
                        "width":       ci.upper - ci.lower,
                        "n_failed":    ci.n_failed,
                        "method_used": ci.method_used,
                        # BCa only. z0 is estimated from the draws and so
                        # moves with the seed; the acceleration comes from
                        # the jackknife and should not.
                        "acceleration":    ci.acceleration,
                        "bias_correction": ci.bias_correction,
                    })
                    done += 1

                # The point estimate does not depend on the seed. If it ever
                # does, a generator is leaking into the estimator and every
                # number in this file is suspect -- so it is checked here
                # rather than assumed.
                if len(points) != 1:
                    print(f"\nFATAL: point estimate varied with the seed for "
                          f"{label}, {method}, B={B}: {sorted(points)}\n"
                          f"This is a bug, not noise. Nothing written.",
                          file=sys.stderr)
                    return 2

                sd = float(np.std(lowers, ddof=1))
                elapsed = time.perf_counter() - t0
                print(f"  {label:<16} {method:<10} B={B:<5} "
                      f"sd(lower)={sd:.5f}  [{elapsed:6.1f}s]  "
                      f"{done}/{total}", flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        for line in provenance():
            fh.write(line + "\n")
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nwrote {len(rows)} rows to {OUT}")
    print(f"total {time.perf_counter() - started:.0f}s")

    # A fallback is not an error, but it means those rows describe a
    # different method than the one requested, and averaging them together
    # would understate the requested method's noise. Said out loud rather
    # than left in a column for someone to notice.
    if fallbacks:
        print("\nFALLBACKS -- these rows did not run the method requested:")
        for (ds, asked, B, ran), n in sorted(fallbacks.items()):
            print(f"  {ds:<16} asked {asked:<10} B={B:<5} ran "
                  f"{ran!s:<10} {n}/{len(SEEDS)} seeds")
        print("  Filter on ci_method == method_used before comparing methods.")
    else:
        print("\nno fallbacks: every interval ran the method requested")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
