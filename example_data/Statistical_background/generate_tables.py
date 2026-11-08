"""
Build every table on the statistical background page as a ratings file.

The page (``guide/statistical_background.md``) works through its examples as
2x2 tables and expected values. This module writes each one out as the
ratings a study would have produced: one row per subject, one column per
rater, read by ``RatingsCleaner`` like any other dataset here.
``reproduce_page.py`` then runs interrater on them and checks the page.

Two kinds of table:

- **Built.** The page gives the cells, so the rows are those cells, in order.
  The prevalence sweep gives expected proportions; each of its tables has
  the fewest slides at which every cell is a whole number, so it is exact
  too.
- **Simulated** (seed 42). The biopsy panel and the abstracts are drawn at
  random, so their coefficients land near the page's values, not on them.

Usage::

    python generate_tables.py           # write the CSVs next to this file
    python generate_tables.py --check   # fail if a committed CSV differs
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction
from math import lcm
from pathlib import Path

import numpy as np
import pandas as pd

_ERROR_PREFIX = "[generate_tables] Error"

HERE = Path(__file__).resolve().parent
RANDOM_STATE = 42
ID_COLUMN = "subject_id"
NOT_RATED = "-"                      # the abstracts' missing marker

# Prevalence sweep: skill fixed, prevalence moves (page: "When the condition is rare")
SWEEP_SENSITIVITY = Fraction(80, 100)
SWEEP_SPECIFICITY = Fraction(95, 100)
SWEEP_PREVALENCES = (50, 30, 20, 10, 5, 2)   # percent

# Biopsy panel (page: "More than two raters")
PANEL_N_BIOPSIES = 1000
PANEL_N_READERS = 32
PANEL_PREVALENCE = 0.10
PANEL_SENSITIVITY = 0.90
PANEL_SPECIFICITY = 0.90

# Abstracts (page: "When not everyone rates everything")
ABSTRACTS_N = 1000
ABSTRACTS_PREVALENCE = 0.12
ABSTRACTS_SEES = {"A": range(0, 1000), "B": range(0, 400), "C": range(300, 500)}
ABSTRACTS_SENSITIVITY = {"A": 0.85, "B": 0.80, "C": 0.75}
ABSTRACTS_SPECIFICITY = {"A": 0.95, "B": 0.93, "C": 0.92}


# Built tables -------------------------------------------------------------- #

def build_two_rater_table(
    cells: dict[tuple[str, str], int], prefix: str
) -> pd.DataFrame:
    """
    One row per subject from the cells of a 2x2 table.

    ``cells`` maps (rater A's label, rater B's label) to a count. Rows come
    out cell by cell, in the order given.
    """
    rows = []
    for (label_a, label_b), count in cells.items():
        if not isinstance(count, int) or count < 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: cell ({label_a}, {label_b}) must be a "
                f"non-negative int, got {count!r}"
            )
        rows.extend([(label_a, label_b)] * count)
    width = len(str(len(rows)))
    return pd.DataFrame({
        ID_COLUMN: [f"{prefix}_{i + 1:0{width}d}" for i in range(len(rows))],
        "A": [a for a, _ in rows],
        "B": [b for _, b in rows],
    })


def build_sweep_cells(prevalence_percent: int) -> dict[tuple[str, str], int]:
    """
    Expected 2x2 cells at one prevalence, for two raters independent given
    the truth, on the fewest slides at which every cell is a whole number.
    """
    p = Fraction(prevalence_percent, 100)
    sens, spec = SWEEP_SENSITIVITY, SWEEP_SPECIFICITY
    shares = {
        ("cancer", "cancer"):   p * sens * sens + (1 - p) * (1 - spec) ** 2,
        ("cancer", "healthy"):  p * sens * (1 - sens) + (1 - p) * (1 - spec) * spec,
        ("healthy", "cancer"):  p * (1 - sens) * sens + (1 - p) * spec * (1 - spec),
        ("healthy", "healthy"): p * (1 - sens) ** 2 + (1 - p) * spec * spec,
    }
    if sum(shares.values()) != 1:
        raise ValueError(f"{_ERROR_PREFIX}: sweep shares do not sum to 1")
    n_slides = lcm(*(share.denominator for share in shares.values()))
    return {key: int(share * n_slides) for key, share in shares.items()}


# Simulated tables ---------------------------------------------------------- #

def simulate_biopsy_panel() -> pd.DataFrame:
    """32 readers on 1,000 biopsies, each independent given the truth."""
    rng = np.random.default_rng(RANDOM_STATE)
    truth = rng.random(PANEL_N_BIOPSIES) < PANEL_PREVALENCE
    draws = rng.random((PANEL_N_BIOPSIES, PANEL_N_READERS))
    positive = np.where(truth[:, None],
                        draws < PANEL_SENSITIVITY,
                        draws < 1 - PANEL_SPECIFICITY)
    frame = pd.DataFrame(np.where(positive, "cancer", "healthy"),
                         columns=[f"reader{r + 1:02d}" for r in range(PANEL_N_READERS)])
    frame.insert(0, ID_COLUMN, [f"biopsy_{i + 1:04d}" for i in range(PANEL_N_BIOPSIES)])
    return frame


def simulate_abstracts() -> pd.DataFrame:
    """
    1,000 abstracts and three reviewers. A screens all, B the first 400, C
    abstracts 300-499. Cells a reviewer did not screen hold ``NOT_RATED``.
    """
    rng = np.random.default_rng(RANDOM_STATE)
    truth = rng.random(ABSTRACTS_N) < ABSTRACTS_PREVALENCE
    draws = rng.random((ABSTRACTS_N, len(ABSTRACTS_SEES)))
    frame = pd.DataFrame({ID_COLUMN: [f"abstract_{i + 1:04d}" for i in range(ABSTRACTS_N)]})
    for column, reviewer in enumerate(ABSTRACTS_SEES):
        include = np.where(truth,
                           draws[:, column] < ABSTRACTS_SENSITIVITY[reviewer],
                           draws[:, column] < 1 - ABSTRACTS_SPECIFICITY[reviewer])
        labels = np.where(include, "include", "exclude").astype(object)
        seen = np.isin(np.arange(ABSTRACTS_N), list(ABSTRACTS_SEES[reviewer]))
        labels[~seen] = NOT_RATED
        frame[reviewer] = labels
    return frame


# All tables ---------------------------------------------------------------- #

def build_all_tables() -> dict[str, pd.DataFrame]:
    """Every table on the page, keyed by file name."""
    tables = {
        "professors_data.csv": build_two_rater_table(
            {("pass", "pass"): 82, ("pass", "fail"): 8,
             ("fail", "pass"): 8, ("fail", "fail"): 2}, "essay"),
        "strict_grader_data.csv": build_two_rater_table(
            {("pass", "pass"): 50, ("pass", "fail"): 0,
             ("fail", "pass"): 40, ("fail", "fail"): 10}, "essay"),
        "pathologists_data.csv": build_two_rater_table(
            {("cancer", "cancer"): 15, ("cancer", "healthy"): 50,
             ("healthy", "cancer"): 50, ("healthy", "healthy"): 885}, "slide"),
        "biopsies_2raters_data.csv": build_two_rater_table(
            {("cancer", "cancer"): 90, ("cancer", "healthy"): 90,
             ("healthy", "cancer"): 90, ("healthy", "healthy"): 730}, "biopsy"),
    }
    for percent in SWEEP_PREVALENCES:
        tables[f"prevalence_{percent:02d}_data.csv"] = build_two_rater_table(
            build_sweep_cells(percent), "slide")
    tables["biopsy_panel_data.csv"] = simulate_biopsy_panel()
    tables["abstracts_data.csv"] = simulate_abstracts()
    return tables


def to_csv_text(frame: pd.DataFrame) -> str:
    """The exact text written, so a check compares what a file would hold."""
    return frame.to_csv(index=False, lineterminator="\n")


def write_tables(folder: Path = HERE) -> list[Path]:
    """Write every table to ``folder``; return the paths written."""
    written = []
    for name, frame in build_all_tables().items():
        path = folder / name
        path.write_text(to_csv_text(frame), encoding="utf-8")
        written.append(path)
    return written


def find_mismatches(folder: Path = HERE) -> list[str]:
    """
    Names of committed CSVs that are missing or differ from what this module
    builds. Empty means the files are reproducible.
    """
    mismatches = []
    for name, frame in build_all_tables().items():
        path = folder / name
        if not path.exists():
            mismatches.append(f"{name}: missing")
        elif path.read_text(encoding="utf-8") != to_csv_text(frame):
            mismatches.append(f"{name}: differs from what generate_tables builds")
    return mismatches


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true",
                        help="compare the committed CSVs instead of writing them")
    args = parser.parse_args(argv)
    if args.check:
        mismatches = find_mismatches()
        for line in mismatches:
            print(f"FAIL {line}")
        print(f"{len(build_all_tables()) - len(mismatches)} of "
              f"{len(build_all_tables())} tables reproduce exactly")
        return 1 if mismatches else 0
    for path in write_tables():
        print(f"wrote {path.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
