"""
Compare interrater against other implementations: the comparison.

Merges ``interrater_values.csv`` with ``python_tools.csv`` and
``r_tools.csv`` into ``results.csv``, one row per value another tool
computed. ``interrater_value`` is NaN where the package refuses that
coefficient on that data; ``package_note`` says why.

Each compared row is within tolerance or not. The tolerance is 5e-6 for
irrCAC, which rounds to 5 decimals, and 1e-10 for the other tools.

Results as of 2.2.0: every compared value is within tolerance except 6 of
irr's ``kripp.alpha`` values, on the 5 datasets that are complete and have
more than 2 raters, by up to 0.0025 (Fleiss 1971). On those datasets the
package and krippendorff agree to 1e-15, and irrCAC agrees within its
rounding. On the incomplete datasets irr agrees. Why irr differs is
uncertain.

Run last. The order is in ``measure.py``'s docstring.

Dependencies: pandas.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
OUT = HERE / "results.csv"

KEYS = ["dataset", "weights", "coef", "quantity"]
TOLERANCE = {"irrCAC": 5e-6}
DEFAULT_TOLERANCE = 1e-10


def main():
    package = pd.read_csv(HERE / "interrater_values.csv").rename(
        columns={"value": "interrater_value", "note": "package_note"}
    )[KEYS + ["interrater_value", "package_note"]]
    tools = pd.concat(
        [pd.read_csv(HERE / "python_tools.csv"), pd.read_csv(HERE / "r_tools.csv")],
        ignore_index=True,
    ).rename(columns={"value": "tool_value"})

    table = tools.merge(package, on=KEYS, how="left")
    table["diff"] = table["tool_value"] - table["interrater_value"]
    table["tolerance"] = table["tool"].map(TOLERANCE).fillna(DEFAULT_TOLERANCE)
    compared = table["interrater_value"].notna()
    table["within_tolerance"] = (table["diff"].abs() <= table["tolerance"]).where(compared)
    table = table[["dataset", "weights", "coef", "quantity", "tool", "tool_version",
                   "tool_value", "interrater_value", "diff", "tolerance",
                   "within_tolerance", "package_note"]]
    table.to_csv(OUT, index=False)

    done = table[compared]
    summary = done.groupby(["quantity", "tool"]).agg(
        compared=("diff", "size"),
        max_abs_diff=("diff", lambda d: d.abs().max()),
        outside=("within_tolerance", lambda w: int((~w.astype(bool)).sum())),
    )
    print(f"{OUT.name}: {len(table)} rows, {len(done)} compared, "
          f"{len(table) - len(done)} computed only by the other tool\n")
    print(summary.to_string())
    outside = done[~done["within_tolerance"].astype(bool)]
    if len(outside):
        print("\nOutside tolerance:")
        print(outside[["dataset", "weights", "coef", "tool", "tool_value",
                       "interrater_value", "diff"]].to_string(index=False))


if __name__ == "__main__":
    main()
