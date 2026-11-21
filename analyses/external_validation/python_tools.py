"""
Compare interrater against other implementations: the Python tools.

Reads ``data/<dataset>.csv`` (written by ``measure.py``) and computes what
statsmodels, scikit-learn and krippendorff compute on it. Writes
``python_tools.csv``, one row per (dataset, weights, coefficient, quantity,
tool), with the tool's version.

What each tool covers:

- statsmodels ``fleiss_kappa``: Fleiss' kappa (``method="fleiss"``) and
  Brennan-Prediger (``method="randolph"``, the unweighted form), on complete
  data only, since it takes a table of counts with the same number of
  raters per subject.
- statsmodels ``cohens_kappa``: Cohen's kappa, unweighted, linear and
  quadratic, and its standard error under the null (``std_kappa0``).
- scikit-learn ``cohen_kappa_score``: Cohen's kappa, unweighted, linear and
  quadratic.
- krippendorff ``alpha``: nominal alpha, and interval alpha, which is
  alpha with quadratic weights on 1..K codes.

Cohen's kappa is computed on the subjects both raters rated, as the
package does. K is read from ``interrater_values.csv``, because a category
nobody chose does not appear in the data. The weightings for a dataset are
the ones ``interrater_values.csv`` has for it (none for Fleiss 1971, whose
categories are unordered).

Run after ``measure.py``. The order is in that file's docstring.

Dependencies: numpy, pandas, statsmodels, scikit-learn, krippendorff. Not
the package.
"""

from __future__ import annotations

from importlib.metadata import version
from pathlib import Path

import krippendorff
import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score
from statsmodels.stats import inter_rater


HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "python_tools.csv"

VERSIONS = {
    "statsmodels": version("statsmodels"),
    "scikit-learn": version("scikit-learn"),
    "krippendorff": version("krippendorff"),
}


def main():
    package = pd.read_csv(HERE / "interrater_values.csv")
    n_categories = package.groupby("dataset")["K"].first()
    weightings_of = package.groupby("dataset")["weights"].unique()
    rows = []

    def add(dataset, weights, coef, quantity, tool, value):
        rows.append(dict(dataset=dataset, weights=weights, coef=coef,
                         quantity=quantity, tool=tool,
                         tool_version=VERSIONS[tool], value=float(value)))

    for dataset, K in n_categories.items():
        codes = pd.read_csv(DATA / f"{dataset}.csv").to_numpy(dtype=float)
        labels = list(range(1, K + 1))
        complete = codes[~np.isnan(codes).any(axis=1)].astype(int)
        weightings = list(weightings_of[dataset])

        # krippendorff takes raters x subjects, NaN for missing.
        add(dataset, "unweighted", "kripp", "estimate", "krippendorff",
            krippendorff.alpha(reliability_data=codes.T, value_domain=labels,
                               level_of_measurement="nominal"))
        if "quadratic" in weightings:
            add(dataset, "quadratic", "kripp", "estimate", "krippendorff",
                krippendorff.alpha(reliability_data=codes.T, value_domain=labels,
                                   level_of_measurement="interval"))

        if len(complete) == len(codes):
            counts = np.stack([(complete == k).sum(axis=1) for k in labels], axis=1)
            add(dataset, "unweighted", "fleiss", "estimate", "statsmodels",
                inter_rater.fleiss_kappa(counts, method="fleiss"))
            add(dataset, "unweighted", "bp", "estimate", "statsmodels",
                inter_rater.fleiss_kappa(counts, method="randolph"))

        if codes.shape[1] == 2:
            a, b = complete[:, 0], complete[:, 1]
            table = np.zeros((K, K))
            np.add.at(table, (a - 1, b - 1), 1)
            for w in weightings:
                wt = None if w == "unweighted" else w
                sm = inter_rater.cohens_kappa(table, wt=wt)
                add(dataset, w, "cohen", "estimate", "statsmodels", sm.kappa)
                add(dataset, w, "cohen", "null_se", "statsmodels", sm.std_kappa0)
                add(dataset, w, "cohen", "estimate", "scikit-learn",
                    cohen_kappa_score(a, b, labels=labels, weights=wt))

    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"{OUT.name}: {len(rows)} rows")


if __name__ == "__main__":
    main()
