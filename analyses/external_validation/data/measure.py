"""
Compare interrater against other implementations: the package's side.

Other packages compute the same coefficients: irrCAC and irr in R;
statsmodels, scikit-learn and krippendorff in Python. Each tool's value is
compared with the package's, on the same data.

This script does two things:

1. Writes each dataset, after cleaning, to ``data/<dataset>.csv``: one
   column per rater, categories coded 1..K in semantics order, missing
   ratings left empty. Every tool reads these files, so all of them see the
   same matrix. The comparison tests the coefficient arithmetic, not the
   cleaning.
2. Writes ``interrater_values.csv``: one row per (dataset, weights,
   coefficient, quantity). The quantity is ``estimate`` (the point
   estimate) or ``null_se`` (the standard error of the ``parametric_null``
   interval, for Fleiss, Scott and Cohen). A coefficient the package
   refuses on that data gets NaN, and ``note`` says why.

Coding categories as 1..K makes the other tools' linear and quadratic
weights the same as the package's, since both are built from the distance
between positions.

Fleiss (1971) publishes only counts per subject. Its rater slots are made
up (``example_data/made_up_ratings.py``), so the coefficients that read
rater identity (Scott, Cohen, Conger, Light) are skipped for it.

Run from the repository root, in this order:

    python analyses/external_validation/measure.py
    python analyses/external_validation/python_tools.py
    Rscript analyses/external_validation/r_tools.R
    python analyses/external_validation/compare.py

Dependencies: the package, numpy, pandas.
"""

from __future__ import annotations

import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from interrater import (
    CleanerConfig, OrderedCategories, RatingsCleaner, RatingsSemantics,
    BrennanPredigerAnalyzer, CohenKappaAnalyzer, CongerKappaAnalyzer,
    FleissKappaAnalyzer, GwetAC1Analyzer, KrippendorffAlphaAnalyzer,
    LightKappaAnalyzer, ScottsAnalyzer,
)
from interrater.analyzers.configs import (
    BrennanPredigerConfig, CohenKappaConfig, CongerKappaConfig,
    FleissKappaConfig, GwetAC1Config, KrippendorffAlphaConfig,
    LightKappaConfig, ScottsConfig,
)


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
EXAMPLE = REPO / "example_data"
DATA_OUT = HERE / "data"
OUT = HERE / "interrater_values.csv"

Z_95 = 1.959963984540054

# (name, folder, data file, semantics yaml, cleaner yaml)
DATASETS = [
    ("gwet", "Gwet_example", "gwet_cont3x3_data.csv",
     "gwet_cont3x3_semantics.yaml", "gwet_cont3x3_cleaner_config.yaml"),
    ("synthetic", "Synthetic_2raters", "synthetic_2raters_data.txt",
     "synthetic_2raters_semantics.yaml", "synthetic_2raters_cleaner_config.yaml"),
    ("sr1_prelim", "SR1_OSF2018", "SR1_prelim_data.csv",
     "SR1_prelim_semantics.yaml", "SR1_prelim_cleaner_config.yaml"),
    ("sr1_main", "SR1_OSF2018", "SR1_main_data.txt",
     "SR1_main_semantics.yaml", "SR1_main_cleaner_config.yaml"),
    ("sr2", "SR2_JSRC2019", "SR2_data.txt",
     "SR2_semantics.yaml", "SR2_cleaner_config.yaml"),
    ("sr3", "SR3_JGIM2019", "SR3_data.txt",
     "SR3_semantics.yaml", "SR3_cleaner_config.yaml"),
    ("sr4", "SR4_FIS_CVD", "SR4_data.txt",
     "SR4_semantics.yaml", "SR4_cleaner_config.yaml"),
    ("sr5", "SR5_SciRep2020", "SR5_data.txt",
     "SR5_semantics.yaml", "SR5_cleaner_config.yaml"),
]

# coefficient -> (analyzer, config, extra config fields)
COEFFICIENTS = {
    "fleiss": (FleissKappaAnalyzer, FleissKappaConfig, {"balanced_data_only": False}),
    "ac1":    (GwetAC1Analyzer, GwetAC1Config, {}),
    "bp":     (BrennanPredigerAnalyzer, BrennanPredigerConfig, {}),
    "kripp":  (KrippendorffAlphaAnalyzer, KrippendorffAlphaConfig, {}),
    "scott":  (ScottsAnalyzer, ScottsConfig, {}),
    "cohen":  (CohenKappaAnalyzer, CohenKappaConfig, {}),
    "conger": (CongerKappaAnalyzer, CongerKappaConfig, {}),
    "light":  (LightKappaAnalyzer, LightKappaConfig, {}),
}
RATER_IDENTITY = {"scott", "cohen", "conger", "light"}
NULL_SE = {"fleiss", "scott", "cohen"}


# -- Loading ---------------------------------------------------------------- #
#
# The same loader as analyses/seed_stability/measure.py: .csv files ship
# with a header, .txt files are tab-separated and headerless.

def load_ratings(folder, data, semantics, cleaner, weighting=None):
    """Load a shipped dataset through the cleaning pipeline.

    With ``weighting`` ('linear' or 'quadratic'), the same categories are
    declared ordered with that weighting, which is how a weighted analysis
    is asked for.
    """
    d = EXAMPLE / folder
    sem = RatingsSemantics.build_from_yaml(str(d / semantics))
    if weighting is not None:
        sem = RatingsSemantics(
            raters=sem.raters,
            categories=OrderedCategories(sem.categories.names, weighting=weighting),
            n_subjects=sem.n_subjects,
            unique_subjects=sem.unique_subjects,
        )
    cfg = CleanerConfig.build_from_yaml(str(d / cleaner))
    path = d / data
    df = pd.read_csv(
        path,
        sep="," if path.suffix == ".csv" else "\t",
        header=0 if path.suffix == ".csv" else None,
        dtype=str,
    )
    return RatingsCleaner(df, sem, config=cfg).ratings


def load_fleiss1971():
    """Fleiss (1971) counts, as ratings in made-up rater slots."""
    spec = importlib.util.spec_from_file_location(
        "made_up_ratings", EXAMPLE / "made_up_ratings.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    counts = pd.read_csv(
        EXAMPLE / "Fleiss1971" / "Fleiss1971_psychiatric_data_nij.csv", index_col=0
    )
    return module.build_made_up_ratings(counts)


# -- Measuring -------------------------------------------------------------- #

def point_estimate(analyzer):
    """The point estimate, computed the way ``estimate()`` computes it."""
    return float(analyzer._estimate_from_input(analyzer.estimator_input))


def null_se(analyzer_cls, config_cls, extra, ratings):
    """Standard error behind the ``parametric_null`` 95% interval."""
    config = config_cls(ci_method="parametric_null", n_resamples=1, **extra)
    ci = analyzer_cls(ratings, config).estimate().confidence_interval
    return (ci.upper - ci.lower) / (2 * Z_95)


def measure(name, ratings_by_weights, rows):
    """Write the dataset's matrix and append the package's values to rows."""
    unweighted = ratings_by_weights["unweighted"]
    index = np.asarray(unweighted.category_index)
    codes = np.where(index < 0, "", (index + 1).astype(str))
    pd.DataFrame(
        codes, columns=[f"r{i + 1}" for i in range(index.shape[1])]
    ).to_csv(DATA_OUT / f"{name}.csv", index=False)

    n_categories = len(unweighted.semantics.categories.names)
    for weights, ratings in ratings_by_weights.items():
        for coef, (analyzer_cls, config_cls, extra) in COEFFICIENTS.items():
            if name == "fleiss1971" and coef in RATER_IDENTITY:
                continue
            row = dict(dataset=name, K=n_categories, n_raters=index.shape[1],
                       weights=weights, coef=coef)
            try:
                analyzer = analyzer_cls(ratings, config_cls(**extra))
                rows.append(row | dict(quantity="estimate",
                                       value=point_estimate(analyzer), note=""))
            except (TypeError, ValueError, NotImplementedError) as e:
                rows.append(row | dict(quantity="estimate", value=np.nan,
                                       note=f"{type(e).__name__}: {e}"[:200]))
                continue
            if coef in NULL_SE:
                try:
                    rows.append(row | dict(
                        quantity="null_se", note="",
                        value=null_se(analyzer_cls, config_cls, extra, ratings),
                    ))
                except (TypeError, ValueError) as e:
                    rows.append(row | dict(quantity="null_se", value=np.nan,
                                           note=f"{type(e).__name__}: {e}"[:200]))


def main():
    warnings.simplefilter("ignore")   # parametric_null warns on every build
    DATA_OUT.mkdir(exist_ok=True)
    rows = []
    measure("fleiss1971", {"unweighted": load_fleiss1971()}, rows)
    for name, *files in DATASETS:
        by_weights = {"unweighted": load_ratings(*files)}
        if len(by_weights["unweighted"].semantics.categories.names) >= 3:
            for weighting in ("linear", "quadratic"):
                by_weights[weighting] = load_ratings(*files, weighting=weighting)
        measure(name, by_weights, rows)
    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)
    computed = table["value"].notna().sum()
    print(f"{OUT.name}: {len(table)} rows, {computed} computed, "
          f"{len(table) - computed} refused")


if __name__ == "__main__":
    main()
