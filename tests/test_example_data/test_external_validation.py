"""
The package against other implementations, on every push.

``analyses/external_validation/results.csv`` holds values computed by
irrCAC and irr (R), and statsmodels, scikit-learn and krippendorff
(Python), on the shipped datasets. Those tools are not installed here.
These tests recompute the package's side and check it against the saved
values, so a change that moves a coefficient or a null standard error fails
here, naming the dataset, coefficient and tool.

The package's values are computed by ``measure.py`` in that folder, so the
analysis and these tests cannot drift apart. The data files the other tools
read are checked too: if cleaning changes what a dataset becomes, the saved
values no longer describe it, and the scripts need rerunning.

Tolerance: 5e-6 for irrCAC, which rounds to 5 decimals; 1e-10 otherwise.

Not checked: irr's ``kripp.alpha`` on complete data with more than 2
raters. There it differs from the package by up to 0.0025, while the package
and krippendorff agree to 1e-15 and irrCAC agrees within its rounding
(``compare.py`` has the details). Those rows are listed in ``IRR_KRIPP``.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


FOLDER = (
    Path(__file__).resolve().parent.parent.parent / "analyses" / "external_validation"
)
RESULTS = pd.read_csv(FOLDER / "results.csv")
COMPARED = RESULTS[RESULTS["interrater_value"].notna()]
KEYS = ["dataset", "weights", "coef", "quantity"]

TOLERANCE = {"irrCAC": 5e-6}
DEFAULT_TOLERANCE = 1e-10

IRR_KRIPP = {
    ("fleiss1971", "unweighted"),
    ("sr1_prelim", "unweighted"),
    ("sr1_prelim", "quadratic"),
    ("sr3", "unweighted"),
    ("sr4", "unweighted"),
    ("sr5", "unweighted"),
}


def _load_measure():
    spec = importlib.util.spec_from_file_location("measure", FOLDER / "measure.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def measure():
    return _load_measure()


@pytest.fixture(scope="module")
def datasets(measure):
    return measure.all_datasets()


@pytest.fixture(scope="module")
def package(measure, datasets):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # parametric_null warns on every build
        rows = []
        for name, by_weights in datasets.items():
            rows += measure.package_rows(name, by_weights)
    return pd.DataFrame(rows).set_index(KEYS)["value"]


def _skipped(row):
    return (row.tool == "irr" and row.coef == "kripp"
            and (row.dataset, row.weights) in IRR_KRIPP)


class TestAgainstOtherTools:

    @pytest.mark.parametrize("dataset", sorted(COMPARED["dataset"].unique()))
    def test_values_match(self, package, dataset):
        mismatches = []
        for row in COMPARED[COMPARED["dataset"] == dataset].itertuples():
            if _skipped(row):
                continue
            ours = package[(row.dataset, row.weights, row.coef, row.quantity)]
            tolerance = TOLERANCE.get(row.tool, DEFAULT_TOLERANCE)
            if not abs(ours - row.tool_value) <= tolerance:
                mismatches.append(
                    f"{row.weights} {row.coef} {row.quantity} vs {row.tool} "
                    f"{row.tool_version}: {ours:.10f} vs {row.tool_value:.10f}"
                )
        assert not mismatches, "\n".join(mismatches)

    def test_refused_rows_still_refused(self, package):
        """A coefficient only another tool computes is still refused here."""
        refused = RESULTS[RESULTS["interrater_value"].isna()]
        now_computed = [
            key for key in refused[KEYS].drop_duplicates().itertuples(index=False)
            if not np.isnan(package[tuple(key)])
        ]
        assert not now_computed, (
            f"now computed by the package, so results.csv is out of date: {now_computed}"
        )

    def test_irr_kripp_still_differs(self):
        """The excluded rows are exactly the ones outside tolerance in results.csv."""
        outside = COMPARED[COMPARED["within_tolerance"] == False]  # noqa: E712
        assert set(zip(outside["tool"], outside["coef"])) == {("irr", "kripp")}
        assert set(zip(outside["dataset"], outside["weights"])) == IRR_KRIPP


class TestDataFiles:

    @pytest.mark.parametrize("dataset", sorted(RESULTS["dataset"].unique()))
    def test_matches_cleaning(self, datasets, dataset):
        index = np.asarray(datasets[dataset]["unweighted"].category_index)
        expected = np.where(index < 0, 0, index + 1)
        saved = pd.read_csv(FOLDER / "data" / f"{dataset}.csv").fillna(0).to_numpy(int)
        np.testing.assert_array_equal(saved, expected)
