"""
The statistical background page and the package agree.

``example_data/Statistical_background/`` holds every table on
``guide/statistical_background.md`` as ratings, the module that builds them,
and the script that runs interrater on them and checks the page. These tests
run both, so the suite fails when the tables, the page's values or the
package drift apart.
"""

import importlib.util
from pathlib import Path

import pytest

FOLDER = Path(__file__).resolve().parents[2] / "example_data" / "Statistical_background"


def _load_module(name):
    """Import a script from the folder without putting it on sys.path."""
    spec = importlib.util.spec_from_file_location(name, FOLDER / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def generate_tables():
    return _load_module("generate_tables")


@pytest.fixture(scope="module")
def reproduce_page():
    return _load_module("reproduce_page")


def test_committed_tables_reproduce(generate_tables):
    assert generate_tables.find_mismatches() == []


def test_page_numbers_reproduce(reproduce_page):
    checks = reproduce_page.run(verbose=False)
    assert checks.failures == []
    assert checks.n_checks > 0


def test_a_wrong_page_value_fails(reproduce_page, monkeypatch):
    """The checks can fail: a page value off in its last digit is caught."""
    professors = dict(reproduce_page.PAGE["professors"], kappa="0.112")
    monkeypatch.setitem(reproduce_page.PAGE, "professors", professors)
    checks = reproduce_page.run(verbose=False)
    assert checks.failures == ["professors: Cohen's kappa: got 0.11111, page 0.112"]
