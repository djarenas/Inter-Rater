"""
Every YAML file shipped in example_data/ loads with the current code.

The files are found by globbing, so a new example is tested without being
listed here. Each file is its own test, so a failure names the file.

The analyzer configs are not tied to one coefficient; they are loaded with
FleissKappaConfig, as the docs do. A weights file is loaded through the
semantics file that names it. A file that matches no rule fails, so a new
kind of YAML gets a rule here rather than going untested.
"""

from pathlib import Path

import pytest

from interrater import CleanerConfig, RatingsSemantics
from interrater.analyzers.configs import (
    FleissKappaConfig,
    PairwisePermutationConfig,
    SubsetPermutationConfig,
)

EXAMPLE_DATA_DIR = Path(__file__).resolve().parents[2] / "example_data"
EXAMPLE_YAMLS = sorted(EXAMPLE_DATA_DIR.glob("*/*.yaml"))


def _load_weights(path):
    """Load the semantics in the same folder that names this weights file."""
    users = [p for p in path.parent.glob("*_semantics.yaml")
             if path.name in p.read_text()]
    assert users, f"{path.name} is not named by any semantics file"
    for semantics in users:
        RatingsSemantics.build_from_yaml(semantics)


LOADERS = {
    "_semantics.yaml": RatingsSemantics.build_from_yaml,
    "_cleaner_config.yaml": CleanerConfig.build_from_yaml,
    "_analyzer_config.yaml": FleissKappaConfig.build_from_yaml,
    "_pairwise_permutation_config.yaml": PairwisePermutationConfig.build_from_yaml,
    "_subset_permutation_config.yaml": SubsetPermutationConfig.build_from_yaml,
    "_weights.yaml": _load_weights,
}


def test_example_yamls_found():
    assert len(EXAMPLE_YAMLS) > 0


@pytest.mark.parametrize("path", EXAMPLE_YAMLS, ids=lambda p: p.name)
def test_example_yaml_loads(path):
    for ending, load in LOADERS.items():
        if path.name.endswith(ending):
            load(path)
            return
    pytest.fail(f"{path.name}: no loader rule; add one to LOADERS")
