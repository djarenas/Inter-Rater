"""
What `from interrater import ...` offers.

Everything an analysis hands back -- results, their intervals and tests, the
profiles they carry -- can be named from the top level, so a user can write
an isinstance check or a type hint without reaching into
`interrater.intermediates`.
"""

import subprocess
import sys

import interrater
from interrater.intermediates.results.analyzer_results import COEFFICIENT_RESULTS


def test_every_exported_name_imports_once():
    assert len(interrater.__all__) == len(set(interrater.__all__))
    missing = [name for name in interrater.__all__ if not hasattr(interrater, name)]
    assert missing == []


def test_what_an_analysis_returns_is_exported():
    """Each coefficient's result class, and the types found on a result."""
    returned = {cls.__name__ for cls in COEFFICIENT_RESULTS.values()} | {
        "PairwisePermutationResult", "SubsetPermutationResult", "MultiKSubsetResult",
        "BootstrapCIResult", "ParametricCIResult", "HypothesisTestResult",
        "BalanceProfile", "RaterProfile", "AgreementProfile", "AgreementParts",
        "NotComputed",
    }
    assert returned - set(interrater.__all__) == set()


def test_import_does_not_load_matplotlib():
    """The plotting names resolve on first use, not at import."""
    code = "import sys, interrater; print('matplotlib' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"
