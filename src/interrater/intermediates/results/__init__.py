"""Everything an analysis produces: estimates, intervals, and test outcomes.

Four kinds of thing live here, and the folder layout says which is which.

``analyzer_results/`` holds the single-coefficient results, all descending from
``PointEstimateResult``. ``permutation_results/`` holds the results of
recomputing a coefficient over subsets of the raters, all descending from
``PermutationResult``. Both of those trees share one root, ``AnalysisResult``,
which sits here rather than in either subfolder because it belongs to neither
alone.

``ci_results/`` holds the intervals themselves, descending from ``CIResult``.
They are not analysis results and are not on that tree: an interval has no
point estimate and no opinion about the data. They are held *by* results,
which is why they live beside them rather than under them.

The remaining module at this level stands outside every tree.
``HypothesisTestResult`` answers a different question -- is the coefficient
distinguishable from zero -- and carries no point estimate of its own to
inherit a contract for.
"""

from .analysis_result import AnalysisResult
from .hypothesis_test_result import HypothesisTestResult

from .ci_results import (
    CIResult,
    BootstrapCIResult,
    ParametricCIResult,
)

from .analyzer_results import (
    PointEstimateResult,
    FleissKappaResult,
    ScottsPiResult,
    GwetsAC1Result,
    BrennanPredigerResult,
    KrippendorffAlphaResult,
    CohenKappaResult,
    CongerKappaResult,
    LightKappaResult,
)

from .permutation_results import (
    PermutationResult,
    PairwisePermutationResult,
    SubsetPermutationResult,
    MultiKSubsetResult,
)

__all__ = [
    "AnalysisResult",
    "CIResult",
    "BootstrapCIResult",
    "ParametricCIResult",
    "PointEstimateResult",
    "FleissKappaResult",
    "ScottsPiResult",
    "GwetsAC1Result",
    "BrennanPredigerResult",
    "KrippendorffAlphaResult",
    "CohenKappaResult",
    "CongerKappaResult",
    "LightKappaResult",
    "PermutationResult",
    "PairwisePermutationResult",
    "SubsetPermutationResult",
    "MultiKSubsetResult",
    "HypothesisTestResult",
]
