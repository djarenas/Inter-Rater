"""Configuration for bootstrap convergence diagnostics (stub).

This file exists so that TYPE_CHECKING imports and isinstance checks in
AnalyzerConfig resolve. The convergence testing feature is not yet
implemented -- every analyzer raises NotImplementedError if
convergence_test is not None. Fields and logic will be added when the
feature is built.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ConvergenceTest:
    """
    I am the configuration for bootstrap convergence diagnostics.

    When set on an AnalyzerConfig, I tell the analyzer to verify that
    n_resamples is large enough for stable CI bounds. Not yet implemented
    -- this is a stub that defines the type so imports and isinstance
    checks work.

    Attributes
    ----------
    (to be defined when convergence testing is implemented)
    """

    def summary(self) -> dict:
        """Return configuration as a dictionary for audit trail."""
        return {"status": "stub -- not yet implemented"}
