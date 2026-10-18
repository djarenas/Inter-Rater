"""Base class for all analysis results."""

from __future__ import annotations

from abc import ABC, abstractmethod

_ERROR_PREFIX = "[AnalysisResult] Error"

class AnalysisResult(ABC):
    """
    I am the base for all inter-rater analysis results.

    Every result must answer "what happened?" through summary(), which is the
    interface for printing, serializing, and feeding to visualizers. A result
    holds the outcome of an analysis -- the coefficient, its confidence
    interval if computed, and the metadata needed to understand and visualize
    it. I am immutable: a result is a fact about a dataset at a moment in time
    and should not change.

    I do not know how to visualize myself. That is the visualizer's job. I just
    hold the data. This separation means a new visualizer can be written
    without touching result classes, and a new result type can be added without
    requiring new visualizers -- as long as it holds the right data.
    """

    # Public API -------------------------------------------------------------  #

    @abstractmethod
    def summary(self) -> dict:
        """
        Return everything I represent, in plain Python types (dict, list, float, str).

        This is the contract for downstream consumers: reporters know how to
        serialize dicts, visualizers dispatch on dict structure, and users can
        always read the output without understanding my internal shape.

        Returns
        -------
        dict
            A plain-type representation of my contents. Keys and structure
            depend on my type.
        """
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: summary() must be implemented by subclass"
        )
