"""
What a run of the package did, written for a person to read.

Reports are terminal. Unlike an intermediate, nothing downstream consumes one:
an intermediate exists because a visualiser or another analyser needs it, while
a report exists to be read, printed, or filed beside the config that produced
it. That is why these live apart from ``intermediates`` -- they are a
presentation layer, a sibling of ``visualizers`` in the same medium as prose.

Every report holds everything it computed and offers two views of it:
``summary()`` in plain types for writing to disk, and ``render(verbosity)``
built from that summary so the two cannot drift.
"""

from .category_report import CategoryReport
from .cleaning_report import CleaningReport
from .drop_report import DropPass, DropReport
from .matching_report import MatchingReport

__all__ = [
    "CleaningReport",
    "MatchingReport",
    "DropReport",
    "DropPass",
    "CategoryReport",
]
