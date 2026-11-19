"""
Build test Ratings from a small table of counts.

Some tests are easiest to read as counts: for each subject, how many raters
chose each category. ``ratings_from_counts`` writes each row out left to
right in category order and names the columns slot_1, slot_2, ... The raters
are made up, so use it only where rater identity does not matter.

Built by hand, not through the cleaner: on unbalanced counts the last slots
hold only a rating or two, and the cleaner's minimum per rater would drop
them and change the counts the test asked for.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from interrater.contracts import EMPTY
from interrater.data_objects.ratings import Ratings
from interrater.semantics.ratings_semantics import RatingsSemantics


def ratings_from_counts(
    counts: np.ndarray, category_names: Sequence[str] | None = None
) -> Ratings:
    """Ratings whose per-subject counts are ``counts`` exactly."""
    counts = np.asarray(counts)
    n_subjects, n_categories = counts.shape
    n_slots = int(counts.sum(axis=1).max())
    index = np.full((n_subjects, n_slots), EMPTY, dtype=np.int8)
    for row, row_counts in enumerate(counts):
        written = np.repeat(np.arange(n_categories), row_counts)
        index[row, : written.size] = written
    names = list(category_names) if category_names is not None else [
        f"column {j}" for j in range(n_categories)
    ]
    semantics = RatingsSemantics.build_from_dict({
        "raters": [f"slot_{i + 1}" for i in range(n_slots)],
        "n_subjects": n_subjects,
        "categories": {"kind": "unordered", "names": names},
    })
    return Ratings(index, semantics)
