"""
Turn a table of counts into a Ratings object with made-up raters.

Some published datasets give only counts: for each subject, how many raters
chose each category (Fleiss 1971 is one). interrater analyses Ratings, so the
examples and tests that start from counts use this helper.

- The raters are made up. Counts do not say who rated what, so each row is
  written out left to right in category order and the columns are named
  slot_1, slot_2, ... They are positions, not people.
- Use the result only with coefficients that treat raters as
  interchangeable: Fleiss, Scott, Gwet's AC1, Brennan-Prediger and
  Krippendorff. Cohen, Conger, Light, rater profiles and rater injectors
  would compute something from the made-up columns, and it would mean
  nothing.
- This file is not part of the package on purpose. interrater has no
  supported route from counts to Ratings.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from interrater import Ratings, RatingsSemantics
from interrater.contracts import EMPTY

_ERROR_PREFIX = "[made_up_ratings] Error"


def build_made_up_ratings(counts: pd.DataFrame) -> Ratings:
    """
    Ratings from a subjects x categories table of counts.

    Parameters
    ----------
    counts : pd.DataFrame
        One row per subject, one column per category, non-negative integer
        counts. The index gives the subject ids and the columns give the
        category names, in order. K is the number of columns, so a category
        nobody chose still needs its column of zeros.

    Returns
    -------
    Ratings
        One column per slot, as many slots as the largest row total. A
        subject with fewer ratings has empty cells at the end of its row.
    """
    if not isinstance(counts, pd.DataFrame):
        raise TypeError(
            f"{_ERROR_PREFIX}: counts must be a pandas DataFrame, got "
            f"{type(counts).__name__}."
        )
    values = counts.to_numpy()
    if values.dtype.kind not in "iu" or (values < 0).any():
        raise ValueError(
            f"{_ERROR_PREFIX}: counts must be non-negative integers."
        )

    n_subjects, n_categories = values.shape
    n_slots = int(values.sum(axis=1).max())
    category_index = np.full((n_subjects, n_slots), EMPTY, dtype=np.int8)
    for row, row_counts in enumerate(values):
        written = np.repeat(np.arange(n_categories), row_counts)
        category_index[row, : written.size] = written

    semantics = RatingsSemantics.build_from_dict({
        "raters": [f"slot_{i + 1}" for i in range(n_slots)],
        "n_subjects": n_subjects,
        "categories": {
            "kind": "unordered",
            "names": [str(name) for name in counts.columns],
        },
    })
    return Ratings(category_index, semantics, subject_ids=list(counts.index))
