"""
Build test Ratings the way users get them: through the cleaner.

A test that writes ``Ratings(codes, semantics)`` by hand can hold a state the
supported path never produces, and then checks behaviour no user can reach.
``clean_ratings`` takes the same codes, writes them out as the labelled table
a user would load, and runs it through ``RatingsCleaner``. The data a test
gets is therefore data the cleaner accepts.

It also refuses to change the data quietly: if the cleaner would have dropped
a rater or a subject, or read a cell differently, it fails and says what
changed. A test that needs such a state on purpose -- a defensive refusal --
builds its Ratings by hand and says so in its docstring.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from interrater.cleaners import CleanerConfig, FrameLayout, RatingsCleaner
from interrater.contracts import EMPTY
from interrater.data_objects.ratings import Ratings
from interrater.semantics.ratings_semantics import RatingsSemantics

_ERROR_PREFIX = "[clean_ratings] Error"

# The id column's header in the table handed to the cleaner
_ID_COLUMN = "subject_id"


def clean_ratings(
    category_index: np.ndarray,
    semantics: RatingsSemantics,
    subject_ids: Sequence[object] | None = None,
) -> Ratings:
    """
    Ratings for these codes, built by the cleaner from a labelled table.

    Parameters
    ----------
    category_index : np.ndarray
        (n_subjects, n_raters) codes: a position in
        ``semantics.categories.names``, or EMPTY for no rating.
    semantics : RatingsSemantics
        The rating task. Its rater names are the table's columns.
    subject_ids : Sequence[object] | None
        Ids for the rows, written to an id column. None leaves the table
        without one, and the cleaner numbers the subjects itself.

    Returns
    -------
    Ratings
        What the cleaner built: the same codes, raters and ids.

    Raises
    ------
    AssertionError
        If the cleaner would have changed anything. The test's data is then
        not data a user can reach through the cleaner.
    """
    codes = np.asarray(category_index)
    names = semantics.categories.names
    raters = [str(rater) for rater in semantics.raters]

    table = pd.DataFrame(
        [[None if code == EMPTY else names[code] for code in row] for row in codes],
        columns=raters,
        dtype=object,
    )
    layout = FrameLayout(use_column_names=True)
    if subject_ids is not None:
        table.insert(0, _ID_COLUMN, [str(i) for i in subject_ids])
        layout = FrameLayout(use_column_names=True, subject_id_col=_ID_COLUMN)

    # unused_categories: keep -- a test's categories are often declared for
    # the task, not chosen by every rater in its small sample.
    config = CleanerConfig(frame_layout=layout, unused_categories="keep")
    ratings = RatingsCleaner(table, semantics, config).ratings

    cleaned_raters = [str(rater) for rater in ratings.semantics.raters]
    if cleaned_raters != raters:
        raise AssertionError(
            f"{_ERROR_PREFIX}: the cleaner kept raters {cleaned_raters} of "
            f"{raters}. Build this Ratings by hand if the test needs it."
        )
    if ratings.category_index.shape != codes.shape:
        raise AssertionError(
            f"{_ERROR_PREFIX}: the cleaner gave shape "
            f"{ratings.category_index.shape}, not {codes.shape}: it dropped "
            f"subjects. Build this Ratings by hand if the test needs it."
        )
    if not np.array_equal(ratings.category_index, codes):
        changed = int((ratings.category_index != codes).sum())
        raise AssertionError(
            f"{_ERROR_PREFIX}: the cleaner read {changed} cell(s) differently."
        )
    if subject_ids is not None and list(ratings.subject_ids) != [str(i) for i in subject_ids]:
        raise AssertionError(f"{_ERROR_PREFIX}: the cleaner changed the subject ids.")
    return ratings
