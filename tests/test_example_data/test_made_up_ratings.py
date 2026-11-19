"""
The counts-to-Ratings helper in example_data/.

``made_up_ratings.py`` lives outside the package, so it is loaded by path.
These tests check that it writes the counts back out exactly. That the
coefficients on Fleiss 1971 come out right is checked in
test_reference_data.py.
"""

import numpy as np
import pandas as pd
import pytest

from interrater.array_utils import category_counts

# Small tables written out here, so no example folder exists only for tests
UNBALANCED = pd.DataFrame(
    [[2, 1, 0], [0, 3, 0], [4, 0, 0], [1, 1, 1], [0, 0, 2], [2, 0, 2]],
    columns=["a", "b", "c"],
)


class TestExpansion:

    def test_counts_come_back_exactly(self, fleiss1971_counts, build_made_up_ratings):
        for counts in (fleiss1971_counts, UNBALANCED):
            ratings = build_made_up_ratings(counts)
            recounted = category_counts(ratings.category_index, ratings.n_categories)
            np.testing.assert_array_equal(recounted, counts.to_numpy())

    def test_names_and_ids_are_kept(self, fleiss1971_counts, build_made_up_ratings):
        ratings = build_made_up_ratings(fleiss1971_counts)
        assert list(ratings.subject_ids) == list(fleiss1971_counts.index)
        assert ratings.semantics.categories.names == tuple(fleiss1971_counts.columns)
        assert ratings.semantics.raters == tuple(f"slot_{i}" for i in range(1, 7))

    def test_short_rows_get_empty_slots(self, build_made_up_ratings):
        ratings = build_made_up_ratings(UNBALANCED)
        per_subject = (ratings.category_index >= 0).sum(axis=1)
        np.testing.assert_array_equal(per_subject, UNBALANCED.sum(axis=1))

    def test_a_subject_rated_once_is_refused(self, build_made_up_ratings):
        counts = pd.DataFrame([[2, 0], [1, 0]], columns=["yes", "no"])
        with pytest.raises(ValueError):
            build_made_up_ratings(counts)
