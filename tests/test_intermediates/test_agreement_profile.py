"""
AgreementProfile: what the data says about agreement, before any coefficient.

Several tests compare against the agreement functions the coefficients use,
so the profile cannot drift from them: its prevalence must be the p_j in
Fleiss', Scott's and Gwet's chance terms, and its rater marginals the ones in
Cohen's and Conger's.
"""

import math

import numpy as np
import pytest

from interrater import EMPTY, NotComputed
from interrater.analyzers.agreement_functions import (
    category_proportions,
    observed_agreement,
    rater_marginals,
)
from interrater.data_objects.ratings import Ratings
from interrater.intermediates.data_profiles import (
    AgreementProfile,
    PrevalenceBiasIndex,
)
from interrater.intermediates.estimator_parts import EstimatorInput, NijMatrix
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics

E = EMPTY


def _ratings(index, raters, cats=("yes", "no", "maybe")):
    sem = RatingsSemantics(raters=list(raters), categories=UnorderedCategories(cats))
    return Ratings(category_index=np.array(index, dtype=np.int8), semantics=sem)


def _profile(index, raters, cats=("yes", "no", "maybe")):
    # From the array, as internal steps build it: Ratings refuses a subject
    # rated once, and these tests need one.
    return AgreementProfile.build_from_category_index(
        np.array(index, dtype=np.int8), tuple(cats), tuple(raters)
    )


# Unbalanced: subject 4 has an extra rater, subject 5 is rated once
UNBALANCED = [
    [0, 0, E],
    [0, 1, E],
    [1, 1, 1],
    [2, 1, E],
    [0, 0, 1],
    [E, 2, E],
]
BALANCED = [[0, 0, 1], [1, 1, 1], [2, 1, 2], [0, 0, 0], [1, 2, 1]]


def _pool_input(index, k=3):
    data = EstimatorInput(category_index=np.array(index, dtype=np.int8), n_categories=k)
    rows = np.flatnonzero((data.category_index != E).sum(axis=1) >= 2)
    return data, rows


class TestMatchesTheCoefficients:

    def test_prevalence_is_the_p_j_the_chance_terms_use(self):
        """Averaged per subject, on subjects rated at least twice."""
        profile = _profile(UNBALANCED, "abc")
        data, rows = _pool_input(UNBALANCED)
        expected = category_proportions(data.take_subjects(rows).nij)
        np.testing.assert_allclose(list(profile.prevalence.values()), expected)

    def test_rater_marginals_are_the_ones_conger_uses(self):
        profile = _profile(UNBALANCED, "abc")
        data, rows = _pool_input(UNBALANCED)
        expected = rater_marginals(data, rows)
        for r, rater in enumerate("abc"):
            np.testing.assert_allclose(
                list(profile.rater_marginals[rater].values()), expected[r]
            )

    def test_specific_agreement_decomposes_po_on_balanced_data(self):
        """Po = Sum_j p_j * s_j when every subject has the same raters."""
        profile = AgreementProfile.build_from_ratings(_ratings(BALANCED, "abc"))
        data, _ = _pool_input(BALANCED)
        decomposed = sum(
            profile.prevalence[c] * profile.specific_agreement[c]
            for c in profile.category_names
        )
        assert decomposed == pytest.approx(observed_agreement(data.nij))

    def test_specific_agreement_for_two_raters(self):
        """2 n_jj / (row_j + col_j), the usual two-rater form."""
        index = [[0, 0], [0, 1], [1, 1], [1, 1], [0, 0], [1, 0]]
        profile = AgreementProfile.build_from_ratings(_ratings(index, "ab", ("yes", "no")))
        # yes: n_jj = 2, row (a says yes) = 3, col (b says yes) = 3
        assert profile.specific_agreement["yes"] == pytest.approx(2 * 2 / (3 + 3))
        # no: n_jj = 2, row = 3, col = 3
        assert profile.specific_agreement["no"] == pytest.approx(2 * 2 / (3 + 3))


class TestPoolAndEdgeCases:

    def test_subjects_rated_once_are_left_out(self):
        profile = _profile(UNBALANCED, "abc")
        assert profile.n_subjects == 5
        assert profile.n_ratings == 12

    def test_unused_category(self):
        profile = AgreementProfile.build_from_ratings(
            _ratings([[0, 0], [1, 1], [0, 1]], "ab")
        )
        assert profile.prevalence["maybe"] == 0.0
        assert math.isnan(profile.specific_agreement["maybe"])
        assert "n/a" in str(profile)


class TestPrevalenceBias:

    def test_two_raters_two_categories(self):
        index = [[0, 0], [0, 1], [1, 1], [1, 0], [0, 0]]
        profile = AgreementProfile.build_from_ratings(_ratings(index, ("Ann", "Bob"), ("yes", "no")))
        expected = PrevalenceBiasIndex.build_from_pair(
            np.array(index), ("yes", "no"), ("Ann", "Bob")
        )
        assert profile.prevalence_bias == expected

    def test_otherwise_says_why_not(self):
        profile = AgreementProfile.build_from_ratings(_ratings(BALANCED, "abc"))
        assert isinstance(profile.prevalence_bias, NotComputed)
        assert "3 raters and 3 categories" in profile.prevalence_bias.reason


class TestFromCounts:

    def test_counts_only(self):
        counts = np.array([[2, 0], [1, 1], [0, 2], [2, 0], [1, 0]])
        profile = AgreementProfile.build_from_nij_matrix(NijMatrix(counts))
        assert profile.category_names == ("column 0", "column 1")
        assert profile.n_subjects == 4
        assert isinstance(profile.rater_marginals, NotComputed)
        assert isinstance(profile.prevalence_bias, PrevalenceBiasIndex)
        assert profile.prevalence_bias.prevalence_index == pytest.approx((2 - 1) / 4)
        assert isinstance(profile.prevalence_bias.bias_index, NotComputed)

    def test_more_than_two_ratings_per_subject(self):
        counts = np.array([[3, 0], [1, 2], [2, 1]])
        profile = AgreementProfile.build_from_nij_matrix(NijMatrix(counts))
        assert isinstance(profile.prevalence_bias, NotComputed)
        assert "up to 3 ratings per subject" in profile.prevalence_bias.reason

    def test_counts_and_ratings_agree(self):
        """The same data as counts or as ratings gives the same numbers."""
        from_ratings = _profile(UNBALANCED, "abc")
        data, _ = _pool_input(UNBALANCED)
        from_counts = AgreementProfile.build_from_nij_matrix(data.nij)
        np.testing.assert_allclose(
            list(from_ratings.prevalence.values()), list(from_counts.prevalence.values())
        )
        np.testing.assert_allclose(
            list(from_ratings.specific_agreement.values()),
            list(from_counts.specific_agreement.values()),
        )
        assert from_ratings.n_subjects == from_counts.n_subjects


class TestSelfChecks:

    def _kwargs(self, **overrides):
        base = dict(
            category_names=("yes", "no"),
            n_subjects=3,
            n_ratings=6,
            prevalence={"yes": 0.5, "no": 0.5},
            specific_agreement={"yes": 0.6, "no": 0.6},
            rater_marginals={"a": {"yes": 0.5, "no": 0.5}},
            prevalence_bias=NotComputed("not two raters"),
        )
        base.update(overrides)
        return base

    def test_valid(self):
        assert AgreementProfile(**self._kwargs()).n_subjects == 3

    def test_prevalence_must_sum_to_one(self):
        with pytest.raises(ValueError, match="sum to"):
            AgreementProfile(**self._kwargs(prevalence={"yes": 0.5, "no": 0.6}))

    def test_rater_row_must_sum_to_one(self):
        with pytest.raises(ValueError, match="rater 'a'"):
            AgreementProfile(**self._kwargs(rater_marginals={"a": {"yes": 0.9, "no": 0.5}}))

    def test_specific_agreement_within_zero_and_one(self):
        with pytest.raises(ValueError, match="outside"):
            AgreementProfile(**self._kwargs(specific_agreement={"yes": 1.2, "no": 0.5}))

    def test_too_few_ratings_refused(self):
        with pytest.raises(ValueError, match="cannot cover"):
            AgreementProfile(**self._kwargs(n_ratings=5))

    def test_prevalence_bias_must_use_these_categories(self):
        wrong = PrevalenceBiasIndex.build_from_pair(
            np.array([[0, 0], [1, 1]]), ("no", "yes"), ("a", "b")
        )
        with pytest.raises(ValueError, match="positive first"):
            AgreementProfile(**self._kwargs(prevalence_bias=wrong))

    def test_missing_rater_marginals_must_say_why(self):
        with pytest.raises(TypeError, match="NotComputed"):
            AgreementProfile(**self._kwargs(rater_marginals=None))


class TestOutputs:

    def test_summary(self):
        profile = AgreementProfile.build_from_ratings(_ratings(BALANCED, "abc"))
        s = profile.summary()
        assert s["n_subjects"] == 5
        assert s["category_names"] == ["yes", "no", "maybe"]
        assert set(s["rater_marginals"]) == {"a", "b", "c"}
        assert s["prevalence_bias"] == {"not_computed": profile.prevalence_bias.reason}

    def test_printed(self):
        text = str(AgreementProfile.build_from_ratings(_ratings(BALANCED, "abc")))
        assert "5 subjects rated at least twice, 15 ratings" in text
        assert "How often each category was used (share of ratings)" in text
        assert "How often raters agreed on each category" in text
        assert "Prevalence and bias: not computed (defined for two raters" in text

    def test_rater_rows_sit_under_overall(self):
        lines = str(AgreementProfile.build_from_ratings(_ratings(BALANCED, "abc"))).splitlines()
        overall = next(i for i, line in enumerate(lines) if line.strip().startswith("overall"))
        assert [lines[overall + n].split()[0] for n in (1, 2, 3)] == ["a", "b", "c"]
