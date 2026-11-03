"""
AgreementParts: the observed and chance terms, and the coefficient they make.
"""

import math

import pytest

from interrater.analyzers.agreement_functions.group_agreement import (
    agreement_corrected_by_chance,
)
from interrater.intermediates.results import AgreementParts


def _agreement(po=0.8, pe=0.5, weighted=False):
    return AgreementParts(observed=po, chance=pe, kind="agreement", weighted=weighted)


def _disagreement(do=0.2, de=0.5):
    return AgreementParts(observed=do, chance=de, kind="disagreement", weighted=False)


class TestCoefficient:

    def test_agreement(self):
        assert _agreement(0.8, 0.5).coefficient() == pytest.approx(0.6)

    def test_disagreement(self):
        """Krippendorff: 1 - Do / De."""
        assert _disagreement(0.2, 0.5).coefficient() == pytest.approx(0.6)

    def test_same_kappa_from_different_parts(self):
        """The reason to keep both terms: 0.40 can mean very different data."""
        skewed = _agreement(0.952, 0.92)
        spread = _agreement(0.70, 0.50)
        assert skewed.coefficient() == pytest.approx(0.40)
        assert spread.coefficient() == pytest.approx(0.40)

    def test_pe_of_one_makes_no_coefficient(self):
        assert math.isnan(_agreement(1.0, 1.0).coefficient())

    def test_de_of_zero_makes_no_coefficient(self):
        assert math.isnan(_disagreement(0.0, 0.0).coefficient())

    def test_pe_a_hair_below_one_makes_no_coefficient(self):
        """Same np.isclose test as agreement_corrected_by_chance."""
        assert math.isnan(_agreement(1.0, 1.0 - 1e-9).coefficient())
        assert math.isnan(agreement_corrected_by_chance(1.0, 1.0 - 1e-9))

    def test_nan_terms_make_no_coefficient(self):
        assert math.isnan(_agreement(float("nan"), float("nan")).coefficient())


class TestReproduces:

    def test_matching_point(self):
        assert _agreement(0.8, 0.5).reproduces(0.6)

    def test_wrong_point(self):
        assert not _agreement(0.8, 0.5).reproduces(0.61)

    def test_nan_matches_nan_only(self):
        undefined = _agreement(1.0, 1.0)
        assert undefined.reproduces(float("nan"))
        assert not undefined.reproduces(0.0)
        assert not _agreement(0.8, 0.5).reproduces(float("nan"))


class TestOutputs:

    def test_agreement_labels(self):
        assert str(_agreement(0.8, 0.5)) == "Po = 0.800, Pe = 0.500"

    def test_disagreement_labels(self):
        assert str(_disagreement(0.2, 0.5)) == "Do = 0.200, De = 0.500"

    def test_weighted_says_so(self):
        assert str(_agreement(0.8, 0.5, weighted=True)) == (
            "weighted Po = 0.800, weighted Pe = 0.500"
        )

    def test_summary(self):
        assert _agreement(0.8, 0.5).summary() == {
            "kind": "agreement", "observed": 0.8, "chance": 0.5, "weighted": False,
        }


class TestSelfChecks:

    def test_unknown_kind_refused(self):
        with pytest.raises(ValueError, match="kind"):
            AgreementParts(observed=0.8, chance=0.5, kind="agree", weighted=False)

    def test_agreement_outside_zero_one_refused(self):
        with pytest.raises(ValueError, match=r"Po must be in \[0, 1\]"):
            _agreement(1.2, 0.5)

    def test_negative_disagreement_refused(self):
        with pytest.raises(ValueError, match="De must be >= 0"):
            _disagreement(0.2, -0.1)

    def test_weighted_must_be_bool(self):
        with pytest.raises(TypeError, match="weighted must be bool"):
            AgreementParts(observed=0.8, chance=0.5, kind="agreement", weighted=1)

    def test_terms_must_be_numbers(self):
        with pytest.raises(TypeError, match="observed must be a number"):
            AgreementParts(observed="0.8", chance=0.5, kind="agreement", weighted=False)
