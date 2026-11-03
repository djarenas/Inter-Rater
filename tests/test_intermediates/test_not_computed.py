"""
NotComputed: a value that does not apply here, with the reason.

It replaces None for fields that can legitimately be absent, so a reader can
tell "does not apply" from "was never filled in".
"""

import dataclasses

import pytest

from interrater import NotComputed

REASON = "bias index needs each rater's ratings; built from category counts only"


class TestNotComputed:

    def test_holds_its_reason(self):
        assert NotComputed(REASON).reason == REASON

    def test_prints_its_reason(self):
        assert str(NotComputed(REASON)) == f"not computed: {REASON}"

    def test_summary_carries_the_reason(self):
        assert NotComputed(REASON).summary() == {"not_computed": REASON}

    def test_empty_reason_refused(self):
        with pytest.raises(ValueError, match="must say why"):
            NotComputed("   ")

    def test_non_string_reason_refused(self):
        with pytest.raises(TypeError, match="reason must be str"):
            NotComputed(None)

    def test_has_no_truth_value(self):
        """`if value:` would mistake it for a number; refuse instead."""
        with pytest.raises(TypeError, match="isinstance"):
            bool(NotComputed(REASON))

    def test_is_frozen(self):
        with pytest.raises(dataclasses.FrozenInstanceError):
            NotComputed(REASON).reason = "other"

    def test_equal_when_reasons_match(self):
        assert NotComputed(REASON) == NotComputed(REASON)
        assert NotComputed(REASON) != NotComputed("another reason")
