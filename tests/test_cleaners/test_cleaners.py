"""
Tests for the cleaners sub-package.

The cleaner takes a messy DataFrame and produces a validated Ratings object
with a CleaningReport receipt for every decision. Tests are organised
bottom-up: the two delegate config objects (FrameLayout, InsufficientDataPolicy)
first, then CleanerConfig that holds them, then the RatingsCleaner pipeline.

Key behaviours:
- Defaults produce a valid object with strict, reproducible settings.
- Type guards catch booleans passed as ints, floats passed as ints, etc.
- from_dict rejects unknown keys so a YAML typo is caught at load time.
- build_from_yaml round-trips through a YAML file on disk.
- Equality compares all fields, not identity.
- CleanerConfig's double-claim check uses folded matching, so a
  case-insensitive config catches "NA" in missing_values and "na" in
  value_map even though the strings differ.
- A path instead of a DataFrame gets a helpful message.
- Column names must match declared raters (when use_column_names=True).
- Encoding folds case, strips whitespace, matches missing markers and
  repairs, and rejects or absorbs unknowns.
- Config is checked against semantics before data is touched.
- Unused categories are caught or kept per policy.
- Raters and subjects below thresholds are dropped per method.
- The report counts everything: nulls, missing, repairs, unknowns, drops.
"""

import numpy as np
import pandas as pd
import pytest

from interrater.contracts import MIN_RATINGS_PER_RATER, MIN_RATINGS_PER_SUBJECT
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.cleaners.cleaner_config import CleanerConfig
from interrater.cleaners.frame_layout import FrameLayout
from interrater.cleaners.insufficient_data_policy import InsufficientDataPolicy
from interrater.cleaners.ratings_cleaner import RatingsCleaner

# Mirror the category names from conftest so tests can construct data inline.
CAT_NAMES: tuple[str, ...] = ("include", "exclude")


# ============================================================================ #
#  FrameLayout                                                                 #
# ============================================================================ #


class TestFrameLayoutDefaults:
    """Default construction gives the strict, trust-the-header setting."""

    def test_defaults(self):
        fl = FrameLayout()
        assert fl.use_column_names is True
        assert fl.subject_id_col is None

    def test_custom(self):
        fl = FrameLayout(use_column_names=False, subject_id_col="abstract_id")
        assert fl.use_column_names is False
        assert fl.subject_id_col == "abstract_id"


class TestFrameLayoutValidation:
    """Type guards on the two fields."""

    def test_use_column_names_int_rejected(self):
        """1 is truthy but not True."""
        with pytest.raises(TypeError, match="boolean"):
            FrameLayout(use_column_names=1)

    def test_subject_id_col_int_rejected(self):
        """A column name must be a string, not an index."""
        with pytest.raises(TypeError, match="column name"):
            FrameLayout(subject_id_col=0)


class TestFrameLayoutFromDict:
    """from_dict round-trips and rejects unknowns."""

    def test_round_trip(self):
        spec = {"use_column_names": False, "subject_id_col": "doc_id"}
        fl = FrameLayout.from_dict(spec)
        assert fl.use_column_names is False
        assert fl.subject_id_col == "doc_id"

    def test_unknown_key(self):
        with pytest.raises(ValueError, match="unknown key"):
            FrameLayout.from_dict({"use_column_names": True, "typo": 1})


class TestFrameLayoutEquality:
    """Equality compares fields, not identity."""

    def test_equal(self):
        assert FrameLayout() == FrameLayout()

    def test_not_equal(self):
        assert FrameLayout() != FrameLayout(use_column_names=False)

    def test_not_comparable(self):
        assert FrameLayout().__eq__("not a layout") is NotImplemented


# ============================================================================ #
#  InsufficientDataPolicy                                                      #
# ============================================================================ #


class TestInsufficientDataPolicyDefaults:
    """Default construction uses the contract floors and raters_first."""

    def test_defaults(self):
        p = InsufficientDataPolicy()
        assert p.min_ratings_per_rater == MIN_RATINGS_PER_RATER
        assert p.min_ratings_per_subject == MIN_RATINGS_PER_SUBJECT
        assert p.method == "raters_first"

    def test_custom(self):
        p = InsufficientDataPolicy(
            min_ratings_per_rater=5,
            min_ratings_per_subject=3,
            method="iterate",
        )
        assert p.min_ratings_per_rater == 5
        assert p.min_ratings_per_subject == 3
        assert p.method == "iterate"


class TestInsufficientDataPolicyValidation:
    """Type guards, floor enforcement, and method validation."""

    def test_bool_rejected(self):
        """True is 1 in Python but not an acceptable integer here."""
        with pytest.raises(TypeError, match="integer"):
            InsufficientDataPolicy(min_ratings_per_rater=True)

    def test_float_rejected(self):
        """A float is not an acceptable integer."""
        with pytest.raises(TypeError, match="integer"):
            InsufficientDataPolicy(min_ratings_per_rater=2.5)

    def test_below_floor_rater(self):
        with pytest.raises(ValueError, match="floor"):
            InsufficientDataPolicy(min_ratings_per_rater=0)

    def test_below_floor_subject(self):
        with pytest.raises(ValueError, match="floor"):
            InsufficientDataPolicy(min_ratings_per_subject=1)

    def test_invalid_method(self):
        with pytest.raises(ValueError, match="method"):
            InsufficientDataPolicy(method="backwards")


class TestInsufficientDataPolicyFromDict:
    """from_dict round-trips and rejects unknowns."""

    def test_round_trip(self):
        spec = {
            "min_ratings_per_rater": 3,
            "min_ratings_per_subject": 4,
            "method": "subjects_first",
        }
        p = InsufficientDataPolicy.from_dict(spec)
        assert p.min_ratings_per_rater == 3
        assert p.method == "subjects_first"

    def test_unknown_key(self):
        with pytest.raises(ValueError, match="unknown key"):
            InsufficientDataPolicy.from_dict({"method": "iterate", "typo": 1})


class TestInsufficientDataPolicyEquality:
    """Equality compares all three fields."""

    def test_equal(self):
        assert InsufficientDataPolicy() == InsufficientDataPolicy()

    def test_not_equal(self):
        assert InsufficientDataPolicy() != InsufficientDataPolicy(
            min_ratings_per_rater=10
        )

    def test_not_comparable(self):
        assert InsufficientDataPolicy().__eq__("nope") is NotImplemented


# ============================================================================ #
#  CleanerConfig                                                               #
# ============================================================================ #


class TestCleanerConfigDefaults:
    """Default construction gives strict, reproducible settings."""

    def test_defaults(self):
        cfg = CleanerConfig()
        assert isinstance(cfg.frame_layout, FrameLayout)
        assert cfg.case_sensitive is False
        assert cfg.strip_whitespace is True
        assert cfg.missing_values == ("",)
        assert cfg.value_map == {}
        assert cfg.unknown_policy == "error"
        assert cfg.unused_categories == "error"
        assert isinstance(cfg.insufficient_data_drops, InsufficientDataPolicy)

    def test_custom(self):
        cfg = CleanerConfig(
            case_sensitive=True,
            strip_whitespace=False,
            missing_values=("NA", "N/A"),
            value_map={"inc": "include"},
            unknown_policy="empty",
            unused_categories="keep",
        )
        assert cfg.case_sensitive is True
        assert cfg.strip_whitespace is False
        assert cfg.missing_values == ("NA", "N/A")
        assert cfg.value_map == {"inc": "include"}
        assert cfg.unknown_policy == "empty"
        assert cfg.unused_categories == "keep"

    def test_value_map_defensive_copy(self):
        """Mutating the original dict after construction does not affect the config."""
        original = {"inc": "include"}
        cfg = CleanerConfig(value_map=original)
        original["typo"] = "exclude"
        assert "typo" not in cfg.value_map


class TestCleanerConfigValidation:
    """Type guards and domain checks on each field."""

    def test_case_sensitive_int_rejected(self):
        with pytest.raises(TypeError, match="boolean"):
            CleanerConfig(case_sensitive=1)

    def test_strip_whitespace_int_rejected(self):
        with pytest.raises(TypeError, match="boolean"):
            CleanerConfig(strip_whitespace=0)

    def test_missing_values_int_rejected(self):
        """Missing values must be strings, not integers."""
        with pytest.raises(TypeError, match="strings"):
            CleanerConfig(missing_values=[99])

    def test_value_map_int_key_rejected(self):
        """Keys must be strings."""
        with pytest.raises(TypeError, match="string to string"):
            CleanerConfig(value_map={1: "include"})

    def test_value_map_empty_key_rejected(self):
        """An empty key belongs in missing_values, not here."""
        with pytest.raises(ValueError, match="empty key"):
            CleanerConfig(value_map={"": "include"})

    def test_value_map_conflicting_keys_rejected(self):
        """'Y' and 'y' match the same cells; they cannot mean different things."""
        with pytest.raises(ValueError, match="'Y' -> 'include' and 'y' -> 'exclude'"):
            CleanerConfig(value_map={"Y": "include", "y": "exclude"})

    def test_value_map_agreeing_keys_accepted(self):
        """The same meaning twice is a harmless duplicate."""
        cfg = CleanerConfig(value_map={"Y": "include", "y": "include"})
        assert len(cfg.value_map) == 2

    def test_value_map_case_keys_distinct_when_case_sensitive(self):
        cfg = CleanerConfig(case_sensitive=True, value_map={"Y": "include", "y": "exclude"})
        assert cfg.value_map["y"] == "exclude"

    def test_unknown_policy_bad(self):
        with pytest.raises(ValueError, match="unknown_policy"):
            CleanerConfig(unknown_policy="ignore")

    def test_unused_categories_bad(self):
        with pytest.raises(ValueError, match="unused_categories"):
            CleanerConfig(unused_categories="ignore")

    def test_double_claim_exact(self):
        """A string in both missing_values and value_map is a contradiction."""
        with pytest.raises(ValueError, match="twice"):
            CleanerConfig(missing_values=("NA",), value_map={"NA": "include"})

    def test_double_claim_folded(self):
        """Case-insensitive matching catches 'na' vs 'NA'."""
        with pytest.raises(ValueError, match="twice"):
            CleanerConfig(
                case_sensitive=False,
                missing_values=("NA",),
                value_map={"na": "include"},
            )

    def test_double_claim_case_sensitive_no_collision(self):
        """When case-sensitive, 'NA' and 'na' are different strings."""
        cfg = CleanerConfig(
            case_sensitive=True,
            missing_values=("NA",),
            value_map={"na": "include"},
        )
        assert "na" in cfg.value_map


class TestCleanerConfigFromDict:
    """from_dict delegates to sub-objects and rejects unknowns."""

    def test_flat(self):
        spec = {"case_sensitive": True, "unknown_policy": "empty"}
        cfg = CleanerConfig.from_dict(spec)
        assert cfg.case_sensitive is True
        assert cfg.unknown_policy == "empty"

    def test_nested(self):
        spec = {
            "frame_layout": {"use_column_names": False},
            "insufficient_data_drops": {"method": "iterate"},
        }
        cfg = CleanerConfig.from_dict(spec)
        assert cfg.frame_layout.use_column_names is False
        assert cfg.insufficient_data_drops.method == "iterate"

    def test_unknown_key(self):
        with pytest.raises(ValueError, match="unknown key"):
            CleanerConfig.from_dict({"typo": True})

    def test_delegate_error_chained(self):
        """A delegate's from_dict error is wrapped with the config prefix."""
        spec = {"frame_layout": {"use_column_names": "yes"}}
        with pytest.raises(TypeError, match=r"\[CleanerConfig\]"):
            CleanerConfig.from_dict(spec)


class TestCleanerConfigYaml:
    """build_from_yaml loads from a YAML file on disk."""

    def test_round_trip(self, tmp_path):
        """A valid YAML file produces the expected config."""
        path = tmp_path / "config.yaml"
        path.write_text(
            "case_sensitive: true\n"
            "unknown_policy: empty\n"
            "frame_layout:\n"
            "  use_column_names: false\n"
        )
        cfg = CleanerConfig.build_from_yaml(path)
        assert cfg.case_sensitive is True
        assert cfg.unknown_policy == "empty"
        assert cfg.frame_layout.use_column_names is False

    def test_non_dict_rejected(self, tmp_path):
        """A YAML file containing a list instead of a mapping is rejected."""
        path = tmp_path / "bad.yaml"
        path.write_text("- item1\n- item2\n")
        with pytest.raises(TypeError, match="mapping"):
            CleanerConfig.build_from_yaml(path)

    def test_missing_file(self):
        """A nonexistent path raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            CleanerConfig.build_from_yaml("/nonexistent/config.yaml")

    def test_to_dict_round_trips(self):
        """Every setting survives to_dict -> from_dict, so a report can carry it."""
        cfg = CleanerConfig(
            frame_layout=FrameLayout(use_column_names=False, subject_id_col=None),
            case_sensitive=True, strip_whitespace=False,
            missing_values=("NA", "o"), value_map={"inc": "include"},
            unknown_policy="empty", unused_categories="keep",
            insufficient_data_drops=InsufficientDataPolicy(
                min_ratings_per_rater=3, min_ratings_per_subject=2, method="iterate",
            ),
        )
        assert CleanerConfig.from_dict(cfg.to_dict()) == cfg

    def test_value_map_conflict_refused_on_load(self, tmp_path):
        """Caught when the YAML is read, before any semantics or data."""
        path = tmp_path / "config.yaml"
        path.write_text("value_map:\n  'Y': include\n  'y ': exclude\n")
        with pytest.raises(ValueError, match="match the same cells"):
            CleanerConfig.build_from_yaml(path)


class TestCleanerConfigEquality:
    """Equality compares all eight fields."""

    def test_equal(self):
        assert CleanerConfig() == CleanerConfig()

    def test_not_equal_flag(self):
        assert CleanerConfig() != CleanerConfig(case_sensitive=True)

    def test_not_equal_delegate(self):
        assert CleanerConfig() != CleanerConfig(
            frame_layout=FrameLayout(use_column_names=False)
        )

    def test_not_comparable(self):
        assert CleanerConfig().__eq__("nope") is NotImplemented


# ============================================================================ #
#  Repr                                                                        #
# ============================================================================ #


class TestRepr:
    """__repr__ returns a string that identifies the object."""

    def test_frame_layout_repr(self):
        fl = FrameLayout(use_column_names=False, subject_id_col="doc_id")
        r = repr(fl)
        assert r.startswith("FrameLayout(")
        assert "use_column_names=False" in r
        assert "doc_id" in r

    def test_insufficient_data_policy_repr(self):
        p = InsufficientDataPolicy(min_ratings_per_rater=5, method="iterate")
        r = repr(p)
        assert r.startswith("InsufficientDataPolicy(")
        assert "5" in r
        assert "iterate" in r

    def test_cleaner_config_repr(self):
        cfg = CleanerConfig(case_sensitive=True, unknown_policy="empty")
        r = repr(cfg)
        assert r.startswith("CleanerConfig(")
        assert "case_sensitive" in r
        assert "'empty'" in r


# ============================================================================ #
#  RatingsCleaner -- type guards                                               #
# ============================================================================ #


class TestTypeGuards:
    """_check_constructor_inputs rejects wrong types early."""

    def test_path_rejected(self, semantics):
        """A path gets a helpful error mentioning pandas.read_csv."""
        with pytest.raises(TypeError, match="read_csv"):
            RatingsCleaner("data.csv", semantics)

    def test_non_dataframe_rejected(self, semantics):
        """A numpy array is not a DataFrame."""
        with pytest.raises(TypeError, match="DataFrame"):
            RatingsCleaner(np.array([[0, 1]]), semantics)

    def test_non_semantics_rejected(self, clean_df):
        """A dict is not a RatingsSemantics."""
        with pytest.raises(TypeError, match="RatingsSemantics"):
            RatingsCleaner(clean_df, {"raters": 2})

    def test_non_config_rejected(self, clean_df, semantics):
        """A dict is not a CleanerConfig."""
        with pytest.raises(TypeError, match="CleanerConfig"):
            RatingsCleaner(clean_df, semantics, config={"case_sensitive": True})

    def test_duplicate_columns_rejected(self, semantics):
        """Two columns with the same name produce an ambiguous selection."""
        df = pd.DataFrame({"ann": [0, 1, 0], "bob": [1, 0, 1]})
        df.columns = ["ann", "ann"]  # force duplicate
        with pytest.raises(ValueError, match="repeated"):
            RatingsCleaner(df, semantics)


# ============================================================================ #
#  Happy path                                                                  #
# ============================================================================ #


class TestHappyPath:
    """A clean frame goes through the whole pipeline without incident."""

    def test_clean_with_column_names(self, clean_df, semantics):
        """Columns named ann, bob match declared raters -- accepted."""
        rc = RatingsCleaner(clean_df, semantics)
        assert rc.ratings.n_subjects == 4
        assert rc.ratings.n_raters == 2
        assert rc.ratings.n_categories == 2

    def test_clean_without_column_names(self, two_cat):
        """use_column_names=False takes rater names by position."""
        df = pd.DataFrame({"R1": ["include", "exclude"], "R2": ["exclude", "include"]})
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(frame_layout=FrameLayout(use_column_names=False))
        rc = RatingsCleaner(df, sem, cfg)
        assert rc.ratings.n_subjects == 2
        # Rater names come from semantics, not columns
        assert list(rc.ratings.semantics.raters) == ["ann", "bob"]

    def test_subject_id_extraction(self, two_cat):
        """subject_id_col lifts a column out as subject ids."""
        df = pd.DataFrame({
            "doc_id": ["A", "B", "C"],
            "ann":    ["include", "exclude", "include"],
            "bob":    ["exclude", "include", "exclude"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(
            frame_layout=FrameLayout(subject_id_col="doc_id"),
        )
        rc = RatingsCleaner(df, sem, cfg)
        np.testing.assert_array_equal(rc.ratings.subject_ids, ["A", "B", "C"])

    def test_none_config_uses_defaults(self, clean_df, semantics):
        """config=None takes strict defaults."""
        rc = RatingsCleaner(clean_df, semantics, config=None)
        assert rc.report.matching.n_null == 0
        assert rc.report.matching.n_missing == 0


# ============================================================================ #
#  Config vs semantics                                                         #
# ============================================================================ #


class TestConfigAgainstSemantics:
    """Checks that catch contradictions between config and semantics."""

    def test_value_map_bad_target(self, clean_df, semantics):
        """A repair into a category that does not exist is a contradiction."""
        cfg = CleanerConfig(value_map={"inc": "INCLUDE"})  # "INCLUDE" != "include"
        with pytest.raises(ValueError, match="does not have"):
            RatingsCleaner(clean_df, semantics, cfg)

    def test_missing_marker_is_category(self, clean_df, semantics):
        """'include' cannot both be a category and mean nothing."""
        cfg = CleanerConfig(missing_values=("include",))
        with pytest.raises(ValueError, match="declared category"):
            RatingsCleaner(clean_df, semantics, cfg)

    @staticmethod
    def _semantics(names):
        cats = UnorderedCategories(names=names, labels={n: n for n in names})
        return RatingsSemantics(raters=2, categories=cats)

    @pytest.mark.parametrize("names, config", [
        (("A", "a", "B"), CleanerConfig()),
        (("A", "A ", "B"), CleanerConfig(case_sensitive=True)),
    ], ids=["case", "trailing-space"])
    def test_names_matched_as_one_refused(self, names, config):
        """
        Folded alike, every answer for either would be coded as one category.

        The frame's columns are wrong on purpose: the names are refused
        before a single cell is read, so the column error never comes up.
        """
        frame = pd.DataFrame({"wrong": ["A"], "columns": ["B"]})
        with pytest.raises(ValueError, match="cannot be told apart"):
            RatingsCleaner(frame, self._semantics(names), config)

    def test_names_kept_apart_when_case_sensitive(self):
        frame = pd.DataFrame({"rater0": ["A", "a", "B"], "rater1": ["a", "A", "B"]})
        ratings = RatingsCleaner(
            frame, self._semantics(("A", "a", "B")), CleanerConfig(case_sensitive=True),
        ).ratings
        assert ratings.category_index.tolist() == [[0, 1], [1, 0], [2, 2]]


class TestFoldCollisions:
    """The helper both checks use: which different strings fold alike."""

    def test_groups_in_given_order(self):
        from interrater.cleaners._string_utils import fold_collisions

        found = fold_collisions(
            ["b", "A", "B", "a", "c"], case_sensitive=False, strip_whitespace=True,
        )
        assert found == [("b", "B"), ("A", "a")]

    def test_repeats_are_not_collisions(self):
        from interrater.cleaners._string_utils import fold_collisions

        assert fold_collisions(["A", "A"], case_sensitive=False, strip_whitespace=True) == []

    def test_respects_the_matching_rule(self):
        from interrater.cleaners._string_utils import fold_collisions

        values = ["A", "a", "B", "B "]
        assert fold_collisions(values, case_sensitive=True, strip_whitespace=False) == []
        assert fold_collisions(values, case_sensitive=True, strip_whitespace=True) == [("B", "B ")]


# ============================================================================ #
#  Column matching                                                             #
# ============================================================================ #


class TestColumnMatching:
    """_interpret_frame matches columns to declared raters."""

    def test_wrong_column_names(self, two_cat):
        """Column names that do not match declared raters are rejected."""
        df = pd.DataFrame({"R1": ["include"], "R2": ["exclude"]})
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        with pytest.raises(ValueError, match="rater columns"):
            RatingsCleaner(df, sem)

    def test_wrong_column_count_positional(self, two_cat):
        """With use_column_names=False, column count must still match."""
        df = pd.DataFrame({"a": ["include"], "b": ["exclude"], "c": ["include"]})
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(frame_layout=FrameLayout(use_column_names=False))
        with pytest.raises(ValueError, match="rater column"):
            RatingsCleaner(df, sem, cfg)

    def test_subject_id_col_missing(self, clean_df, semantics):
        """A subject_id_col that does not exist in the frame is an error."""
        cfg = CleanerConfig(
            frame_layout=FrameLayout(subject_id_col="nonexistent"),
        )
        with pytest.raises(KeyError, match="nonexistent"):
            RatingsCleaner(clean_df, semantics, cfg)

    def test_n_subjects_mismatch(self, two_cat):
        """When semantics declares n_subjects, the row count must match."""
        df = pd.DataFrame({"ann": ["include", "exclude"], "bob": ["exclude", "include"]})
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat, n_subjects=10)
        with pytest.raises(ValueError, match="row"):
            RatingsCleaner(df, sem)


# ============================================================================ #
#  Encoding                                                                    #
# ============================================================================ #


class TestEncoding:
    """_encode normalises strings and crosses the boundary to integers."""

    def test_case_folding(self, two_cat):
        """Case-insensitive matching catches INCLUDE as include."""
        df = pd.DataFrame({"ann": ["INCLUDE", "exclude"], "bob": ["exclude", "INCLUDE"]})
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(case_sensitive=False)
        rc = RatingsCleaner(df, sem, cfg)
        assert rc.ratings.n_subjects == 2

    def test_whitespace_stripping(self, two_cat):
        """Whitespace around a category name is stripped before matching."""
        df = pd.DataFrame(
            {"ann": ["  include ", "exclude"], "bob": ["exclude", " include"]},
        )
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(strip_whitespace=True)
        rc = RatingsCleaner(df, sem, cfg)
        assert rc.ratings.n_subjects == 2

    def test_null_always_missing(self, two_cat):
        """A null cell is always missing, regardless of missing_values."""
        df = pd.DataFrame({
            "ann": ["include", None, "exclude"],
            "bob": ["exclude", "include", "include"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        rc = RatingsCleaner(df, sem)
        assert rc.report.matching.n_null == 1

    def test_missing_marker_matched(self, two_cat):
        """A string in missing_values becomes EMPTY."""
        df = pd.DataFrame({
            "ann": ["include", "NA", "exclude"],
            "bob": ["exclude", "include", "include"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(missing_values=("NA",))
        rc = RatingsCleaner(df, sem, cfg)
        assert rc.report.matching.missing_matched["NA"] == 1

    def test_repair_applied(self, two_cat):
        """A value_map entry transforms a typo into the right category."""
        df = pd.DataFrame({"ann": ["inc", "exclude"], "bob": ["exclude", "include"]})
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(value_map={"inc": "include"})
        rc = RatingsCleaner(df, sem, cfg)
        assert rc.report.matching.repairs["inc"] == 1
        # "inc" was repaired to "include" which is index 0
        assert rc.ratings.category_index[0, 0] == 0

    def test_unknown_error(self, two_cat):
        """unknown_policy='error' stops on the first unknown value."""
        df = pd.DataFrame({"ann": ["???", "exclude"], "bob": ["exclude", "include"]})
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(unknown_policy="error")
        with pytest.raises(ValueError, match="neither a declared category"):
            RatingsCleaner(df, sem, cfg)

    def test_unknown_empty(self, two_cat):
        """unknown_policy='empty' treats unknown values as no rating."""
        df = pd.DataFrame({
            "ann": ["???", "exclude", "include"],
            "bob": ["exclude", "include", "exclude"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(unknown_policy="empty")
        rc = RatingsCleaner(df, sem, cfg)
        assert rc.report.matching.n_unknown == 1

    def test_unknown_empty_counts_as_no_rating_everywhere(self, two_cat):
        """
        A null, a missing marker and an unknown treated as empty are all no
        rating. Before 2.1.0 the headline counted 2 here and the per-rater
        lines 3.
        """
        df = pd.DataFrame({
            "ann": ["???", "NA", "include", "exclude"],
            "bob": [None, "include", "exclude", "include"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(unknown_policy="empty", missing_values=("NA",))
        m = RatingsCleaner(df, sem, cfg).report.matching
        assert (m.n_null, m.n_matched_missing, m.n_unknown) == (1, 1, 1)
        assert m.n_missing == 3
        assert sum(m.per_rater_missing.values()) == 3
        assert sum(k * n for k, n in m.subject_missing_distribution.items()) == 3
        assert m.summary()["n_missing"] == 3

    def test_case_sensitive_encoding(self, two_cat):
        """case_sensitive=True: 'Include' does not match 'include'."""
        df = pd.DataFrame({
            "ann": ["Include", "exclude", "include"],
            "bob": ["exclude", "include", "exclude"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(case_sensitive=True, unknown_policy="error")
        with pytest.raises(ValueError, match="neither a declared category"):
            RatingsCleaner(df, sem, cfg)

    def test_no_strip_whitespace(self, two_cat):
        """strip_whitespace=False: ' include ' stays as-is and is unknown."""
        df = pd.DataFrame({
            "ann": [" include ", "exclude", "include"],
            "bob": ["exclude",  "include", "exclude"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(strip_whitespace=False, unknown_policy="error")
        with pytest.raises(ValueError, match="neither a declared category"):
            RatingsCleaner(df, sem, cfg)


# ============================================================================ #
#  Unused categories                                                           #
# ============================================================================ #


class TestUnusedCategories:
    """_resolve_unused_categories enforces the policy."""

    def _three_cat_frame(self):
        """Helper: 2 subjects, 2 raters, only 'mild' and 'severe' used."""
        cats = UnorderedCategories(("mild", "moderate", "severe"))
        df = pd.DataFrame({
            "ann": ["mild", "severe"],
            "bob": ["severe", "mild"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=cats)
        return df, sem

    def test_unused_error(self):
        """unused_categories='error' stops when a category is never chosen."""
        df, sem = self._three_cat_frame()
        cfg = CleanerConfig(unused_categories="error")
        with pytest.raises(ValueError, match="never chosen"):
            RatingsCleaner(df, sem, cfg)

    def test_unused_keep(self):
        """unused_categories='keep' preserves K even when a category is empty."""
        df, sem = self._three_cat_frame()
        cfg = CleanerConfig(unused_categories="keep")
        rc = RatingsCleaner(df, sem, cfg)
        assert rc.ratings.n_categories == 3
        assert rc.report.categories.unused == ("moderate",)
        assert rc.report.categories.action == "kept"


# ============================================================================ #
#  Drops                                                                       #
# ============================================================================ #


class TestDrops:
    """_drop removes raters and subjects below thresholds."""

    def _thin_rater_frame(self):
        """
        Helper: 3 raters, 3 subjects. Cara has only 1 non-null rating.

            subject  | ann     | bob     | cara
            ---------+---------+---------+---------
            sub0     | include | exclude | None
            sub1     | exclude | include | None
            sub2     | include | exclude | include   <- cara's only rating
        """
        cats = UnorderedCategories(("include", "exclude"))
        df = pd.DataFrame({
            "ann":  ["include", "exclude", "include"],
            "bob":  ["exclude", "include", "exclude"],
            "cara": [None,      None,      "include"],
        })
        sem = RatingsSemantics(raters=["ann", "bob", "cara"], categories=cats)
        return df, sem

    def test_drop_rater_first(self):
        """raters_first drops cara (1 rating < min=2) before checking subjects."""
        df, sem = self._thin_rater_frame()
        cfg = CleanerConfig(unknown_policy="empty")
        rc = RatingsCleaner(df, sem, cfg)
        # Cara dropped, 2 raters survive
        assert rc.ratings.n_raters == 2
        assert "cara" in rc.report.drops.dropped_raters
        assert rc.report.semantics_changed is True

    def test_report_records_the_rule_it_applied(self):
        """The drop report carries the configured minimums, not defaults of its own."""
        df, sem = self._thin_rater_frame()
        policy = InsufficientDataPolicy(min_ratings_per_rater=3, min_ratings_per_subject=2)
        drops = RatingsCleaner(
            df, sem, CleanerConfig(insufficient_data_drops=policy),
        ).report.drops
        assert (drops.min_ratings_per_rater, drops.min_ratings_per_subject) == (3, 2)
        assert "at least 3 ratings per rater and 2 per subject" in drops.render()

    def test_drop_subjects_first(self):
        """subjects_first sweeps subjects before raters."""
        df, sem = self._thin_rater_frame()
        cfg = CleanerConfig(
            unknown_policy="empty",
            insufficient_data_drops=InsufficientDataPolicy(method="subjects_first"),
        )
        rc = RatingsCleaner(df, sem, cfg)
        # Under subjects_first with min_per_subject=2: every subject has
        # at least 2 ratings from ann and bob, so no subjects dropped.
        # Then raters: cara has 1 rating < min=2, so cara is dropped.
        assert rc.ratings.n_raters == 2
        assert "cara" in rc.report.drops.dropped_raters

    def test_nothing_dropped_on_clean_data(self, clean_df, semantics):
        """Clean data with everyone above thresholds drops nothing."""
        rc = RatingsCleaner(clean_df, semantics)
        assert rc.report.drops.dropped_nothing is True
        assert rc.report.n_raters_lost == 0
        assert rc.report.n_subjects_lost == 0

    def test_drop_iterate(self, two_cat):
        """
        iterate cascades across three passes.

            sub0: ann=include  bob=exclude  cara=None
            sub1: ann=exclude  bob=include  cara=None
            sub2: ann=include  bob=None     cara=include
            sub3: ann=None     bob=None     cara=include

        Pass 1: no raters below 2, sub3 has 1 rating -> drop sub3.
        Pass 2: cara now has 1 rating -> drop cara; sub2 now has 1 -> drop sub2.
        Pass 3: ann=2, bob=2, sub0=2, sub1=2 -> converged.
        """
        df = pd.DataFrame({
            "ann":  ["include", "exclude", "include", None],
            "bob":  ["exclude", "include", None,      None],
            "cara": [None,      None,      "include", "include"],
        })
        sem = RatingsSemantics(
            raters=["ann", "bob", "cara"], categories=two_cat,
        )
        cfg = CleanerConfig(
            insufficient_data_drops=InsufficientDataPolicy(method="iterate"),
        )
        rc = RatingsCleaner(df, sem, cfg)

        assert rc.ratings.n_subjects == 2
        assert rc.ratings.n_raters == 2

        d = rc.report.drops
        assert d.method == "iterate"
        assert d.converged is True
        assert d.n_passes == 3
        # Pass 1: 0 raters, 1 subject (sub3)
        assert d.passes[0].n_raters == 0
        assert d.passes[0].n_subjects == 1
        # Pass 2: 1 rater (cara), 1 subject (sub2)
        assert d.passes[1].n_raters == 1
        assert d.passes[1].n_subjects == 1
        # Pass 3: nothing -> converged
        assert d.passes[2].dropped_nothing is True
        # Counts at removal, not in the original array: cara had 2 ratings to
        # start with but went in pass 2 holding 1, after sub3 was removed.
        # Before 2.1.0 this reported 2, against a minimum of 2.
        assert d.dropped_raters == {"cara": 1}
        # sub3 went in pass 1 with 1 rating; sub2 in pass 2, after cara, with
        # 1. The original counts would say {1: 1, 2: 1}.
        assert d.dropped_subject_distribution == {1: 2}
        assert d.n_subjects_dropped == 2
        assert set(d.dropped_subject_ids) == {"subject2", "subject3"}

        # Survived: 2 of 4 input rows
        assert rc.report.survived == 0.5
        assert rc.report.n_subjects_lost == 2
        assert rc.report.n_raters_lost == 1
        assert rc.report.semantics_changed is True

    # Each case survives cleaning and cascades: before 2.1.0 every one of them
    # reported a casualty holding at least the minimum.
    @pytest.mark.parametrize("method, seed", [
        ("raters_first", 8), ("raters_first", 10), ("raters_first", 11),
        ("subjects_first", 1), ("subjects_first", 2), ("subjects_first", 6),
        ("subjects_first", 7), ("subjects_first", 8),
        ("iterate", 1), ("iterate", 6), ("iterate", 7),
    ])
    def test_every_casualty_held_less_than_the_minimum(self, method, seed):
        """
        What the report says a casualty held is what it was judged on, so it
        is always below the minimum. Reporting counts from before the sweep
        broke this whenever one axis's removals thinned the other.
        """
        rng = np.random.default_rng(seed)
        values = np.where(rng.random((40, 6)) < 0.4, None,
                          rng.choice(["include", "exclude"], (40, 6)))
        df = pd.DataFrame(values, columns=[f"r{i}" for i in range(6)])
        df = df[df.notna().sum(axis=1) > 0]
        sem = RatingsSemantics(
            raters=list(df.columns),
            categories=UnorderedCategories(("include", "exclude")),
        )
        policy = InsufficientDataPolicy(
            min_ratings_per_rater=20, min_ratings_per_subject=3, method=method,
        )
        cfg = CleanerConfig(insufficient_data_drops=policy, unused_categories="keep")
        d = RatingsCleaner(df, sem, cfg).report.drops
        assert all(n < 20 for n in d.dropped_raters.values()), d.dropped_raters
        assert all(k < 3 for k in d.dropped_subject_distribution), d.dropped_subject_distribution

    def test_drop_iterate_nothing_dropped(self, clean_df, semantics):
        """iterate on clean data: one pass that drops nothing, converged."""
        cfg = CleanerConfig(
            insufficient_data_drops=InsufficientDataPolicy(method="iterate"),
        )
        rc = RatingsCleaner(clean_df, semantics, cfg)

        assert rc.report.drops.dropped_nothing is True
        assert rc.report.drops.converged is True
        assert rc.report.drops.n_passes == 1
        assert rc.report.drops.passes[0].dropped_nothing is True

    def test_drop_cascade_to_empty(self, two_cat):
        """A frame where every rater is below threshold cascades to nothing."""
        df = pd.DataFrame({
            "ann": ["include", None],
            "bob": [None,      "exclude"],
        })
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        # ann has 1 rating, bob has 1 -> both dropped -> 0 subjects
        # -> all categories unused -> "never chosen" error
        with pytest.raises(ValueError, match="never chosen"):
            RatingsCleaner(df, sem)

    def test_drop_report_details(self):
        """Drop report records rater counts, pass details, and convergence."""
        df, sem = self._thin_rater_frame()
        rc = RatingsCleaner(df, sem)
        d = rc.report.drops

        assert d.method == "raters_first"
        assert d.converged is True
        assert d.n_passes == 1
        assert d.passes[0].number == 1
        assert d.passes[0].n_raters == 1
        assert d.passes[0].n_subjects == 0
        # Cara had 1 rating in the original array
        assert d.dropped_raters == {"cara": 1}
        assert d.dropped_subject_ids == ()
        assert d.n_raters_dropped == 1
        assert d.n_subjects_dropped == 0
        assert d.dropped_nothing is False


# ============================================================================ #
#  Report correctness                                                          #
# ============================================================================ #


class TestReport:
    """The CleaningReport is the receipt for every decision."""

    def test_shape_in(self, clean_df, semantics):
        """shape_in reflects the frame as handed in."""
        rc = RatingsCleaner(clean_df, semantics)
        assert rc.report.shape_in == (4, 2)

    def test_shape_out(self, clean_df, semantics):
        """shape_out reflects the Ratings that came out."""
        rc = RatingsCleaner(clean_df, semantics)
        assert rc.report.shape_out == (4, 2)

    def test_survived(self, clean_df, semantics):
        """All rows survived on clean data."""
        rc = RatingsCleaner(clean_df, semantics)
        assert rc.report.survived == 1.0

    def test_semantics_unchanged_on_clean(self, clean_df, semantics):
        """No raters or categories dropped means semantics did not change."""
        rc = RatingsCleaner(clean_df, semantics)
        assert rc.report.semantics_changed is False

    def test_matching_counts_clean(self, clean_df, semantics):
        """On clean data: 8 cells, 0 null, 0 missing, 0 repaired."""
        rc = RatingsCleaner(clean_df, semantics)
        m = rc.report.matching
        assert m.n_cells == 8
        assert m.n_null == 0
        assert m.n_missing == 0
        assert m.n_repaired == 0
        assert m.n_unknown == 0

    def test_matching_counts_messy(self, two_cat, messy_df):
        """
        MESSY_DF has 8 cells: 1 null (bob row 3), 1 missing marker (ann "NA"),
        1 repair (ann "inc"), 0 unknowns (we configure repair + missing).

            ann: include, NA, inc, exclude
            bob: exclude, exclude, include, None
        """
        sem = RatingsSemantics(raters=["ann", "bob"], categories=two_cat)
        cfg = CleanerConfig(
            missing_values=("NA",),
            value_map={"inc": "include"},
        )
        rc = RatingsCleaner(messy_df, sem, cfg)
        m = rc.report.matching
        assert m.n_cells == 8
        assert m.n_null == 1
        assert m.missing_matched["NA"] == 1
        assert m.repairs["inc"] == 1
        assert m.n_unknown == 0

    def test_category_counts(self, clean_df, semantics):
        """
        CLEAN_DF: ann=[include, exclude, include, exclude],
                  bob=[exclude, exclude, include, include].
        include: 4 times, exclude: 4 times.
        """
        rc = RatingsCleaner(clean_df, semantics)
        counts = rc.report.categories.counts
        assert counts["include"] == 4
        assert counts["exclude"] == 4

    def test_multiple_problems_per_rater(self, two_cat):
        """
        Mixed problems in one frame: per-rater breakdowns are correct.

            ann:  include, NA,      inc,     exclude, include, exclude
            bob:  exclude, exclude, include, None,    ???,     include
            cara: include, include, exclude, include, exclude, exclude
        """
        df = pd.DataFrame({
            "ann":  ["include", "NA",      "inc",     "exclude", "include", "exclude"],
            "bob":  ["exclude", "exclude", "include", None,      "???",     "include"],
            "cara": ["include", "include", "exclude", "include", "exclude", "exclude"],
        })
        sem = RatingsSemantics(
            raters=["ann", "bob", "cara"], categories=two_cat,
        )
        cfg = CleanerConfig(
            missing_values=("NA",),
            value_map={"inc": "include"},
            unknown_policy="empty",
        )
        rc = RatingsCleaner(df, sem, cfg)
        m = rc.report.matching

        assert m.n_cells == 18
        assert m.n_null == 1                     # bob sub3
        assert m.missing_matched["NA"] == 1      # ann sub1
        assert m.repairs["inc"] == 1             # ann sub2
        assert m.unknown_values["???"] == 1      # bob sub4

        # Per-rater null: only bob had a null
        assert m.per_rater_null["ann"] == 0
        assert m.per_rater_null["bob"] == 1
        assert m.per_rater_null["cara"] == 0

        # Per-rater unknown: only bob had an unknown
        assert m.per_rater_unknown["ann"] == 0
        assert m.per_rater_unknown["bob"] == 1
        assert m.per_rater_unknown["cara"] == 0

        # Per-rater missing (total EMPTY cells per column):
        # ann: 1 (NA), bob: 2 (null + unknown), cara: 0
        assert m.per_rater_missing["ann"] == 1
        assert m.per_rater_missing["bob"] == 2
        assert m.per_rater_missing["cara"] == 0

        # Subject missing distribution: 3 subjects with 0 missing, 3 with 1
        assert m.subject_missing_distribution == {0: 3, 1: 3}
