"""
Tests for MatchingReport, DropPass, DropReport, CategoryReport, and
CleaningReport.

Tested bottom-up: the three sub-reports first (they are delegates of
CleaningReport), then CleaningReport itself.

Key behaviours:
- Frozen: every report is immutable after construction.
- Properties derive from stored counts, never from data.
- summary() returns plain types only, so the record can be serialised.
- render() produces a human-readable string; verbosity controls detail,
  never content.
- __repr__ delegates to render().
"""

import pytest

from interrater.reports.matching_report import MatchingReport
from interrater.reports.drop_report import DropPass, DropReport
from interrater.reports.category_report import CategoryReport
from interrater.reports.cleaning_report import CleaningReport


# ============================================================================ #
#  DropPass                                                                    #
# ============================================================================ #


class TestDropPass:
    """DropPass is a single sweep of the drop loop."""

    def test_construction(self):
        p = DropPass(number=1, n_raters=2, n_subjects=3)
        assert p.number == 1
        assert p.n_raters == 2
        assert p.n_subjects == 3

    def test_dropped_nothing_true(self):
        p = DropPass(number=1, n_raters=0, n_subjects=0)
        assert p.dropped_nothing is True

    def test_dropped_nothing_false_raters(self):
        p = DropPass(number=1, n_raters=1, n_subjects=0)
        assert p.dropped_nothing is False

    def test_dropped_nothing_false_subjects(self):
        p = DropPass(number=1, n_raters=0, n_subjects=1)
        assert p.dropped_nothing is False

    def test_frozen(self):
        p = DropPass(number=1, n_raters=0, n_subjects=0)
        with pytest.raises(AttributeError):
            p.number = 2


# ============================================================================ #
#  DropReport                                                                  #
# ============================================================================ #


class TestDropReportDefaults:
    """Default construction means nothing was dropped."""

    def test_defaults(self, drops_nothing):
        r = drops_nothing
        assert r.method == "raters_first"
        assert r.converged is True
        assert r.passes == ()
        assert r.dropped_raters == {}
        assert r.dropped_subject_ids == ()
        assert r.dropped_subject_distribution == {}


class TestDropReportProperties:
    """Derived properties count what the stored data holds."""

    def test_n_passes(self, drops_some):
        assert drops_some.n_passes == 1

    def test_n_raters_dropped(self, drops_some):
        assert drops_some.n_raters_dropped == 1

    def test_n_subjects_dropped(self, drops_some):
        assert drops_some.n_subjects_dropped == 2

    def test_dropped_nothing_true(self, drops_nothing):
        assert drops_nothing.dropped_nothing is True

    def test_dropped_nothing_false(self, drops_some):
        assert drops_some.dropped_nothing is False


class TestDropReportSummary:
    """summary() returns plain types with the right keys."""

    def test_keys(self, drops_some):
        s = drops_some.summary()
        expected = {
            "method", "min_ratings_per_rater", "min_ratings_per_subject",
            "converged", "n_passes", "passes",
            "n_raters_dropped", "n_subjects_dropped",
            "dropped_raters", "dropped_subject_ids",
            "dropped_subject_distribution",
        }
        assert set(s.keys()) == expected

    def test_passes_are_dicts(self, drops_some):
        s = drops_some.summary()
        assert isinstance(s["passes"], list)
        assert all(isinstance(p, dict) for p in s["passes"])
        assert s["passes"][0] == {"pass": 1, "raters": 1, "subjects": 2}

    def test_dropped_subject_ids_are_list(self, drops_some):
        """summary() converts tuples to lists for JSON compatibility."""
        s = drops_some.summary()
        assert isinstance(s["dropped_subject_ids"], list)

    def test_nothing_dropped_summary(self, drops_nothing):
        s = drops_nothing.summary()
        assert s["n_raters_dropped"] == 0
        assert s["n_subjects_dropped"] == 0
        assert s["passes"] == []


class TestDropReportRender:
    """render() produces readable text at two verbosity levels."""

    def test_nothing_dropped_shortcut(self, drops_nothing):
        """Says what was checked, and that nobody of either kind fell below it."""
        text = drops_nothing.render()
        assert "rule       : at least 2 ratings per rater and 2 per subject (raters_first)" in text
        assert "no raters or subjects dropped" in text

    def test_rule_shown_when_something_dropped(self, drops_some):
        text = drops_some.render()
        assert "rule       : at least 2 ratings per rater and 2 per subject (" in text

    def test_rule_in_summary(self, drops_some):
        s = drops_some.summary()
        assert (s["min_ratings_per_rater"], s["min_ratings_per_subject"]) == (2, 2)

    def test_summary_shows_counts(self, drops_some):
        text = drops_some.render("summary")
        assert "1 removed" in text        # raters
        assert "2 removed" in text         # subjects
        assert "carol" in text

    def test_full_shows_pass_detail(self):
        """Full verbosity shows per-pass breakdown when there are multiple."""
        r = DropReport(
            method="iterate",
            min_ratings_per_rater=2, min_ratings_per_subject=2,
            converged=True,
            passes=(DropPass(1, 1, 0), DropPass(2, 0, 3), DropPass(3, 0, 0)),
            dropped_raters={"carol": 1},
            dropped_subject_ids=("a", "b", "c"),
            dropped_subject_distribution={1: 3},
        )
        text = r.render("full")
        assert "pass 1" in text
        assert "pass 2" in text

    def test_full_shows_subject_ids(self, drops_some):
        text = drops_some.render("full")
        assert "doc_E" in text
        assert "doc_F" in text

    def test_non_converged_text(self):
        r = DropReport(
            method="iterate",
            min_ratings_per_rater=2, min_ratings_per_subject=2,
            converged=False,
            passes=(DropPass(1, 0, 5),),
            dropped_subject_ids=tuple(f"s{i}" for i in range(5)),
            dropped_subject_distribution={1: 5},
        )
        text = r.render()
        assert "STOPPED AT THE PASS LIMIT" in text

    def test_subject_ids_with_ints(self):
        """Integer subject IDs render without error (widened type)."""
        r = DropReport(
            method="raters_first",
            min_ratings_per_rater=2, min_ratings_per_subject=2,
            passes=(DropPass(1, 0, 2),),
            dropped_subject_ids=(0, 1),
            dropped_subject_distribution={1: 2},
        )
        text = r.render("full")
        assert "0" in text
        assert "1" in text

    def test_more_than_three_raters_truncated(self):
        """Only the first 3 raters are shown, with a +N more suffix."""
        raters = {f"r{i}": i + 1 for i in range(5)}
        r = DropReport(
            method="raters_first",
            min_ratings_per_rater=2, min_ratings_per_subject=2,
            passes=(DropPass(1, 5, 0),),
            dropped_raters=raters,
        )
        text = r.render()
        assert "+2 more" in text


# ============================================================================ #
#  MatchingReport                                                              #
# ============================================================================ #


class TestMatchingReportProperties:
    """Derived properties sum the stored dicts."""

    def test_n_missing_clean(self, matching_clean):
        assert matching_clean.n_missing == 0

    def test_n_missing_messy(self, matching_messy):
        # 1 null + 1 matched "NA" + 1 unknown treated as empty
        assert matching_messy.n_missing == 3

    def test_n_repaired(self, matching_messy):
        assert matching_messy.n_repaired == 1

    def test_n_unknown(self, matching_messy):
        assert matching_messy.n_unknown == 1

    def test_n_repaired_zero(self, matching_clean):
        assert matching_clean.n_repaired == 0

    def test_n_unknown_zero(self, matching_clean):
        assert matching_clean.n_unknown == 0


class TestMatchingReportSummary:
    """summary() returns plain types with the right keys."""

    def test_keys(self, matching_messy):
        s = matching_messy.summary()
        expected = {
            "n_cells", "n_null", "n_matched_missing", "n_missing", "n_repaired", "n_unknown",
            "missing_matched", "repairs", "unknown_values",
            "per_rater_null", "per_rater_missing", "per_rater_unknown",
            "subject_missing_distribution",
        }
        assert set(s.keys()) == expected

    def test_derived_values_in_summary(self, matching_messy):
        s = matching_messy.summary()
        assert s["n_missing"] == 3
        assert s["n_repaired"] == 1
        assert s["n_unknown"] == 1

    def test_clean_summary(self, matching_clean):
        s = matching_clean.summary()
        assert s["n_missing"] == 0
        assert s["repairs"] == {}


class TestMatchingReportRender:
    """render() produces readable text at two verbosity levels."""

    def test_summary_shows_cell_count(self, matching_clean):
        text = matching_clean.render()
        assert "8" in text

    def test_repairs_shown_when_present(self, matching_messy):
        text = matching_messy.render()
        assert "repaired" in text
        assert "'inc'" in text

    def test_repairs_absent_when_none_configured(self, matching_clean):
        text = matching_clean.render()
        assert "repaired" not in text

    def test_unknowns_shown(self, matching_messy):
        text = matching_messy.render()
        assert "unknown" in text
        assert "'???'" in text

    def test_full_shows_per_rater(self, matching_messy):
        text = matching_messy.render("full")
        assert "ann" in text
        assert "bob" in text

    def test_full_shows_subject_distribution(self, matching_messy):
        text = matching_messy.render("full")
        assert ("no rating per subject (before drops): 2 subjects with 0, "
                "1 subject with 1, 1 subject with 2") in text

    def test_division_by_zero_guard(self):
        """Zero cells does not raise."""
        r = MatchingReport(n_cells=0, n_null=0)
        text = r.render()
        assert "0" in text

    def test_headline_lists_parts_that_add_up(self, matching_messy):
        """'No rating' is null + marker + unknown, and says so."""
        text = matching_messy.render()
        assert ("no rating  : 3 (37.5%): 1 arrived null, 1 matched a missing "
                "marker, 1 unknown treated as empty") in text

    def test_headline_is_just_zero_when_nothing_missing(self, matching_clean):
        lines = matching_clean.render().splitlines()
        assert "    no rating  : 0" in lines

    def test_headline_omits_unknown_part_when_none(self):
        r = MatchingReport(n_cells=4, n_null=1, missing_matched={"NA": 1})
        assert "unknown" not in r.render()

    def test_per_rater_unknown_is_part_of_the_count(self, matching_messy):
        # 8 cells over 2 raters: 4 rows each.
        assert "ann: 2 rated, 2 with no rating, 1 of them unknown" in matching_messy.render("full")


class TestMatchingReportAddsUp:
    """
    "No rating" means one thing everywhere: a cell that ended up empty. Before
    2.1.0 the headline left out unknowns treated as empty while the per-rater
    lines counted them, so one report gave two totals. The report now refuses
    to exist when its breakdowns disagree with its totals.
    """

    def _kwargs(self, **overrides):
        kwargs = dict(
            n_cells=8, n_null=1, missing_matched={"NA": 1},
            unknown_values={"???": 1},
            per_rater_null={"bob": 1},
            per_rater_missing={"ann": 2, "bob": 1},
            per_rater_unknown={"ann": 1},
            subject_missing_distribution={0: 2, 1: 1, 2: 1},
        )
        kwargs.update(overrides)
        return kwargs

    def test_consistent_report_accepted(self):
        assert MatchingReport(**self._kwargs()).n_missing == 3

    @pytest.mark.parametrize("overrides, name", [
        (dict(per_rater_missing={"ann": 1, "bob": 1}), "per_rater_missing"),
        (dict(per_rater_null={"bob": 2}), "per_rater_null"),
        (dict(per_rater_unknown={"ann": 0}), "per_rater_unknown"),
        (dict(subject_missing_distribution={0: 3, 2: 1}), "subject_missing_distribution"),
    ])
    def test_breakdown_that_disagrees_is_refused(self, overrides, name):
        with pytest.raises(ValueError, match=f"{name} adds up to"):
            MatchingReport(**self._kwargs(**overrides))

    def test_more_empty_cells_than_cells_refused(self):
        with pytest.raises(ValueError, match="out of 2 cells"):
            MatchingReport(n_cells=2, n_null=3)

    def test_breakdowns_not_recorded_are_not_checked(self):
        """A report with only totals is incomplete, not inconsistent."""
        assert MatchingReport(n_cells=8, n_null=2).n_missing == 2


# ============================================================================ #
#  CategoryReport                                                              #
# ============================================================================ #


class TestCategoryReportRare:
    """The rare property filters categories by threshold."""

    def test_no_rare(self, categories_all_used):
        assert categories_all_used.rare == {}

    def test_rare_detected(self):
        r = CategoryReport(
            counts={"include": 200, "exclude": 3},
            policy="error",
            rare_threshold=10,
        )
        assert r.rare == {"exclude": 3}

    def test_zero_not_rare(self):
        """Zero-count categories are unused, not rare (0 < n is False)."""
        r = CategoryReport(
            counts={"include": 200, "exclude": 0},
            policy="keep",
            unused=("exclude",),
        )
        assert r.rare == {}

    def test_custom_threshold(self):
        r = CategoryReport(
            counts={"include": 200, "exclude": 50},
            policy="error",
            rare_threshold=100,
        )
        assert r.rare == {"exclude": 50}


class TestCategoryReportSummary:
    """summary() returns plain types with the right keys."""

    def test_keys(self, categories_all_used):
        s = categories_all_used.summary()
        expected = {
            "counts", "unused", "rare", "rare_threshold", "policy", "action",
        }
        assert set(s.keys()) == expected

    def test_unused_is_list(self, categories_with_unused):
        """summary() converts tuple to list for JSON compatibility."""
        s = categories_with_unused.summary()
        assert isinstance(s["unused"], list)
        assert s["unused"] == ["maybe"]


class TestCategoryReportRender:
    """render() produces readable text at two verbosity levels."""

    def test_summary_shows_percentages(self, categories_all_used):
        text = categories_all_used.render()
        assert "include" in text
        assert "50%" in text

    def test_full_shows_counts(self, categories_all_used):
        text = categories_all_used.render("full")
        assert "100" in text

    def test_unused_shown(self, categories_with_unused):
        text = categories_with_unused.render()
        assert "unused" in text
        assert "maybe" in text
        assert "kept" in text

    def test_rare_shown(self):
        r = CategoryReport(
            counts={"include": 200, "exclude": 3},
            policy="error",
            rare_threshold=10,
        )
        text = r.render()
        assert "exclude 1% (rare)" in text
        assert "estimates involving it are imprecise" in text

    @pytest.mark.parametrize("verbosity", ["summary", "full"])
    def test_no_line_looks_like_a_category(self, verbosity):
        """'rare' and 'unused' mark a row; they never sit where a name would."""
        r = CategoryReport(
            counts={"include": 200, "exclude": 3, "maybe": 0},
            policy="keep", unused=("maybe",), rare_threshold=10,
        )
        names = [line.strip().split()[0] for line in r.render(verbosity).splitlines()[1:] if line.strip()]
        assert not {"rare", "unused", "rare:", "unused:"} & set(names)

    def test_full_marks_rows(self):
        r = CategoryReport(
            counts={"include": 200, "exclude": 3, "maybe": 0},
            policy="keep", unused=("maybe",), rare_threshold=10,
        )
        lines = r.render("full").splitlines()
        assert next(l for l in lines if l.strip().startswith("exclude")).endswith("(rare)")
        assert next(l for l in lines if l.strip().startswith("maybe")).endswith("(unused, kept)")
        include = next(l for l in lines if l.strip().startswith("include"))
        assert "(rare)" not in include and "(unused" not in include

    def test_all_zero_counts_no_division_error(self):
        """All zero counts should not raise (total or 1 guard)."""
        r = CategoryReport(counts={"a": 0, "b": 0}, policy="keep")
        text = r.render()
        assert "Categories" in text


# ============================================================================ #
#  CleaningReport                                                              #
# ============================================================================ #


class TestCleaningReportProperties:
    """Derived properties combine stored shapes and sub-reports."""

    def test_n_subjects_lost_none(self, report_stable):
        assert report_stable.n_subjects_lost == 0

    def test_n_subjects_lost_some(self, report_changed):
        # shape_in (6, 3), shape_out (4, 2)
        assert report_changed.n_subjects_lost == 2

    def test_n_raters_lost(self, report_changed):
        assert report_changed.n_raters_lost == 1

    def test_survived_full(self, report_stable):
        assert report_stable.survived == 1.0

    def test_survived_partial(self, report_changed):
        # 4 / 6
        assert abs(report_changed.survived - 4 / 6) < 1e-12

    def test_survived_zero_rows(self):
        """Zero input rows does not raise."""
        from interrater.semantics.categories import UnorderedCategories
        from interrater.semantics.ratings_semantics import RatingsSemantics
        from interrater.cleaners.cleaner_config import CleanerConfig

        sem = RatingsSemantics(
            raters=["ann", "bob"], categories=UnorderedCategories(("a", "b")),
        )
        r = CleaningReport(
            config=CleanerConfig(),
            declared_semantics=sem,
            effective_semantics=sem,
            shape_in=(0, 2),
            shape_out=(0, 2),
            matching=MatchingReport(n_cells=0, n_null=0),
            drops=DropReport(method="raters_first", min_ratings_per_rater=2, min_ratings_per_subject=2),
            categories=CategoryReport(counts={"a": 0, "b": 0}, policy="keep"),
        )
        assert r.survived == 0.0

    def test_semantics_changed_false(self, report_stable):
        assert report_stable.semantics_changed is False

    def test_semantics_changed_true(self, report_changed):
        assert report_changed.semantics_changed is True

    def test_k_changed_false(self, report_stable):
        assert report_stable.k_changed is False

    def test_k_changed_true(self):
        """k_changed is True when a category was dropped."""
        from interrater.semantics.categories import UnorderedCategories
        from interrater.semantics.ratings_semantics import RatingsSemantics
        from interrater.cleaners.cleaner_config import CleanerConfig

        declared = RatingsSemantics(
            raters=["ann", "bob"],
            categories=UnorderedCategories(("a", "b", "c")),
        )
        effective = RatingsSemantics(
            raters=["ann", "bob"],
            categories=UnorderedCategories(("a", "b")),
        )
        r = CleaningReport(
            config=CleanerConfig(),
            declared_semantics=declared,
            effective_semantics=effective,
            shape_in=(4, 2),
            shape_out=(4, 2),
            matching=MatchingReport(n_cells=8, n_null=0),
            drops=DropReport(method="raters_first", min_ratings_per_rater=2, min_ratings_per_subject=2),
            categories=CategoryReport(
                counts={"a": 4, "b": 4, "c": 0},
                policy="drop",
                action="dropped",
                unused=("c",),
            ),
        )
        assert r.k_changed is True


class TestCleaningReportSummary:
    """summary() returns plain types, delegating to sub-reports."""

    def test_keys(self, report_stable):
        s = report_stable.summary()
        expected = {
            "shape_in", "shape_out", "n_subjects_lost", "n_raters_lost",
            "survived", "semantics_changed",
            "declared_semantics", "effective_semantics",
            "matching", "drops", "categories",
            "_package_version", "cleaner_config",
        }
        assert set(s.keys()) == expected

    def test_record_carries_its_config_and_version(self, report_stable):
        """The saved record says how the data was cleaned, and by which release."""
        from interrater import __version__
        from interrater.cleaners.cleaner_config import CleanerConfig

        s = report_stable.summary()
        assert s["_package_version"] == __version__
        assert CleanerConfig.from_dict(s["cleaner_config"]) == report_stable.config

    def test_record_survives_yaml(self, report_changed, tmp_path):
        import yaml

        path = tmp_path / "report.yaml"
        report_changed.to_yaml(path)
        assert yaml.safe_load(path.read_text())["cleaner_config"] == report_changed.config.to_dict()

    def test_shape_is_list(self, report_stable):
        """Tuples become lists for JSON compatibility."""
        s = report_stable.summary()
        assert isinstance(s["shape_in"], list)
        assert isinstance(s["shape_out"], list)

    def test_semantics_are_dicts(self, report_stable):
        s = report_stable.summary()
        assert isinstance(s["declared_semantics"], dict)
        assert isinstance(s["effective_semantics"], dict)

    def test_sub_reports_are_dicts(self, report_stable):
        s = report_stable.summary()
        assert isinstance(s["matching"], dict)
        assert isinstance(s["drops"], dict)
        assert isinstance(s["categories"], dict)


class TestCleaningReportRender:
    """render() produces readable text; verbosity is validated."""

    def test_invalid_verbosity_raises(self, report_stable):
        with pytest.raises(ValueError, match="verbosity"):
            report_stable.render("garbage")

    def test_summary_header(self, report_stable):
        text = report_stable.render("summary")
        assert "CleaningReport" in text
        assert "100.0% of rows kept" in text

    def test_full_contains_all_sections(self, report_stable):
        text = report_stable.render("full")
        assert "Matching" in text
        assert "Drops" in text
        assert "Categories" in text

    def test_semantics_section_absent_when_unchanged(self, report_stable):
        text = report_stable.render("full")
        assert "Semantics" not in text

    def test_semantics_section_present_when_changed(self, report_changed):
        text = report_changed.render()
        assert "Semantics" in text
        assert "declared" in text
        assert "effective" in text

    def test_dropped_rater_named_in_semantics(self, report_changed):
        text = report_changed.render()
        assert "carol" in text

    def test_k_changed_warning(self):
        """When K changes, the render warns about coefficient impact."""
        from interrater.semantics.categories import UnorderedCategories
        from interrater.semantics.ratings_semantics import RatingsSemantics
        from interrater.cleaners.cleaner_config import CleanerConfig

        declared = RatingsSemantics(
            raters=["ann", "bob"],
            categories=UnorderedCategories(("a", "b", "c")),
        )
        effective = RatingsSemantics(
            raters=["ann", "bob"],
            categories=UnorderedCategories(("a", "b")),
        )
        r = CleaningReport(
            config=CleanerConfig(),
            declared_semantics=declared,
            effective_semantics=effective,
            shape_in=(4, 2),
            shape_out=(4, 2),
            matching=MatchingReport(n_cells=8, n_null=0),
            drops=DropReport(method="raters_first", min_ratings_per_rater=2, min_ratings_per_subject=2),
            categories=CategoryReport(
                counts={"a": 4, "b": 4, "c": 0},
                policy="drop",
                action="dropped",
                unused=("c",),
            ),
        )
        text = r.render()
        assert "CATEGORIES" in text
        assert "Kappa is unaffected" in text


class TestCleaningReportToYaml:
    """to_yaml writes a file that yaml.safe_load can read back."""

    def test_round_trip(self, report_stable, tmp_path):
        import yaml

        path = tmp_path / "report.yaml"
        report_stable.to_yaml(path)
        loaded = yaml.safe_load(path.read_text())
        assert isinstance(loaded, dict)
        assert loaded["shape_in"] == [4, 2]
        assert loaded["survived"] == 1.0
        assert "matching" in loaded
        assert "drops" in loaded
        assert "categories" in loaded


class TestCleaningReportFrozen:
    """CleaningReport is immutable."""

    def test_frozen(self, report_stable):
        with pytest.raises(AttributeError):
            report_stable.shape_in = (10, 10)


# ============================================================================ #
#  Repr                                                                        #
# ============================================================================ #


class TestRepr:
    """__repr__ delegates to render() for all four report classes."""

    def test_matching_report(self, matching_messy):
        r = repr(matching_messy)
        assert r.startswith("Matching")

    def test_drop_report_nothing(self, drops_nothing):
        r = repr(drops_nothing)
        assert r.startswith("Drops")

    def test_drop_report_some(self, drops_some):
        r = repr(drops_some)
        assert r.startswith("Drops")

    def test_category_report(self, categories_all_used):
        r = repr(categories_all_used)
        assert r.startswith("Categories")

    def test_cleaning_report(self, report_stable):
        r = repr(report_stable)
        assert r.startswith("CleaningReport")


# ============================================================================ #
#  Wording                                                                     #
# ============================================================================ #


class TestReportWording:
    """Plain plurals, pass lines that say what they removed, labels and totals."""

    @pytest.mark.parametrize("verbosity", ["summary", "full"])
    def test_no_bracketed_plurals(self, report_changed, verbosity):
        assert "(s)" not in report_changed.render(verbosity)
        assert "(es)" not in report_changed.render(verbosity)

    def test_pass_lines_say_what_was_removed(self):
        r = DropReport(
            method="iterate", min_ratings_per_rater=2, min_ratings_per_subject=2,
            passes=(DropPass(1, 1, 0), DropPass(2, 0, 0)),
            dropped_raters={"carol": 1},
        )
        text = r.render("full")
        assert "pass 1: removed 1 rater and 0 subjects" in text
        assert "pass 2: nothing removed" in text
        assert "-0" not in text

    def test_held_counts_read_as_sentences(self, drops_some):
        assert "2 subjects holding 1 rating" in drops_some.render()

    def test_categories_header_counts_cleaned_ratings(self, categories_all_used):
        assert categories_all_used.render().splitlines()[0] == (
            "Categories (cleaned data: 200 ratings)"
        )

    def test_label_shown_when_it_adds_something(self):
        r = CategoryReport(
            counts={"1": 62, "0": 1, "-1": 259}, policy="error",
            labels={"1": "include", "0": "uncertain", "-1": "Exclude"},
        )
        text = r.render("full")
        assert "1 (include)" in text and "-1 (Exclude)" in text

    def test_label_differing_only_in_case_not_repeated(self):
        r = CategoryReport(
            counts={"include": 5, "exclude": 5}, policy="error",
            labels={"include": "Include", "exclude": "Exclude"},
        )
        assert "(Include)" not in r.render("full")
