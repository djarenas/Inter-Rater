"""The only way in: a messy frame becomes a validated Ratings, with a receipt."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import numpy as np

from ..contracts import EMPTY, INDEX_DTYPE, LABEL_DTYPE
from ..data_objects.ratings import Ratings
from ..semantics.categories import Categories
from ..semantics.ratings_semantics import RatingsSemantics
from . import _drop_utils as _d
from . import _string_utils as _s
from ..reports.category_report import CategoryReport
from ..reports.cleaning_report import CleaningReport
from ..reports.drop_report import DropPass, DropReport
from ..reports.matching_report import MatchingReport
from .cleaner_config import CleanerConfig
from .insufficient_data_policy import MAX_DROP_PASSES

_ERROR_PREFIX = "[RatingsCleaner] Error"


class RatingsCleaner:
    """
    I turn a messy pandas dataframe into a Ratings, and hand you a receipt for what I did. The 
    receipt is the CleaningReport.

    I am the preferred entrance into Ratings. The ``Ratings`` constructor is a 
    strict version that would not produce a report nor allow you to supply 
    detailed choices.

    I check the config against the semantics before touching a single row. A
    repair pointing at a category that does not exist, or a missing marker that
    is also a valid response, is a contradiction between two documents and you
    should hear about it before you hear about row 400.

    The data is first cleaned in string space. Normalising, matching and
    repairing all happen on text. Once every cell is either a valid category
    name or nothing, everything is encoded into int8 so that the Ratings object
    holds data that can be easily and efficiently used by estimators.

    All repair is accounted for. Everything I change is counted, attributed to the
    rule that caused it, and broken down per rater -- because four hundred
    absent cells spread evenly is by design and four hundred in one column
    is somebody who stopped halfway, and only the breakdown tells them apart.

    Parameters
    ----------
    frame:
        The data, as strings. Reading the file is your job: by the time a frame
        exists, whether it had a header has already been decided.
    semantics:
        Who the raters are and what they could say.
    config:
        Every decision I am allowed to make. Defaults are strict.
    """

    # Attributes -------------------------------------------------------------  #

    _frame:               pd.DataFrame
    _declared_semantics:  RatingsSemantics
    _effective_semantics: RatingsSemantics
    _config:              CleanerConfig

    ratings: Ratings
    report:  CleaningReport


    # Constructor ------------------------------------------------------------  #

    def __init__(
        self,
        frame:     pd.DataFrame,
        semantics: RatingsSemantics,
        config:    CleanerConfig | None = None,
    ) -> None:
        self._check_constructor_inputs(frame, semantics, config)

        self._frame              = frame
        self._declared_semantics = semantics
        self._config             = config or CleanerConfig()

        self.ratings, self.report = self._clean()

    # Helper methods ---------------------------------------------------------  #

    @staticmethod
    def _check_constructor_inputs(
        frame:     pd.DataFrame,
        semantics: RatingsSemantics,
        config:    CleanerConfig | None,
    ) -> None:
        """
        Check that what was handed in is the right kind of thing.

        Runs before anything is stored, so nothing invalid is ever held. The
        checks here are about types and about the one structural property of a
        frame I cannot work around; everything about whether the *contents*
        make sense happens later, once there is something to look at.

        Two cases are worth naming rather than leaving to a generic message.

        Reading the file is the caller's job, so a path arrives here more often
        than any other wrong type. It gets told what to do about it.

        A pandas frame is allowed to have two columns with the same name, and
        then ``frame[name]`` hands back a frame rather than a column. Nothing
        downstream survives that, and the failure it produces -- an ambiguous
        truth value, raised from inside a numpy comparison -- says nothing
        about the cause. With ``use_column_names`` on, the mismatch against the
        declared raters would catch it; with it off, the labels are never
        looked at and this is the only place it can be found.
        """
        if isinstance(frame, (str, Path)):
            raise TypeError(
                f"{_ERROR_PREFIX}: frame must be a pandas DataFrame, got the "
                f"path {str(frame)!r}. Reading the file should occur before this, so that "
                f"how it was read -- which separator, whether it had a header, "
                f"what counted as a null -- stays a decision you can see. Try "
                f"pandas.read_csv first."
            )

        if not isinstance(frame, pd.DataFrame):
            raise TypeError(
                f"{_ERROR_PREFIX}: frame must be a pandas DataFrame, got "
                f"{type(frame).__name__}."
            )

        duplicated = frame.columns[frame.columns.duplicated()].tolist()
        if duplicated:
            raise ValueError(
                f"{_ERROR_PREFIX}: the frame has {len(duplicated)} repeated "
                f"column name(s) (e.g. {duplicated[:3]}). Selecting one of "
                f"them returns two columns rather than one, which nothing "
                f"downstream can read. Rename them before cleaning."
            )

        if not isinstance(semantics, RatingsSemantics):
            raise TypeError(
                f"{_ERROR_PREFIX}: semantics must be a RatingsSemantics, got "
                f"{type(semantics).__name__}. Build one directly or with "
                f"RatingsSemantics.build_from_yaml."
            )

        if config is not None and not isinstance(config, CleanerConfig):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be a CleanerConfig or None, "
                f"got {type(config).__name__}. None takes the strict "
                f"defaults; build one directly or with "
                f"CleanerConfig.build_from_yaml."
            )

    def _clean(self) -> tuple[Ratings, CleaningReport]:
        """
        Run the whole pipeline. Returns the Ratings and the account of it.

        The report is returned beside the Ratings, never attached to it: a
        Ratings describes its contents, not its history, and two datasets
        holding the same values are the same data however they were made.
        """
        self._check_config_against_semantics()

        subject_ids, columns = self._interpret_frame()
        codes, matching      = self._encode(columns)
        codes, subject_ids, raters, drops = self._drop(codes, subject_ids)
        codes, cats, categories           = self._resolve_unused_categories(codes)

        self._effective_semantics = self._effective_semantics_for(raters, cats)

        # Ratings names itself, but not that the array it rejected was one I
        # produced. That is the part worth knowing: the failure is in what the
        # cleaning left behind, not in the file.
        try:
            ratings = Ratings(codes, self._effective_semantics, subject_ids)
        except (TypeError, ValueError) as error:
            raise type(error)(f"{_ERROR_PREFIX}: {error}") from error

        report = CleaningReport(
            config              = self._config,
            declared_semantics  = self._declared_semantics,
            effective_semantics = self._effective_semantics,
            shape_in            = (len(self._frame), len(self._frame.columns)),
            shape_out           = (ratings.n_subjects, ratings.n_raters),
            matching            = matching,
            drops               = drops,
            categories          = categories,
        )
        return ratings, report

    def _check_config_against_semantics(self) -> None:
        """
        Checks that the choices in the cleaner configuration are 
        consistent with the semantics of the Ratings. For example, 
        value_map may be inconsistent in trying to map strings to a
        category that does not exist in semantics.
        """
        self._check_names_distinct_after_folding()

        names  = self._declared_semantics.categories.names
        folded = {self._fold(n): n for n in names}

        unknown = sorted(
            f"{written!r} -> {meant!r}"
            for written, meant in self._config.value_map.items()
            if meant not in names
        )
        if unknown:
            raise ValueError(
                f"{_ERROR_PREFIX}: the config repairs {len(unknown)} string(s) into "
                f"categories this rating task does not have (e.g. {unknown[:3]}). "
                f"Declared categories are {list(names)}."
            )

        collides = sorted(
            marker for marker in self._config.missing_values
            if self._fold(marker) in folded and marker != ""
        )
        if collides:
            raise ValueError(
                f"{_ERROR_PREFIX}: the config treats {collides[:3]} as meaning no rating, "
                f"but {[folded[self._fold(c)] for c in collides[:3]]} is a "
                f"declared category. A response a rater was offered cannot "
                f"also be a marker for absence."
            )

    def _check_names_distinct_after_folding(self) -> None:
        """
        Refuse category names the matching rule cannot tell apart.

        Matched as one, every answer for either would be coded as the same
        category, the other would look unused, and K would still count both.
        The semantics alone cannot decide this: 'A' and 'a' are distinct
        under case_sensitive=True. So it is checked here, the first point
        where both the names and the matching rule are known, before any
        cell is read.
        """
        collisions = _s.fold_collisions(
            self._declared_semantics.categories.names,
            case_sensitive   = self._config.case_sensitive,
            strip_whitespace = self._config.strip_whitespace,
        )
        if collisions:
            shown = ", ".join(" and ".join(repr(n) for n in g) for g in collisions[:3])
            rule = _s.describe_folding(
                case_sensitive   = self._config.case_sensitive,
                strip_whitespace = self._config.strip_whitespace,
            )
            raise ValueError(
                f"{_ERROR_PREFIX}: the categories {shown} cannot be told apart. "
                f"{rule} The cleaner would code them as one category. Rename "
                f"one, or set case_sensitive=True / strip_whitespace=False if "
                f"the difference is meant."
            )

    def _interpret_frame(self) -> tuple[np.ndarray, list]:
        """
        Work out which column is which, and lift out the subject ids.

        The subject id column, if one was named, is taken out first so the
        rest are rater columns. Whether their labels have to match the
        declared raters is the layout's call; when they do not, names come
        from the semantics by position.

        ``n_subjects`` is checked here, against the frame as handed in. It was
        an expectation about the input, so checking it after dropping would
        turn every legitimate drop into a failure.
        """
        layout  = self._config.frame_layout
        frame   = self._frame
        columns = list(frame.columns)

        subject_ids = None
        if layout.subject_id_col is not None:
            if layout.subject_id_col not in columns:
                raise KeyError(
                    f"{_ERROR_PREFIX}: subject_id_col {layout.subject_id_col!r} is not a "
                    f"column of this frame. Columns are {columns}."
                )
            subject_ids = np.asarray(
                frame[layout.subject_id_col].astype(str), dtype=LABEL_DTYPE
            )
            columns = [c for c in columns if c != layout.subject_id_col]

        declared = self._declared_semantics.n_subjects
        if declared is not None and len(frame) != declared:
            raise ValueError(
                f"{_ERROR_PREFIX}: the frame has {len(frame):,} row(s) but the semantics "
                f"expect {declared:,}. That expectation is about the data "
                f"handed in, not what survives cleaning -- either this is not "
                f"the file the rating task describes, or set n_subjects to None."
            )

        expected = list(self._declared_semantics.raters)
        if layout.use_column_names:
            if columns != expected:
                raise ValueError(
                    f"{_ERROR_PREFIX}: the frame's rater columns are {columns} but the "
                    f"semantics declare {expected}. Columns must match by name "
                    f"and in order. A column nobody declared is a mistake to "
                    f"look into, not a rater to invent -- set "
                    f"use_column_names=False to take names by position instead."
                )
        elif len(columns) != len(expected):
            raise ValueError(
                f"{_ERROR_PREFIX}: the frame has {len(columns)} rater column(s) but the "
                f"semantics declare {len(expected)}: {expected}."
            )

        # Ratings are matched as text, and a decimal column has lost how its
        # values were written: '1' became 1.0 and '1.50' became 1.5. pandas
        # does this to a numeric column with a blank cell. Refused outright
        # rather than matched, since a match on such a column is luck.
        decimal = [c for c in columns if frame[c].dtype.kind == "f"]
        if decimal:
            raise TypeError(
                f"{_ERROR_PREFIX}: rater column(s) {decimal} were read as decimal "
                f"numbers ({frame[decimal[0]].dtype}), usually because a numeric "
                f"column has blank cells. Decimals lose how a rating was written "
                f"('1' becomes 1.0, '1.50' becomes 1.5), so they cannot be "
                f"matched to categories. Read the file as text, e.g. "
                f"pd.read_csv(path, dtype=str)."
            )

        return subject_ids, columns

    def _encode(self, columns: list) -> tuple[np.ndarray, MatchingReport]:
        """
        Cross the boundary: every cell becomes a category position or EMPTY.

        One lookup does all the matching -- categories, missing markers and
        repairs alike -- so a string cannot match one rule under normalisation
        and a different one later. Counting happens per distinct token rather
        than per cell, which is why the report can attribute a cell to the
        rule that caught it instead of only totalling outcomes.

        Anything that matches nothing is an unknown, and what happens then is
        ``unknown_policy``: stop, or record it and treat it as no rating.
        """
        cfg    = self._config
        names  = self._declared_semantics.categories.names
        raters = list(self._declared_semantics.raters)

        lookup, origin = _s.build_lookup(
            names, cfg.missing_values, cfg.value_map, EMPTY,
            case_sensitive   = cfg.case_sensitive,
            strip_whitespace = cfg.strip_whitespace,
        )
        # Folded tokens are what we match on; reports name what was declared.
        from_missing = {self._fold(m): m for m in cfg.missing_values}
        from_repair  = {self._fold(k): k for k in cfg.value_map}

        n_rows = len(self._frame)
        codes  = np.empty((n_rows, len(columns)), dtype=INDEX_DTYPE)

        n_null          = 0
        missing_matched = {m: 0 for m in cfg.missing_values}
        repairs         = {k: 0 for k in cfg.value_map}
        unknown_values: dict[str, int] = {}
        per_rater_null    = {}
        per_rater_missing = {}
        per_rater_unknown = {}

        for position, column in enumerate(columns):
            rater  = raters[position]
            values = self._frame[column]
            null   = values.isna().to_numpy()

            folded = np.array(
                [
                    "" if is_null else self._fold(str(raw))
                    for raw, is_null in zip(values.to_numpy(), null)
                ],
                dtype=object,
            )

            column_codes = np.full(n_rows, EMPTY, dtype=INDEX_DTYPE)
            unknown_here = 0

            for token, count in zip(*np.unique(folded[~null], return_counts=True)):
                where = (folded == token) & ~null
                if token in lookup:
                    column_codes[where] = lookup[token]
                    kind = origin[token]
                    if kind == "missing":
                        missing_matched[from_missing[token]] += int(count)
                    elif kind == "repair":
                        repairs[from_repair[token]] += int(count)
                else:
                    unknown_here += int(count)
                    unknown_values[token] = unknown_values.get(token, 0) + int(count)

            if unknown_here and cfg.unknown_policy == "error":
                examples = sorted(unknown_values, key=lambda t: -unknown_values[t])[:3]
                raise ValueError(
                    f"{_ERROR_PREFIX}: rater {rater!r} gave {unknown_here:,} response(s) "
                    f"that are neither a declared category nor a marker for "
                    f"absence (e.g. {examples}). Declared categories are "
                    f"{list(names)}. Add a value_map entry to repair them, "
                    f"list them in missing_values, or set "
                    f"unknown_policy='empty' to treat them as no rating."
                )

            codes[:, position]        = column_codes
            per_rater_null[rater]     = int(null.sum())
            per_rater_unknown[rater]  = unknown_here
            per_rater_missing[rater]  = int((column_codes == EMPTY).sum())
            n_null                   += int(null.sum())

        report = MatchingReport(
            n_cells         = n_rows * len(columns),
            n_null          = n_null,
            missing_matched = missing_matched,
            repairs         = repairs,
            unknown_values  = unknown_values,
            per_rater_null    = per_rater_null,
            per_rater_missing = per_rater_missing,
            per_rater_unknown = per_rater_unknown,
            subject_missing_distribution = _d.distribution(
                (codes == EMPTY).sum(axis=1)
            ),
        )
        return codes, report

    def _drop(
        self, codes: np.ndarray, subject_ids: np.ndarray | None
    ) -> tuple[np.ndarray, np.ndarray, list[str], DropReport]:
        """
        Remove raters and subjects holding too little data to be worth keeping.

        Which order, and whether to repeat, is the policy's call -- dropping
        is not symmetric, so a sweep one way gives a different answer from a
        sweep the other. Counts are taken before the drop so the report can
        say how much each casualty held when it went: the count it was judged
        on, in the pass that removed it, not its count before any dropping.

        Subject ids are generated here if the frame had none, because a row
        must be nameable before it can be reported as dropped.
        """
        policy = self._config.insufficient_data_drops
        raters = list(self._declared_semantics.raters)

        if subject_ids is None:
            subject_ids = np.array(
                [f"subject{i}" for i in range(codes.shape[0])], dtype=LABEL_DTYPE
            )

        if policy.method == "iterate":
            # raters_first is arbitrary here: the maximal set satisfying both
            # thresholds is unique, so sweep order affects the path but not the
            # converged endpoint.
            (keep_s, keep_r, per_pass, converged,
             subject_held, rater_held) = _d.sweep_until_settled(
                codes, policy.min_ratings_per_rater,
                policy.min_ratings_per_subject, True, MAX_DROP_PASSES,
            )
        else:
            keep_s, keep_r, subject_held, rater_held = _d.sweep_once(
                codes, policy.min_ratings_per_rater,
                policy.min_ratings_per_subject,
                raters_first = policy.method == "raters_first",
            )
            per_pass  = [(int((~keep_r).sum()), int((~keep_s).sum()))]
            converged = True

        dropped_raters = {
            raters[i]: int(rater_held[i])
            for i in np.flatnonzero(~keep_r)
        }
        dropped_ids = tuple(subject_ids[~keep_s].tolist())

        report = DropReport(
            method    = policy.method,
            min_ratings_per_rater   = policy.min_ratings_per_rater,
            min_ratings_per_subject = policy.min_ratings_per_subject,
            converged = converged,
            passes    = tuple(
                DropPass(i + 1, r, s) for i, (r, s) in enumerate(per_pass)
            ),
            dropped_raters      = dropped_raters,
            dropped_subject_ids = dropped_ids,
            dropped_subject_distribution = _d.distribution(
                subject_held[~keep_s]
            ),
        )
        return (
            codes[np.ix_(keep_s, keep_r)],
            subject_ids[keep_s],
            [r for r, keep in zip(raters, keep_r) if keep],
            report,
        )

    def _resolve_unused_categories(
        self, codes: np.ndarray
    ) -> tuple[np.ndarray, Categories, CategoryReport]:
        """
        Decide what happens to categories nobody chose.

        Under ``error`` I stop. Under ``keep`` the scheme is untouched: K is a
        fact about what raters could have said, not about what they did say.
        I never change the declared categories.

        I return categories, never semantics. Building the declaration the
        surviving data describes happens in one place, so there is a single
        answer to what the effective rating task is.
        """
        cats   = self._declared_semantics.categories
        names  = cats.names
        counts = {
            name: int((codes == index).sum())
            for index, name in enumerate(names)
        }
        unused = tuple(name for name, n in counts.items() if n == 0)
        policy = self._config.unused_categories

        if unused and policy == "error":
            raise ValueError(
                f"{_ERROR_PREFIX}: {len(unused)} declared categor(ies) were never chosen "
                f"by anybody: {list(unused)}. That is a fact about the cohort "
                f"as often as it is a mistake, but it is also the shape a "
                f"mistake takes. Set unused_categories='keep' to accept it, or "
                f"remove them from your semantics file if they should not be "
                f"in the scheme."
            )

        report = CategoryReport(
            counts = counts,
            policy = policy,
            action = "kept",
            unused = unused,
            labels = dict(self._declared_semantics.categories.labels),
        )
        return codes, cats, report

    def _effective_semantics_for(
        self, raters: list[str], categories: Categories
    ) -> RatingsSemantics:
        """
        The declaration the surviving data actually describes.

        A rater who fell below a threshold is gone, and a Ratings must carry a
        declaration matching what it holds -- so when anything was dropped this
        is not the object handed in. That is not me overruling the user: it is
        the consequence of a threshold they set. But it does mean
        ``ratings.semantics`` differs from theirs, and the report says so in
        the one place they will read.

        ``n_subjects`` is deliberately dropped. It was an expectation about
        the frame handed in, already checked before any row was touched, and
        carrying it forward would make it an assertion about what survived --
        which it never was.

        Built once, here, so there is a single answer to what the effective
        rating task is.
        """
        if (
            list(raters) == list(self._declared_semantics.raters)
            and categories == self._declared_semantics.categories
        ):
            return self._declared_semantics

        try:
            return RatingsSemantics(
                raters          = raters,
                categories      = categories,
                n_subjects      = None,
                unique_subjects = self._declared_semantics.unique_subjects,
            )
        except (TypeError, ValueError) as error:
            raise type(error)(
                f"{_ERROR_PREFIX}: what survived cleaning is not a rating "
                f"task that can exist: {error}"
            ) from error

    def _fold(self, value: str) -> str:
        return _s.fold(
            value,
            case_sensitive   = self._config.case_sensitive,
            strip_whitespace = self._config.strip_whitespace,
        )