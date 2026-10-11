"""How to read the frame, as opposed to what the rating task means."""

from __future__ import annotations


class FrameLayout:
    """
    I describe the file, not the rating task.

    Everything else a cleaner is configured with -- which strings mean nothing,
    what to do with an undeclared response, how thin the data may get -- is a
    decision about the research and travels unchanged to a second export. I am
    the part that does not: which column holds the subject ids, whether the
    header row can be trusted. A different spreadsheet of the same rating task needs
    a different me.

    I take a frame, never a path. Reading the file is the caller's job, and by
    the time a frame reaches a cleaner the question of whether the file had a
    header has already been answered -- pandas assigns integer labels when it
    did not, so nothing downstream could tell the difference anyway. What I
    decide is narrower: whether to trust the labels that are there.

    Attributes
    ----------
    use_column_names:
        True means the frame's column labels must equal the declared raters, in
        order, and anything else is an error. A column nobody declared is a
        mistake to look into, not a rater to invent. False ignores the labels
        entirely and takes rater names from the semantics by position.
    subject_id_col:
        Which column holds subject ids, or None to generate them positionally.

        Independent of ``use_column_names``. That flag governs whether the
        *rater* labels have to match the declared raters; this one names a
        single column to lift out. A frame always has column labels, so a
        subject id column is findable either way -- and a file whose rater
        columns are called R1..R4 while the raters are ann and bob is exactly
        the case where you want the labels ignored and the id column kept.

        The only requirement is that the column exists, which the cleaner
        checks when it has the frame in front of it.
    """

    # Attributes -------------------------------------------------------------  #

    use_column_names: bool
    subject_id_col:   str | None

    # Constructor ------------------------------------------------------------  #

    def __init__(
        self,
        use_column_names: bool = True,
        subject_id_col:   str | None = None,
    ) -> None:
        self.use_column_names = use_column_names
        self.subject_id_col   = subject_id_col
        self._validate()

    def _validate(self) -> None:
        name = type(self).__name__

        if not isinstance(self.use_column_names, bool):
            raise TypeError(
                f"{name}.use_column_names must be a boolean, got "
                f"{type(self.use_column_names).__name__}."
            )

        if self.subject_id_col is None:
            return

        if not isinstance(self.subject_id_col, str):
            raise TypeError(
                f"{name}.subject_id_col must be a column name or None, got "
                f"{type(self.subject_id_col).__name__}."
            )

    # Classmethods -----------------------------------------------------------  #

    @classmethod
    def from_dict(cls, spec: dict) -> FrameLayout:
        unknown = sorted(set(spec) - {"use_column_names", "subject_id_col"})
        if unknown:
            raise ValueError(
                f"frame_layout has unknown key(s) {unknown}. Known: "
                f"use_column_names, subject_id_col."
            )
        return cls(**spec)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FrameLayout):
            return NotImplemented
        return (
            self.use_column_names == other.use_column_names
            and self.subject_id_col == other.subject_id_col
        )

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return (
            f"FrameLayout(use_column_names={self.use_column_names}, "
            f"subject_id_col={self.subject_id_col!r})"
        )