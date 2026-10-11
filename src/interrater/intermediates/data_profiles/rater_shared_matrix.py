"""How many subjects each pair of raters both rated."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...contracts import MIN_OVERLAP_HARD

_ERROR_PREFIX = "[RaterSharedMatrix] Error"


@dataclass(frozen=True, eq=False)
class RaterSharedMatrix:
    """
    I count, for every pair of raters, how many subjects both of them rated.

    Every pairwise coefficient rests on one of my off-diagonal cells. Two
    raters who co-rated eight subjects produce a kappa, and it is worth
    almost nothing; two who co-rated four hundred produce one worth reading.
    Nothing in the coefficient itself says which you are looking at, which is
    why I am worth keeping beside it.

    My diagonal is each rater's own total -- the number of subjects that rater
    rated at all, since a rater trivially shares every one of them with
    themselves. That makes it redundant with ``n_per_rater``, but it is the
    natural value for the cell and it is useful to have to hand: no pair can
    share more subjects than either of its members rated, which is a check
    worth being able to make.

    Every question I answer is about pairs, so I never report a cell twice and
    never report the diagonal. Pairs come back as ``(a, b, value)`` with
    ``a < b``, so a list of them has one entry per pair rather than one per
    ordered pair.

    I hold indices, not names. A bootstrapper resampling rows has no use for
    rater names and should not have to carry them; whoever renders me has the
    semantics to hand and can supply them. That is also why my shape does not
    depend on where I came from: a two-rater dataset gives me a 2x2 matrix
    whose only pair is (0, 1), and every method behaves the same as it would
    on twenty raters.

    Attributes
    ----------
    shared:
        (n_raters, n_raters) integer, symmetric. Cell (a, b) is how many
        subjects both a and b rated.

    Examples
    --------
    >>> import numpy as np
    >>> m = RaterSharedMatrix(np.array([[10, 4, 1], [4, 8, 2], [1, 2, 9]]))
    >>> m.thinnest_pair()
    (0, 2, 1)
    >>> m.below(2)
    [(0, 2, 1)]
    """

    # Attributes -------------------------------------------------------------  #

    shared: np.ndarray

    # Constructor ------------------------------------------------------------  #

    def __post_init__(self) -> None:
        object.__setattr__(self, "shared", np.asarray(self.shared))
        self._validate()

    def _validate(self) -> None:
        shared = self.shared

        if shared.ndim != 2 or shared.shape[0] != shared.shape[1]:
            raise ValueError(
                f"{_ERROR_PREFIX}: shared must be a square matrix, one row "
                f"and column per rater, got shape {shared.shape}."
            )

        if shared.shape[0] < 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: shared must cover at least 2 raters, got "
                f"{shared.shape[0]}. With one rater there are no pairs, and "
                f"every question I answer is about a pair."
            )

        if shared.dtype.kind not in "iu":
            raise TypeError(
                f"{_ERROR_PREFIX}: shared must hold integers -- these are "
                f"counts of subjects -- got dtype {shared.dtype}."
            )

        if not np.array_equal(shared, shared.T):
            raise ValueError(
                f"{_ERROR_PREFIX}: shared must be symmetric. How many "
                f"subjects two raters both rated cannot depend on which of "
                f"them is named first."
            )

        if shared.min() < 0:
            raise ValueError(
                f"{_ERROR_PREFIX}: shared holds a negative count "
                f"({shared.min()}). These are numbers of subjects."
            )

        diagonal = np.diag(shared)
        exceeds  = shared > diagonal[:, None]
        if exceeds.any():
            row, col = np.argwhere(exceeds)[0]
            raise ValueError(
                f"{_ERROR_PREFIX}: raters {row} and {col} are recorded as "
                f"sharing {shared[row, col]:,} subject(s), but rater {row} "
                f"only rated {diagonal[row]:,}. A pair cannot share more "
                f"subjects than either of its members saw."
            )

    # Public API -------------------------------------------------------------  #

    @property
    def n_raters(self) -> int:
        return int(self.shared.shape[0])

    @property
    def n_pairs(self) -> int:
        """How many unordered pairs there are: n(n-1)/2."""
        n = self.n_raters
        return n * (n - 1) // 2

    @property
    def per_rater(self) -> np.ndarray:
        """My diagonal: how many subjects each rater rated."""
        return np.diag(self.shared).copy()

    def value_for(self, rater_a: int, rater_b: int) -> int:
        """How many subjects this pair both rated."""
        if rater_a == rater_b:
            raise ValueError(
                f"{_ERROR_PREFIX}: {rater_a} is not a pair with itself. Use "
                f"per_rater for how many subjects one rater rated."
            )
        return int(self.shared[rater_a, rater_b])

    def _upper_triangle(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        My off-diagonal cells once, as parallel arrays of a, b and value.

        Every pair question starts here. ``triu_indices`` with ``k=1`` skips
        the diagonal and takes each pair once, in row-major order -- which is
        already sorted by rater index, so a stable sort on the values keeps
        that as the tie-break without sorting on it explicitly.
        """
        rows, cols = np.triu_indices(self.n_raters, k=1)
        return rows, cols, self.shared[rows, cols]

    def sorted_pairs(self) -> list[tuple[int, int, int]]:
        """
        Every pair as ``(a, b, shared)``, thinnest first.

        One entry per pair, never per ordered pair, and never the diagonal.
        Ties break on rater index so the order is stable between runs.

        Note that this materialises one tuple per pair, and pairs grow as the
        square of raters: four raters make six, but a thousand crowd workers
        make half a million. If you only want the worst of them, ``below`` and
        ``thinnest_pair`` answer without building the list.
        """
        rows, cols, values = self._upper_triangle()
        order = np.argsort(values, kind="stable")
        return [
            (int(rows[i]), int(cols[i]), int(values[i])) for i in order
        ]

    def thinnest_pair(self) -> tuple[int, int, int]:
        """
        The pair sharing the fewest subjects, as ``(a, b, shared)``.

        The weakest link: whatever coefficient this pair produces is the least
        trustworthy one in the set.

        Found by a single pass rather than by sorting, because the ordering of
        everything else is not wanted here. That matters at crowdsourcing
        scale, where sorting half a million pairs to read the first is most of
        a second for no reason.
        """
        rows, cols, values = self._upper_triangle()
        i = int(np.argmin(values))
        return int(rows[i]), int(cols[i]), int(values[i])

    def below(self, threshold: int = MIN_OVERLAP_HARD) -> list[tuple[int, int, int]]:
        """
        Pairs sharing fewer than ``threshold`` subjects, thinnest first.

        At the default a coefficient cannot be computed for these at all:
        with fewer than two co-rated subjects there is nothing to estimate
        expected agreement from. Above that floor a pair may still be
        untrustworthy, and how far above depends on prevalence -- so a higher
        threshold here is a judgement, not a rule.

        Filtered before sorting, so the cost follows the number of thin pairs
        rather than the number of pairs.
        """
        rows, cols, values = self._upper_triangle()
        keep = values < threshold
        rows, cols, values = rows[keep], cols[keep], values[keep]
        order = np.argsort(values, kind="stable")
        return [
            (int(rows[i]), int(cols[i]), int(values[i])) for i in order
        ]

    def summary(self) -> dict:
        """Everything I hold, in plain types."""
        return {
            "n_raters":   self.n_raters,
            "n_pairs":    self.n_pairs,
            "per_rater":  self.per_rater.tolist(),
            "shared":     self.shared.tolist(),
            "pairs":      [list(pair) for pair in self.sorted_pairs()],
        }

    # Display ----------------------------------------------------------------  #

    _STYLES = ("grid", "pairs")

    _MAX_GRID_RATERS = 12
    """Past this a grid stops fitting a terminal and the list is clearer."""

    def render(
        self,
        rater_names: object = None,
        style:       str = "grid",
        max_pairs:   int = 10,
        header:      bool = True,
    ) -> str:
        """
        Me, as a person would want to read me.

        Two views, because they answer different questions. The grid shows a
        whole row at once, which is how you notice a rater who is thin against
        *everybody* rather than thin against one partner -- a difference the
        sorted list hides. The list shows the weak pairs first and stays
        readable at any number of raters, which the grid does not.

        Names are a parameter rather than something I hold: whoever has the
        semantics can pass them, and everything below that has no use for
        them. Without names, raters are reported by index.

        Past a dozen raters the grid stops fitting a terminal, so it falls
        back to the list rather than producing something unreadable.

        ``header`` names me at the top, which is right when I am printed on my
        own and redundant when a parent has already said what I am.
        """
        if style not in self._STYLES:
            raise ValueError(
                f"{_ERROR_PREFIX}: style must be one of {list(self._STYLES)}, "
                f"got {style!r}."
            )

        labels = (
            [str(name) for name in rater_names]
            if rater_names is not None
            else [str(i) for i in range(self.n_raters)]
        )
        if len(labels) != self.n_raters:
            raise ValueError(
                f"{_ERROR_PREFIX}: given {len(labels)} rater name(s) for "
                f"{self.n_raters} rater(s)."
            )

        if style == "grid" and self.n_raters <= self._MAX_GRID_RATERS:
            return self._render_grid(labels, header)

        note = (
            f"a grid of {self.n_raters} raters does not fit; showing the "
            f"thinnest pairs instead"
            if style == "grid" else None
        )
        return self._render_pairs(labels, max_pairs, note, header)

    # Helper methods ---------------------------------------------------------  #

    def _render_grid(self, labels: list[str], header: bool = True) -> str:
        """
        The upper triangle, with each rater's own total down the diagonal.

        Only the upper half is filled. The lower half repeats it, and leaving
        it blank makes the symmetry visible rather than something you check by
        eye. The diagonal is set apart in brackets because it is a rater's own
        total, not a pair -- printing it as a bare number invites reading it as
        one.
        """
        cells = [
            [
                f"[{self.shared[r, c]:,}]" if r == c
                else f"{self.shared[r, c]:,}" if c > r
                else ""
                for c in range(self.n_raters)
            ]
            for r in range(self.n_raters)
        ]
        name_width = max(len(name) for name in labels)
        widths = [
            max([len(labels[c])] + [len(row[c]) for row in cells])
            for c in range(self.n_raters)
        ]

        column_header = " " * (name_width + 6) + "  ".join(
            labels[c].rjust(widths[c]) for c in range(self.n_raters)
        )
        lines = ["RaterSharedMatrix"] if header else []
        lines += [
            f"    {self.n_raters} rater(s), {self.n_pairs} pair(s); "
            f"[n] is that rater's own total",
            "",
            column_header,
        ]
        for r in range(self.n_raters):
            body = "  ".join(cells[r][c].rjust(widths[c]) for c in range(self.n_raters))
            lines.append(f"    {labels[r].ljust(name_width)}  {body}")
        return "\n".join(lines)

    def _render_pairs(
        self,
        labels:    list[str],
        max_pairs: int,
        note:      str | None = None,
        header:    bool = True,
    ) -> str:
        """
        Pairs thinnest first, truncated. Readable at any number of raters.

        When a grid was asked for and refused, ``note`` says so. Falling back
        silently would leave somebody wondering why they got a different shape
        than they asked for.
        """
        pairs = self.sorted_pairs()
        shown = pairs[:max_pairs]
        pad   = "    "
        lines = ["RaterSharedMatrix"] if header else []
        lines.append(f"{pad}{self.n_raters} rater(s), {self.n_pairs} pair(s)")
        if note:
            lines.append(f"{pad}{note}")
        for a, b, value in shown:
            lines.append(f"{pad}{labels[a]} + {labels[b]}: {value:,}")
        if len(pairs) > len(shown):
            lines.append(f"{pad}... {len(pairs) - len(shown)} more pair(s)")
        return "\n".join(lines)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, RaterSharedMatrix):
            return NotImplemented
        return np.array_equal(self.shared, other.shared)

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        thin = self.thinnest_pair()
        return (
            f"RaterSharedMatrix({self.n_raters} raters, {self.n_pairs} pairs, "
            f"thinnest {thin[0]}+{thin[1]}={thin[2]:,})"
        )
