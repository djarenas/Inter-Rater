"""Prevalence and bias indices for two raters and two categories."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...contracts import EMPTY
from ..not_computed import NotComputed

_ERROR_PREFIX = "[PrevalenceBiasIndex] Error"

_COUNTS_ONLY = (
    "needs each rater's ratings; this was built from category counts only, "
    "which do not say which rater gave which answer"
)


@dataclass(frozen=True)
class PrevalenceBiasIndex:
    """
    I hold the prevalence and bias indices of Byrt, Bishop and Carlin (1993)
    for two raters and two categories, with their orientation.

    Kappa can be low while raters agree on almost every subject, when one
    category dominates (prevalence) or when one rater says "positive" more
    often than the other (bias). These two numbers show which.

    The signs depend on which category is "positive" and which rater comes
    first, so I carry both and say them in words:

    - "positive" is the first declared category;
    - prevalence_index = (a - d) / N, where a is the subjects both raters
      called positive and d the subjects both called negative. Above zero:
      more agreement on the positive category than on the negative one.
    - bias_index = P(rater_b says positive) - P(rater_a says positive)
      = (b - c) / N, as in R's epiR. Above zero: rater_b says positive more
      often. It needs to know which rater gave which answer, so from
      category counts alone it is NotComputed.

    I store the four cells and work the indices out from them, so they can
    never disagree with the counts.

    Attributes
    ----------
    positive_category, negative_category : str
        The first and second declared categories.
    rater_a, rater_b : str | NotComputed
        The two raters, in order; NotComputed from counts only.
    n_both_positive : int
        a: subjects both raters called positive.
    n_both_negative : int
        d: subjects both raters called negative.
    n_disagree : int
        b + c: subjects the raters disagreed on.
    n_only_b_positive : int | NotComputed
        b: rater_b positive, rater_a negative.
    n_only_a_positive : int | NotComputed
        c: rater_a positive, rater_b negative.
    """

    positive_category: str
    negative_category: str
    rater_a: str | NotComputed
    rater_b: str | NotComputed
    n_both_positive: int
    n_both_negative: int
    n_disagree: int
    n_only_b_positive: int | NotComputed
    n_only_a_positive: int | NotComputed

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        for name in ("positive_category", "negative_category"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must be a non-empty str, "
                    f"got {value!r}"
                )
        if self.positive_category == self.negative_category:
            raise ValueError(
                f"{_ERROR_PREFIX}: positive and negative category are both "
                f"{self.positive_category!r}"
            )

        for name in ("n_both_positive", "n_both_negative", "n_disagree"):
            _check_count(name, getattr(self, name))
        if self.n_subjects < 1:
            raise ValueError(f"{_ERROR_PREFIX}: no subjects")

        # Rater identity is either fully known or fully absent
        known = [
            not isinstance(getattr(self, n), NotComputed)
            for n in ("rater_a", "rater_b", "n_only_b_positive", "n_only_a_positive")
        ]
        if any(known) and not all(known):
            raise ValueError(
                f"{_ERROR_PREFIX}: rater_a, rater_b and the two one-sided "
                f"counts must all be known or all be NotComputed"
            )
        if all(known):
            for name in ("rater_a", "rater_b"):
                value = getattr(self, name)
                if not isinstance(value, str) or not value:
                    raise TypeError(
                        f"{_ERROR_PREFIX}: {name} must be a non-empty str, "
                        f"got {value!r}"
                    )
            if self.rater_a == self.rater_b:
                raise ValueError(
                    f"{_ERROR_PREFIX}: rater_a and rater_b are both "
                    f"{self.rater_a!r}"
                )
            _check_count("n_only_b_positive", self.n_only_b_positive)
            _check_count("n_only_a_positive", self.n_only_a_positive)
            if self.n_only_b_positive + self.n_only_a_positive != self.n_disagree:
                raise ValueError(
                    f"{_ERROR_PREFIX}: the one-sided counts "
                    f"({self.n_only_b_positive} + {self.n_only_a_positive}) do "
                    f"not add up to n_disagree ({self.n_disagree})"
                )

    @classmethod
    def from_pair(
        cls,
        category_index: np.ndarray,
        category_names: tuple[str, str],
        rater_names: tuple[str, str],
    ) -> PrevalenceBiasIndex:
        """
        Build from two raters' category positions on the subjects both rated.

        Parameters
        ----------
        category_index : np.ndarray
            Shape (n_subjects, 2): column 0 is rater_a, column 1 rater_b,
            each cell 0 (positive) or 1 (negative). Rows where either is
            EMPTY are left out.
        category_names, rater_names : tuple[str, str]
            In declared order.
        """
        index = np.asarray(category_index)
        if index.ndim != 2 or index.shape[1] != 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: needs exactly 2 rater columns, got shape "
                f"{index.shape}"
            )
        shared = index[(index != EMPTY).all(axis=1)]
        a_pos = shared[:, 0] == 0
        b_pos = shared[:, 1] == 0
        return cls(
            positive_category=category_names[0],
            negative_category=category_names[1],
            rater_a=rater_names[0],
            rater_b=rater_names[1],
            n_both_positive=int((a_pos & b_pos).sum()),
            n_both_negative=int((~a_pos & ~b_pos).sum()),
            n_disagree=int((a_pos != b_pos).sum()),
            n_only_b_positive=int((b_pos & ~a_pos).sum()),
            n_only_a_positive=int((a_pos & ~b_pos).sum()),
        )

    @classmethod
    def from_counts(
        cls,
        counts: np.ndarray,
        category_names: tuple[str, str],
    ) -> PrevalenceBiasIndex:
        """
        Build from per-subject category counts, two ratings per subject.

        Subjects with fewer than two ratings are left out, as everywhere in
        the package. Without rater identity the bias index is NotComputed.

        Parameters
        ----------
        counts : np.ndarray
            Shape (n_subjects, 2): how many ratings each subject got in the
            positive and the negative category.
        category_names : tuple[str, str]
            In declared order.

        Raises
        ------
        ValueError
            If any subject has more than two ratings: these are not two
            raters.
        """
        counts = np.asarray(counts)
        if counts.ndim != 2 or counts.shape[1] != 2:
            raise ValueError(
                f"{_ERROR_PREFIX}: needs exactly 2 categories, got shape "
                f"{counts.shape}"
            )
        per_subject = counts.sum(axis=1)
        if (per_subject > 2).any():
            raise ValueError(
                f"{_ERROR_PREFIX}: some subjects have more than two ratings; "
                f"the indices are defined for two raters"
            )
        pair = counts[per_subject == 2]
        reason = NotComputed(_COUNTS_ONLY)
        return cls(
            positive_category=category_names[0],
            negative_category=category_names[1],
            rater_a=NotComputed("rater names are not known from category counts"),
            rater_b=NotComputed("rater names are not known from category counts"),
            n_both_positive=int((pair[:, 0] == 2).sum()),
            n_both_negative=int((pair[:, 1] == 2).sum()),
            n_disagree=int((pair[:, 0] == 1).sum()),
            n_only_b_positive=reason,
            n_only_a_positive=reason,
        )

    # Public API ------------------------------------------------------------- #

    @property
    def n_subjects(self) -> int:
        """N: subjects both raters rated."""
        return self.n_both_positive + self.n_both_negative + self.n_disagree

    @property
    def prevalence_index(self) -> float:
        """(a - d) / N. Above zero: more agreement on the positive category."""
        return (self.n_both_positive - self.n_both_negative) / self.n_subjects

    @property
    def bias_index(self) -> float | NotComputed:
        """(b - c) / N. Above zero: rater_b says positive more often."""
        if isinstance(self.n_only_b_positive, NotComputed):
            return NotComputed(f"bias index {_COUNTS_ONLY}")
        return (self.n_only_b_positive - self.n_only_a_positive) / self.n_subjects

    def __str__(self) -> str:
        pos, neg = repr(self.positive_category), repr(self.negative_category)
        pi = self.prevalence_index
        if pi > 0:
            pi_words = f"more agreement on {pos} than on {neg}"
        elif pi < 0:
            pi_words = f"more agreement on {neg} than on {pos}"
        else:
            pi_words = f"as much agreement on {pos} as on {neg}"

        bi = self.bias_index
        if isinstance(bi, NotComputed):
            bi_line = f"BI {bi}"
            who = ""
        else:
            if bi > 0:
                bi_words = f"{self.rater_b} says {pos} more often than {self.rater_a}"
            elif bi < 0:
                bi_words = f"{self.rater_a} says {pos} more often than {self.rater_b}"
            else:
                bi_words = f"{self.rater_a} and {self.rater_b} say {pos} equally often"
            bi_line = f"BI = {bi:+.3f}: {bi_words}"
            who = f", {self.rater_a} vs {self.rater_b}"
        return (
            f"Prevalence and bias (positive = {pos}{who}, "
            f"{self.n_subjects:,} subjects)\n"
            f"    PI = {pi:+.3f}: {pi_words}\n"
            f"    {bi_line}"
        )

    def summary(self) -> dict:
        """
        Return my record in plain types.

        Returns
        -------
        dict
            The categories, the raters, the four cells, and both indices.
            Anything not computed appears as ``{"not_computed": reason}``.
        """
        def plain(value):
            return value.summary() if isinstance(value, NotComputed) else value

        return {
            "positive_category": self.positive_category,
            "negative_category": self.negative_category,
            "rater_a": plain(self.rater_a),
            "rater_b": plain(self.rater_b),
            "n_subjects": self.n_subjects,
            "n_both_positive": self.n_both_positive,
            "n_both_negative": self.n_both_negative,
            "n_disagree": self.n_disagree,
            "n_only_b_positive": plain(self.n_only_b_positive),
            "n_only_a_positive": plain(self.n_only_a_positive),
            "prevalence_index": self.prevalence_index,
            "bias_index": plain(self.bias_index),
        }


def _check_count(name: str, value: object) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(
            f"{_ERROR_PREFIX}: {name} must be int, got {type(value).__name__}"
        )
    if value < 0:
        raise ValueError(f"{_ERROR_PREFIX}: {name} must be >= 0, got {value}")
