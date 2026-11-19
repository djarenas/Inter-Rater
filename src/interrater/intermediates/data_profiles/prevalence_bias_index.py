"""Prevalence and bias indices for two raters and two categories."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...contracts import EMPTY

_ERROR_PREFIX = "[PrevalenceBiasIndex] Error"


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
      often.

    I store the four cells and work the indices out from them, so they can
    never disagree with the counts.

    Attributes
    ----------
    positive_category, negative_category : str
        The first and second declared categories.
    rater_a, rater_b : str
        The two raters, in order.
    n_both_positive : int
        a: subjects both raters called positive.
    n_both_negative : int
        d: subjects both raters called negative.
    n_only_b_positive : int
        b: rater_b positive, rater_a negative.
    n_only_a_positive : int
        c: rater_a positive, rater_b negative.
    """

    positive_category: str
    negative_category: str
    rater_a: str
    rater_b: str
    n_both_positive: int
    n_both_negative: int
    n_only_b_positive: int
    n_only_a_positive: int

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        for name in ("positive_category", "negative_category", "rater_a", "rater_b"):
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

        if self.rater_a == self.rater_b:
            raise ValueError(
                f"{_ERROR_PREFIX}: rater_a and rater_b are both {self.rater_a!r}"
            )

        for name in ("n_both_positive", "n_both_negative",
                     "n_only_b_positive", "n_only_a_positive"):
            _check_count(name, getattr(self, name))
        if self.n_subjects < 1:
            raise ValueError(f"{_ERROR_PREFIX}: no subjects")

    @classmethod
    def build_from_pair(
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
            n_only_b_positive=int((b_pos & ~a_pos).sum()),
            n_only_a_positive=int((a_pos & ~b_pos).sum()),
        )

    # Public API ------------------------------------------------------------- #

    @property
    def n_disagree(self) -> int:
        """b + c: subjects the raters disagreed on."""
        return self.n_only_b_positive + self.n_only_a_positive

    @property
    def n_subjects(self) -> int:
        """N: subjects both raters rated."""
        return self.n_both_positive + self.n_both_negative + self.n_disagree

    @property
    def prevalence_index(self) -> float:
        """(a - d) / N. Above zero: more agreement on the positive category."""
        return (self.n_both_positive - self.n_both_negative) / self.n_subjects

    @property
    def bias_index(self) -> float:
        """(b - c) / N. Above zero: rater_b says positive more often."""
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
        if bi > 0:
            bi_words = f"{self.rater_b} says {pos} more often than {self.rater_a}"
        elif bi < 0:
            bi_words = f"{self.rater_a} says {pos} more often than {self.rater_b}"
        else:
            bi_words = f"{self.rater_a} and {self.rater_b} say {pos} equally often"
        return (
            f"Prevalence and bias (positive = {pos}, {self.rater_a} vs "
            f"{self.rater_b}, {self.n_subjects:,} subjects)\n"
            f"    PI = {pi:+.3f}: {pi_words}\n"
            f"    BI = {bi:+.3f}: {bi_words}"
        )

    def summary(self) -> dict:
        """
        Return my record in plain types.

        Returns
        -------
        dict
            The categories, the raters, the four cells, and both indices.
        """
        return {
            "positive_category": self.positive_category,
            "negative_category": self.negative_category,
            "rater_a": self.rater_a,
            "rater_b": self.rater_b,
            "n_subjects": self.n_subjects,
            "n_both_positive": self.n_both_positive,
            "n_both_negative": self.n_both_negative,
            "n_disagree": self.n_disagree,
            "n_only_b_positive": self.n_only_b_positive,
            "n_only_a_positive": self.n_only_a_positive,
            "prevalence_index": self.prevalence_index,
            "bias_index": self.bias_index,
        }


def _check_count(name: str, value: object) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(
            f"{_ERROR_PREFIX}: {name} must be int, got {type(value).__name__}"
        )
    if value < 0:
        raise ValueError(f"{_ERROR_PREFIX}: {name} must be >= 0, got {value}")
