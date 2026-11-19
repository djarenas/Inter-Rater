"""What the data says about agreement, whatever coefficient is computed on it."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ...contracts import EMPTY
from ..not_computed import NotComputed
from .prevalence_bias_index import PrevalenceBiasIndex

if TYPE_CHECKING:
    from ...data_objects.ratings import Ratings

_ERROR_PREFIX = "[AgreementProfile] Error"
_TOL = 1e-9


@dataclass(frozen=True)
class AgreementProfile:
    """
    I describe agreement in the data itself, before any coefficient: how
    often each category was used, how often raters agreed on each one, how
    each rater used the categories. In the case of two raters and two
    categories, I also hold information on the prevalence and bias indices.

    This information helps a reader understand the context of agreement
    coefficients. For example, kappa falls when one category dominates, or when raters differ in how
    often they use a category, even if they agree on almost every subject.

    I am computed from the subjects every coefficient uses: those rated at
    least twice. A subject rated once forms no pair, so it enters neither
    agreement nor chance.

    Attributes
    ----------
    category_names : tuple[str, ...]
        In declared order. From a nij matrix, which carries no names, they
        are "column 0", "column 1", ...
    n_subjects : int
        Subjects rated at least twice: the pool the coefficients use.
    n_ratings : int
        Ratings on those subjects.
    prevalence : dict[str, float]
        p_j for each category: the share of each subject's ratings in that
        category, averaged over subjects. The same p_j Fleiss', Scott's and
        Gwet's chance terms use; on unbalanced data it differs from a plain
        share of all ratings (see the agreement-functions README).
    specific_agreement : dict[str, float]
        For each category j, the share of rater pairs involving j that
        agree on j: Sum_i n_ij(n_ij - 1) / Sum_i n_ij(n_i - 1). NaN for a
        category nobody used. Unweighted, even for weighted coefficients.
        On balanced data, Po = Sum_j p_j * specific_agreement_j; on
        unbalanced data that identity does not hold, because Po averages
        per subject and this pools pairs.
    rater_marginals : dict[str, dict[str, float]]
        For each rater, the share of their ratings (on these subjects) in
        each category.
    prevalence_bias : PrevalenceBiasIndex | NotComputed
        Only for two raters and two categories; otherwise NotComputed with
        the reason.
    """

    category_names: tuple
    n_subjects: int
    n_ratings: int
    prevalence: dict
    specific_agreement: dict
    rater_marginals: dict
    prevalence_bias: PrevalenceBiasIndex | NotComputed

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        names = self.category_names
        if not isinstance(names, tuple) or len(names) < 2:
            raise TypeError(
                f"{_ERROR_PREFIX}: category_names must be a tuple of at least "
                f"2 names, got {names!r}"
            )
        if len(set(names)) != len(names):
            raise ValueError(f"{_ERROR_PREFIX}: category_names repeat: {names}")

        for name in ("n_subjects", "n_ratings"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: {name} must be an int >= 0, got {value!r}"
                )
        if self.n_ratings < 2 * self.n_subjects:
            raise ValueError(
                f"{_ERROR_PREFIX}: {self.n_ratings} ratings cannot cover "
                f"{self.n_subjects} subjects rated at least twice"
            )

        _check_shares("prevalence", self.prevalence, names, self.n_subjects > 0)

        if list(self.specific_agreement) != list(names):
            raise ValueError(
                f"{_ERROR_PREFIX}: specific_agreement must have one entry per "
                f"category, in declared order"
            )
        for cat, value in self.specific_agreement.items():
            if not (math.isnan(value) or -_TOL <= value <= 1 + _TOL):
                raise ValueError(
                    f"{_ERROR_PREFIX}: specific agreement for {cat!r} is "
                    f"{value}, outside [0, 1]"
                )

        if isinstance(self.rater_marginals, dict) and self.rater_marginals:
            for rater, shares in self.rater_marginals.items():
                # A rater who rated none of these subjects has no shares
                if all(math.isnan(v) for v in shares.values()):
                    if list(shares) != list(names):
                        raise ValueError(
                            f"{_ERROR_PREFIX}: rater {rater!r} must have one "
                            f"share per category"
                        )
                    continue
                _check_shares(f"rater {rater!r}", shares, names, True)
        else:
            raise TypeError(
                f"{_ERROR_PREFIX}: rater_marginals must be a non-empty dict, "
                f"got {type(self.rater_marginals).__name__}"
            )

        if isinstance(self.prevalence_bias, PrevalenceBiasIndex):
            if len(names) != 2 or (
                self.prevalence_bias.positive_category,
                self.prevalence_bias.negative_category,
            ) != names:
                raise ValueError(
                    f"{_ERROR_PREFIX}: prevalence_bias must use these two "
                    f"categories, positive first"
                )
        elif not isinstance(self.prevalence_bias, NotComputed):
            raise TypeError(
                f"{_ERROR_PREFIX}: prevalence_bias must be PrevalenceBiasIndex "
                f"or NotComputed, got {type(self.prevalence_bias).__name__}"
            )

    @classmethod
    def build_from_ratings(cls, ratings: Ratings) -> AgreementProfile:
        """Build from ratings, where each rater's answers are known."""
        return cls.build_from_category_index(
            ratings.category_index,
            tuple(ratings.semantics.categories.names),
            tuple(ratings.semantics.raters),
        )

    @classmethod
    def build_from_category_index(
        cls,
        category_index: np.ndarray,
        category_names: tuple,
        rater_names: tuple,
    ) -> AgreementProfile:
        """
        Build from a (n_subjects, n_raters) array of category positions.

        Parameters
        ----------
        category_index : np.ndarray
            EMPTY where a rater did not rate a subject.
        category_names, rater_names : tuple
            In declared order.
        """
        index = np.asarray(category_index)
        k = len(category_names)
        present = index != EMPTY
        pool = present.sum(axis=1) >= 2
        index = index[pool]
        present = present[pool]

        counts = np.zeros((index.shape[0], k), dtype=np.int64)
        for j in range(k):
            counts[:, j] = (index == j).sum(axis=1)

        marginals = {}
        for r, rater in enumerate(rater_names):
            column = index[present[:, r], r].astype(np.int64)
            if column.size == 0:
                marginals[str(rater)] = {c: float("nan") for c in category_names}
            else:
                shares = np.bincount(column, minlength=k) / column.size
                marginals[str(rater)] = dict(zip(category_names, map(float, shares)))

        if len(rater_names) == 2 and k == 2:
            prevalence_bias = PrevalenceBiasIndex.build_from_pair(
                index, tuple(category_names), tuple(map(str, rater_names))
            )
        else:
            prevalence_bias = NotComputed(
                f"defined for two raters and two categories; this data has "
                f"{len(rater_names)} raters and {k} categories"
            )
        names = tuple(category_names)
        counts = counts.astype(np.float64)
        n_i = counts.sum(axis=1)
        n_subjects = int(counts.shape[0])
        if n_subjects:
            prevalence = (counts / n_i[:, None]).mean(axis=0)
        else:
            prevalence = np.full(len(names), np.nan)
        agree = (counts * (counts - 1)).sum(axis=0)
        involved = (counts * (n_i[:, None] - 1)).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            specific = np.where(involved > 0, agree / involved, np.nan)
        return cls(
            category_names=names,
            n_subjects=n_subjects,
            n_ratings=int(n_i.sum()),
            prevalence=dict(zip(names, map(float, prevalence))),
            specific_agreement=dict(zip(names, map(float, specific))),
            rater_marginals=marginals,
            prevalence_bias=prevalence_bias,
        )

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        names = self.category_names
        widths = [max(len(c), 5) for c in names]

        def row(label: str, values, label_width: int) -> str:
            cells = "  ".join(
                f"{_cell(v):>{w}}" for v, w in zip(values, widths)
            )
            return f"    {label:<{label_width}}  {cells}"

        raters = list(self.rater_marginals)
        label_width = max([len("overall")] + [len(r) for r in raters])
        header = row("", names, label_width)

        lines = [
            f"Agreement profile ({self.n_subjects:,} subjects rated at least "
            f"twice, {self.n_ratings:,} ratings)",
            "",
            "  How often each category was used (share of ratings)",
            header,
            row("overall", [self.prevalence[c] for c in names], label_width),
        ]
        for rater, shares in self.rater_marginals.items():
            lines.append(row(rater, [shares[c] for c in names], label_width))

        lines += [
            "",
            "  How often raters agreed on each category "
            "(share of rater pairs involving it that agree on it)",
            header,
            row("", [self.specific_agreement[c] for c in names], label_width),
            "",
        ]
        if isinstance(self.prevalence_bias, NotComputed):
            lines.append(
                f"  Prevalence and bias: not computed ({self.prevalence_bias.reason})"
            )
        else:
            lines.extend("  " + line for line in str(self.prevalence_bias).splitlines())
        return "\n".join(lines)

    def summary(self) -> dict:
        """
        Return my record in plain types. Anything not computed appears as
        ``{"not_computed": reason}``.
        """
        def plain(value):
            return value.summary() if hasattr(value, "summary") else value

        return {
            "category_names": list(self.category_names),
            "n_subjects": self.n_subjects,
            "n_ratings": self.n_ratings,
            "prevalence": dict(self.prevalence),
            "specific_agreement": dict(self.specific_agreement),
            "rater_marginals": {r: dict(s) for r, s in self.rater_marginals.items()},
            "prevalence_bias": plain(self.prevalence_bias),
        }


def _check_shares(label: str, shares: dict, names: tuple, must_sum: bool) -> None:
    if not isinstance(shares, dict) or list(shares) != list(names):
        raise ValueError(
            f"{_ERROR_PREFIX}: {label} must have one share per category, in "
            f"declared order"
        )
    values = np.array(list(shares.values()), dtype=np.float64)
    if not must_sum:
        if not np.isnan(values).all():
            raise ValueError(
                f"{_ERROR_PREFIX}: {label} must be NaN when no subject was "
                f"rated twice"
            )
        return
    if np.isnan(values).any() or (values < -_TOL).any():
        raise ValueError(f"{_ERROR_PREFIX}: {label} has a missing or negative share")
    if abs(values.sum() - 1.0) > 1e-6:
        raise ValueError(
            f"{_ERROR_PREFIX}: {label} shares sum to {values.sum():.6f}, not 1"
        )


def _cell(value) -> str:
    """A table cell: a name as is, a share to three places, NaN as n/a."""
    if isinstance(value, str):
        return value
    return "n/a" if math.isnan(value) else f"{value:.3f}"
