"""Base class for the results that describe a coefficient on a panel of raters."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

if TYPE_CHECKING:
    from ..data_profiles.balance_profile import BalanceProfile
    from ..data_profiles.rater_profile import RaterProfile

_ERROR_PREFIX = "[PanelResult] Error"

# The weightings a result can record; "identity" is unweighted
_WEIGHTS_TYPES = ("identity", "linear", "quadratic", "custom")


@dataclass(frozen=True, kw_only=True)
class PanelResult(ABC):
    """
    I am the base for the results an analysis returns: a coefficient's
    result (PointEstimateResult) and a permutation analysis's result
    (PermutationResult). Both describe a coefficient computed on a panel of
    raters, so I hold what they share: the panel's profiles, the interval
    method, and the weighting, and I validate them once.

    Every result must answer "what happened?" through summary(), which is the
    interface for printing, serializing, and feeding to visualizers. I am
    immutable: a result is a fact about a dataset at a moment in time and
    should not change.

    I do not know how to visualize myself. That is the visualizer's job. I just
    hold the data. This separation means a new visualizer can be written
    without touching result classes, and a new result type can be added without
    requiring new visualizers -- as long as it holds the right data.

    My fields are keyword-only, so subclasses can keep their own fields
    positional.

    Attributes
    ----------
    balance_profile : BalanceProfile
        How many ratings each subject has.
    rater_profile : RaterProfile
        Who the raters are, how much each rated, how much each pair
        overlapped. Describes the same subjects as balance_profile.
    ci_method : str
        The interval method: one of my class's ``_CI_METHODS``.
    weights_type : str, default "identity"
        The weighting: "identity" (unweighted), "linear", "quadratic" or
        "custom". Not "identity" only where ``_accepts_weights()`` is True.
    weight_matrix : list[list[float]] | None, default None
        The weight matrix that was used, as nested lists; None exactly when
        weights_type is "identity".
    """

    # The interval methods my class allows. Fleiss, Cohen and Scott add
    # "parametric_null".
    _CI_METHODS: ClassVar[frozenset[str]] = frozenset({"percentile", "bca"})

    balance_profile: BalanceProfile
    rater_profile: RaterProfile
    ci_method: str
    weights_type: str = "identity"
    weight_matrix: list[list[float]] | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate the profiles, the interval method and the weighting."""
        self._check_profiles()
        self._check_ci_method()
        self._check_weights()

    def _check_profiles(self) -> None:
        """Both profiles are there, of the right type, about the same data."""
        from ..data_profiles.balance_profile import BalanceProfile
        from ..data_profiles.rater_profile import RaterProfile

        if self.balance_profile is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: balance_profile must not be None"
            )
        if not isinstance(self.balance_profile, BalanceProfile):
            raise TypeError(
                f"{_ERROR_PREFIX}: balance_profile must be BalanceProfile, "
                f"got {type(self.balance_profile).__name__}"
            )
        if not isinstance(self.rater_profile, RaterProfile):
            raise TypeError(
                f"{_ERROR_PREFIX}: rater_profile must be RaterProfile, got "
                f"{type(self.rater_profile).__name__}"
            )
        if self.rater_profile.n_subjects != self.balance_profile.n_subjects:
            raise ValueError(
                f"{_ERROR_PREFIX}: rater_profile has "
                f"{self.rater_profile.n_subjects} subjects but "
                f"balance_profile has {self.balance_profile.n_subjects}"
            )

    def _check_ci_method(self) -> None:
        """The interval method is one my class allows."""
        if not isinstance(self.ci_method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: ci_method must be str, "
                f"got {type(self.ci_method).__name__}"
            )
        if self.ci_method == "parametric_null" and (
            self.ci_method not in self._CI_METHODS
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: the package does not implement a parametric "
                f"null for this result, so ci_method must be one of "
                f"{sorted(self._CI_METHODS)}, got 'parametric_null'"
            )
        if self.ci_method not in self._CI_METHODS:
            raise ValueError(
                f"{_ERROR_PREFIX}: ci_method must be one of "
                f"{sorted(self._CI_METHODS)}, got '{self.ci_method}'"
            )

    def _check_weights(self) -> None:
        """The weighting is stated, matches its matrix, and is allowed."""
        if not isinstance(self.weights_type, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: weights_type must be str, "
                f"got {type(self.weights_type).__name__}"
            )
        if self.weights_type not in _WEIGHTS_TYPES:
            raise ValueError(
                f"{_ERROR_PREFIX}: weights_type must be one of "
                f"{list(_WEIGHTS_TYPES)}, got '{self.weights_type}'"
            )
        if (self.weights_type == "identity") != (self.weight_matrix is None):
            raise ValueError(
                f"{_ERROR_PREFIX}: weight_matrix must be None exactly when "
                f"weights_type is 'identity' (unweighted); got weights_type "
                f"'{self.weights_type}' with "
                f"{'no matrix' if self.weight_matrix is None else 'a matrix'}"
            )
        if self.weight_matrix is not None:
            matrix = np.asarray(self.weight_matrix, dtype=float)
            if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
                raise ValueError(
                    f"{_ERROR_PREFIX}: weight_matrix must be square, "
                    f"got shape {matrix.shape}"
                )
        if self.weights_type != "identity" and not self._accepts_weights():
            raise ValueError(
                f"{_ERROR_PREFIX}: this coefficient has no weighted form, so "
                f"weights_type must be 'identity', got '{self.weights_type}'"
            )

    def _accepts_weights(self) -> bool:
        """Whether my coefficient has a weighted form. Subclasses say."""
        # Just an extra precaution to ensure proper inheritance of the
        # abstract method.
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: _accepts_weights() must be implemented by "
            f"subclass"
        )

    # Public API ------------------------------------------------------------- #

    @abstractmethod
    def summary(self) -> dict:
        """
        Return everything I represent, in plain Python types (dict, list,
        float, str).

        This is the contract for downstream consumers: reporters know how to
        serialize dicts, visualizers dispatch on dict structure, and users can
        always read the output without understanding my internal shape.

        Returns
        -------
        dict
            A plain-type representation of my contents. Keys and structure
            depend on my type.
        """
        # Just an extra precaution to ensure proper inheritance of the
        # abstract method.
        raise NotImplementedError(
            f"{_ERROR_PREFIX}: summary() must be implemented by subclass"
        )

    # Helpers for subclasses ------------------------------------------------- #

    def _panel_summary(self) -> dict:
        """The shared part of summary(): profiles, interval method, weights."""
        return {
            "ci_method": self.ci_method,
            "weights_type": self.weights_type,
            "weight_matrix": self.weight_matrix,
            "balance": self.balance_profile.summary(),
            "raters": self.rater_profile.summary(),
        }

    def _subjects_phrase(self) -> str:
        """'1,617 subjects'."""
        return f"{self.balance_profile.n_subjects:,} subjects"

    def _raters_phrase(self) -> str:
        """'4 raters'."""
        return f"{self.rater_profile.n_raters} raters"
