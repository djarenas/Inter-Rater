"""Result of a hypothesis test for inter-rater agreement coefficients."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_ERROR_PREFIX = "[HypothesisTestResult] Error"

# The coefficients with a hypothesis test, by the names the pairwise and
# subset configs already use.
TESTED_COEFFICIENTS: frozenset[str] = frozenset([
    "fleiss_kappa",
    "scotts_pi",
    "gwet_ac1",
    "brennan_prediger",
    "krippendorff_alpha",
    "cohen_kappa",
    "conger_kappa",
    "light_kappa",
])


def compute_permutation_p_value(
    n_at_or_above: int, n_at_or_below: int, n_usable: int, alternative: str
) -> float:
    """
    Compute a permutation p-value.

    Parameters
    ----------
    n_at_or_above : int
        Permutations whose coefficient was at or above the observed one.
    n_at_or_below : int
        Permutations whose coefficient was at or below the observed one.
    n_usable : int
        Permutations whose coefficient was defined (not NaN).
    alternative : str
        "greater", "less" or "two-sided".

    Returns
    -------
    float
        "greater": (n_at_or_above + 1) / (n_usable + 1)
        "less": (n_at_or_below + 1) / (n_usable + 1)
        "two-sided": twice the smaller of the two, capped at 1

    The +1 counts the observed data as one of the arrangements, so the
    p-value is never zero. For permutations drawn at random, as here, it is
    valid but conservative (Phipson & Smyth 2010, Section 6), and it is
    what scipy.stats.permutation_test computes.

    Two-sided is twice the smaller one-sided p-value, capped at 1: Fisher's
    rule, which gives equal weight to departures in either direction
    (Kulinskaya 2008, Section 1), and scipy.stats.permutation_test's
    default. It reads each tail directly, so it does not assume the null
    distribution is centred on 0. When the null distribution is
    asymmetric, ``alternative="greater"`` or ``"less"`` gives the one-sided
    p-value instead (Gibbons & Pratt 1975).

    References
    ----------
    Gibbons, J. D., & Pratt, J. W. (1975). P-values: interpretation and
        methodology. The American Statistician, 29(1), 20-25.
    Kulinskaya, E. (2008). On two-sided p-values for non-symmetric
        distributions. arXiv:0810.2124.
    Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never
        be zero: calculating exact P-values when permutations are randomly
        drawn. Statistical Applications in Genetics and Molecular Biology,
        9(1), Article 39.
    """
    greater = (n_at_or_above + 1) / (n_usable + 1)
    less = (n_at_or_below + 1) / (n_usable + 1)
    if alternative == "greater":
        return greater
    if alternative == "less":
        return less
    if alternative == "two-sided":
        return min(1.0, 2.0 * min(greater, less))
    raise ValueError(
        f"{_ERROR_PREFIX}: alternative must be 'two-sided', 'greater' or "
        f"'less', got {alternative!r}"
    )


@dataclass(frozen=True)
class HypothesisTestResult:
    """
    I hold the result of a hypothesis test against H₀: coefficient = 0.

    I store the test outcome and enough metadata to interpret it.
    Different methods produce different statistics: z_statistic
    is only meaningful for parametric_null, while null_distribution
    is only meaningful for permutation.

    Attributes
    ----------
    p_value : float
        The p-value of the test, for the alternative recorded in
        ``alternative``. Must be in [0, 1].
    alpha : float
        The significance level used. Must be in (0, 1).
    reject_null : bool
        Whether H₀ is rejected at the given alpha. Must be consistent
        with p_value and alpha.
    method : str
        The method used: "permutation" or "parametric_null".
    point : float
        The observed point estimate being tested. Must not be NaN: an
        undefined coefficient has no agreement to test.
    coefficient : str
        Which coefficient was tested, e.g. 'fleiss_kappa'. Every
        coefficient returns this same class, so without it a result read
        on its own -- printed, summarised or plotted -- could not say what
        it is a test of. One of TESTED_COEFFICIENTS.
    z_statistic : float | None
        Z-score (parametric_null only, None for permutation).
    null_variance : float | None
        The null variance used (parametric_null only). Must be >= 0.
    null_distribution : np.ndarray | None
        The permutation distribution (permutation only, None for parametric).
    n_permutations : int | None
        Number of permutations performed (permutation only). Must be >= 1.
    alternative : str | None
        "two-sided", "greater" or "less". None if not recorded.
    n_usable : int | None
        Permutations whose coefficient was defined (permutation only).
        Draws that came out NaN are skipped, and this says how many were
        left. At most n_permutations.
    n_at_or_above : int | None
        Usable permutations at or above the observed point.
    n_at_or_below : int | None
        Usable permutations at or below the observed point.
    scheme : str | None
        The permutation scheme: "within_rater_shuffle" or
        "pool_and_redistribute" (permutation only). None if not recorded.
    random_state : int | None
        The seed of the permutations (permutation only). None if not
        recorded.

    When the three counts are recorded, p_value must follow from them, as
    in Phipson & Smyth (2010), counting the observed data as one of the
    possible arrangements::

        greater:   (n_at_or_above + 1) / (n_usable + 1)
        less:      (n_at_or_below + 1) / (n_usable + 1)
        two-sided: min(1, 2 * min(greater, less))

    so anyone reading the result can recompute the p-value from it.
    """

    # Required fields
    p_value: float
    alpha: float
    reject_null: bool
    method: str
    point: float
    coefficient: str

    # Parametric-null fields (None for permutation)
    z_statistic: float | None = None
    null_variance: float | None = None

    # Permutation fields (None for parametric)
    null_distribution: np.ndarray | None = None
    n_permutations: int | None = None
    scheme: str | None = None
    random_state: int | None = None

    # How the p-value was computed
    alternative: str | None = None
    n_usable: int | None = None
    n_at_or_above: int | None = None
    n_at_or_below: int | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate hypothesis test result attributes."""

        # Refuse an undefined coefficient first. Every comparison with NaN
        # is False, so a permutation test on it counts no null value as
        # extreme and reports p = 0 -- a rejection about a number that
        # does not exist.
        if isinstance(self.point, (int, float)) and np.isnan(self.point):
            raise ValueError(
                f"{_ERROR_PREFIX}: the coefficient is undefined (NaN) on "
                f"these ratings, so there is no agreement to test. This "
                f"happens when, for example, every rating is in the same "
                f"category."
            )

        # Validate p_value
        if not isinstance(self.p_value, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: p_value must be float, "
                f"got {type(self.p_value).__name__}"
            )
        if not (0.0 <= self.p_value <= 1.0):
            raise ValueError(
                f"{_ERROR_PREFIX}: p_value must be in [0, 1], "
                f"got {self.p_value}"
            )

        # Validate alpha
        if not isinstance(self.alpha, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: alpha must be float, "
                f"got {type(self.alpha).__name__}"
            )
        if not (0.0 < self.alpha < 1.0):
            raise ValueError(
                f"{_ERROR_PREFIX}: alpha must be in (0, 1), "
                f"got {self.alpha}"
            )

        # Validate reject_null
        if not isinstance(self.reject_null, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: reject_null must be bool, "
                f"got {type(self.reject_null).__name__}"
            )

        # Validate reject_null consistency with p_value and alpha
        # (boundary case p_value == alpha accepts either decision)
        if self.p_value < self.alpha and not self.reject_null:
            raise ValueError(
                f"{_ERROR_PREFIX}: reject_null is False but "
                f"p_value ({self.p_value}) < alpha ({self.alpha})"
            )
        if self.p_value > self.alpha and self.reject_null:
            raise ValueError(
                f"{_ERROR_PREFIX}: reject_null is True but "
                f"p_value ({self.p_value}) > alpha ({self.alpha})"
            )

        # Validate method
        if not isinstance(self.method, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: method must be str, "
                f"got {type(self.method).__name__}"
            )
        if self.method not in {"permutation", "parametric_null"}:
            raise ValueError(
                f"{_ERROR_PREFIX}: method must be 'permutation' or "
                f"'parametric_null', got '{self.method}'"
            )

        # Validate point
        if not isinstance(self.point, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: point must be float, "
                f"got {type(self.point).__name__}"
            )

        # Validate coefficient
        if not isinstance(self.coefficient, str):
            raise TypeError(
                f"{_ERROR_PREFIX}: coefficient must be str, "
                f"got {type(self.coefficient).__name__}"
            )
        if self.coefficient not in TESTED_COEFFICIENTS:
            raise ValueError(
                f"{_ERROR_PREFIX}: coefficient must be one of "
                f"{sorted(TESTED_COEFFICIENTS)}, got {self.coefficient!r}"
            )

        # Cross-validate method-specific fields
        if self.method == "parametric_null":
            if self.z_statistic is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: z_statistic must not be None "
                    f"for method='parametric_null'"
                )
            if self.null_variance is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: null_variance must not be None "
                    f"for method='parametric_null'"
                )
        if self.method == "permutation":
            if self.n_permutations is None:
                raise ValueError(
                    f"{_ERROR_PREFIX}: n_permutations must not be None "
                    f"for method='permutation'"
                )

        # Validate optional field values when present
        if self.z_statistic is not None:
            if not isinstance(self.z_statistic, (int, float)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: z_statistic must be float, "
                    f"got {type(self.z_statistic).__name__}"
                )
        if self.null_variance is not None:
            if not isinstance(self.null_variance, (int, float)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: null_variance must be float, "
                    f"got {type(self.null_variance).__name__}"
                )
            if self.null_variance < 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: null_variance must be >= 0, "
                    f"got {self.null_variance}"
                )
        if self.null_distribution is not None:
            if not isinstance(self.null_distribution, np.ndarray):
                raise TypeError(
                    f"{_ERROR_PREFIX}: null_distribution must be np.ndarray, "
                    f"got {type(self.null_distribution).__name__}"
                )
        if self.n_permutations is not None:
            if not isinstance(self.n_permutations, int):
                raise TypeError(
                    f"{_ERROR_PREFIX}: n_permutations must be int, "
                    f"got {type(self.n_permutations).__name__}"
                )
            if self.n_permutations < 1:
                raise ValueError(
                    f"{_ERROR_PREFIX}: n_permutations must be >= 1, "
                    f"got {self.n_permutations}"
                )

        self._validate_counts()
        self._validate_scheme_and_seed()

    def _validate_scheme_and_seed(self) -> None:
        """Validate scheme and random_state; both belong to permutation only."""
        if self.method == "parametric_null" and (
            self.scheme is not None or self.random_state is not None
        ):
            raise ValueError(
                f"{_ERROR_PREFIX}: scheme and random_state describe permutations "
                f"and must be None for method='parametric_null'"
            )
        if self.scheme is not None and self.scheme not in {
            "within_rater_shuffle", "pool_and_redistribute"
        }:
            raise ValueError(
                f"{_ERROR_PREFIX}: scheme must be 'within_rater_shuffle' or "
                f"'pool_and_redistribute', got {self.scheme!r}"
            )
        if self.random_state is not None and (
            isinstance(self.random_state, bool)
            or not isinstance(self.random_state, (int, np.integer))
        ):
            raise TypeError(
                f"{_ERROR_PREFIX}: random_state must be int or None, "
                f"got {type(self.random_state).__name__}"
            )

    def _validate_counts(self) -> None:
        """Validate alternative and the counts, and that p_value follows from them."""
        if self.alternative is not None and self.alternative not in {
            "two-sided", "greater", "less"
        }:
            raise ValueError(
                f"{_ERROR_PREFIX}: alternative must be 'two-sided', "
                f"'greater' or 'less', got {self.alternative!r}"
            )

        counts = {
            "n_usable": self.n_usable,
            "n_at_or_above": self.n_at_or_above,
            "n_at_or_below": self.n_at_or_below,
        }
        given = {k: v for k, v in counts.items() if v is not None}
        if not given:
            return
        if len(given) < 3:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_usable, n_at_or_above and "
                f"n_at_or_below are recorded together; got only "
                f"{sorted(given)}"
            )
        for name, value in counts.items():
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise TypeError(
                    f"{_ERROR_PREFIX}: {name} must be int, "
                    f"got {type(value).__name__}"
                )
        if self.n_usable < 1:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_usable must be >= 1, got {self.n_usable}"
            )
        if self.n_permutations is not None and self.n_usable > self.n_permutations:
            raise ValueError(
                f"{_ERROR_PREFIX}: n_usable ({self.n_usable}) cannot exceed "
                f"n_permutations ({self.n_permutations})"
            )
        for name in ("n_at_or_above", "n_at_or_below"):
            if not 0 <= counts[name] <= self.n_usable:
                raise ValueError(
                    f"{_ERROR_PREFIX}: {name} must be in [0, n_usable="
                    f"{self.n_usable}], got {counts[name]}"
                )
        if self.alternative is None:
            raise ValueError(
                f"{_ERROR_PREFIX}: the counts are recorded, so alternative "
                f"must be too: the p-value depends on it"
            )

        expected = compute_permutation_p_value(
            self.n_at_or_above, self.n_at_or_below, self.n_usable,
            self.alternative,
        )
        if not np.isclose(self.p_value, expected, rtol=0, atol=1e-12):
            raise ValueError(
                f"{_ERROR_PREFIX}: p_value ({self.p_value}) does not follow "
                f"from the counts: alternative={self.alternative!r} gives "
                f"{expected}"
            )

    # Public API ------------------------------------------------------------- #

    def __str__(self) -> str:
        """
        Human-readable summary of the hypothesis test result.

        Examples
        --------
        Parametric:
            H₀: Fleiss' κ = 0, z = 3.42, two-sided p = 0.0006, reject at α = 0.05
        Permutation:
            H₀: Fleiss' κ = 0, two-sided p = 0.0030 (5000 permutations,
            within_rater_shuffle, seed 42), reject at α = 0.05. Two-sided p
            is twice the smaller one-sided p (Fisher's rule; Kulinskaya 2008)
        """
        from .analyzer_results import result_class_for

        # The coefficient as every other result prints it, e.g. "Fleiss' κ"
        name = result_class_for(self.coefficient).NAME

        # Format p-value: use scientific notation for very small values
        if self.p_value < 0.000001:
            p_str = "p < 0.000001"
        else:
            p_str = f"p = {self.p_value:.6f}"

        if self.alternative is not None:
            p_str = f"{self.alternative} {p_str}"

        # Decision
        decision = "reject" if self.reject_null else "fail to reject"

        if self.method == "parametric_null" and self.z_statistic is not None:
            return (
                f"H₀: {name} = 0, z = {self.z_statistic:.2f}, "
                f"{p_str}, {decision} at α = {self.alpha}"
            )
        elif self.method == "permutation" and self.n_permutations is not None:
            if self.n_usable is not None and self.n_usable < self.n_permutations:
                draws = (
                    f"{self.n_usable} of {self.n_permutations} permutations "
                    f"usable"
                )
            else:
                draws = f"{self.n_permutations} permutations"
            if self.scheme is not None:
                draws += f", {self.scheme}"
            if self.random_state is not None:
                draws += f", seed {self.random_state}"
            note = (
                ". Two-sided p is twice the smaller one-sided p (Fisher's "
                "rule; Kulinskaya 2008)"
                if self.alternative == "two-sided" else ""
            )
            return (
                f"H₀: {name} = 0, {p_str} ({draws}), "
                f"{decision} at α = {self.alpha}{note}"
            )
        else:
            return (
                f"H₀: {name} = 0, {p_str}, method={self.method}, "
                f"{decision} at α = {self.alpha}"
            )

    def summary(self) -> dict:
        """
        Return a plain-type summary of the hypothesis test result.

        Returns
        -------
        dict
            Keys: p_value, alpha, reject_null, method, point,
            and method-specific fields.
        """
        result = {
            "p_value": self.p_value,
            "alpha": self.alpha,
            "reject_null": self.reject_null,
            "method": self.method,
            "point": self.point,
            "coefficient": self.coefficient,
        }

        if self.z_statistic is not None:
            result["z_statistic"] = self.z_statistic
        if self.null_variance is not None:
            result["null_variance"] = self.null_variance
        if self.n_permutations is not None:
            result["n_permutations"] = self.n_permutations
        for name in ("alternative", "n_usable", "n_at_or_above", "n_at_or_below"):
            value = getattr(self, name)
            if value is not None:
                result[name] = value if name == "alternative" else int(value)
        if self.scheme is not None:
            result["scheme"] = self.scheme
        if self.random_state is not None:
            result["random_state"] = int(self.random_state)

        return result
