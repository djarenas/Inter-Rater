"""Configuration for within-subject rater injection."""

from __future__ import annotations

from dataclasses import dataclass

from .injector_config import InjectorConfig

_ERROR_PREFIX = "[WithinSubjectInjectorConfig] Error"


@dataclass(frozen=True)
class WithinSubjectInjectorConfig(InjectorConfig):
    """
    I configure a rater who copies real answers at each subject.

    At every subject, the injected rater randomly selects one of the real
    raters who rated that subject and copies their answer. This is the
    most structure-preserving strategy: where the panel agreed, the new
    rater agrees too; where they disagreed, the new rater picks a side.

    With noise_rate > 0, I introduce manufactured disagreement: at each
    subject, with probability noise_rate the injected rater picks a
    category that none of the real raters used for that subject. If all
    categories were already used, the rater falls back to copying
    normally (no unused category available).

    This simulates "a rater who tried" -- someone reading from the same
    signal as the panel, just noisily.

    Statistically, this is a within-subject bootstrap of rater identity,
    with optional injected noise.

    Attributes
    ----------
    noise_rate : float
        Probability of picking an unused category instead of copying a
        real answer. 0.0 means pure copying (default). 1.0 means always
        pick an unused category when possible.
    n_new_raters : int
        Inherited from InjectorConfig. How many synthetic raters to
        inject.
    random_state : int or None
        Inherited from InjectorConfig. Seed for reproducibility.
    """

    noise_rate: float = 0.0

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate within-subject-specific fields."""
        try:
            super().__post_init__()
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # noise_rate: float in [0.0, 1.0], no bools
        if isinstance(self.noise_rate, bool):
            raise TypeError(
                f"{_ERROR_PREFIX}: noise_rate must be a number, "
                f"got bool"
            )
        if not isinstance(self.noise_rate, (int, float)):
            raise TypeError(
                f"{_ERROR_PREFIX}: noise_rate must be a number, "
                f"got {type(self.noise_rate).__name__}"
            )
        if not (0.0 <= self.noise_rate <= 1.0):
            raise ValueError(
                f"{_ERROR_PREFIX}: noise_rate must be in [0.0, 1.0], "
                f"got {self.noise_rate}"
            )

    # Public API ------------------------------------------------------------- #

    def summary(self) -> dict:
        """
        Return configuration as a dictionary for documentation and audit.

        Returns
        -------
        dict
            Base fields plus noise_rate.
        """
        result = super().summary()
        result["noise_rate"] = self.noise_rate
        return result
