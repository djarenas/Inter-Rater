"""Configuration for distribution-based rater injection."""

from __future__ import annotations

from dataclasses import dataclass

from .injector_config import InjectorConfig

_ERROR_PREFIX = "[DistributionInjectorConfig] Error"


@dataclass(frozen=True)
class DistributionInjectorConfig(InjectorConfig):
    """
    I configure a rater drawn from a categorical distribution.

    The injected rater has no relationship to the real data -- every
    rating is an independent draw from the prevalence vector. This is
    the null model: "what does agreement look like if we add someone
    who guesses at these base rates?"

    When prevalences is None, the injector draws uniformly across all
    categories at generate time. When provided, prevalences must be
    valid probabilities (non-negative, sum to 1.0). The injector
    validates that len(prevalences) == n_categories.

    For binary categories, this reduces to Bernoulli draws.

    Attributes
    ----------
    prevalences : tuple of float, or None
        Probability of each category, in category order. None means
        uniform (equal probability for all categories, computed at
        generate time). When provided, must sum to 1.0 (within 1e-6)
        and have length == n_categories -- the length check happens
        at generate time in the injector.
    n_new_raters : int
        Inherited from InjectorConfig.
    random_state : int or None
        Inherited from InjectorConfig.
    """

    prevalences: tuple[float, ...] | None = None

    # Constructor ------------------------------------------------------------ #

    def __post_init__(self) -> None:
        """Validate distribution-specific fields."""
        try:
            super().__post_init__()
        except (TypeError, ValueError) as exc:
            raise type(exc)(f"{_ERROR_PREFIX}: {exc}") from exc

        # prevalences: None is valid (uniform at generate time)
        if self.prevalences is not None:
            if not isinstance(self.prevalences, tuple):
                raise TypeError(
                    f"{_ERROR_PREFIX}: prevalences must be a tuple "
                    f"or None, got {type(self.prevalences).__name__}"
                )

            if len(self.prevalences) == 0:
                raise ValueError(
                    f"{_ERROR_PREFIX}: prevalences cannot be empty. "
                    f"Provide one probability per category, or None "
                    f"for uniform."
                )

            for i, p in enumerate(self.prevalences):
                if isinstance(p, bool):
                    raise TypeError(
                        f"{_ERROR_PREFIX}: prevalences[{i}] must be "
                        f"numeric, got bool"
                    )
                if not isinstance(p, (int, float)):
                    raise TypeError(
                        f"{_ERROR_PREFIX}: prevalences[{i}] must be "
                        f"numeric, got {type(p).__name__}"
                    )
                if p < 0:
                    raise ValueError(
                        f"{_ERROR_PREFIX}: prevalences[{i}] must be "
                        f">= 0, got {p}"
                    )

            total = sum(self.prevalences)
            if abs(total - 1.0) > 1e-6:
                raise ValueError(
                    f"{_ERROR_PREFIX}: prevalences must sum to 1.0, "
                    f"got {total:.6f}"
                )

    # Public API ------------------------------------------------------------- #

    @classmethod
    def build_from_dict(
        cls, data_dict: dict
    ) -> DistributionInjectorConfig:
        """
        Build a config from a plain dictionary.

        I convert prevalences from a list (as YAML provides) to a
        tuple, then delegate to the base class.

        Parameters
        ----------
        data_dict : dict
            Configuration as plain Python types.

        Returns
        -------
        DistributionInjectorConfig
            Instantiated config object with validation applied.
        """
        # Hydrate: list -> tuple for prevalences
        if "prevalences" in data_dict:
            val = data_dict["prevalences"]
            if isinstance(val, list):
                data_dict["prevalences"] = tuple(val)

        return super().build_from_dict(data_dict)

    def summary(self) -> dict:
        """
        Return configuration as a dictionary for documentation and audit.

        Returns
        -------
        dict
            Base fields plus prevalences.
        """
        result = super().summary()
        result["prevalences"] = self.prevalences
        return result
