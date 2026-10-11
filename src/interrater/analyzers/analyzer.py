"""Base analyzer for inter-rater agreement coefficients."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..data_objects.ratings import Ratings
    from ..intermediates import AnalysisResult
    from ..intermediates.analysis_results import HypothesisTestResult
    from ..intermediates.data_profiles import BalanceProfile
    from .configs.analyzer_config import AnalyzerConfig

_ERROR_PREFIX = "[Analyzer] Error"


class Analyzer(ABC):
    """
    I am the abstract base for all inter-rater agreement analyzers.

    Every analyzer takes Ratings and a Config, validates both, computes a
    BalanceProfile that describes the data structure, and provides an
    estimate() method that returns an AnalysisResult.

    I am responsible for:
    - Storing valid Ratings and Config objects
    - Applying the injection pipeline (if config.injection is set),
      which adds synthetic raters to Ratings before any analysis
    - Computing BalanceProfile (relevant for most coefficients)
    - Providing the public estimate() interface
    - Delegating to subclasses for coefficient-specific logic

    Subclasses override __init__ to add their own validation
    (e.g., GroupAnalyzer checks balanced_data_only, WeightableGroupAnalyzer
    builds Weights, FleissKappaAnalyzer validates ratings type).

    Some subclasses may provide alternate constructors (classmethods) that accept
    NijMatrix directly instead of Ratings.

    Attributes
    ----------
    ratings : Ratings | None
        The validated ratings data. Type varies by analyzer:
        FleissKappaAnalyzer may require BalancedRatings if config.balanced_data_only=True.
        None if initialized from NijMatrix via alternate constructor.
    config : AnalyzerConfig
        The analysis configuration (FleissKappaConfig, etc.).
        Defines bootstrap resamples, confidence level, CI method, weights, etc.
    balance_profile : BalanceProfile
        Computed from ratings or nij_matrix. Describes how subjects were rated 
        (constant n per subject, distribution of subjects by rater count, fill rate, etc.). 
        Stored here for access by all analyzers and included in results.
    """

    # Attribute declaration

    # Populated at this level
    ratings: Ratings | None
    config: AnalyzerConfig
    balance_profile: BalanceProfile

    def __init__(
        self, ratings: Ratings | None, config: AnalyzerConfig
    ) -> None:
        """
        Initialize analyzer: validate inputs and compute BalanceProfile.

        Parameters
        ----------
        ratings : Ratings | None
            The ratings data to analyze. Can be None if subclass provides
            alternate constructor that uses NijMatrix directly.
        config : AnalyzerConfig
            The analysis configuration.

        Raises
        ------
        TypeError
            If ratings is not a Ratings object (when provided) or config is not AnalyzerConfig.
        ValueError
            If config is None.
        """
        from ..data_objects.ratings import Ratings as RatingsBase
        from .configs import AnalyzerConfig as AnalyzerConfigClass

        # Validate ratings type if provided
        if ratings is not None and not isinstance(ratings, RatingsBase):
            raise TypeError(
                f"{_ERROR_PREFIX}: ratings must be a Ratings object, "
                f"got {type(ratings).__name__}"
            )

        # Validate config is not None
        if config is None:
            raise ValueError(f"{_ERROR_PREFIX}: config cannot be None")

        # Validate config type
        if not isinstance(config, AnalyzerConfigClass):
            raise TypeError(
                f"{_ERROR_PREFIX}: config must be an AnalyzerConfig, "
                f"got {type(config).__name__}"
            )

        # Store inputs
        self.ratings = ratings
        self.config = config

        # Apply injection pipeline (before BalanceProfile)
        self.apply_injection()

        # Compute BalanceProfile from ratings if available
        from ..intermediates.data_profiles import BalanceProfile

        if self.ratings is not None:
            self.balance_profile = BalanceProfile.from_ratings(self.ratings)
        else:
            # BalanceProfile will be computed by subclass from nij_matrix
            self.balance_profile = None  # type: ignore

    # Public API ------------------------------------------------------------- #

    @abstractmethod
    def estimate(self) -> AnalysisResult:
        """
        Compute the agreement coefficient and return the result.

        Returns
        -------
        AnalysisResult
            Subclass-specific result (FleissKappaResult, CohenKappaResult, etc.)
            holding point estimate, confidence interval, and relevant metadata.

        Raises
        ------
        NotImplementedError
            If called on the abstract base class.
        ValueError
            If ratings or config are invalid for this analyzer.
        """
        raise NotImplementedError(
            f"[Analyzer] Error: estimate() must be implemented by subclass, "
            f"got {type(self).__name__}"
        )

    @abstractmethod
    def hypothesis_test(self, method: str = "permutation") -> HypothesisTestResult:
        """
        Test the null hypothesis that the agreement coefficient equals zero.

        Each leaf class implements this with the methods available for that
        coefficient. Common methods:
        - "permutation": permute ratings to build null distribution (requires Ratings)
        - "parametric_null": use analytical null variance formula (coefficient-specific)

        Parameters
        ----------
        method : str, default "permutation"
            The method to use for the hypothesis test.

        Returns
        -------
        HypothesisTestResult
            Test statistic, p-value, and metadata.

        Raises
        ------
        NotImplementedError
            If called on the abstract base class, or if the requested method
            is not implemented for this coefficient.
        ValueError
            If method="permutation" but ratings is None.
        """
        raise NotImplementedError(
            f"[Analyzer] Error: hypothesis_test() must be implemented by subclass, "
            f"got {type(self).__name__}"
        )

    def apply_injection(self) -> None:
        """
        I apply the injection pipeline from config to the ratings.

        I iterate through each injection step in ``self.config.injection``,
        instantiate the corresponding injector, and call ``inject()`` to
        produce a new Ratings with synthetic raters appended. Each step
        chains -- the output Ratings of one step becomes the input to the
        next. The new Ratings validates itself through its own constructor.

        I update ``self.ratings`` in place. If ``self.config.injection``
        is None or ``self.ratings`` is None (NijMatrix path), I do nothing.

        Called automatically in ``__init__`` before BalanceProfile is built,
        so all downstream computation sees the injected data.
        """
        if self.config.injection is None or self.ratings is None:
            return

        from ..injectors.registry import get_injector_class

        for step in self.config.injection:
            try:
                injector_cls = get_injector_class(step)
                injector = injector_cls(self.ratings, step)
                self.ratings = injector.inject()
            except Exception as e:
                raise type(e)(
                    f"[Analyzer] Error in apply_injection: "
                    f"injection step {type(step).__name__} failed: {e}"
                ) from e

    # Private helpers -------------------------------------------------------- #

    def _check_permutation_available(self) -> None:
        """
        Guard: permutation hypothesis test requires Ratings.

        Permutation tests shuffle individual rating assignments to build the
        null distribution. This requires the original Ratings object -- it
        cannot be done from an NijMatrix because rater identity is lost.

        Raises
        ------
        ValueError
            If self.ratings is None (initialized from NijMatrix).
        """
        if self.ratings is None:
            raise ValueError(
                f"[Analyzer] Error: method='permutation' requires Ratings data. "
                f"Use method='parametric_null' when initialized from NijMatrix."
            )
