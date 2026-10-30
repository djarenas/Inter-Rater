"""
Shared fixtures for test_example_data/.

Integration tests that load the shipped example datasets and run analyzers
end-to-end.  Each fixture returns either a Ratings (via RatingsCleaner) or
a NijMatrix (via pd.read_csv + NijMatrix), depending on the data format.

Construction recipe for ratings pathway:
    YAML -> RatingsSemantics.build_from_yaml
    YAML -> CleanerConfig.build_from_yaml
    csv/txt -> pd.read_csv -> RatingsCleaner(df, semantics, config) -> .ratings

Construction recipe for nij pathway:
    csv/txt -> pd.read_csv -> np.ndarray -> NijMatrix
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from interrater import (
    CleanerConfig,
    NijMatrix,
    OrderedCategories,
    RatingsCleaner,
    RatingsSemantics,
)
from interrater.analyzers.configs import (
    AnalyzerConfig,
    FleissKappaConfig,
    GwetAC1Config,
    BrennanPredigerConfig,
    KrippendorffAlphaConfig,
    ScottsConfig,
    CohenKappaConfig,
    CongerKappaConfig,
    LightKappaConfig,
    PairwisePermutationConfig,
    SubsetPermutationConfig,
    PermutationBootstrapMeanConfig,
)


# -- Paths ------------------------------------------------------------------ #

EXAMPLE_DATA_DIR = Path(__file__).resolve().parent.parent.parent / "example_data"

FLEISS_DIR     = EXAMPLE_DATA_DIR / "Fleiss1971"
GWET_DIR       = EXAMPLE_DATA_DIR / "Gwet_example"
SYNTHETIC_DIR  = EXAMPLE_DATA_DIR / "Synthetic_2raters"
SR1_DIR        = EXAMPLE_DATA_DIR / "SR1_OSF2018"
SR2_DIR        = EXAMPLE_DATA_DIR / "SR2_JSRC2019"
SR3_DIR        = EXAMPLE_DATA_DIR / "SR3_JGIM2019"
SR4_DIR        = EXAMPLE_DATA_DIR / "SR4_FIS_CVD"
SR5_DIR        = EXAMPLE_DATA_DIR / "SR5_CocCardio"
UNBALANCED_DIR = EXAMPLE_DATA_DIR / "Unbalanced_nij"


# -- Resample count --------------------------------------------------------- #

# Efron & Tibshirani (1986) section 9, p. 72: the rough minimum for a useful
# percentile interval. Every config in this file uses it, including the
# permutation ones, so that no interval these tests build ever reports
# below_recommended_resamples -- which would be correct of it, and noise here.
_INTEGRATION_RESAMPLES = 250


# -- Helpers ---------------------------------------------------------------- #

def _load_ratings(data_path: Path, semantics_path: Path,
                  cleaner_config_path: Path,
                  weighting: str | None = None) -> "Ratings":
    """
    Load a ratings dataset via the RatingsCleaner pipeline.

    With ``weighting`` ('linear' or 'quadratic'), the same file is read under
    ordered categories instead: same names, same order, declared ordinal.
    That is how a weighted analysis is asked for, since weights are declared
    in the semantics and not in the analyzer config.
    """
    semantics = RatingsSemantics.build_from_yaml(str(semantics_path))
    if weighting is not None:
        semantics = RatingsSemantics(
            raters          = semantics.raters,
            categories      = OrderedCategories(
                semantics.categories.names, weighting=weighting
            ),
            n_subjects      = semantics.n_subjects,
            unique_subjects = semantics.unique_subjects,
        )
    config = CleanerConfig.build_from_yaml(str(cleaner_config_path))

    # Determine separator: .csv -> comma, .txt -> tab
    sep = "," if data_path.suffix == ".csv" else "\t"

    # Determine whether the file has a header row.
    # This is a property of the file format, NOT of use_column_names.
    # use_column_names controls rater-label matching inside RatingsCleaner,
    # not whether the CSV/TXT has a header.  All shipped .csv files have
    # headers; all .txt files are headerless.
    has_header = data_path.suffix == ".csv"

    df = pd.read_csv(
        data_path,
        sep=sep,
        header=0 if has_header else None,
        dtype=str,
    )

    rc = RatingsCleaner(df, semantics, config=config)
    return rc.ratings


def _load_nij(data_path: Path, *, index_col: int | str | None = None,
              sep: str | None = None) -> NijMatrix:
    """
    Load an nij (counts) file into a NijMatrix.

    Parameters
    ----------
    data_path : Path
        CSV or TXT file with subject rows x category columns.
    index_col : int | str | None
        Column to use as the row index (subject IDs).  Pass 0 for the
        Fleiss CSV (unnamed first column), ``"Subject"`` for Gwet nij,
        or ``None`` for headerless TXT files.
    sep : str | None
        Separator; inferred from extension if None.
    """
    if sep is None:
        sep = "," if data_path.suffix == ".csv" else "\t"

    has_header = data_path.suffix == ".csv"

    df = pd.read_csv(
        data_path,
        sep=sep,
        header=0 if has_header else None,
        index_col=index_col,
    )

    return NijMatrix(df.values.astype(np.int32))


def _integration_config(config_cls, **overrides):
    """
    Build an analyzer config for integration testing.

    Percentile at B=250, which is two deliberate choices.

    **Percentile, because these tests assert that the point estimate falls
    inside its interval.** Under percentile that is a sound assertion. Under
    BCa it is not an invariant at all -- strong bias correction can put the
    interval entirely to one side of the point legitimately, which is why
    ``PointEstimateResult`` stopped requiring it. Asserting it on a BCa
    interval is asserting something that is allowed to be false.

    **B=250, because that is the minimum Efron & Tibshirani (1986) section 9
    give for a useful percentile interval.** Below it, every interval these
    tests build would report ``below_recommended_resamples`` and say so in
    its own description -- correct of it, and noise here.

    This file previously used BCa at B=50, which was under the published
    minimum of 1000 for that method and paid for a jackknife besides: BCa
    needs one estimator call per subject regardless of B, so on the larger
    fixtures the jackknife cost far more than the bootstrap. Percentile at
    250 runs the four largest fixtures in 0.2s against 1.8s for that.

    BCa is not left untested end-to-end -- see ``bca_fleiss_config`` below.
    Point estimates are deterministic and do not depend on any of this.
    """
    defaults = dict(
        n_resamples=_INTEGRATION_RESAMPLES,
        random_state=42,
        ci_method="percentile",
    )
    defaults.update(overrides)
    return config_cls(**defaults)


def _bca_config(config_cls, **overrides):
    """
    Build a BCa config at the resample count the literature asks for.

    B=1000 is Efron & Tibshirani's rough minimum for BCa. Used by a single
    end-to-end test, on a small fixture: the jackknife costs one estimator
    call per subject, so this belongs on 30 subjects rather than 1617.
    """
    defaults = dict(
        n_resamples=1000,
        random_state=42,
        ci_method="bca",
    )
    defaults.update(overrides)
    return config_cls(**defaults)


# -- NijMatrix fixtures ----------------------------------------------------- #

@pytest.fixture(scope="module")
def fleiss1971_nij():
    """NijMatrix from Fleiss (1971) psychiatric diagnosis data (30x5, 6 raters)."""
    return _load_nij(
        FLEISS_DIR / "Fleiss1971_psychiatric_data_nij.csv",
        index_col=0,
    )


@pytest.fixture(scope="module")
def gwet_nij():
    """NijMatrix from Gwet cont3x3 data (100x3, 2 raters)."""
    return _load_nij(
        GWET_DIR / "gwet_cont3x3_data_nij.csv",
        index_col="Subject",
    )


@pytest.fixture(scope="module")
def unbalanced_nij():
    """NijMatrix from unbalanced synthetic data (variable raters)."""
    return _load_nij(
        UNBALANCED_DIR / "unbalanced_data_nij.txt",
    )


# -- Ratings fixtures ------------------------------------------------------- #

@pytest.fixture(scope="module")
def gwet_ratings():
    """Ratings from Gwet cont3x3 data (100 subjects, 2 raters, K=3)."""
    return _load_ratings(
        GWET_DIR / "gwet_cont3x3_data.csv",
        GWET_DIR / "gwet_cont3x3_semantics.yaml",
        GWET_DIR / "gwet_cont3x3_cleaner_config.yaml",
    )


@pytest.fixture(scope="module")
def synthetic_ratings():
    """Ratings from synthetic 2-rater data (100 subjects, 2 raters, K=3)."""
    return _load_ratings(
        SYNTHETIC_DIR / "synthetic_2raters_data.txt",
        SYNTHETIC_DIR / "synthetic_2raters_semantics.yaml",
        SYNTHETIC_DIR / "synthetic_2raters_cleaner_config.yaml",
    )


@pytest.fixture(scope="module")
def sr1_prelim_ratings():
    """Ratings from SR1 preliminary data (199 subjects, 4 raters, K=3, balanced)."""
    return _load_ratings(
        SR1_DIR / "SR1_prelim_data.csv",
        SR1_DIR / "SR1_prelim_semantics.yaml",
        SR1_DIR / "SR1_prelim_cleaner_config.yaml",
    )


@pytest.fixture(scope="module")
def sr1_main_ratings():
    """Ratings from SR1 main data (1617 subjects, 4 raters, K=3, unbalanced)."""
    return _load_ratings(
        SR1_DIR / "SR1_main_data.txt",
        SR1_DIR / "SR1_main_semantics.yaml",
        SR1_DIR / "SR1_main_cleaner_config.yaml",
    )


@pytest.fixture(scope="module")
def sr2_ratings():
    """Ratings from SR2 data (81 subjects, 4 raters, K=3, nearly balanced)."""
    return _load_ratings(
        SR2_DIR / "SR2_data.txt",
        SR2_DIR / "SR2_semantics.yaml",
        SR2_DIR / "SR2_cleaner_config.yaml",
    )


@pytest.fixture(scope="module")
def sr3_ratings():
    """Ratings from SR3 FIS Psych data (770 subjects, 3 raters, K=2, balanced)."""
    return _load_ratings(
        SR3_DIR / "SR_FIS_Psych_data.txt",
        SR3_DIR / "SR_FIS_Psych_semantics.yaml",
        SR3_DIR / "SR_FIS_Psych_cleaner_config.yaml",
    )


@pytest.fixture(scope="module")
def sr4_ratings():
    """Ratings from SR4 FIS CVD data (787 subjects, 4 raters, K=2, balanced)."""
    return _load_ratings(
        SR4_DIR / "SR_FIS_CVD_data.txt",
        SR4_DIR / "SR_FIS_CVD_semantics.yaml",
        SR4_DIR / "SR_FIS_CVD_cleaner_config.yaml",
    )


@pytest.fixture(scope="module")
def sr5_ratings():
    """Ratings from SR5 Cocaine Cardio data (881 subjects, 3 raters, K=2, balanced)."""
    return _load_ratings(
        SR5_DIR / "SR_cocaine_cardio_data.txt",
        SR5_DIR / "SR_cocaine_cardio_semantics.yaml",
        SR5_DIR / "SR_cocaine_cardio_cleaner_config.yaml",
    )


# -- Config fixtures -------------------------------------------------------- #

@pytest.fixture(scope="session")
def fast_fleiss_config():
    """FleissKappaConfig for integration CI checks."""
    return _integration_config(FleissKappaConfig)


@pytest.fixture(scope="session")
def fast_gwet_config():
    """GwetAC1Config for integration CI checks."""
    return _integration_config(GwetAC1Config)


@pytest.fixture(scope="session")
def fast_bp_config():
    """BrennanPredigerConfig for integration CI checks."""
    return _integration_config(BrennanPredigerConfig)


@pytest.fixture(scope="session")
def fast_krippendorff_config():
    """KrippendorffAlphaConfig for integration CI checks."""
    return _integration_config(KrippendorffAlphaConfig)


@pytest.fixture(scope="session")
def fast_cohen_config():
    """CohenKappaConfig for integration CI checks."""
    return _integration_config(CohenKappaConfig)


@pytest.fixture(scope="session")
def fast_scotts_config():
    """ScottsConfig for integration CI checks."""
    return _integration_config(ScottsConfig)


@pytest.fixture(scope="session")
def fast_conger_config():
    """CongerKappaConfig for integration CI checks."""
    return _integration_config(CongerKappaConfig)


@pytest.fixture(scope="session")
def fast_light_config():
    """LightKappaConfig for integration CI checks."""
    return _integration_config(LightKappaConfig)


@pytest.fixture(scope="session")
def subset_config():
    """Build a SubsetPermutationConfig for integration testing.

    Returns a callable taking the subset size k, since the tests that use
    it compare behaviour across k rather than fixing one.

    Uses Fleiss' kappa with balanced_data_only=False, because the dataset
    this exists for is SR1 main: partially crossed, which is the only
    shipped data where a subset's measurable pool and its complete cases
    differ. On balanced data the two coincide and the test proves nothing.
    """
    def _make(k: int):
        return SubsetPermutationConfig(
            coefficient_to_compute="fleiss_kappa",
            config_of_coefficient=_integration_config(
                FleissKappaConfig, balanced_data_only=False,
            ),
            subset_k=k,
        )

    return _make


@pytest.fixture(scope="session")
def pairwise_config():
    """Build a PairwisePermutationConfig for integration testing.

    Returns a callable, because the mean-bootstrap method is what varies
    across these tests and one fixture per method would be three
    near-identical fixtures.

    Both resample counts -- the per-pair bootstrap and the mean bootstrap --
    use _INTEGRATION_RESAMPLES, the same number every other config here
    uses. The per-pair intervals are the reason that matters: a pairwise
    run builds one interval per rater pair, so an under-resampled setting
    is not one noisy interval but C(n,2) of them.

    Usage:
        config = pairwise_config()               # no mean bootstrap
        config = pairwise_config("correlated")
    """
    def _make(mean_method: str | None = None):
        kwargs = dict(
            coefficient_to_compute="cohen_kappa",
            config_of_coefficient=CohenKappaConfig(
                n_resamples=_INTEGRATION_RESAMPLES,
                random_state=42,
            ),
        )
        if mean_method is not None:
            kwargs["bootstrap_mean_config"] = PermutationBootstrapMeanConfig(method=mean_method)
        return PairwisePermutationConfig(**kwargs)

    return _make


_ORDERED_DATASETS = {
    "sr1_prelim": (
        SR1_DIR / "SR1_prelim_data.csv",
        SR1_DIR / "SR1_prelim_semantics.yaml",
        SR1_DIR / "SR1_prelim_cleaner_config.yaml",
    ),
    "gwet": (
        GWET_DIR / "gwet_cont3x3_data.csv",
        GWET_DIR / "gwet_cont3x3_semantics.yaml",
        GWET_DIR / "gwet_cont3x3_cleaner_config.yaml",
    ),
    "synthetic": (
        SYNTHETIC_DIR / "synthetic_2raters_data.txt",
        SYNTHETIC_DIR / "synthetic_2raters_semantics.yaml",
        SYNTHETIC_DIR / "synthetic_2raters_cleaner_config.yaml",
    ),
}


@pytest.fixture(scope="module")
def ordered_ratings():
    """Load a dataset with its categories declared ordered.

    Returns a callable rather than ratings, because the weighted tests need
    several (dataset, weighting) combinations. Weights are declared in the
    semantics, so the weighted analysis runs the ordinary config on ratings
    whose categories are ordered; the unordered fixtures (sr1_prelim_ratings,
    gwet_ratings, synthetic_ratings) give the unweighted numbers.

    Usage:
        ratings = ordered_ratings("sr1_prelim", "quadratic")
    """
    cache: dict = {}

    def _make(dataset: str, weighting: str):
        key = (dataset, weighting)
        if key not in cache:
            cache[key] = _load_ratings(*_ORDERED_DATASETS[dataset], weighting=weighting)
        return cache[key]

    return _make


@pytest.fixture(scope="session")
def fleiss_unbalanced_config():
    """FleissKappaConfig that accepts unbalanced data.

    FleissKappaConfig defaults to balanced_data_only=True and rejects data
    where subjects were rated different numbers of times. SR1 main is such
    a dataset, and the approximation path needs exercising. Everything else
    matches _integration_config, so the resample count and method stay in
    one place rather than being written out again at the call site.
    """
    return _integration_config(FleissKappaConfig, balanced_data_only=False)


@pytest.fixture(scope="session")
def bca_fleiss_config():
    """FleissKappaConfig running real BCa at B=1000.

    The one config in this file that exercises the BCa path end-to-end.
    Paired with a small fixture on purpose -- see _bca_config.
    """
    return _bca_config(FleissKappaConfig)
