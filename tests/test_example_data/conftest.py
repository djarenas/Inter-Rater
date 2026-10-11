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


# -- Helpers ---------------------------------------------------------------- #

def _load_ratings(data_path: Path, semantics_path: Path,
                  cleaner_config_path: Path) -> "Ratings":
    """Load a ratings dataset via the RatingsCleaner pipeline."""
    semantics = RatingsSemantics.build_from_yaml(str(semantics_path))
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


def _fast_config(config_cls, **overrides):
    """
    Build a fast analyzer config for integration testing.

    Uses n_resamples=50 and random_state=42 for CI structural checks.
    Point estimates are deterministic and do not depend on these values.
    """
    defaults = dict(
        n_resamples=50,
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
    """FleissKappaConfig with n_resamples=50 for fast CI checks."""
    return _fast_config(FleissKappaConfig)


@pytest.fixture(scope="session")
def fast_gwet_config():
    """GwetAC1Config with n_resamples=50 for fast CI checks."""
    return _fast_config(GwetAC1Config)


@pytest.fixture(scope="session")
def fast_bp_config():
    """BrennanPredigerConfig with n_resamples=50 for fast CI checks."""
    return _fast_config(BrennanPredigerConfig)


@pytest.fixture(scope="session")
def fast_krippendorff_config():
    """KrippendorffAlphaConfig with n_resamples=50 for fast CI checks."""
    return _fast_config(KrippendorffAlphaConfig)


@pytest.fixture(scope="session")
def fast_cohen_config():
    """CohenKappaConfig with n_resamples=50 for fast CI checks."""
    return _fast_config(CohenKappaConfig)


@pytest.fixture(scope="session")
def fast_scotts_config():
    """ScottsConfig with n_resamples=50 for fast CI checks."""
    return _fast_config(ScottsConfig)


@pytest.fixture(scope="session")
def fast_conger_config():
    """CongerKappaConfig with n_resamples=50 for fast CI checks."""
    return _fast_config(CongerKappaConfig)


@pytest.fixture(scope="session")
def fast_light_config():
    """LightKappaConfig with n_resamples=50 for fast CI checks."""
    return _fast_config(LightKappaConfig)
