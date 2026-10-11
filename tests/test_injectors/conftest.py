"""
Shared fixtures for injector tests.

Two groups: configs (standalone, no data needed) and ratings (for injector
tests that need real data to inject into).

Config fixtures:
- One valid instance per config type, with non-default values where possible.

Ratings fixtures:
- SIMPLE_RATINGS: 4 subjects, 3 raters, 2 categories. Fully crossed (no
  EMPTY cells). Small enough that pool counts, majority/minority, and
  within-subject copying can be verified by inspection.
- THREE_CAT_RATINGS: 4 subjects, 3 raters, 3 categories. Used for
  contrarian testing (unused categories at some subjects).

All are module-level constants so no test can mutate what another sees.
Fixtures return fresh copies where mutation matters.
"""

import numpy as np
import pytest

from interrater.contracts import INDEX_DTYPE
from interrater.semantics.categories import UnorderedCategories
from interrater.semantics.ratings_semantics import RatingsSemantics
from interrater.data_objects.ratings import Ratings
from interrater.injectors.configs import (
    InjectorConfig,
    DistributionInjectorConfig,
    WithinSubjectInjectorConfig,
    PooledInjectorConfig,
    RaterBootstrapInjectorConfig,
    ExtremalInjectorConfig,
)


# -- Config instances -------------------------------------------------------- #

BASE_CONFIG = InjectorConfig(n_new_raters=2, random_state=42)

DISTRIBUTION_CONFIG = DistributionInjectorConfig(
    prevalences=(0.6, 0.4), n_new_raters=2, random_state=42,
)

DISTRIBUTION_CONFIG_UNIFORM = DistributionInjectorConfig(
    n_new_raters=1, random_state=42,
)

WITHIN_SUBJECT_CONFIG = WithinSubjectInjectorConfig(
    noise_rate=0.0, n_new_raters=2, random_state=42,
)

WITHIN_SUBJECT_CONFIG_NOISY = WithinSubjectInjectorConfig(
    noise_rate=0.5, n_new_raters=2, random_state=42,
)

POOLED_CONFIG = PooledInjectorConfig(
    smoothing=0.0, n_new_raters=2, random_state=42,
)

POOLED_CONFIG_SMOOTH = PooledInjectorConfig(
    smoothing=0.5, n_new_raters=2, random_state=42,
)

POOLED_CONFIG_EXCLUDE = PooledInjectorConfig(
    exclude_raters=(2,), n_new_raters=1, random_state=42,
)

BOOTSTRAP_CONFIG = RaterBootstrapInjectorConfig(
    source_raters=(0, 1), n_new_raters=2, random_state=42,
)

EXTREMAL_CONFIG_MAJORITY = ExtremalInjectorConfig(
    strategy="majority", n_new_raters=1, random_state=42,
)

EXTREMAL_CONFIG_MINORITY = ExtremalInjectorConfig(
    strategy="minority", n_new_raters=1, random_state=42,
)

EXTREMAL_CONFIG_CONTRARIAN = ExtremalInjectorConfig(
    strategy="contrarian", n_new_raters=1, random_state=42,
)

EXTREMAL_CONFIG_MEDIAN = ExtremalInjectorConfig(
    strategy="median", n_new_raters=1, random_state=42,
)


# -- Ratings datasets -------------------------------------------------------- #

# SIMPLE: 4 subjects, 3 raters, 2 categories (yes=0, no=1).
# Fully crossed, no EMPTYs.
#
#   subject | alice | bob | carol | majority | minority
#   --------+-------+-----+-------+----------+---------
#   s0      |   0   |  0  |   1   |    0     |    1
#   s1      |   1   |  1  |   0   |    1     |    0
#   s2      |   0   |  0  |   0   |    0     |    0 (unanimous)
#   s3      |   1   |  0  |   1   |    1     |    0
#
# Pool (all 12 cells): 0 appears 6 times, 1 appears 6 times.

TWO_CAT = UnorderedCategories(("yes", "no"))

SIMPLE_SEMANTICS = RatingsSemantics(
    raters=["alice", "bob", "carol"],
    categories=TWO_CAT,
)

SIMPLE_INDEX = np.array(
    [[0, 0, 1],
     [1, 1, 0],
     [0, 0, 0],
     [1, 0, 1]],
    dtype=INDEX_DTYPE,
)

SIMPLE_SUBJECT_IDS = ("s0", "s1", "s2", "s3")


# THREE_CAT: 4 subjects, 3 raters, 3 categories (a=0, b=1, c=2).
# Designed so that at some subjects not all categories are used,
# giving contrarian something to pick.
#
#   subject | alice | bob | carol | used      | unused
#   --------+-------+-----+-------+-----------+-------
#   s0      |   0   |  0  |   1   | {0, 1}    | {2}
#   s1      |   1   |  1  |   1   | {1}       | {0, 2}
#   s2      |   0   |  1  |   2   | {0, 1, 2} | {}
#   s3      |   2   |  2  |   0   | {0, 2}    | {1}

THREE_CAT = UnorderedCategories(("a", "b", "c"))

THREE_CAT_SEMANTICS = RatingsSemantics(
    raters=["alice", "bob", "carol"],
    categories=THREE_CAT,
)

THREE_CAT_INDEX = np.array(
    [[0, 0, 1],
     [1, 1, 1],
     [0, 1, 2],
     [2, 2, 0]],
    dtype=INDEX_DTYPE,
)


# -- Fixtures ---------------------------------------------------------------- #


@pytest.fixture
def base_config():
    return BASE_CONFIG


@pytest.fixture
def distribution_config():
    return DISTRIBUTION_CONFIG


@pytest.fixture
def distribution_config_uniform():
    return DISTRIBUTION_CONFIG_UNIFORM


@pytest.fixture
def within_subject_config():
    return WITHIN_SUBJECT_CONFIG


@pytest.fixture
def within_subject_config_noisy():
    return WITHIN_SUBJECT_CONFIG_NOISY


@pytest.fixture
def pooled_config():
    return POOLED_CONFIG


@pytest.fixture
def pooled_config_smooth():
    return POOLED_CONFIG_SMOOTH


@pytest.fixture
def pooled_config_exclude():
    return POOLED_CONFIG_EXCLUDE


@pytest.fixture
def bootstrap_config():
    return BOOTSTRAP_CONFIG


@pytest.fixture
def extremal_config_majority():
    return EXTREMAL_CONFIG_MAJORITY


@pytest.fixture
def extremal_config_minority():
    return EXTREMAL_CONFIG_MINORITY


@pytest.fixture
def extremal_config_contrarian():
    return EXTREMAL_CONFIG_CONTRARIAN


@pytest.fixture
def extremal_config_median():
    return EXTREMAL_CONFIG_MEDIAN


@pytest.fixture
def simple_ratings():
    """4 subjects, 3 raters, 2 categories, fully crossed."""
    return Ratings(
        SIMPLE_INDEX.copy(), SIMPLE_SEMANTICS,
        subject_ids=SIMPLE_SUBJECT_IDS,
    )


@pytest.fixture
def three_cat_ratings():
    """4 subjects, 3 raters, 3 categories, some unused per subject."""
    return Ratings(
        THREE_CAT_INDEX.copy(), THREE_CAT_SEMANTICS,
        subject_ids=SIMPLE_SUBJECT_IDS,
    )
