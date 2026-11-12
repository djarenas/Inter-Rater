"""
A misspelt key is refused the same way everywhere.

Every ``build_from_dict`` in the package, and ``build_categories``, sends its
keys through ``yaml_utils.refuse_unknown_keys``. So a typo gets a ValueError
that names the class, the nearest valid key and every valid key, whichever
file it was in.
"""

import pytest

from interrater import (
    CleanerConfig,
    FleissKappaConfig,
    RatingsSemantics,
    SubsetPermutationConfig,
)
from interrater.analyzers.configs.permutation_bootstrap_mean_config import (
    PermutationBootstrapMeanConfig,
)
from interrater.cleaners.frame_layout import FrameLayout
from interrater.cleaners.insufficient_data_policy import InsufficientDataPolicy
from interrater.injectors.configs import ExtremalInjectorConfig
from interrater.semantics.categories import build_categories

_SEMANTICS = {"raters": ["a", "b"], "n_subjects": 10,
              "categories": {"kind": "unordered", "names": ["x", "y"]}}

# (owner in the message, call, misspelt key, key it is nearest to)
CASES = [
    ("FleissKappaConfig",
     lambda: FleissKappaConfig.build_from_dict({"n_resample": 100}),
     "n_resample", "n_resamples"),
    ("SubsetPermutationConfig",
     lambda: SubsetPermutationConfig.build_from_dict({"subset_kk": 2}),
     "subset_kk", "subset_k"),
    ("HypothesisTestConfig",
     lambda: FleissKappaConfig.build_from_dict({"hypothesis_test": {"n_permutation": 100}}),
     "n_permutation", "n_permutations"),
    ("PermutationBootstrapMeanConfig",
     lambda: PermutationBootstrapMeanConfig.build_from_dict({"metod": "correlated"}),
     "metod", "method"),
    ("ExtremalInjectorConfig",
     lambda: ExtremalInjectorConfig.build_from_dict({"stratgy": "majority"}),
     "stratgy", "strategy"),
    ("RatingsSemantics",
     lambda: RatingsSemantics.build_from_dict({**_SEMANTICS, "rater": ["a", "b"]}),
     "rater", "raters"),
    ("UnorderedCategories",
     lambda: build_categories({"kind": "unordered", "names": ["x", "y"], "label": ["X", "Y"]}),
     "label", "labels"),
    ("CleanerConfig",
     lambda: CleanerConfig.build_from_dict({"missing_value": ["-"]}),
     "missing_value", "missing_values"),
    ("CleanerConfig",
     lambda: CleanerConfig.build_from_dict({"frame_layuot": {}}),
     "frame_layuot", "frame_layout"),
    ("FrameLayout",
     lambda: FrameLayout.build_from_dict({"subject_id_column": "id"}),
     "subject_id_column", "subject_id_col"),
    ("InsufficientDataPolicy",
     lambda: InsufficientDataPolicy.build_from_dict({"min_rating_per_rater": 2}),
     "min_rating_per_rater", "min_ratings_per_rater"),
]


@pytest.mark.parametrize("owner, call, typo, meant", CASES,
                         ids=[f"{case[0]}:{case[2]}" for case in CASES])
def test_misspelt_key_refused_the_same_way(owner, call, typo, meant):
    with pytest.raises(ValueError) as err:
        call()
    message = str(err.value)
    assert message.startswith(f"[{owner}] Error: unknown key(s) [{typo!r}]")
    assert f"{typo!r} -> did you mean {meant!r}?" in message
    assert "Valid keys are: [" in message
    assert f"{meant!r}" in message.split("Valid keys are:")[1]
