"""YAML key normalization for configuration objects.

I hold the canonical-key alias registry and a function to normalize
user-facing YAML keys to their canonical field names. Every config's
``build_from_dict`` calls ``normalize_yaml_key`` before constructing
the dataclass, so users can write the short form in YAML and the code
always sees the canonical field name.

Adding a new alias is a one-line edit to YAML_KEY_ALIASES.
"""

from __future__ import annotations

# ------------------------------------------------- YAML Key Aliases #

YAML_KEY_ALIASES: dict[str, list[str]] = {
    "hypothesis_test_config": ["hypothesis_test"],
    "bootstrap_mean_config": ["bootstrap_mean"],
}
"""Map canonical field name to accepted alternate spellings.

The canonical key is always accepted. The aliases are the extra
spellings a user may write in YAML or pass in a dict. If multiple
keys from the same alias group are present, normalize_yaml_key
raises ValueError.
"""

_ERROR_PREFIX = "[yaml_utils] Error"


# ---------------------------------------------------- Normalization #

def normalize_yaml_key(data_dict: dict, canonical_key: str) -> None:
    """
    Resolve an alias to the canonical key in place.

    I look up the aliases for ``canonical_key`` in YAML_KEY_ALIASES.
    If an alias is present in ``data_dict``, I rename it to the
    canonical key. If both an alias and the canonical key are present,
    I raise ValueError.

    Parameters
    ----------
    data_dict : dict
        The raw configuration dictionary. Mutated in place.
    canonical_key : str
        The canonical field name (must be a key in YAML_KEY_ALIASES).

    Raises
    ------
    KeyError
        If ``canonical_key`` is not registered in YAML_KEY_ALIASES.
    ValueError
        If multiple keys from the same alias group are present.
    """
    if canonical_key not in YAML_KEY_ALIASES:
        raise KeyError(
            f"{_ERROR_PREFIX}: '{canonical_key}' is not registered "
            f"in YAML_KEY_ALIASES"
        )

    aliases = YAML_KEY_ALIASES[canonical_key]
    all_keys = [canonical_key] + aliases
    present = [k for k in all_keys if k in data_dict]

    if len(present) > 1:
        raise ValueError(
            f"{_ERROR_PREFIX}: multiple keys for the same field: "
            f"{present}. Pass one of {all_keys}, not multiple."
        )

    if len(present) == 1 and present[0] != canonical_key:
        data_dict[canonical_key] = data_dict.pop(present[0])
