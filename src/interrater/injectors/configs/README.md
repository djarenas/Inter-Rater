# Configs

Frozen dataclasses that configure rater injection strategies. Each config
pairs 1:1 with an injector class and validates its own fields at
construction time. The injector reads the config but never writes it.

## Inheritance

`InjectorConfig` is the base. It holds the two fields every injector
needs (`n_new_raters`, `random_state`) and provides the shared API
(`build_from_dict`, `build_from_yaml`, `summary`). Five leaf configs
inherit it and add strategy-specific fields.

## Config classes

| Config                         | Injector                  | Added field(s)                     |
|--------------------------------|---------------------------|------------------------------------|
| `InjectorConfig`               | (base -- not used alone)  | `n_new_raters`, `random_state`     |
| `DistributionInjectorConfig`   | `DistributionInjector`    | `prevalences`                      |
| `WithinSubjectInjectorConfig`  | `WithinSubjectInjector`   | `noise_rate`                       |
| `PooledInjectorConfig`         | `PooledInjector`          | `smoothing`, `exclude_raters`      |
| `RaterBootstrapInjectorConfig` | `RaterBootstrapInjector`  | `source_raters`                    |
| `ExtremalInjectorConfig`       | `ExtremalInjector`        | `strategy`                         |

## Shared API

Every config inherits three class/instance methods from `InjectorConfig`:

- **`build_from_dict(data_dict)`** -- class method. Builds a config from
  a plain dictionary. Leaf configs override to hydrate YAML types
  (list -> tuple) before delegating to `super()`.
- **`build_from_yaml(path)`** -- class method. Loads a YAML file,
  validates that all keys are recognized fields, then delegates to
  `build_from_dict`.
- **`summary()`** -- instance method. Returns the config as a plain
  dictionary for audit trails. Leaf configs call `super().summary()`
  and extend with their own fields. Includes `_package_version` for
  reproducibility.

## Validation split

Configs validate what they can validate alone -- no data needed:

- Probabilities are non-negative and sum to 1.0
- Indices are non-negative integers, no duplicates
- Rates are in [0.0, 1.0]
- Strategy strings are from the valid set
- Bools are rejected where ints or floats are expected

Injectors validate what requires the `Ratings` data:

- Rater indices are in bounds for this dataset
- At least one rater remains after exclusion
- The pool is not empty after filtering
- Prevalences length matches n_categories

This split keeps configs self-contained and testable without real data.

## Dependencies

- None outside the package. Configs import only from each other
  (leaf -> base) and the standard library.
- `build_from_yaml` imports `yaml` inside the method body to keep the
  top-level import clean.
- `summary` imports `__version__` from the package root inside the
  method body (three levels up) for the audit trail.
