# Reports

What a run of the package did, written for a person to read. Reports are
terminal -- nothing downstream consumes one. They are a presentation layer,
a sibling of `visualizers` in the medium of prose.

Every report holds everything it computed and offers two views of it:
`summary()` in plain types for writing to disk, and `render(verbosity)`
built from that summary so the two cannot drift.

## Report classes

| Class             | What it records                                          |
|-------------------|----------------------------------------------------------|
| `MatchingReport`  | What happened to the strings: nulls, missing markers, repairs, unknowns. Per rater and per subject. |
| `DropReport`      | What was removed for thinness: which raters, which subjects, in how many passes, and whether the loop converged. |
| `DropPass`        | One sweep of the drop loop -- how many raters and subjects fell in that pass. |
| `CategoryReport`  | Which declared categories were used, which were unused, which are rare (used but too seldom to estimate). |
| `CleaningReport`  | The composite: all three sub-reports, plus the declared and effective semantics, shapes in and out, and the config that produced it. |

## CleaningReport

The top-level receipt. It holds:

- **declared_semantics** -- the rating task as written, before any row was read.
- **effective_semantics** -- the rating task the surviving data actually
  describes. Identical to the declared one unless a rater or a category
  was dropped.
- **config** -- the `CleanerConfig` that was applied (imported under
  `TYPE_CHECKING` to avoid a circular dependency with `cleaners/`).
- **matching, drops, categories** -- the three sub-reports.

Key properties: `semantics_changed`, `k_changed`, `survived`,
`n_subjects_lost`, `n_raters_lost`.

`summary()` returns the full record in plain types, including the cleaner
config that produced it (`cleaner_config`, in the shape
`CleanerConfig.from_dict` reads) and the package version
(`_package_version`). `to_yaml(path)` writes it to disk, so the record alone
says how the data was cleaned and by which release.

`render(verbosity)` validates verbosity (`"summary"` or `"full"`) and
delegates to each sub-report's `render()`. When semantics changed, an
extra section names what moved and why it matters for the arithmetic.

## Dependencies

- `semantics` -- `RatingsSemantics`, `DistanceCategories` (for the
  weight-shift warning when K changes on an ordinal scale)
- `cleaners` -- `CleanerConfig` (annotation only, under `TYPE_CHECKING`)
