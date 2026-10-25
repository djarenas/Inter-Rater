# Cleaners

The only way in: a messy pandas DataFrame becomes a validated `Ratings`,
with a receipt for everything that changed. The receipt is the
`CleaningReport`.

## Configuration layer

Every decision the cleaner makes is declared before the data is seen.
Three objects hold those decisions, bottom-up:

| Class                    | What it declares                                        |
|--------------------------|---------------------------------------------------------|
| `FrameLayout`            | Whether the frame has column headers and which column (if any) holds subject IDs. |
| `InsufficientDataPolicy` | Minimum ratings per rater and per subject, and the sweep order for dropping thin data (`raters_first`, `subjects_first`, or `iterate`). |
| `CleanerConfig`          | Everything else -- case sensitivity, whitespace stripping, missing-value markers, value repairs, policies for unknown responses and unused categories, plus the two delegates above. |

`CleanerConfig.from_dict` and `CleanerConfig.build_from_yaml` build a
config from a plain dictionary or a YAML file, so the entire cleaning
specification can live in a versioned file alongside the analysis.

## RatingsCleaner

The public entry point. Takes a DataFrame, a `RatingsSemantics`, and a
`CleanerConfig`, and produces two things: a `Ratings` and a
`CleaningReport`.

The pipeline runs in four stages:

1. **Config-vs-semantics check** -- repairs pointing at nonexistent
   categories, missing markers that collide with valid responses, and
   other contradictions between the two documents are caught before a
   single row is touched.
2. **String cleaning** -- normalisation (case folding, whitespace
   stripping), missing-value detection, and value-map repairs, all in
   text space.
3. **Encoding** -- every cell that survived cleaning is mapped to its
   integer category code; everything else becomes `EMPTY`.
4. **Dropping** -- raters and subjects below the minimum rating
   thresholds are swept out, one or more passes depending on the
   declared method.

All repair is accounted for. Everything that changes is counted,
attributed to the rule that caused it, and broken down per rater.

## Private modules

| Module            | What it does                                            |
|-------------------|---------------------------------------------------------|
| `_string_utils.py`| `fold`, `fold_all`, `fold_collisions`, `describe_folding`, `build_lookup` -- string normalisation, finding different strings that fold alike, and lookup-table construction. Free functions over plain strings; the config and the cleaner call the same code so matching is guaranteed consistent. |
| `_drop_utils.py`  | `rating_counts`, `sweep_once`, `sweep_until_settled`, `distribution` -- row/column dropping over the code array. Free functions over numpy; the cleaner supplies thresholds and order, gets back boolean masks. |

## Data flow

```
DataFrame (messy strings)
    │
    ├──► CleanerConfig ◄── FrameLayout
    │         │             InsufficientDataPolicy
    │         ▼
    ├──► config-vs-semantics check (_string_utils.fold_collisions, build_lookup)
    │         │
    │         ▼
    ├──► string cleaning (_string_utils.fold, fold_all)
    │         │
    │         ▼
    ├──► encoding (category names → int8 codes)
    │         │
    │         ▼
    └──► dropping (_drop_utils.sweep_once / sweep_until_settled)
              │
              ▼
        Ratings + CleaningReport
```

## Dependencies

- `semantics` -- `RatingsSemantics` and `Categories` for declaring the task
- `data_objects` -- `Ratings` as the output
- `reports` -- `CleaningReport`, `MatchingReport`, `CategoryReport`, `DropReport`
- `contracts` -- `EMPTY`, `INDEX_DTYPE`, `LABEL_DTYPE`
