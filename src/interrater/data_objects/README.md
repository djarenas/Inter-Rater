# Data Objects

The validated arrays that hold ratings data. If a data object exists, it
is complete -- its constructor raises rather than repairs.

## The Hierarchy

| Class             | Guarantee                                        | What It Buys You                       |
|-------------------|--------------------------------------------------|----------------------------------------|
| `Ratings`         | Valid values, shape matches semantics, no idle raters, every subject rated at least twice | Accepts any rating layout              |
| `BalancedRatings` | + every subject rated by the same *number* of raters | Exact Fleiss' kappa (no approximation) |

`BalancedRatings` narrows the contract. Fleiss' kappa with
`balanced_data_only: true` converts its input to `BalancedRatings`, so
unbalanced data is refused when that object is built, with the reason,
before anything is computed.

## What Ratings Holds

Three things, and nothing else:

- **category_index** -- `(n_subjects, n_raters)` int8 array. Each cell is
  a position in `semantics.categories.names`, or `EMPTY` (-1). Read-only
  after construction.
- **subject_ids** -- `(n_subjects,)` object array. One label per row,
  generated as `subject0, subject1, ...` when the data had none.
- **semantics** -- a `RatingsSemantics` declaring who the raters are, what
  the categories mean, and what the data was expected to look like.

Ratings stores nothing about itself. How skewed the categories are,
whether every subject was rated the same number of times -- those are
facts *about* it, computed when asked (`check_my_balance()`,
`to_nij_matrix()`), not stored. Keeping them out means it stays cheap
to build, which matters when a bootstrapper is making thousands of
arrays a second.

## Validation

All validation lives in `_ratings_validation.py`, separate from the
classes. Each check takes the concrete class name so that a
`BalancedRatings` failing a base check says `BalancedRatings`, not `Ratings`.

The validation chain:
1. dtype (must be integer, not float -- NaN means EMPTY was not used)
2. shape (2-D, at least MIN_SUBJECTS rows)
3. rater count matches semantics
4. subject count matches semantics (if declared)
5. values in range for the declared categories (also checked before the
   cast to int8, which would wrap 256 to 0)
6. subject ids match row count
7. subject ids unique (unless semantics say otherwise)
8. no idle raters (every declared rater rated at least one subject)
9. every subject rated at least twice (MIN_RATINGS_PER_SUBJECT non-EMPTY cells per row)

BalancedRatings adds: constant raters per subject.

## Dependencies

- `semantics` (Categories, RatingsSemantics)
- `contracts.py` (EMPTY, INDEX_DTYPE, MIN_SUBJECTS, MIN_RATINGS_PER_SUBJECT)
- `array_utils.py` (balance and crossing checks)
