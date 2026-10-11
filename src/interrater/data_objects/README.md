# Data Objects

The validated arrays that hold ratings data. If a data object exists, it
is complete -- its constructor raises rather than repairs.

## The Hierarchy

| Class             | Guarantee                                        | What It Buys You                       |
|-------------------|--------------------------------------------------|----------------------------------------|
| `Ratings`         | Valid values, shape matches semantics, no idle raters | Accepts any rating layout              |
| `BalancedRatings` | + every subject rated by the same *number* of raters | Exact Fleiss' kappa (no approximation) |
| `PairedRatings`   | + exactly two raters, no gaps                    | Cohen's kappa, individual marginals    |

Each subclass narrows the contract. An analyzer that takes `PairedRatings`
cannot receive three raters by accident -- the precondition lives in the
type, not a runtime check buried in the computation.

## What Ratings Holds

Three things, and nothing else:

- **category_index** -- `(n_subjects, n_raters)` int8 array. Each cell is
  a position in `semantics.categories.names`, or `EMPTY` (-1). Read-only
  after construction.
- **subject_ids** -- `(n_subjects,)` object array. One label per row,
  generated as `subject0, subject1, ...` when the data had none.
- **semantics** -- a `RatingsSemantics` declaring who the raters are, what
  the categories mean, and what the data was expected to look like.

Ratings computes nothing about itself. How skewed the categories are,
whether every subject was rated the same number of times -- those are
facts *about* it, not part of it. Keeping them out means it stays cheap
to build, which matters when a bootstrapper is making thousands of
arrays a second.

## Validation

All validation lives in `_ratings_validation.py`, separate from the
classes. Each check takes the concrete class name so that a
`PairedRatings` failing a base check says `PairedRatings`, not `Ratings`.

The validation chain:
1. dtype (must be integer, not float -- NaN means EMPTY was not used)
2. shape (2-D, at least MIN_SUBJECTS rows)
3. rater count matches semantics
4. subject count matches semantics (if declared)
5. values in range for the declared categories
6. subject ids match row count
7. subject ids unique (unless semantics say otherwise)
8. no idle raters (every declared rater rated at least one subject)
9. no unrated subjects (every row has at least one non-EMPTY cell)

BalancedRatings adds: constant raters per subject.
PairedRatings adds: exactly two raters, fully crossed (no gaps).

## Dependencies

- `semantics` (Categories, RatingsSemantics)
- `contracts.py` (EMPTY, INDEX_DTYPE, MIN_SUBJECTS)
- `array_utils.py` (balance and crossing checks)
