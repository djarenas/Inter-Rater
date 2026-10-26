# Semantics

What a rating task looks like before the data is seen. Semantics declare
the raters, the categories, and the structural expectations -- everything
a `Ratings` object validates itself against.

## Categories

Categories own K -- the number of responses a rater was allowed to give.
K is between 2 and 127: ratings store each category's position as int8,
so a larger scale is refused when it is declared.
Three subclasses handle different measurement scales:

| Subclass              | Scale     | Distances                         |
|-----------------------|-----------|-----------------------------------|
| `UnorderedCategories` | Nominal   | None -- every disagreement is equal |
| `OrderedCategories`   | Ordinal   | Derived from position (linear or quadratic) |
| `WeightedCategories`  | Custom    | Supplied as an explicit weight matrix |

The base class (`Categories`) is abstract. `DistanceCategories` is the
intermediate ABC shared by `OrderedCategories` and `WeightedCategories`
-- both produce a weight matrix, just from different sources.

`build_categories()` dispatches a plain dict to the right subclass,
which is the path YAML declarations take.

## RatingsSemantics

Sits between Categories and Ratings in the dependency chain. Holds:

- **raters** -- names (list of strings) or a count (integer, which
  generates `rater0`, `rater1`, ...).
- **categories** -- a `Categories` instance.
- **n_subjects** -- optional expected subject count.
- **unique_subjects** -- whether subject IDs must be unique (default True).

Type guards catch the traps Python's type system leaves open: a string
is a `Sequence`, a bool is an `int`. If a semantics can be built, a
`Ratings` that matches it can be built too.

## Dependencies

None beyond `contracts.py`. This is the bottom of the dependency chain.
