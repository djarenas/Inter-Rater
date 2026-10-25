# Injectors

Synthetic rater injection for sensitivity analysis. Add a fake rater to
your panel and see what happens to agreement.

Each injector takes a `Ratings` object and an `InjectorConfig`, generates
one or more synthetic rater columns, and returns a **new** `Ratings` with
those columns appended. The original data is never modified.

## The Spectrum

Injectors are ordered by how much real structure they preserve:

| Intuitive Name | Class                    | What It Simulates                                                                 |
|----------------|--------------------------|-----------------------------------------------------------------------------------|
| **Echo**       | `WithinSubjectInjector`  | Copies a real answer at each subject. Correlated with the panel, tunable noise.   |
| **Clone**      | `RaterBootstrapInjector` | Resamples from one or more real raters' columns. Copies biases, scrambles judgments.|
| **Chance**     | `PooledInjector`         | Draws from the global pool of all ratings. Right base rates, no subject signal.    |
| **Chance**     | `DistributionInjector`   | Draws from a fixed categorical distribution. Pure null model.                      |
| **Conformist** | `ExtremalInjector`       | Always picks the majority category at each subject. Goes with the crowd.           |
| **Dissenter**  | `ExtremalInjector`       | Picks the least common used category. Unpopular option someone did choose.         |
| **Contrarian** | `ExtremalInjector`       | Picks a category no one used. True contrarian. Falls back to dissenter if none.    |
| **Median**     | `ExtremalInjector`       | Picks the middle-prevalence category. Useful with many raters.                     |

The further down the table, the less the injected rater "knows" about
the real subjects -- until the extremal strategies, which know everything
but act deterministically.

## Contrarian vs Dissenter

The `"minority"` strategy picks the **least common category that was
actually used** at each subject -- the unpopular choice someone did
make. This is a dissenter who picks a side, just the least popular one.

The `"contrarian"` strategy picks a category that **no one** used at
that subject -- the zero-prevalence option. The true lowest prevalence
is zero, and only the contrarian goes there. When all categories were
used (no unused category available), the contrarian falls back to
minority behavior.

In short: the dissenter picks the road less traveled; the contrarian
picks the road no one took.

## Config and Injector

Every injector has a paired config:

- **Config** validates its own fields (probabilities sum to 1, indices
  are non-negative, strategy is a valid string).
- **Injector** validates compatibility with the `Ratings` at generate
  time (indices in bounds, enough raters to compute a mode, etc.).

Configs are frozen dataclasses. Injectors are stateless beyond their
`ratings` and `config` -- call `inject()` and get a new `Ratings` back.

## Chaining

Injections can be chained. Each call returns a new `Ratings`, and
injected raters are named `injected_rater0`, `injected_rater1`, etc.
When chaining, naming picks up where the previous injection left off,
so names never collide.

Seeds don't collide either. Each step's generator is seeded with the pair
`[random_state, step_position]`, the step's own seed and its place in the
pipeline (0, 1, 2, ...), which numpy mixes into one seed. Two steps given
the same `random_state` therefore produce different raters, and a whole
pipeline still reruns identically. Reordering steps changes the results,
since order is part of the pipeline. `random_state=None` means fresh
randomness on every run.

A later step draws from every rater present when it runs, including
those injected by earlier steps. A pooled or bootstrap step after a
distribution step learns partly from the distribution step's raters.
