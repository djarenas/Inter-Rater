## Why this is built out of objects, and what the flexibility buys

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

A coefficient is a formula, and a formula is a function. So it is fair to ask
what the types are for, and why we spent so much time on a OOP design.

The answer is not the coefficients. It is everything you might want to do
*around* them. Once you have a number, the interesting questions are about how
it would change:

| question | operation |
| --- | --- |
| How precise is this estimate? | resample subjects with replacement |
| What is its bias? | leave out one subject at a time |
| Was one rater pulling it? | leave out one rater at a time |
| Which pairs of raters disagree? | every pair in turn |
| Does any group of three hold together? | every subset of that size |
| Would other raters have agreed? | resample raters with replacement |
| What if the panel had one more rater? | inject a synthetic rater of a declared kind |

Every one of those crosses with every coefficient, and with whatever form your
data arrived in. Most of the combinations are fine. Some are not, and the ones
that are not fail *quietly*. This package throws an error on configurations that 
make unreasonable requests.

**Resampling subjects and resampling raters are not symmetric.** A subject
drawn twice is two subjects that happen to be identical, which is what a
bootstrap is for. A rater drawn twice is one rater agreeing perfectly with
themselves, which inflates every coefficient: on three raters, 77% of draws
contain a duplicate and 11% return exactly 1.0.

**Dropping a rater changes what the data describes.** Dropping a subject does
not. One narrows the rating task; the other narrows the sample.

**And what you can do depends on what you are holding.** A matrix of category
counts is enough for any coefficient that treats raters as interchangeable, and
resampling its rows is a perfectly good subject bootstrap. But it has no
columns: rater identity was discarded when the counts were made, so leaving one
rater out is not a harder computation, it is an impossible one. Nothing about a
bare array says which of those you have.

That is the case for types. Not that they make the arithmetic cleaner (they do
not) but that the number of *legal* combinations here is large, the number of
illegal ones is larger, and several of the illegal ones return a plausible
number rather than an error. A function will compute leave-one-rater-out on
counts that have no raters in them. An object that never had rater identity
cannot be asked for it.

It also means the awkward cases stay possible. Subsets of three raters,
weighted variants, a coefficient nobody has written yet: these compose with the
existing resamplers rather than requiring a new function for each pairing. That
is harder to arrange with free functions, which may be part of why these
analyses are not routine in packages that offer them.