## What if one of your raters were someone else?

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

A coefficient of 0.68 computed from three reviewers is a statement about those
three. It does not say whether a fourth reviewer who rubber-stamped the
majority would have raised it, or whether one careless addition would have
destroyed it, and those are the scenarios a researcher may ask about: *how
much does this number depend on these particular people?* Resampling subjects,
dropping a rater, or taking every pair all leave the panel as it was. As far
as we know, in other agreement packages the panel is the one thing the
analysis can never change.

`interrater` lets you add a rater who does not exist. You declare what kind of
rater they are and get back a new `Ratings` with a synthetic column appended.
The original is never modified, and the new columns are named
`injected_rater0`, `injected_rater1`, so an injected panel cannot be mistaken
for a real one.

### Eight kinds of rater

Ordered by how much the synthetic rater knows about your subjects, from
copying a real answer to ignoring the subjects entirely:

| kind | class | what it simulates |
| --- | --- | --- |
| echo | `WithinSubjectInjector` | copies a real answer at each subject, with tunable noise: a reviewer who follows the panel |
| clone | `RaterBootstrapInjector` | resamples one real rater's column: keeps their base rates, scrambles their judgments |
| chance, pooled rates | `PooledInjector` | draws from the pool of all observed ratings: correct base rates, no subject signal |
| chance, declared rates | `DistributionInjector` | draws from a distribution you declare: a pure null rater |
| conformist | `ExtremalInjector` | always picks the majority category: goes with the crowd |
| dissenter | `ExtremalInjector` | picks the least common category anyone actually used |
| contrarian | `ExtremalInjector` | picks a category nobody used on that subject |
| median | `ExtremalInjector` | picks the middle-prevalence category |

The dissenter picks the road less travelled; the contrarian picks the road no
one took.

Every injector is seeded (`random_state`, 42 by default, as everywhere in the
package; pass your own, or `None` for fresh randomness). The sampling
injectors are random by nature, and even the extremal strategies break ties at
random when a subject has no single majority category.

```python
from interrater.injectors import ExtremalInjector
from interrater.injectors.configs import ExtremalInjectorConfig

augmented = ExtremalInjector(
    ratings,
    ExtremalInjectorConfig(n_new_raters=1, strategy="majority", random_state=42),
).inject()

FleissKappaAnalyzer(augmented, config).estimate()
```

Injections chain, so a panel can be stressed by several hypothetical reviewers
at once without name collisions. Each step is seeded with its own
`random_state` and its position in the chain, so a pipeline reruns
identically and reordering it changes the result.

### What the eight kinds are for

They are reference points, not portraits of real reviewers. Two of them are
impossible under blind rating: a conformist and a contrarian both need to see
the other answers before giving their own. Others a real rater could
produce: someone guessing at the base rates, someone guessing uniformly, a
clone's scrambled judgments with intact habits, and, where ratings are not
blind, a noisy follower. What the set does is bound the movement. The
conformist is the most flattering single rater a panel can gain, the
contrarian the most damaging, and a real added rater lands between them.

Run the eight on one of the screening rounds shipped with the package
(`example_data/SR3_JGIM2019`: 770 abstracts, three real reviewers, Fleiss'
kappa 0.675). One injected rater each, `random_state=42`:

| added rater | Fleiss' κ | change |
| --- | --- | --- |
| conformist | 0.755 | +0.080 |
| echo, no noise | 0.732 | +0.057 |
| dissenter | 0.679 | +0.004 |
| echo, 20% noise | 0.397 | −0.278 |
| chance, pooled rates | 0.324 | −0.351 |
| clone of reviewer 0 | 0.306 | −0.369 |
| uniform guesser | 0.128 | −0.547 |
| contrarian | −0.226 | −0.901 |

The table is ordered by how much the added rater knows about the subjects,
and the coefficient follows that order. Raters correlated with the panel push
it up; raters independent of the panel push it down; a rater who
systematically opposes the panel pushes it below zero. Even the upper bound
is worth knowing: a fourth reviewer who contributes no judgment and simply
agrees with the other three raises the coefficient by 0.08, so a panel
coefficient cannot be read as a count of independent readers.

The size of the drops has a simple account. With four raters there are six
pairs, three of them real. A fourth rater who carries no information about
the subjects brings the three pairs that include them to roughly chance
agreement, so the mean over pairs falls to about half its value: 0.5 × 0.675 is
0.34, and the pooled-rates row gives 0.324. The clone of reviewer 0 is the
same thing with one reviewer's base rates instead of the panel's, and lands at
0.306.

The uniform guesser falls further because it also changes what chance means.
The real reviewers mark about nine abstracts in ten as exclude, so two of them
agree by chance about 82% of the time; a rater who marks half of everything
exclude agrees with a real reviewer only about half the time. Observed
agreement falls from 0.943 to 0.721, and although the chance term also falls,
from 0.824 to 0.680, because the pooled prevalence has moved toward the
middle, it falls by less than the agreement does, and the coefficient ends at
0.128.

The echo with 20% noise sits between the correlated and the independent
rows, which is what it is: a reviewer who follows the panel four times in five
and guesses the fifth.

A question these kinds could serve, and this page does not answer, is the
sensitivity of a panel: given a panel of a certain size and agreement, how
different does an added rater have to be before the pairwise table shows
them? Injecting a rater of a declared kind and looking at the pairs is the
experiment; the answer would be a property of the panel rather than of any
one rater, and would need to be studied across panels rather than read off one.

The injectors are in the architecture early on purpose. Once a `Ratings`
cannot be modified and every analyzer accepts any `Ratings`, a synthetic rater
is one more class, and the simulations that use it, the sensitivity study
above and whatever else turns out to be fun, cost nothing to add later and a
great deal to retrofit.

### What it is not

This is sensitivity analysis, not data generation. An injected coefficient
describes a hypothetical panel and belongs in a robustness paragraph, never in
the line where you report your study's agreement. The naming and the
immutability exist so that the distinction survives being written down: the
column is called `injected_rater0`, and the original `Ratings` is still there.

One connection worth making explicit. Resampling raters with replacement is
refused as a confidence interval method, because a rater drawn twice agrees
perfectly with themselves and the interval narrows around an agreement that
was never observed. The same mechanic is offered here as the clone injector,
where the duplication is the declared intent: you asked for a rater with
reviewer 0's habits, and the result is labelled as such.
