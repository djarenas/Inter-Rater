## What if one of your raters were someone else?

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

Every operation above perturbs the data: drop a subject, drop a rater, take the
pairs. The panel itself stays fixed, and in every other agreement package it
stays fixed permanently: the raters you had are the raters the analysis knows
about.

But the question a reviewer actually asks is about the panel. *How much does
this kappa depend on these particular people?* A coefficient of 0.68 computed
from three reviewers is a statement about those three. It does not say whether a
fourth reviewer who rubber-stamped the majority would have raised it, or whether
one careless addition would have destroyed it, and those are the scenarios
worth knowing before you report the number.

Version 2 answers that by letting you add a rater who does not exist. You
declare what kind of rater they are, and you get back a new `Ratings` object
with a synthetic column appended. The original is never modified, and the new
columns are named `injected_rater0`, `injected_rater1`, so an injected panel can
never be mistaken for a real one.

There are eight kinds, which can be ordered by how much the fake rater knows
about your actual subjects:

| kind | class | what it simulates |
| --- | --- | --- |
| **Echo** | `WithinSubjectInjector` | copies a real answer at each subject, with tunable noise: a reviewer who follows the panel |
| **Clone** | `RaterBootstrapInjector` | resamples one real rater's column: keeps their base rates, scrambles their judgments |
| **Chance** | `PooledInjector` | draws from the pool of all observed ratings: correct base rates, no subject signal |
| **Chance** | `DistributionInjector` | draws from a distribution you declare: a pure null rater |
| **Conformist** | `ExtremalInjector` | always picks the majority category: goes with the crowd |
| **Dissenter** | `ExtremalInjector` | picks the least common category anyone actually used |
| **Contrarian** | `ExtremalInjector` | picks a category nobody used on that subject: goes where nobody went |
| **Median** | `ExtremalInjector` | picks the middle-prevalence category |

The dissenter picks the road less travelled; the contrarian picks the road no
one took.

Every injector takes a `random_state`. Pass one: the sampling injectors are
random by nature, and even the extremal strategies break ties at random when a
subject has no single majority category.

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
at once without name collisions.

### What this tells you

Run the spectrum on one of the screening rounds here (770 abstracts, three real
reviewers, Fleiss' kappa 0.675) and the panel's dependence on its own
composition becomes visible (one injected rater each, `random_state=42`):

| added rater | Fleiss' κ | change |
| --- | --- | --- |
| conformist | 0.755 | **+0.080** |
| echo, no noise | 0.732 | +0.057 |
| dissenter | 0.679 | +0.004 |
| echo, 20% noise | 0.397 | −0.278 |
| chance, pooled rates | 0.324 | −0.351 |
| clone of reviewer 0 | 0.306 | −0.369 |
| uniform guesser | 0.128 | −0.547 |
| contrarian | −0.226 | −0.901 |

The first row is the one worth sitting with. A fourth reviewer who contributes
no independent judgment at all, who simply agrees with whatever the other three
decided, *raises* the reported reliability by 0.08. Agreement statistics cannot
distinguish a panel of careful independent readers from a panel with a
rubber stamp in it, and the direction of the error is the flattering one.

The ordering of that table is not arbitrary either. Raters correlated with the
panel push the coefficient up; raters independent of it push the coefficient
down. That is the same phenomenon as the shared-subject dependence in the
pairwise mean, and the next section treats them together.

### What it is not

This is sensitivity analysis, not data generation. An injected coefficient
describes a hypothetical panel and belongs in a robustness paragraph, never in
the line where you report your study's agreement. The point of the naming and
the immutability is that the distinction survives being written down.

One connection worth making explicit: resampling raters with replacement is
refused here as a confidence interval method, for the reason given above: a
rater drawn twice agrees perfectly with themselves and inflates everything. The
same mechanic is offered as the `Clone` injector, because there the duplication
is the declared intent rather than an accident of the procedure. The difference
between a bug and a feature is whether you asked for it.
