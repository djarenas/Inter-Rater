# interrater

Agreement between raters, for any rating task where there is no ground truth to
check against. For example: systematic review screening, content coding, annotation
quality, evaluation panels, etc.

Version 2 is a complete rewrite in object-oriented Python; see
[What version 2 adds](#what-version-2-adds) below.

## Is this package good for you?

If you are looking for a kappa you can defend later, one that arrives with a
pipeline you can re-run from config files, result files that carry the details
behind it, and plots you can produce directly from those results, this package
is for you.

Why? Because although there is no shortage of ways to compute a kappa, and
several good libraries do it well in a handful of lines of R or Python, what is
harder to come by is reproducible pipelines, result/run details, and knowing
whether the number means what you think it means.

For example, the original Fleiss' kappa assumes every subject was rated the same
number of times. If yours were not, the formula does not apply. A package function
may still run, and nothing afterwards would tell you. Other silent issues may be
that not all your subjects were rated by all raters (hereafter referred to as
unbalanced). A category nobody chose changes coefficients whose chance term counts
categories, and leaves the others untouched. This may not matter for Fleiss' kappa,
but it matters for Brennan-Prediger, Gwet's AC1, and Bennett's S, where the number
of categories goes into the kappa. Not knowing how your package function handled
unchosen categories would be bad.

Also what happens when a rating was not something you expected (e.g. '0' instead of
'o' or 'O'). Even these data cleaning or management issues are hard to handle with
scripts in a way that is easily reproducible and archivable.

What is also missing in existing packages is detail beyond two raters, such as
whether one rater disagreed with the others more than usual, along with
extendability. Dropping a rater who quit halfway changes what your data
describes. A pair of raters who share eight subjects produces a coefficient
that is perfectly valid and perfectly useless.

The `interrater` package is built for the case where you want those things checked rather
than assumed. Every analysis runs through one pipeline: a cleaner that declares
its decisions in a config file and reports what it did, a validated data object
that cannot exist in an invalid state, and analysers that refuse data their
coefficient does not apply to. The intent is that a wrong number is hard to
produce, and that a right one arrives with an account of how it was reached.

The trade is deliberate. There is no shortcut that takes a bare array and
returns a coefficient, because a shortcut is exactly where the guarantees are
lost. If you want the two-line version, other packages do it well and this one
is not trying to compete with them.

## Install

```bash
pip install interrater
```

Python 3.10 or later. Depends on numpy, scipy, pandas, pyyaml and matplotlib.

## Quickstart

An analysis needs three things: your ratings, a declaration of what they mean,
and a declaration of how to read them. The second and third are files rather
than arguments, so that a reviewer can be shown them.

This example runs on data that ships with the repository: the preliminary
screening round of a real systematic review, 199 abstracts rated by four
reviewers. Clone the repository, or copy
[`example_data/SR1_OSF2018/`](https://github.com/djarenas/Inter-Rater/tree/master/example_data/SR1_OSF2018/), and run it as it stands.

`SR1_prelim_semantics.yaml` declares who rated and what the categories mean:

```yaml
raters: 4
categories:
    kind: unordered
    names: ['1', '0', '-1']      # include / unclear / exclude
n_subjects: 199                  # asserted, not assumed
unique_subjects: true
```

`SR1_prelim_cleaner_config.yaml` declares how to read the file and what to do
with the mess in it:

```yaml
frame_layout:
    use_column_names: true
    subject_id_col: "(internal) id"
case_sensitive: false
strip_whitespace: true
missing_values: ['o']            # this rater did not rate this abstract
value_map: {}                    # nothing needs repairing
unknown_policy: 'error'          # an unexpected string stops the run
unused_categories: 'error'
insufficient_data_drops:
    min_ratings_per_rater:   2
    min_ratings_per_subject: 2
    method:                  iterate
```

Then:

```python
from pathlib import Path
import pandas as pd
from interrater import (
    RatingsSemantics, CleanerConfig, RatingsCleaner,
    FleissKappaAnalyzer, FleissKappaConfig,
)

data = Path(https://github.com/djarenas/Inter-Rater/tree/master/example_data/SR1_OSF2018/)

semantics      = RatingsSemantics.build_from_yaml(data / "SR1_prelim_semantics.yaml")
cleaner_config = CleanerConfig.build_from_yaml(data / "SR1_prelim_cleaner_config.yaml")

frame   = pd.read_csv(data / "SR1_prelim_data.csv", dtype=str)
cleaner = RatingsCleaner(frame, semantics, config=cleaner_config)
ratings = cleaner.ratings          # validated; cannot exist in an invalid state

print(cleaner.report.render())     # what cleaning did, before any statistics

config = FleissKappaConfig(n_resamples=4000, random_state=42, ci_method="bca")
result = FleissKappaAnalyzer(ratings, config).estimate()
print(f"Fleiss' kappa = {result.point:.3f}  95% CI "
      f"[{result.confidence_interval.lower:.3f}, {result.confidence_interval.upper:.3f}]")
```

```
CleaningReport
    in         : 199 row(s) x 5 column(s)
    out        : 199 subject(s) x 4 rater(s) (100.0% of rows kept)

Matching
    cells      : 796
    no rating  : 0 (0.0%), of which 0 arrived null
    matched    : 'o' 0

Drops
    nothing dropped (iterate)

Categories
    used       : 1 13%, 0 4%, -1 83%

Fleiss' kappa = 0.493  95% CI [0.383, 0.606]
```

The report is the part worth noticing. Before any coefficient appears you are
told that all 199 rows survived, that 796 cells were read and none were
missing, and that 83% of ratings were "exclude". That last line is not
decoration: under skew that heavy, Fleiss' kappa is held down by its own chance
term, which is why the same data gives 0.828 under Gwet's AC1. You would want
to know that before quoting 0.493, and it arrives without being asked for.

It is more ceremony than a one-line kappa, and that is the trade this package
makes deliberately. The next section says why.

### Then

Every pair of reviewers, with the mean bootstrapped across them correctly:

```python
from interrater import (PairwisePermutationAnalysis, PairwisePermutationConfig,
                        CohenKappaConfig, PermutationBootstrapMeanConfig)

pairwise = PairwisePermutationAnalysis(ratings, PairwisePermutationConfig(
    coefficient_to_compute="cohen_kappa",
    config_of_coefficient=CohenKappaConfig(n_resamples=4000, random_state=42, ci_method="bca"),
    bootstrap_mean_config=PermutationBootstrapMeanConfig(
        method="correlated", n_bootstraps=4000, random_state=42, ci_method="bca"),
)).run()

pairwise.pairwise_values      # (0,1) 0.539  (0,2) 0.453  (0,3) 0.425
                              # (1,2) 0.553  (1,3) 0.522  (2,3) 0.481
pairwise.mean_over_pairs      # 0.496
print(pairwise.bootstrap_ci_of_mean)
                              # [0.387, 0.608] (95% CI, bca, n=4000)
```

And what a fifth reviewer would have done to the number:

```python
from interrater import ExtremalInjector, ExtremalInjectorConfig

augmented = ExtremalInjector(
    ratings,
    ExtremalInjectorConfig(n_new_raters=1, strategy="majority", random_state=42),
).inject()

FleissKappaAnalyzer(augmented, config).estimate().point   # 0.563, up from 0.493
```

A reviewer who contributes no independent judgment, who simply agrees with
whatever the other four decided, raises this panel's reported reliability by
0.07. (The seed matters: with four reviewers and three categories a majority
can tie, and ties are broken at random.) See
[What if one of your raters were someone else?](#what-if-one-of-your-raters-were-someone-else)

Worked configurations for ten datasets, five of them real systematic reviews,
are in [`example_data/`](example_data/), each with the reference values it
should reproduce.

## Background

`interrater` is the successor to
[Inter-Rater v1](https://github.com/djarenas/Inter-Rater), described in
[arXiv:1809.05731](https://arxiv.org/abs/1809.05731) (2018) and archived at
[Zenodo](https://doi.org/10.5281/zenodo.1227660). Version 1 has been used in
published systematic reviews and is cited in peer-reviewed research. Version 2
is a complete rewrite in object-oriented Python, undertaken to fix maintenance
issues and limitations that emerged from real-world use. Version 1 remains
available at git tag `v1.5`.

## Why pairwise, and not only a single number

Version 1 computed a coefficient for every pair of raters as well as one for
the panel, and version 2 keeps that. Two reasons, and the second only became
clear while building version 2.

**A group coefficient cannot tell you that one rater was the problem.** In one
of the screening rounds behind this package, the pairwise kappas run from 0.58
to 0.82 and one reviewer accounts for the whole spread; Fleiss' kappa for the
same data is 0.71, with no hint that anything is uneven. In another systematic review,
three pairs span a 0.05 variation. Those are different situations calling for different
responses (a conversation with one reviewer, or nothing at all) and the
single number is identical in kind for both.

**On unbalanced data, pairwise arithmetic is better defined.** Fleiss' kappa
assumes every subject was rated the same number of times, and he said plainly
that varying numbers were out of scope. Applying it anyway means choosing among
generalisations he did not give. A pair, by contrast, is always analysed
on the subjects both raters saw, so the count is exactly two for every subject
by construction. There is no generalisation to choose and nothing for two
implementations to disagree about.

**Version 2 makes this argument visual.** The pairwise permutation plot shows
every C(n,2) pair grouped by rater, with bootstrap confidence intervals on each
pair and a shaded band for the uncertainty on the mean. A reviewer can see at a
glance whether one rater is pulling the group coefficient, whether pairs cluster
or spread, and how the uncertainty on any single pair compares to the
uncertainty on the average, all things a single number hides.

**The uncertainty on the mean is not the obvious calculation.** Every pair is
computed on the same subjects, so the six pairwise kappas of a four-rater panel
are not six independent numbers; an abstract that both confuses rater 0 and
confuses rater 2 moves several pairs at once. Averaging them and treating the
result as a mean of independent estimates can understate the interval, badly. This was a 
mistake I made in Version 1. In contrast, Version 2
resamples subjects once per draw and recomputes every pair on the
same resampled subjects, which keeps that dependence in the interval. For more information see: [Correlated and independent error](#correlated-and-independent-error).

**And agreement structure is not the only thing you can perturb.** A panel is
also a thing you can add to: a reviewer who agrees with whatever the others
said, a contrarian, someone guessing from a declared distribution. Asking what
such a rater would have done to your coefficient is a harder analysis to
arrange than any single formula, and arranging it is most of why version 2 is
built the way it is: a tight object-oriented design that composes, even though
the ramp-up is steeper than a function call. See
[What if one of your raters were someone else?](#what-if-one-of-your-raters-were-someone-else)

## Why this is built out of objects, and what the flexibility buys

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

## What if one of your raters were someone else?

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

## Two ways to mean "by chance"

Every chance-corrected coefficient must decide what to do when two raters use
the categories at different rates, when one says "include" twice as often as
the other. There are two defensible answers, and which one a coefficient picks
is most of what separates the familiar names from each other.

|  | pooled | individual |
| --- | --- | --- |
| two raters | Scott's pi | Cohen's kappa |
| many raters | Fleiss' kappa | Conger's kappa |

**Pooled** treats the raters as one distribution and counts a difference in
their rates as part of the disagreement. **Individual** lets each rater keep
their own rate and subtracts the difference out before measuring, treating a
systematic difference in threshold as a known feature of the raters rather than
as disagreement.

Version 2 reports both, because the **gap between them is informative in its
own right.** It is zero when two raters use the categories at the same rate and
grows as their rates diverge, so it measures how differently two raters draw
the line, as distinct from how erratic they are. Those are different problems
with different remedies: a conversation about the criteria, or retraining. A
single coefficient cannot tell them apart, and neither can a single pair of
them without the comparison.

## Why resampling, and not a formula

Getting a standard error for kappa right took the literature three papers and a
decade. Fleiss, Cohen and Everitt corrected Cohen's in 1969; Fleiss published
one for many raters in 1971; Fleiss, Nee and Landis corrected *that* in 1979,
checking themselves by Monte Carlo. Every coefficient needs the work done
again, which is part of why Krippendorff's alpha is conventionally
bootstrapped: the analytic route is hard enough that nobody insists on it.

Resampling gives one procedure for all of them, so adding a coefficient costs
nothing. But avoiding algebra is the least of it, and two better reasons come
first.

**The published formulas answer a different question.** Both Fleiss variances
are derived under the hypothesis of no agreement beyond chance. They tell you
whether a kappa could have arisen by accident. Putting an interval around a
kappa that is *not* zero is a different question, and the answer is not the
same number.

**They also assume more.** An analytic variance comes with a sampling model
attached: in 1979, that each subject's ratings are multinomial and that
subjects are independent. If your subjects are clustered, or a rater's mistakes
run together across subjects, the formula is wrong and nothing says so.
Resampling subjects asks only that subjects be exchangeable.

None of which makes the null test unnecessary. "Could this have happened by
chance" is a real question and version 2 answers it, under its own name,
returning its own number, so that the two cannot be mistaken for each other.

## Correlated and independent error

Almost every default in agreement statistics assumes that errors are
independent: that two raters disagree for unrelated reasons, that one abstract
being hard tells you nothing about the next, that two estimates computed from
the same data can be averaged as though they came from separate studies. When
that assumption fails, it fails in a consistent direction: the numbers look
more certain, or more reliable, than they are. Four cases come up in practice,
and version 2 takes a position on each.

**Pairs that share subjects.** The C(n,2) pairwise coefficients are all
computed on the same abstracts, so their errors are correlated by construction.
Bootstrapping each pair separately and assembling a mean from those draws
ignores it; resampling subjects once and recomputing every pair on the same
resample preserves it. The gap is not subtle:

| round | pairs | correlated | independent | understated by |
| --- | --- | --- | --- | --- |
| SR1 preliminary, 199 abstracts | 6 | 0.222 | 0.117 | 47% |
| SR4, 787 abstracts | 6 | 0.089 | 0.053 | 40% |
| SR2, 81 abstracts | 6 | 0.335 | 0.200 | 40% |
| SR3, 770 abstracts | 3 | 0.137 | 0.106 | 23% |

(width of the 95% interval on the mean pairwise Cohen's kappa, 4000 resamples)

Three raters give three pairs and a 23% understatement; four raters give six
pairs and roughly twice that. The more pairs you average, the more dependence
there is to ignore. `correlated` is the default for exactly this reason, and
`independent` is kept available because the comparison is the argument.

**Raters who share a bias.** Two reviewers trained together, reading the same
ambiguous abstract, will often be wrong in the same direction. Every
chance-corrected coefficient reads that as agreement, because agreement is all
it can see. No coefficient can separate "these readers are careful" from "these
readers share a blind spot." But the injectors can show you what the second
one does to your number, which is what the conformist and echo rows of the
table in the previous section are for.

**Subjects that cluster.** Abstracts from the same trial, the same author
group, or the same journal are not independent draws. An analytic variance
derived under independent multinomial subjects is simply wrong on such data and
says nothing about it. The subject bootstrap asks less (only that subjects be
exchangeable), but clustering strains that too, and the honest answer is that
the resampling unit should be the cluster. Version 2 does not yet offer a
clustered bootstrap; the limitation is stated here rather than hidden.

**Errors that run together within a rater.** A reviewer who grows stricter over
a long screening round, or more tired, produces mistakes that are correlated
across subjects rather than scattered among them. This is the assumption
violation that analytic variance formulas are least able to detect, and the
case where leaving a rater out one at a time is more informative than any
interval.

The through-line is that correlated error is not an exotic failure. It is the
normal condition of a screening panel, and most of what distinguishes this
package from a two-line kappa is machinery for noticing it.

## What version 2 adds

- **Confidence intervals that estimate what they claim to.** Bootstrap
  intervals from resampled subjects, with the null test kept separate and
  separately named. See above for why.

- **An auditable cleaning layer.** Which strings meant nothing, which were
  repaired and into what, how thin the data was allowed to get: declared in a
  versioned config file and reported afterwards, rather than reconstructed from
  memory when a reviewer asks. If cleaning changed what the data describes
  (a rater dropped, a category removed), the report says so and hands back a
  corrected declaration.

- **Preconditions in the type system.** Exact Fleiss' kappa requires a constant
  number of raters per subject; that is a property of the data, so it is
  checked before the coefficient is computed rather than assumed by it.

- **Synthetic raters, in eight declared kinds.** Add a reviewer who does not
  exist (a conformist, a contrarian, an echo of the panel with tunable noise, a
  guesser drawing from prevalences you declare) and see what your coefficient
  would have been. The original data is never modified and injected columns are
  named as such. This does not appear to be available in other agreement
  packages, and it is the clearest way to ask how much a reported kappa depends
  on the particular people who produced it. See above.

- **Intervals that account for shared subjects.** The mean over rater pairs is
  bootstrapped by resampling subjects once and recomputing every pair, because
  the pairs are not independent of each other. Treating them as independent
  understates the interval by 23 to 47% on real screening data.

- **Support for unbalanced designs.** Where each subject is rated by a subset
  of the rater pool rather than all of them, which is what most rating tasks
  can actually afford.

- **Pairwise permutation with visualization.** Compute any agreement
  coefficient for every pair of raters, bootstrap each pair and the mean, and
  plot the full agreement structure: which pairs agree, which do not, and
  how much of the spread one rater accounts for. The plot supports seven
  coefficients (Cohen's kappa, Scott's pi, Conger's kappa, Gwet's AC1,
  Fleiss' kappa, Brennan-Prediger, Krippendorff's alpha), so the same data
  can be viewed under different chance models side by side.

## What is coming next

Version 2.0 is the pipeline: cleaning, validated data objects, eight
coefficients, bootstrap intervals, permutation tests, pairwise and subset
analyses, simulated raters. What follows is mostly about driving that pipeline
with less ceremony, and widening what it can compute. This is a statement of
direction, not a schedule.

- **Declarative configuration for plots.** Version 2.0 ships plotting as
  standalone functions: pass a result, get a figure. That is fine once and
  awkward for a figure in a paper, because the choices that shape it (sort
  order, which intervals are drawn, labels) live in the call rather than in a
  file. Plotting moves to the config-plus-plotter pattern already used here for
  cleaning and analysis, and used throughout
  [eventus](https://github.com/djarenas/eventus), so that a figure is
  reproducible from a declared file. The standalone functions stay, with the
  config objects built on top.

- **One YAML story for every config.** A prerequisite for the above, and worth
  having on its own. Cleaning, semantics and the top-level analyzer config load
  from YAML today; the leaf configs do not. When they all do, an entire
  analysis (cleaning, coefficient, resampling, plotting) is declarable in
  files and recoverable from them.

- **Gwet's AC2 and extended weighting schemes.** AC2 is AC1 with weights, so it
  waits on extending the weight set beyond identity, linear and quadratic to
  ordinal, interval and ratio, which also widens Krippendorff's alpha to the
  variants its literature defines. Loading a custom weight matrix from file
  arrives with it. Further coefficients follow once that foundation is in
  place.

- **Subset permutation, properly visualised.** Version 2.0 plots subsets sorted
  by value and as a strip. Neither yet shows what the analysis is for: how
  agreement degrades as raters are removed, and whose absence accounts for it.
  A panel across subset sizes is planned, along with a smaller fix that is
  overdue: results carry rater indices, so plots read `(0, 1)` where they
  should read `(Smith, Jones)`.

- **Orchestrators, for computing every coefficient at once.** Eight
  coefficients today means eight analyzers and eight configs, with the
  bookkeeping left to you. An orchestrator takes one dataset and one resampling
  budget and returns the set, computed on the same resampled subjects so the
  numbers are comparable.

  A function that hands back eight coefficients is also a function that invites
  reporting the largest one, and the position of this package is that the
  comparison across chance models is informative only if the set is reported as
  a set: the gap between Fleiss' kappa and AC1 on skewed data is a fact about
  the data, not a menu. So the design is meant to make that the easy path: a
  primary coefficient declared in the config before anything is computed and
  recorded in the result, output in a fixed order with no ranking and no
  sort-by-value convenience, and each coefficient labelled with the chance
  model it assumes so a reader can see why two numbers differ. Alongside it, a
  utility for bootstrapping the difference between two estimates directly,
  because the tempting shortcut, reading non-overlapping 95% intervals as a
  significant difference, is closer to a p < 0.005 test than a p < 0.05 one
  (Schenker & Gentleman, 2001).

- **Rater influence by replacement.** For each rater: drop them, inject a
  synthetic replacement (a contrarian, a draw from a declared distribution, a
  resample of the rest of the panel) and re-run. One coefficient per "what if
  this rater were someone else" scenario, which distinguishes a rater who was
  already noise from one who was holding the panel together. The injectors ship
  in 2.0; the loop around them does not.

- **Fewer lines, without losing the account of what was decided.** Building an
  analysis today means naming every config explicitly. That is deliberate where
  the choice matters, and it is more ceremony than it should be where the choice
  is obvious: a first coefficient on clean balanced data does not need four
  objects. 2.1 will explore defaults for config values, under two conditions:
  anything assumed rather than declared says so in a warning, and the result
  records which defaults it used, so a figure remains reproducible from what the
  run reports. Where no default can be chosen safely, it stays an error rather
  than a guess. Note that this reduces the ceremony around the pipeline, not the
  pipeline itself: there is still no path that takes a bare array and hands back
  a coefficient, because that is the interface where the guarantees are lost.
  More to come.

Also planned, less immediately: diagnostics for agreement by category and rater
confusion, a visualiser for the balance profile, convergence checking cheap
enough to leave on by default, pairwise variants of the individual-marginal
coefficients, a clustered bootstrap for subjects that are not independent
draws, and generating ratings from a declared true coefficient, so that an
estimator can be checked against an answer known in advance, an interval can be
checked for actually covering it, and a review can be powered before it starts
rather than regretted after.

One thing is deliberately **not** planned: rater bootstrap will not be offered
as a confidence interval method, for the reason given above. It remains
available as an *injector*, where duplicating a rater is the declared intent
rather than an accident of the procedure.

Suggestions and bug reports are welcome at
[the issue tracker](https://github.com/djarenas/Inter-Rater/issues). If you are
waiting on something above, saying so is the most useful way to move it up.

## How to cite

If you use this software, please cite both the software and the paper
describing the method it implements. Machine-readable metadata is in
Citation details are in [`CITATION.cff`](https://github.com/djarenas/Inter-Rater/blob/master/CITATION.cff) (on GitHub, the "Cite this repository" button in the sidebar formats it for you).

Arenas, D. J. (2018). Software for analysis of inter-rater reliability by
permutating pairs of multiple raters.
[arXiv:1809.05731](https://arxiv.org/abs/1809.05731).
Version 1 archived at [Zenodo](https://doi.org/10.5281/zenodo.1227660).

## License

GNU General Public License v3 or later. Same license as version 1.
