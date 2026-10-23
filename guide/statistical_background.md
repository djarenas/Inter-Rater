## Two ways to mean "by chance"

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

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
