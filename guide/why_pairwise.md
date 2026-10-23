## Why pairwise, and not only a single number

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

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
same resampled subjects, which keeps that dependence in the interval. For more information see: [Correlated and independent error](statistical_background.md#correlated-and-independent-error).

**And agreement structure is not the only thing you can perturb.** A panel is
also a thing you can add to: a reviewer who agrees with whatever the others
said, a contrarian, someone guessing from a declared distribution. Asking what
such a rater would have done to your coefficient is a harder analysis to
arrange than any single formula, and arranging it is most of why version 2 is
built the way it is: a tight object-oriented design that composes, even though
the ramp-up is steeper than a function call. See
[What if one of your raters were someone else?](rater_substitution.md)