## Why every pair

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

### Why you may be interested in the pairwise/subset permutations

For a panel of raters, `interrater` offers the option of calculating 
each coefficient (Cohen, Fleiss, AC1 and so on) as permutations of pairs
or subsets. For pairwise permutations, the result would be, for example, a coefficient for every pair of raters, the mean of those pairs, and a
confidence interval on that mean. 

The reason you may be interested is that in practice raters may not be "identical" and 
we need ways to visualize or describe it that are intuitive. For example, you may have a systematic
review where some of the abstract readers have more experience in the exclusion process or more experience in the topic. Furthermore, despite good efforts, some readers may spend more time per abstract than others. Among pathologists, some may be more experienced than others.
A group coefficient does not give you information on one rater that may have
disagreed more with others than the rest. Other examples could be having ratings where some are 
LLMs and others are people. 

The group coefficient and the pairwise set
answer different questions, and sometimes the reviewer's question is the
second one.

### The pairs locate the anomaly

Four reviewers screened 787 abstracts for a systematic review, each abstract
marked include or exclude. Fleiss' kappa for the panel is 0.708, with a 95%
interval of [0.663, 0.750]. The six pairwise Cohen's kappas:

| raters | κ | 95% CI |
| --- | --- | --- |
| rater0, rater1 | 0.824 | [0.768, 0.873] |
| rater0, rater2 | 0.781 | [0.719, 0.831] |
| rater1, rater2 | 0.752 | [0.688, 0.806] |
| rater0, rater3 | 0.687 | [0.611, 0.753] |
| rater2, rater3 | 0.595 | [0.519, 0.668] |
| rater1, rater3 | 0.583 | [0.507, 0.660] |

The three pairs without rater 3 are the three highest. The three pairs with
rater 3 are the three lowest. One reviewer accounts for the whole spread, and
the panel coefficient of 0.708 would have come out the same from a panel in
which all six pairs sat near 0.70. 

A second review, 770 abstracts and three reviewers, gives Fleiss' kappa of
0.675 and pairwise kappas of 0.649, 0.681 and 0.697, a span of 0.05. Here the
single number seems to be the whole story and the pairs confirm it.

Those two situations call for different responses: a conversation with one
reviewer, or nothing at all. The group coefficient is the same kind of number
in both. The pairwise permutation plot is built to make the difference
visible at a glance: every pair grouped by rater, an interval on each, and a
shaded band for the interval on the mean, so that one rater pulling the panel
shows up as a cluster of low bars under one name. Rater 3 in the first
example is not necessarily the wrong one: maybe rater 3 had more insight
about the project and could share it with the others, or maybe rater 3 just
read faster. The point is that this printout and its plot find an anomaly
that a group coefficient and its interval would have hidden.

### The pairs are well defined when the panel is not

Fleiss' kappa assumes every subject was rated the same number of times, and
Fleiss said plainly that varying numbers were out of scope. Applying it to an
unbalanced panel means choosing among generalisations he did not give, and
two implementations can choose differently. A pair has no such choice. It is
computed on the subjects both raters saw, so the number of ratings per subject
is exactly two by construction, and Cohen's formula applies as written.

That comes with a cost. On an unbalanced panel, different pairs are computed
on different subject sets, and the mean of the pairs is a mean across them.
Light defined that mean for complete data, where every pair shares every
subject. On incomplete data the package takes the unweighted mean of the
pairs, each on the subjects it shares. Weighting the pairs by how many
subjects they share is a choice we have not taken, and we do not know of a
standard form for it.

### The interval on the mean

Every pair is computed on the same subjects, so the six pairwise kappas of a
four-rater panel are not six independent numbers. An abstract that both
confuses rater 0 and confuses rater 2 moves several pairs at once. Averaging
the pairs and treating the result as a mean of independent estimates
understates the interval. On the 787-abstract panel above, the independent
calculation gives [0.677, 0.729] for the mean; the correct interval is
[0.658, 0.747], and the independent one is 58% of its width. This was a
mistake I made in version 1.

Version 2 resamples subjects once per draw and recomputes every pair on the
same resampled subjects, which keeps the dependence in the interval. The
interval on the mean of pairs then lands where it should: on the same panel,
[0.658, 0.747] against Fleiss' own bootstrap interval of [0.663, 0.750]. For
the arithmetic see the permutation analysis README in the source
(`analyzers/permutation_analysis/README.md`).

On an unbalanced panel the same resampling does one more thing. A pair with
few shared subjects has a wider sampling distribution across draws than a
pair with many, and the unweighted mean gives every pair the same weight, so
the thin pair contributes more of the variance of the mean than its share of
the data. The interval on the mean is widened by the thinnest pair. That is
the right behavior, and it is what the independent-pairs error of version 1
could not see.

### What the package computes

The pairwise analysis takes any two-rater coefficient; the examples above use
Cohen's kappa, whose unweighted mean over all pairs is Light's kappa. Seven
coefficients can be used pairwise, and six in subset permutation (Scott's pi
is two raters by definition).

The interval on the mean is computed only when asked for
(`bootstrap_mean_config`), as the correlated or the independent bootstrap.
Asking for the independent one with a BCa interval is refused at
configuration, because BCa's acceleration needs a jackknife over the same
units the bootstrap resamples, and the independent mean has none.

A panel is also something you can add to. See
[What if one of your raters were someone else?](rater_substitution.md).
