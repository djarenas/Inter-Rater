# Statistical background

*Part of the [interrater docs](https://djarenas.github.io/packages/interrater/interrater/) · `pip install interrater` · [Source on GitHub](https://github.com/djarenas/Inter-Rater)*

## Percent agreement, and agreement by chance

The simplest measure of agreement is the proportion of subjects two raters put
in the same category. Call it p_o. It is easy to compute, easy to explain, but there
are scenarios where it could be badly misleading.

One scenario involves two raters examining the same subjects but with completely different criteria. For example, suppose two professors grade the same hundred essays; each passes 90% of them, a proportion that may be set by
something outside the rater, whether a department curve, a referral
capacity, or a screening protocol's expected yield. Furthermore, in this
scenario nobody told these two what to grade on, and they
chose differently: one reads for creativity, the other for grammar and
structure. 

Here is an example of how they could score:

|  | B: pass | B: fail |  |
| --- | --- | --- | --- |
| **A: pass** | 82 | 8 | 90 |
| **A: fail** | 8 | 2 | 10 |
|  | 90 | 10 | 100 |

p_o = (82 + 2) / 100 = **0.84**

Although nothing connected their judgments, they scored eighty-four percent agreement, which may appear like good agreement, but it isn't. If two graders each pass 90% and their criteria are unrelated, then on any given essay the chance they both pass it is 0.9 × 0.9 = 0.81, and the chance
they both fail it is 0.1 × 0.1 = 0.01. Two strangers would land on **0.82** agreement.

The professors managed 0.84. Subtract the 0.82 that chance already explains
and 0.02 is left: out of a hundred essays, two that chance does not account
for. That subtraction is the move every coefficient below makes; what they
disagree about is how to work out the 0.82.

The example is built to be extreme, deliberately. Making the two sets of
criteria independent is the only way to know in advance roughly where the
answer should land, and a coefficient that calls 0.84 strong agreement here is
not measuring agreement. Real raters are rarely that independent, but they are
rarely identical either, and whatever part of their criteria fails to overlap
contributes agreement of exactly this kind. The grammar/creativity professors
are the far end of a line that every real panel falls somewhere along.

One limit on everything below: these coefficients treat categories as
unordered, so two raters either match or they do not. On a ranked scale a
disagreement of one step and a disagreement of four count the same. Weighted
coefficients exist for that case and are outside this page.

## Cohen's kappa

Cohen's kappa ([Cohen, 1960](#ref-cohen1960)) is one of the many important agreement
statistics that take into account how often raters would agree by chance. 

**Agreement by chance**. 
The agreement by chance, p_e, is what two raters produce when each applies their own proportions independently of the other. It was 0.82 in the
previous example. For any number of
categories it is the calculation from the last section: multiply the two
raters' rates for a category, and sum over the categories.

p_e = Σ_j p_j(A) · p_j(B)

p_e is the **chance term**. Every coefficient on this page is the same
subtraction with a different chance term, and which one a coefficient builds is
most of what separates the familiar names from each other.

Calling that subtraction "agreement corrected for chance" is a simplification,
and it has been pointed out as one. Kappa also admits a reading as a
reliability coefficient ([Kraemer, Periyakoil, & Noda, 2002](#ref-kraemer2002)),
as a parameter of a symmetry model whose assumptions can be tested
([Schuster, 2001](#ref-schuster2001)), and as the mixing weight between
independence and perfect agreement ([Agresti, 1989](#ref-agresti1989)). Those
accounts, and the conditions they put on when one coefficient is enough, are
outside this tutorial.

Cohen's kappa normalizes for it:

κ = (p_o − p_e) / (1 − p_e)

The numerator is the agreement chance does not account for: the two essays.
The denominator is how much such agreement was available in the first place:
1 − p_e is the headroom above chance. So κ reports the share of that headroom
the raters actually took. One means they took all of it, zero means none of it,
and negative values mean they did worse than two strangers.

For the grammar/creativity professors:

κ = (0.84 − 0.82) / (1 − 0.82) = 0.02 / 0.18 = **0.111**

Eighty-four percent agreement, and eleven percent of the possible gains.

## Pooled or individual rates

Both professors so far have used the categories at the same percentages (each passed
90%), which is why one chance term has been enough. When two raters use the
categories at *different* proportions there are two defensible ways to build it.

For example, suppose two different professors grade essays for strength of arguments, and 
their department does not set pass/fail rates. They rank the same essays alike, but one professor tends to be a tougher grader than the other:

|  | B: pass | B: fail |  |
| --- | --- | --- | --- |
| **A: pass** | 50 | 0 | 50 |
| **A: fail** | 40 | 10 | 50 |
|  | 90 | 10 | 100 |

A passes 50, B passes 90, and they agree on 60. Note the zero: there is no
essay A passed that B failed. B passes everyone A passes, and forty more.

**Individual rates.** 
This choice to calculate chance is what Cohen's kappa, discussed above, uses. This term uses each rater's own rates. In the example above, the chance that they both pass is
0.50 × 0.90 and the chance they both fail is 0.50 × 0.10, so p_e = 0.50. Cohen's kappa
is therefore:

κ = (0.60 − 0.50) / 0.50 = **0.200**

**Pooled rates.** 
Scott's pi ([Scott, 1955](#ref-scott1955)) is the two-rater kappa that uses the pooled rate. This choice treats
the two raters as one distribution. Between them they
passed 70% of the time, so p_e = 0.70 × 0.70 + 0.30 × 0.30 = 0.58. 

p_e = Σ_j p̄_j²,  where p̄_j = (p_j(A) + p_j(B)) / 2

Scott's pi is then:

π = (0.60 − 0.58) / 0.42 = **0.048**

Same hundred essays, same sixty agreements, and the two differ by a factor of
four.

The gap between the two coefficients is driven by
the difference in their chance terms, for two raters ([Rae, 1988](#ref-rae1988)):

p_e(pooled) − p_e(individual) = ¼ Σ_j (p_j(A) − p_j(B))²

Here the two raters pass at 0.50 and 0.90, so the right side is
¼ [(0.50 − 0.90)² + (0.50 − 0.10)²] = ¼ × 0.32 = **0.08**, which is exactly
the 0.58 − 0.50 between the two chance terms above.

The further apart the two raters' rates, the wider the gap between the coefficients.

This is the second of the two paradoxes named by
[Feinstein and Cicchetti (1990)](#ref-feinstein1990): κ "will be higher with an
asymmetrical rather than symmetrical imbalance in marginal totals" (p. 543).
Here, A splits 50/50 and B 90/10, and Cohen's kappa rewards the mismatch.

One last thing to note: **The gap between the two coefficients sees threshold differences and is
blind to criteria differences.** The grammar/creativity professors both passed 90%, so pooled and individual coincide and both report 0.111.

**Krippendorff's alpha** ([Krippendorff, 2013](#ref-krippendorff2013)) belongs to the pooled family.
Details in a future section.

## When the condition is rare

Consider a new example. Two pathologists read a thousand microscope slides from
a surveillance programme, looking for evidence of cancer. They agreed on the
criteria in advance and both apply them as well as they can. Because it is
surveillance rather than referral, most of the tissue comes from healthy people
and cancer is rare.

|  | B: cancer | B: healthy |  |
| --- | --- | --- | --- |
| **A: cancer** | 15 | 50 | 65 |
| **A: healthy** | 50 | 885 | 935 |
|  | 65 | 935 | 1000 |

They agree on 900 slides. The hundred they disagree on split evenly, fifty each
way, and both call cancer 6.5% of the time, so there is no threshold
difference here, and pooled and individual chance terms coincide. This is not
the problem from the last section.

That 90% is not spread evenly, though. Of the 65 slides A called cancer, B
confirmed **15**, or 23.1%. On the healthy slides the two agree 94.7% of the
time. Nearly all of the agreement is in the category that nearly all of the
tissue falls into.

The chance term is 0.065² + 0.935² = 0.878. Chance alone would put these two on
the same side of 878 slides. They managed 900. They are twenty-two slides ahead
of chance, out of a possible hundred and twenty-two:

κ = (0.900 − 0.878) / 0.122 = **0.177**

Ninety percent agreement, and a coefficient that reads as barely better than
guessing. The headroom, 1 − p_e, collapses to a small value: 0.122. This is the other paradox
named by [Feinstein and Cicchetti (1990)](#ref-feinstein1990) and
[Cicchetti and Feinstein (1990)](#ref-cicchetti1990).
When one category holds 93.5% of the calls,
two raters could land together by chance, and little is left for skill to claim.

Some coefficients use different terms for chance and may be useful in scenarios like this.

**Bennett's S** ([Bennett, Alpert, & Goldstein, 1954](#ref-bennett1954)) does not use the
raters' rates. With q categories, chance is
1/q: two raters picking blindly between cancer and healthy would agree half
the time, whatever the prevalence.

p_e = 1/q = 0.500,  S = (0.900 − 0.500) / 0.500 = **0.800**

The 1/q term has a cost, and Scott named it in the paper that introduced pi:
categories nobody uses still move it. "Given a two-category sex dimension and
a Po of 60 per cent, the S calculated from this equation would be 0.20. But a
whimsical researcher might add two more categories … thereby increasing S to
0.47, though the two additional categories are not used at all"
([Scott, 1955](#ref-scott1955), p. 322). The raters behaved identically in both
versions. Only the coding sheet changed.

**Gwet's AC1** ([Gwet, 2008](#ref-gwet2008)) takes a middle route: chance agreement should be high when one
category dominates, but not as high as the raters' own rates imply, because
those rates are partly produced by the rarity itself.

p_e = (1/(q − 1)) Σ_j p̄_j(1 − p̄_j) = 2 × 0.065 × 0.935 = 0.122

AC1 = (0.900 − 0.122) / 0.878 = **0.886**

Gwet writes the pooled proportion p̄_j as π_k. This page uses p̄_j, as in the
Scott's pi section, to avoid confusion with Scott's π.

One set of slides and three different coefficients 0.177, 0.800, 0.886. This is not an artifact of one table. Hold the pathologists' skill fixed (each catching 80% of cancers and correctly clearing 95% of healthy slides) and move only the share of slides that are cancerous:

| cancer prevalence | p_o | κ | S | AC1 |
| --- | --- | --- | --- | --- |
| 50% | .7925 | .5754 | .5850 | .5941 |
| 30% | .8375 | .5925 | .6750 | .7297 |
| 20% | .8600 | .5625 | .7200 | .7941 |
| 10% | .8825 | .4629 | .7650 | .8496 |
| 5% | .8938 | .3346 | .7875 | .8736 |
| 2% | .9005 | .1814 | .8010 | .8867 |

Nothing about the pathologists changes down that column. The bottom row is the
pair above at full size (2% prevalence puts their call rate at exactly 6.5%),
so its near-match to the worked table is the same case, not a second one. At a
balanced caseload the three nearly coincide; the fork opens as the condition
gets rarer. Kappa
is the exception to its own direction: it edges up from .5754 to .5925 between
50% and 30% before it begins to fall, so the column rises before it turns.

To finish the section, it is interesting to compare the coefficients using the example of 
the grammar/creativity professors, the two who agreed about nothing:

|  | κ | S | AC1 |
| --- | --- | --- | --- |
| grammar/creativity professors | 0.111 | **0.680** | **0.805** |

## Why interpretability is not trivial

S and AC1 are generous for the same reason in both of the examples in the last section. Both professors
passed 90%, so a chance term that ignores the rates declines to charge them for
it; the proportion of healthy tissues is meant to be small in a surveillance example.
However, the pathology case is an example of an important concept. The 6.5% was truly
decided by the pathologists, reading the tissue in front of them. In the
grading case it is not. That 90% came from a department curve, "you must pass
around 90%," and was never theirs to set. **Nothing in the data distinguishes
the two situations, and it must be up to the researcher to understand which
one they are in.**

This is the distinction [Brennan and Prediger (1981)](#ref-brennan1981)
draw. A margin is "'fixed' whenever the marginal proportions are known to the
assigner before classifying the objects into categories," and "'free'
whenever the marginal proportions are not known to the assigner beforehand"
(p. 690). Their argument for treating the two differently is about what the
raters earned. Two judges "who independently, and with no a priori knowledge,
produce similar marginal distributions must obtain a much higher agreement
rate to obtain a given value of kappa, than two judges who produce radically
different marginals. Cohen's kappa gives the former judges no credit for
producing agreement in marginals. Indeed, they are in a sense penalized"
(p. 692).

Their prescription follows: "when marginals are fixed, coefficient kappa is
found to be appropriate. When either or both of the marginals are free to
vary, however, it is suggested that the 'chance' term in kappa be replaced by
1/n, where n is the number of categories" (p. 687). Their *n* is the category
count, written *q* here.

The two cases above fall on opposite sides of that line. The pathologists'
6.5% was free, and they arrived at it independently; the department's 90% was
fixed before an essay was read. **Which side a study is on is a fact about
what the raters were told, not a fact about the ratings**, so no amount of
looking at the table will recover it.

None of which means Cohen's kappa is malfunctioning ([Kraemer, 1992](#ref-kraemer1992);
[Kraemer, Periyakoil, & Noda, 2002](#ref-kraemer2002)). Cohen's kappa reports how much of their agreement
survives what the caseload hands them for free. AC1 reports whether they behaved
consistently. Both are true, and the question is which one you asked. It is yet
another example of the coefficient you need depending on your situation. 

Lastly, in the example where one whimsical researcher hurt Bennett's S by adding two categories that no one
used; it could be argued that if those two categories were indeed reasonable hypotheses of how some subjects would have been rated, then it offers some information too. Perhaps better subjects could have
been shown to the raters.

## More than two raters

As an example of extending two raters to multiple, first consider a panel reading a thousand biopsies for two categories: cancer and healthy. A tenth of them are cancerous. Every rater
catches 90% of the cancers and correctly clears 90% of the rest. Imagine the data looks like:

|  | B: cancer | B: healthy |  |
| --- | --- | --- | --- |
| **A: cancer** | 90 | 90 | 180 |
| **A: healthy** | 90 | 730 | 820 |
|  | 180 | 820 | 1000 |

p_o = 0.820; the pooled and
individual agreement-by-chance are equal because of the equal category rates:
p_e = 0.18² + 0.82² = 0.7048, and Cohen's kappa and Scott's pi are 0.390.

Fleiss' kappa ([Fleiss, 1971](#ref-fleiss1971)) is Scott's pi for more than
two raters ([Rae, 1988](#ref-rae1988)). For N subjects, k raters, the percent agreement is:

P_o = (1/N) Σ_i [ Σ_j n_ij (n_ij − 1) ] / [ k (k − 1) ]

where the n_ij matrix specifies how many raters assign subject i to category j.

Fleiss' kappa uses pooled values, p̄ the proportion of all ratings falling in the first category,

P_e = p̄² + (1 − p̄)²

κ_F = (P_o − P_e) / (1 − P_e)

To the example above, we can add raters with specificity and sensitivity that are both 90%: 

The theoretical columns are what such a panel gives on average. The observed
columns are one simulated panel of 1,000 biopsies (seed 42; the data and the
script that computes them are in `example_data/Statistical_background/`).

| raters | pairs | theoretical P_o | theoretical P_e | theoretical κ_F | theoretical all agree | observed κ_F [95% CI] | observed all agree |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | 1 | 0.820 | 0.7048 | 0.390 | 82.0% | 0.375 [0.299, 0.444] | 82.6% |
| 4 | 6 | 0.820 | 0.7048 | 0.390 | 65.6% | 0.376 [0.315, 0.427] | 65.6% |
| 6 | 15 | 0.820 | 0.7048 | 0.390 | 53.1% | 0.387 [0.337, 0.435] | 54.6% |
| 16 | 120 | 0.820 | 0.7048 | 0.390 | 18.5% | 0.382 [0.336, 0.422] | 19.1% |
| 32 | 496 | 0.820 | 0.7048 | 0.390 | 3.4% | 0.383 [0.336, 0.422] | 3.8% |

In theory κ is 0.390 at every panel size. P_o averages over the C(k, 2) pairs
([Fleiss, 1971](#ref-fleiss1971), p. 379), and when
raters are interchangeable each pair added agrees at the same rate as the pairs
already counted. Every rater has the same category rates, so P_e does not
change either. This holds because the raters are alike; a panel whose members
differ will not behave this way. 
The simulated panel stays between 0.375 and 0.387, with 0.390 inside every
interval; its observed P_o ran from 0.821 to 0.826 and P_e from 0.714 to 0.722.
Unanimity falls from 82.0% to 3.4% across the same range. It falls even though
every pair agrees at the same proportion throughout.

Lastly, Conger's kappa ([Conger, 1980](#ref-conger1980)) extends Cohen's kappa as Fleiss' extends Scott's pi, with the chance term built from each rater's own rates instead of the pooled rates; [Light (1971)](#ref-light1971) gives the alternative of averaging the pairwise coefficients.

Later sections will cover Fleiss against Conger when raters use the categories
at different rates, how more raters change the standard error and the floor
beneath it, and rater generalizability.

## When not everyone rates everything

Every table so far has been complete: each rater saw every subject. Many real
designs are not, and not by accident. A systematic review may screen thousands
of abstracts with a handful of reviewers, where double-screening everything
would cost more than the review is worth. The usual compromise is that one
reviewer sees everything and a sample is screened twice, so agreement is
measured on the part that overlaps.

Take a thousand abstracts and three reviewers. A screens all of them, B screens
the first four hundred, C screens a two-hundred block that overlaps B on a
hundred. The result:

| | abstracts |
| --- | --- |
| rated once | 500 |
| rated twice | 400 |
| rated three times | 100 |
| **A and B both saw** | **400** |
| **A and C both saw** | **200** |
| **B and C both saw** | **100** |

The first thing to notice is that **half the corpus cannot contribute to any
agreement coefficient at all.** A subject rated once forms no pair, so there is
nothing for it to agree or disagree with. Five hundred abstracts were screened
and are simply absent from every coefficient. That is not a flaw in any
coefficient; it is what measuring agreement means.

**Krippendorff's alpha** ([Krippendorff, 2013](#ref-krippendorff2013)) was
built for this case, and the reason is structural rather than a feature bolted
on. Alpha is not the kappa family's (P_o − P_e)/(1 − P_e) rearranged: it works
in **pairs of ratings**, not subjects, comparing observed disagreement to the
disagreement expected if the same ratings were shuffled. A subject rated once
contributes no pairs and therefore contributes nothing, which needs no special
handling because nothing about the formula assumed every subject would
contribute equally. The kappa family averages over **subjects**, and a subject
with one rating has no agreement to average.

That is the whole of the difference. Alpha is not more tolerant of missing
data; it is counting something that missing data does not disturb.

## Reproducing these numbers

Every table on this page is in
[`example_data/Statistical_background/`](https://github.com/djarenas/Inter-Rater/tree/master/example_data/Statistical_background/)
as ratings that `interrater` reads like any other dataset. The professors,
for example:

```python
from pathlib import Path
import pandas as pd
from interrater import (RatingsSemantics, CleanerConfig, RatingsCleaner,
                        CohenKappaAnalyzer, CohenKappaConfig)

folder = Path("example_data/Statistical_background")
semantics = RatingsSemantics.build_from_yaml(folder / "pass_fail_semantics.yaml")
cleaner_config = CleanerConfig.build_from_yaml(folder / "complete_cleaner_config.yaml")
frame = pd.read_csv(folder / "professors_data.csv", dtype=str)
ratings = RatingsCleaner(frame, semantics, config=cleaner_config).ratings

result = CohenKappaAnalyzer(ratings, CohenKappaConfig(random_state=42)).estimate()
print(result)
```

```
Cohen's κ = 0.111 [-0.102, 0.389] (95% CI, percentile, n=1000)
  Po = 0.840, Pe = 0.820; κ_max = 1.000 (these raters' marginals)
  100 subjects; unweighted
```

The interval includes zero: on a hundred essays, these two are not
distinguishable from strangers.

To check the whole page at once, from that folder:

```
python generate_tables.py --check   # the CSVs are what the page describes
python reproduce_page.py            # interrater gives the page's numbers
```

## References

<a id="ref-agresti1989"></a>
Agresti, A. (1989). An agreement model with kappa as parameter. *Statistics &
Probability Letters*, 7(4), 271–273.
<https://users.stat.ufl.edu/~aa/articles/SPL_1989.pdf>

<a id="ref-bennett1954"></a>
Bennett, E. M., Alpert, R., & Goldstein, A. C. (1954). Communications through
limited-response questioning. *Public Opinion Quarterly*, 18(3), 303–308.
<https://www.jstor.org/stable/2745989>

<a id="ref-brennan1981"></a>
Brennan, R. L., & Prediger, D. J. (1981). Coefficient kappa: Some uses,
misuses, and alternatives. *Educational and Psychological Measurement*, 41(3),
687–699.
<https://doi.org/10.1177/001316448104100307>

<a id="ref-cicchetti1990"></a>
Cicchetti, D. V., & Feinstein, A. R. (1990). High agreement but low kappa: II.
Resolving the paradoxes. *Journal of Clinical Epidemiology*, 43(6), 551–558.
<https://doi.org/10.1016/0895-4356(90)90159-M>

<a id="ref-cohen1960"></a>
Cohen, J. (1960). A coefficient of agreement for nominal scales. *Educational
and Psychological Measurement*, 20(1), 37–46.
<https://doi.org/10.1177/001316446002000104>

<a id="ref-conger1980"></a>
Conger, A. J. (1980). Integration and generalization of kappas for multiple
raters. *Psychological Bulletin*, 88(2), 322–328.
<https://doi.org/10.1037/0033-2909.88.2.322>

<a id="ref-feinstein1990"></a>
Feinstein, A. R., & Cicchetti, D. V. (1990). High agreement but low kappa: I.
The problems of two paradoxes. *Journal of Clinical Epidemiology*, 43(6),
543–549.
<https://doi.org/10.1016/0895-4356(90)90158-L>

<a id="ref-fleiss1971"></a>
Fleiss, J. L. (1971). Measuring nominal scale agreement among many raters.
*Psychological Bulletin*, 76(5), 378–382.
<https://doi.org/10.1037/h0031619>

<a id="ref-gwet2008"></a>
Gwet, K. L. (2008). Computing inter-rater reliability and its variance in the
presence of high agreement. *British Journal of Mathematical and Statistical
Psychology*, 61, 29–48.
<https://doi.org/10.1348/000711006X126600>

<a id="ref-kraemer1992"></a>
Kraemer, H. C. (1992). Measurement of reliability for categorical data in
medical research. *Statistical Methods in Medical Research*, 1(2), 183–199.
<https://doi.org/10.1177/096228029200100204>

<a id="ref-kraemer2002"></a>
Kraemer, H. C., Periyakoil, V. S., & Noda, A. (2002). Kappa coefficients in
medical research. *Statistics in Medicine*, 21(14), 2109–2129.
<https://doi.org/10.1002/sim.1180>

<a id="ref-krippendorff2013"></a>
Krippendorff, K. (2013). *Computing Krippendorff's alpha-reliability.*
Technical report, University of Pennsylvania.
<https://www.asc.upenn.edu/sites/default/files/2021-03/Computing%20Krippendorff%27s%20Alpha-Reliability.pdf>

<a id="ref-light1971"></a>
Light, R. J. (1971). Measures of response agreement for qualitative data: Some
generalizations and alternatives. *Psychological Bulletin*, 76(5), 365–377.
<https://doi.org/10.1037/h0031643>

<a id="ref-rae1988"></a>
Rae, G. (1988). The equivalence of multiple rater kappa statistics and
intraclass correlation coefficients. *Educational and Psychological
Measurement*, 48(2), 367–374.
<https://doi.org/10.1177/0013164488482009>

<a id="ref-schuster2001"></a>
Schuster, C. (2001). Kappa as a parameter of a symmetry model for rater
agreement. *Journal of Educational and Behavioral Statistics*, 26(3), 331–342.
<https://doi.org/10.3102/10769986026003331>

<a id="ref-scott1955"></a>
Scott, W. A. (1955). Reliability of content analysis: The case of nominal
scale coding. *Public Opinion Quarterly*, 19(3), 321–325.
<https://doi.org/10.1086/266577>
