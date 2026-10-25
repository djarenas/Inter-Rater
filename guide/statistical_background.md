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

Although nothing connected their judgments, they scored eighty-four per cent agreement, which may appear like good agreement, but it isn't. If two graders each pass 90% and their criteria are unrelated, then on any given essay the chance they both pass it is 0.9 × 0.9 = 0.81, and the chance
they both fail it is 0.1 × 0.1 = 0.01. Two strangers would land on **0.82** agreement.

The professors managed 0.84. Out of a hundred essays, two of them represent
agreement that chance does not already explain.

The example is built to be extreme, deliberately. Making the two sets of
criteria independent is the only way to know in advance roughly where the
answer should land, and a coefficient that calls 0.84 strong agreement here is
not measuring agreement. Real raters are rarely that independent — but they are
rarely identical either, and whatever part of their criteria fails to overlap
contributes agreement of exactly this kind. The grammar/creativity professors
are the far end of a line that every real panel falls somewhere along.

## Cohen's kappa

Cohen's kappa is one of the many important agreement
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

Cohen's kappa normalizes for it:

κ = (p_o − p_e) / (1 − p_e)

The numerator is the agreement chance does not account for — the two essays.
The denominator is how much such agreement was available in the first place:
1 − p_e is the headroom above chance. So κ reports the share of that headroom
the raters actually took. One means they took all of it, zero means none of it,
and negative values mean they did worse than two strangers.

For the grammar/creativity professors:

κ = (0.84 − 0.82) / (1 − 0.82) = 0.02 / 0.18 = **0.111**

Eighty-four per cent agreement, and eleven per cent of the possible gains.

Note: The denominator is worth a second look, because it is small here. Both
professors pass 90%, so chance alone delivers 0.82 and only 0.18 is left to
earn. Suppose they had shared their criteria after all and agreed on 95 of the
100 essays:

κ = (0.95 − 0.82) / 0.18 = **0.722**

Even agreeing on 98 would give 0.889. High rates leave little room, and Cohen's kappa
measures what happens inside that room. That becomes the whole story when the
thing being rated is rare, as discussed below.

## Pooled or individual rates

Both professors so far have used the categories at the same percentages — each passed
90% — which is why one chance term has been enough. When two raters use the
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
Scott's Pi is the two-rater kappa that uses the pooled rate. This choice treats
the two raters as one distribution. Between them they
passed 70% of the time, so p_e = 0.70 × 0.70 + 0.30 × 0.30 = 0.58. 

p_e = Σ_j p̄_j²,  where p̄_j = (p_j(A) + p_j(B)) / 2

Scott's pi is then:

π = (0.60 − 0.58) / 0.42 = **0.048**

Same hundred essays, same sixty agreements, and the two differ by a factor of
four.

The gap between the two coefficients is driven by
the difference in their chance terms, for two raters:

p_e(pooled) − p_e(individual) = ¼ Σ_j (p_j(A) − p_j(B))²

The further apart the two raters' rates, the wider the gap between the coefficients. 

The grammar/creativity professors both passed 90%, so pooled and individual coincide and both report 0.111. **The gap between the two coefficients sees threshold differences and is
blind to criteria differences.** 

**Krippendorff's alpha** belongs to the pooled family. Details in a future section.

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
way, and both call cancer 6.5% of the time — so there is no threshold
difference here, and pooled and individual chance terms coincide. This is not
the problem from the last section.

The chance term is 0.065² + 0.935² = 0.878. Chance alone would put these two on
the same side of 878 slides. They managed 900. They are twenty-two slides ahead
of chance, out of a possible hundred and twenty-two:

κ = (0.900 − 0.878) / 0.122 = **0.177**

Ninety per cent agreement, and a coefficient that reads as barely better than
guessing. The headroom collapsed. When one category holds 93.5% of the calls,
two raters land together most of the time without trying, and little is left
for skill to claim.

Some coefficients use different terms for chance and may be useful in scenarios like this.

**Bennett's S** refuses to use the raters' rates. With q categories, chance is
1/q — two raters picking blindly between cancer and healthy would agree half
the time, whatever the prevalence.

p_e = 1/q = 0.500,  S = (0.900 − 0.500) / 0.500 = **0.800**

**Gwet's AC1** takes a middle route: chance agreement should be high when one
category dominates, but not as high as the raters' own rates imply, because
those rates are partly produced by the rarity itself.

p_e = (1/(q − 1)) Σ_k π_k(1 − π_k) = 2 × 0.065 × 0.935 = 0.122

AC1 = (0.900 − 0.122) / 0.878 = **0.886**

One set of slides and three different coefficients 0.177, 0.800, 0.886. This is not an artifact of one table. Hold the pathologists' skill fixed — each catching 80% of cancers and correctly clearing 95% of healthy slides — and move only the share of slides that are cancerous:

| cancer prevalence | p_o | κ | S | AC1 |
| --- | --- | --- | --- | --- |
| 50% | .792 | .575 | .585 | .594 |
| 30% | .837 | .592 | .675 | .730 |
| 20% | .860 | .562 | .720 | .794 |
| 10% | .883 | .463 | .765 | .850 |
| 5% | .894 | .335 | .787 | .874 |
| 2% | .900 | .181 | .801 | .887 |

Nothing about the pathologists changes down that column. At a balanced caseload
the three nearly coincide; the fork opens as the condition gets rarer.

Before reaching for the ones that go up, look at what they do to the
grammar/creativity professors — the two who agreed about nothing:

|  | κ | S | AC1 |
| --- | --- | --- | --- |
| grammar/creativity professors | 0.111 | **0.680** | **0.805** |

S and AC1 are generous for the same reason in both cases. Both professors
passed 90%, so a chance term that ignores the rates declines to charge them for
it. In the pathology case that reticence is defensible: the 6.5% was truly
decided by the pathologists, reading the tissue in front of them. In the
grading case it is not. That 90% came from a department curve, "you must pass
around 90%," and was never theirs to set. **Nothing in the data distinguishes
the two situations.** 

None of which means Cohen's kappa is malfunctioning. Cohen's kappa reports how much of their agreement
survives what the caseload hands them for free. AC1 reports whether they behaved
consistently. Both are true, and the question is which one you asked — yet
another example of the coefficient you need depending on your situation. 

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
p_e = 0.18² + 0.82² = 0.7048, and Cohen's kappa and Scott's Pi are 0.390.

Fleiss' kappa is Scott's pi for more than two raters. For N subjects, k raters, the percent agreement is:

P_o = (1/N) Σ_i [ Σ_j n_ij (n_ij − 1) ] / [ k (k − 1) ]

where the n_ij matrix specifies how many raters assign subject i to category j.

Fleiss' kappa uses pooled values, p̄ the proportion of all ratings falling in the first category,

P_e = p̄² + (1 − p̄)²

κ_F = (P_o − P_e) / (1 − P_e)

To the example above, we can add raters with specificity and sensitivity that are both 90%: 

| raters | pairs | P_o | P_e | κ_F | all agree |
| --- | --- | --- | --- | --- | --- |
| 2 | 1 | 0.820 | 0.7048 | 0.390 | 82.0% |
| 4 | 6 | 0.820 | 0.7048 | 0.390 | 65.6% |
| 6 | 15 | 0.820 | 0.7048 | 0.390 | 53.1% |
| 16 | 120 | 0.820 | 0.7048 | 0.390 | 18.5% |
| 32 | 496 | 0.820 | 0.7048 | 0.390 | 3.4% |

κ is 0.390 at every panel size. P_o averages over the C(k, 2) pairs, and when
raters are interchangeable each pair added agrees at the same rate as the pairs
already counted. Every rater has the same category rates, so P_e does not
change either. This holds because the raters are alike; a panel whose members
differ will not behave this way. 
Unanimity falls from 82.0% to 3.4% across the same range. It falls even though
every pair agrees at the same proportion throughout.

Other sections will talk about how more raters
change the standard error and generalizability. Lastly, Conger's kappa extends Cohen's kappa as Fleiss' extends Scott's pi, with the chance term built from each rater's own rates instead of the pooled rates.