# interrater

Agreement between raters, for any rating task where there is no ground truth to
check against. For example: systematic review screening, content coding, annotation
quality, evaluation panels, etc.

- **Install:** `pip install interrater` · [PyPI: interrater](https://pypi.org/project/interrater/)
- **Source:** [GitHub: interrater](https://github.com/djarenas/Inter-Rater)
- **Paper:** [arXiv:1809.05731](https://arxiv.org/abs/1809.05731) (describes the original version)

Version 2 is a complete rewrite in object-oriented Python; see
[What version 2 adds](https://djarenas.github.io/packages/interrater/roadmap/).

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

# Download example_data/ from GitHub first (link above),
# and change this path to where you saved it.
data = Path("example_data/SR1_OSF2018")

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
makes deliberately. [Why objects?](https://djarenas.github.io/packages/interrater/why_objects/)

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
[What if one of your raters were someone else?](https://djarenas.github.io/packages/interrater/rater_substitution/)

Worked configurations for ten datasets, five of them real systematic reviews,
are in [`example_data/`](https://github.com/djarenas/Inter-Rater/tree/master/example_data/), each with the reference values it
should reproduce.

## Background

`interrater` is the successor to
[Inter-Rater v1](https://github.com/djarenas/Inter-Rater/tree/v1.5), described in
[arXiv:1809.05731](https://arxiv.org/abs/1809.05731) (2018) and archived at
[Zenodo](https://doi.org/10.5281/zenodo.1227660). Version 1 has been used in
published systematic reviews and is cited in peer-reviewed research. Version 2
is a complete rewrite in object-oriented Python, undertaken to fix maintenance
issues and limitations that emerged from real-world use. Version 1 remains
available at git tag `v1.5`.

See [Why pairwise?](https://djarenas.github.io/packages/interrater/why_pairwise/) for why the package reports every pair of raters.

See [Why objects?](https://djarenas.github.io/packages/interrater/why_objects/) for why and what of the OOP design choices.

See [Statistical background](https://djarenas.github.io/packages/interrater/statistical_background/) for sections on chance-definitions; resampling choice; correlated and independent errors.

See [What if one of your raters were someone else?](https://djarenas.github.io/packages/interrater/rater_substitution/) for rater injection.

See [roadmap](https://djarenas.github.io/packages/interrater/roadmap/) for info on what V2 added and what is coming in the future.


## How to cite

If you use this software, please cite both the software and the paper
describing the method it implements. Citation details are in [`CITATION.cff`](https://github.com/djarenas/Inter-Rater/blob/master/CITATION.cff) (on GitHub, the "Cite this repository" button in the sidebar formats it for you).

Arenas, D. J. (2018). Software for analysis of inter-rater reliability by
permutating pairs of multiple raters.
[arXiv:1809.05731](https://arxiv.org/abs/1809.05731).
Version 1 archived at [Zenodo](https://doi.org/10.5281/zenodo.1227660).

## License

GNU General Public License v3 or later. Same license as version 1.
