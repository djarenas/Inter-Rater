# AI usage disclosure

Generative AI was used in developing versions 2.0.0, 2.1.0 and 2.2.0 of
`interrater`. The tool was Claude (Anthropic), used as a coding and editing
assistant. Several models were used at different points, including Claude
Fable 5.2 and Claude Opus 5.5.

## Where it was applied

Source code, tests and documentation.

## Nature of the assistance

- **Code.** The object-oriented design, the inheritance hierarchies and the
  use of configuration files are the author's ideas. The coding work took two
  forms. In some cases the author wrote a class
  and the AI edited or refactored it. In others the author specified the
  design of a class (its attributes, methods, helper methods and constructors)
  and the AI implemented it to that specification. The AI also ran
  consistency checks across the package before release.
- **Tests.** This is where the AI was used most heavily. The AI wrote most of
  the test code at the author's direction. The reference values the tests
  check against come from published sources and other packages, not from the
  AI.
- **Documentation.** The work took three forms. In some cases the author
  wrote a draft and the AI edited it. In others the author gave the content
  as bullet points and the AI wrote it up. In others the AI wrote a first
  draft at the author's direction and the author revised it. On the guide
  pages the AI also helped locate and format references,
  check quotations against the source papers, and check arithmetic.
- **Working documents.** The AI helped maintain the project's intermediate
  documents, such as a literature review built up one reference at a time as
  the author obtained the papers. Each entry records what the paper says,
  with quotations checked against the source.

## What is the author's

The author made the architectural, design and statistical decisions. These
include the design of the package (declared semantics, validated objects,
inheritance hierarchies, configuration files, results that record how they
were made), which coefficients and methods to
implement, the pairwise and subset permutation analyses, and how each
coefficient behaves on unbalanced or missing data. The pairwise approach
dates to the original release, whose title names it: "Inter-Rater: Software
for analysis of inter-rater reliability by permutating pairs of multiple
users" (2018, https://doi.org/10.48550/arXiv.1809.05731). Version 1.4
followed in 2020 (https://doi.org/10.5281/zenodo.3630832).

## Human oversight

To the best of the author's knowledge, each change was proposed first and
made only after the author approved it.
The work was iterative: the author checked the code, sent it back with
corrections, and checked it again, and several rounds of this were common
throughout. The author reviewed, edited and validated the AI-generated code
and documentation, and checked every class. As in any software, errors may
remain; reports are welcome (see `CONTRIBUTING.md`).

## Verification

Correctness does not rest on the AI. The test suite (2,297 tests as of 2.1.0)
checks the coefficients and intervals against published worked examples and
reference values, such as Efron & Tibshirani (1986) for the BCa interval and
irrCAC for point estimates, and reproduces the numbers on the statistical
background page. The suite runs on every push.

The AI also wrote the scripts in `analyses/external_validation/`, which
compare the package's point estimates and null standard errors with five
other implementations: irrCAC and irr (R), and statsmodels, scikit-learn and
krippendorff (Python), on the shipped datasets. The AI installed the tools
and ran them; the values come from those packages, not from the AI.

## Not used for

Correspondence with journal editors or reviewers.

## License

The package is distributed under GPL-3.0-or-later. See `LICENSE` for the
terms, including the disclaimer of warranty.
