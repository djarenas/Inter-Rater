# Contributing to interrater

Thank you for taking an interest. Bug reports, questions, documentation fixes
and code are all welcome.

**A note from the maintainer.** I was slow to respond during the v1 years, and
I promise to do better. If something is unclear or broken, open an issue and
let's figure it out together. As of 2026, I am also happy to help with
academic work: if you are preparing a paper and need a hand running your
ratings through the package, ask.

> **WARNING: do not post confidential or identifiable data in a public
> issue.** Issues on GitHub are visible to everyone. Describe the data
> instead, or write to interrater.python@gmail.com, and we will find a safe
> way to share it.

## Questions and bug reports

Open an issue at <https://github.com/djarenas/Inter-Rater/issues>.

For a bug, please include:

- the `interrater` version (`pip show interrater`) and your Python version
- the smallest example that shows the problem: the semantics, the config, and
  a few rows of data are usually enough
- the full error message, or the number you got and the number you expected

If a result disagrees with a published value or with another package, say
which one and where the reference value comes from. Those reports are the most
useful kind.

## Proposing a change

1. Open an issue first for anything larger than a typo, so the approach can be
   agreed before you write code.
2. Fork the repository and create a branch from `main`.
3. Make the change, with tests (see below).
4. Add an entry to `CHANGELOG.md` under an "Unreleased" heading. The file
   follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). A change
   that alters a number or breaks existing code or configuration goes under
   "Breaking" and says what changed.
5. Open a pull request against `main` and describe what the change does and
   why.

Pull requests are reviewed and merged by the maintainer, Daniel J. Arenas.

## Setting up

The package needs Python 3.10 or later.

```
git clone https://github.com/djarenas/Inter-Rater.git
cd Inter-Rater
pip install -e ".[dev]"
```

## Running the tests

```
pytest -q
```

The same suite runs on GitHub for every push and pull request, on Python 3.10
with the oldest dependency versions `pyproject.toml` allows and on 3.11 to
3.13 with the newest. A pull request needs all four to pass.

`tests/README.md` describes how the suite is organised. It is built bottom-up,
so a failure points at the layer that broke.

## What a change needs

These are the rules the package is built on. A contribution that follows them
is much easier to accept.

- **A new coefficient or interval method is checked against a published
  value.** Cite the source in the test.
- **Preconditions are checked when an object is built.** Data or settings a
  coefficient is not defined for are refused with a message that names the
  class, the field and the alternatives.
- **A setting that does not apply is refused with a reason.** It is never
  silently ignored.
- **A run that cannot produce a number says so.** Return a `NotComputed` with
  the reason. Do not return zero or NaN.
- **The point estimate and its resamples use the same function and the same
  subjects.**
- **`Ratings` never changes.** An operation that would change it returns a new
  one.
- **Every new behaviour has a test, and every new refusal has a test that
  triggers it.**

The reasons behind these rules are in `guide/why_objects.md`.

## Documentation

The guide pages live in `guide/`. Numbers quoted on the statistical background
page are reproduced by
`example_data/Statistical_background/reproduce_page.py`; if you change a
number on the page, update that script so it still passes.

## Credit and authorship

Every contribution is credited, in the changelog and in the list of
contributors. People who make a substantial contribution, such as
implementing and validating a method, a significant feature, or sustained
maintenance, may be invited to be listed as authors of the package in
`CITATION.cff`. Authorship of papers about the package is decided paper by
paper.

Authorship is not automatic. It is agreed in advance, ideally when a
substantial contribution starts, not when a paper is being written. If you
are planning something large and want to discuss credit, please raise it
then. Requests for authorship that are not tied to a contribution will be
closed without reply.

## License

By contributing you agree that your contribution is licensed under
GPL-3.0-or-later, the license of the package.
