# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versioning follows [semver](https://semver.org/).

## [2.1.0] - 2026-10-18

### Breaking

Three of these break code that ran on 2.0.0. The minor version does not warn
you, so this section does.

- `result.confidence_interval` is a `CIResult`, not a tuple. Use `.lower` and
  `.upper`. `to_dict()` still emits `ci_lower` and `ci_upper`.
- `bca_ci()` takes a fourth argument, `jackknife_estimates`.
- `percentile_ci()` and `bca_ci()` return a `BootstrapCIResult`, not a tuple.
- `scipy_bca_ci()` is removed. See `analyzers/README.md`.
- `interrater.intermediates.analysis_results` is now
  `interrater.intermediates.results`. Top-level imports are unchanged.
- `intermediates/ci_methods.py` is removed; the module is
  `analyzers/ci_methods.py`.
- Subset permutation intervals resample subjects rated by at least two of the
  k raters, not only those all k rated. Numbers change on partially crossed
  data; balanced data is unaffected. See
  `intermediates/results/permutation_results/README.md`.

### Added

- `BootstrapCIResult` and `ParametricCIResult`: bounds plus `method_used`,
  `n_usable`, `degraded` and `summary()`.
- `recommended_resamples` — 250 for percentile, 1000 for BCa — and
  `below_recommended_resamples`, compared against the draws that survived.
- A `"ci"` key in `to_dict()`, holding the interval's `summary()`.
- A `README.md` in each source subpackage.

### Changed

- `ci_method="parametric_null"` warns twice, from the config and from the
  interval, and both say the width is scaled as though the estimate were zero.
- The resample-count notice no longer reads "B below 1000 recommended", which
  parsed as a recommendation to stay below 1000.
- Cohen's and Conger's `balanced_data_only` warnings state what was set.

### Fixed

- A subset with fewer than two measurable subjects returns
  `method_used = None` rather than a tuple of NaNs.

## [2.0.0] - 2026-10-11

Complete rewrite in object-oriented Python. This changelog begins at 2.1.0;
for what 2.0.0 contained, see
[What version 2 adds](README.md#what-version-2-adds).

## [1.5]

Final release of [Inter-Rater v1](https://github.com/djarenas/Inter-Rater),
described in [arXiv:1809.05731](https://arxiv.org/abs/1809.05731) and archived
at [Zenodo](https://doi.org/10.5281/zenodo.1227660).
