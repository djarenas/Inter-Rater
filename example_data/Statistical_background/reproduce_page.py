"""
Run interrater on every table of the statistical background page and check
the page's numbers.

``generate_tables.py`` writes the tables; this script loads each one through
``RatingsCleaner``, computes the coefficients the page quotes, and compares
them with the values printed on the page (held in ``PAGE`` below, so a
change to the page is a change here too).

- **Built tables** must match the page at the precision it prints.
- **The simulated biopsy panel** must land within error: the page's kappa
  inside interrater's 95% interval, and P_o, P_e and unanimity within
  ``PANEL_TOLERANCE`` of the page's expected values.

``guide/statistical_background_numbers.py`` checks the same page with the
standard library alone, so the two are independent.

Usage::

    python reproduce_page.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from interrater import (
    BrennanPredigerAnalyzer, BrennanPredigerConfig,
    CleanerConfig,
    CohenKappaAnalyzer, CohenKappaConfig,
    FleissKappaAnalyzer, FleissKappaConfig,
    GwetAC1Analyzer, GwetAC1Config,
    RatingsCleaner, RatingsSemantics,
    ScottsAnalyzer, ScottsConfig,
)
from interrater.intermediates.data_profiles import BalanceProfile, RaterProfile

_ERROR_PREFIX = "[reproduce_page] Error"

HERE = Path(__file__).resolve().parent
RANDOM_STATE = 42

# Only the points are checked on built tables, so a short bootstrap will do;
# 250 is the smallest count the results call enough.
N_RESAMPLES_POINT = 250
# The panel's check uses the interval, so it gets a full bootstrap.
N_RESAMPLES_PANEL = 1000
# How far a simulated proportion may sit from its expected value. With 1,000
# biopsies the standard error of P_o or unanimity is about 0.01 (P_e a
# little more), so 0.02 is about two standard errors.
PANEL_TOLERANCE = 0.02

PANEL_SIZES = (2, 4, 6, 16, 32)

# The page's values, by section. Strings keep the precision the page prints.
PAGE = {
    "professors": {"p_o": "0.84", "p_e": "0.82", "kappa": "0.111",
                   "S": "0.680", "AC1": "0.805"},
    "strict_grader": {"p_o": "0.60", "p_e_individual": "0.50", "kappa": "0.200",
                      "p_e_pooled": "0.58", "pi": "0.048", "gap": "0.08"},
    "scott_quote": {"S_2": "0.20", "S_4": "0.47"},
    "pathologists": {"p_o": "0.900", "p_e": "0.878", "headroom": "0.122",
                     "kappa": "0.177", "S_p_e": "0.500", "S": "0.800",
                     "AC1_p_e": "0.122", "AC1": "0.886",
                     "confirmed": 15, "confirmed_share": "0.231",
                     "healthy_agree": "0.947"},
    # prevalence percent: (p_o, kappa, S, AC1)
    "sweep": {50: (".7925", ".5754", ".5850", ".5941"),
              30: (".8375", ".5925", ".6750", ".7297"),
              20: (".8600", ".5625", ".7200", ".7941"),
              10: (".8825", ".4629", ".7650", ".8496"),
              5:  (".8938", ".3346", ".7875", ".8736"),
              2:  (".9005", ".1814", ".8010", ".8867")},
    "sweep_call_rate_at_2": "0.065",
    "biopsies_2raters": {"p_o": "0.820", "p_e": "0.7048", "kappa": "0.390"},
    # raters: (P_o, P_e, kappa_F, all agree)
    "panel": {2:  ("0.820", "0.7048", "0.390", "0.820"),
              4:  ("0.820", "0.7048", "0.390", "0.656"),
              6:  ("0.820", "0.7048", "0.390", "0.531"),
              16: ("0.820", "0.7048", "0.390", "0.185"),
              32: ("0.820", "0.7048", "0.390", "0.034")},
    "abstracts": {"rated_once": 500, "rated_twice": 400, "rated_three": 100,
                  "A_B": 400, "A_C": 200, "B_C": 100},
}


# Checks -------------------------------------------------------------------- #

class PageChecks:
    """
    I compare what interrater computed with what the page prints, and keep
    every comparison that failed.
    """

    failures: list[str]
    n_checks: int

    def __init__(self, verbose: bool = True) -> None:
        self.failures = []
        self.n_checks = 0
        self._verbose = verbose

    def check_rounded(self, label: str, got: float, page_value: str) -> None:
        """``got`` rounds to ``page_value`` at the decimals the page prints."""
        decimals = len(page_value.split(".")[1]) if "." in page_value else 0
        ok = abs(got - float(page_value)) <= 0.5 * 10 ** -decimals + 1e-12
        self._record(label, ok, f"{got:.{decimals + 2}f}", f"page {page_value}")

    def check_within(self, label: str, got: float, page_value: str,
                     tolerance: float) -> None:
        """``got`` is within ``tolerance`` of ``page_value``."""
        ok = abs(got - float(page_value)) <= tolerance
        self._record(label, ok, f"{got:.4f}", f"page {page_value} +/- {tolerance}")

    def check_inside(self, label: str, page_value: str,
                     lower: float, upper: float) -> None:
        """``page_value`` lies inside [lower, upper]."""
        ok = lower <= float(page_value) <= upper
        self._record(label, ok, f"[{lower:.3f}, {upper:.3f}]", f"page {page_value}")

    def check_equal(self, label: str, got: object, want: object) -> None:
        """``got == want`` exactly."""
        self._record(label, got == want, f"{got}", f"want {want}")

    def _record(self, label: str, ok: bool, got: str, want: str) -> None:
        self.n_checks += 1
        if not ok:
            self.failures.append(f"{label}: got {got}, {want}")
        if self._verbose:
            print(f"   {'ok  ' if ok else 'FAIL'} {label:<44} {got:>18}   {want}")


# Loading ------------------------------------------------------------------- #

def load_ratings(data: str, semantics: str | RatingsSemantics,
                 cleaner: str = "complete_cleaner_config.yaml",
                 columns: list[str] | None = None):
    """
    Read one table through the cleaner. ``columns`` keeps only those rater
    columns (with the id), for a smaller panel. Returns (ratings, report).
    """
    if isinstance(semantics, str):
        semantics = RatingsSemantics.build_from_yaml(HERE / semantics)
    frame = pd.read_csv(HERE / data, dtype=str)
    if columns is not None:
        frame = frame[["subject_id", *columns]]
    config = CleanerConfig.build_from_yaml(HERE / cleaner)
    cleaner_run = RatingsCleaner(frame, semantics, config=config)
    return cleaner_run.ratings, cleaner_run.report


def estimate(analyzer, config_class, ratings, **settings):
    """The analyzer's result at a fixed seed."""
    config = config_class(n_resamples=settings.pop("n_resamples", N_RESAMPLES_POINT),
                          random_state=RANDOM_STATE, **settings)
    return analyzer(ratings, config).estimate()


def crosstab(ratings) -> np.ndarray:
    """The 2x2 table of two raters' codes, in category order."""
    codes = ratings.category_index
    q = ratings.n_categories
    table = np.zeros((q, q), dtype=int)
    for a, b in codes:
        table[a, b] += 1
    return table


# Sections ------------------------------------------------------------------ #

def check_professors(checks: PageChecks) -> None:
    page = PAGE["professors"]
    ratings, _ = load_ratings("professors_data.csv", "pass_fail_semantics.yaml")
    cohen = estimate(CohenKappaAnalyzer, CohenKappaConfig, ratings)
    checks.check_rounded("professors: p_o", cohen.parts.observed, page["p_o"])
    checks.check_rounded("professors: p_e", cohen.parts.chance, page["p_e"])
    checks.check_rounded("professors: Cohen's kappa", cohen.point, page["kappa"])
    bp = estimate(BrennanPredigerAnalyzer, BrennanPredigerConfig, ratings)
    checks.check_rounded("professors: Bennett's S", bp.point, page["S"])
    ac1 = estimate(GwetAC1Analyzer, GwetAC1Config, ratings)
    checks.check_rounded("professors: Gwet's AC1", ac1.point, page["AC1"])


def check_strict_grader(checks: PageChecks) -> None:
    page = PAGE["strict_grader"]
    ratings, _ = load_ratings("strict_grader_data.csv", "pass_fail_semantics.yaml")
    cohen = estimate(CohenKappaAnalyzer, CohenKappaConfig, ratings)
    scott = estimate(ScottsAnalyzer, ScottsConfig, ratings)
    checks.check_rounded("strict grader: p_o", cohen.parts.observed, page["p_o"])
    checks.check_rounded("strict grader: p_e individual", cohen.parts.chance,
                         page["p_e_individual"])
    checks.check_rounded("strict grader: Cohen's kappa", cohen.point, page["kappa"])
    checks.check_rounded("strict grader: p_e pooled", scott.parts.chance,
                         page["p_e_pooled"])
    checks.check_rounded("strict grader: Scott's pi", scott.point, page["pi"])
    checks.check_rounded("strict grader: pooled - individual",
                         scott.parts.chance - cohen.parts.chance, page["gap"])

    # Scott's whimsical researcher: same essays, two unused categories added
    quote = PAGE["scott_quote"]
    s_2 = estimate(BrennanPredigerAnalyzer, BrennanPredigerConfig, ratings)
    ratings_4, _ = load_ratings("strict_grader_data.csv",
                                "pass_fail_4categories_semantics.yaml",
                                cleaner="unused_categories_cleaner_config.yaml")
    s_4 = estimate(BrennanPredigerAnalyzer, BrennanPredigerConfig, ratings_4)
    checks.check_rounded("Scott's quote: S with 2 categories", s_2.point, quote["S_2"])
    checks.check_rounded("Scott's quote: S with 4 categories", s_4.point, quote["S_4"])


def check_pathologists(checks: PageChecks) -> None:
    page = PAGE["pathologists"]
    ratings, _ = load_ratings("pathologists_data.csv", "slides_semantics.yaml")
    cohen = estimate(CohenKappaAnalyzer, CohenKappaConfig, ratings)
    checks.check_rounded("pathologists: p_o", cohen.parts.observed, page["p_o"])
    checks.check_rounded("pathologists: p_e", cohen.parts.chance, page["p_e"])
    checks.check_rounded("pathologists: headroom 1 - p_e", 1 - cohen.parts.chance,
                         page["headroom"])
    checks.check_rounded("pathologists: Cohen's kappa", cohen.point, page["kappa"])
    bp = estimate(BrennanPredigerAnalyzer, BrennanPredigerConfig, ratings)
    checks.check_rounded("pathologists: Bennett's S p_e", bp.parts.chance, page["S_p_e"])
    checks.check_rounded("pathologists: Bennett's S", bp.point, page["S"])
    ac1 = estimate(GwetAC1Analyzer, GwetAC1Config, ratings)
    checks.check_rounded("pathologists: AC1 p_e", ac1.parts.chance, page["AC1_p_e"])
    checks.check_rounded("pathologists: Gwet's AC1", ac1.point, page["AC1"])

    table = crosstab(ratings)                 # rows A, columns B; 0 = cancer
    checks.check_equal("pathologists: A's cancers B confirmed", int(table[0, 0]),
                       page["confirmed"])
    checks.check_rounded("pathologists: share confirmed",
                         table[0, 0] / table[0].sum(), page["confirmed_share"])
    checks.check_rounded("pathologists: agreement on A's healthy",
                         table[1, 1] / table[1].sum(), page["healthy_agree"])


def check_sweep(checks: PageChecks) -> None:
    for percent, (p_o, kappa, s, ac1) in PAGE["sweep"].items():
        ratings, _ = load_ratings(f"prevalence_{percent:02d}_data.csv",
                                  "slides_semantics.yaml")
        cohen = estimate(CohenKappaAnalyzer, CohenKappaConfig, ratings)
        bp = estimate(BrennanPredigerAnalyzer, BrennanPredigerConfig, ratings)
        gwet = estimate(GwetAC1Analyzer, GwetAC1Config, ratings)
        checks.check_rounded(f"sweep {percent}%: p_o", cohen.parts.observed, p_o)
        checks.check_rounded(f"sweep {percent}%: kappa", cohen.point, kappa)
        checks.check_rounded(f"sweep {percent}%: S", bp.point, s)
        checks.check_rounded(f"sweep {percent}%: AC1", gwet.point, ac1)
        if percent == 2:
            calls = crosstab(ratings)[0].sum() / ratings.n_subjects
            checks.check_rounded("sweep 2%: A's cancer call rate", calls,
                                 PAGE["sweep_call_rate_at_2"])


def check_biopsies(checks: PageChecks) -> None:
    page = PAGE["biopsies_2raters"]
    ratings, _ = load_ratings("biopsies_2raters_data.csv", "slides_semantics.yaml")
    cohen = estimate(CohenKappaAnalyzer, CohenKappaConfig, ratings)
    scott = estimate(ScottsAnalyzer, ScottsConfig, ratings)
    checks.check_rounded("biopsies, 2 raters: p_o", cohen.parts.observed, page["p_o"])
    checks.check_rounded("biopsies, 2 raters: p_e", cohen.parts.chance, page["p_e"])
    checks.check_rounded("biopsies, 2 raters: Cohen's kappa", cohen.point, page["kappa"])
    checks.check_rounded("biopsies, 2 raters: Scott's pi", scott.point, page["kappa"])


def check_panel(checks: PageChecks) -> None:
    full = RatingsSemantics.build_from_yaml(HERE / "biopsy_panel_semantics.yaml")
    for k in PANEL_SIZES:
        p_o, p_e, kappa, unanimity = PAGE["panel"][k]
        readers = list(full.raters[:k])
        semantics = RatingsSemantics(raters=readers, categories=full.categories,
                                     n_subjects=full.n_subjects)
        ratings, _ = load_ratings("biopsy_panel_data.csv", semantics, columns=readers)
        fleiss = estimate(FleissKappaAnalyzer, FleissKappaConfig, ratings,
                          n_resamples=N_RESAMPLES_PANEL)
        interval = fleiss.confidence_interval
        counts = ratings.to_nij_matrix().matrix
        agree_all = float(np.mean(counts.max(axis=1) == k))
        checks.check_within(f"panel k={k}: P_o", fleiss.parts.observed, p_o,
                            PANEL_TOLERANCE)
        checks.check_within(f"panel k={k}: P_e", fleiss.parts.chance, p_e,
                            PANEL_TOLERANCE)
        checks.check_inside(f"panel k={k}: page kappa in 95% CI", kappa,
                            interval.lower, interval.upper)
        checks.check_within(f"panel k={k}: all agree", agree_all, unanimity,
                            PANEL_TOLERANCE)
        if k == 2:
            # The page: Fleiss' kappa is Scott's pi for two raters
            scott = estimate(ScottsAnalyzer, ScottsConfig, ratings)
            checks.check_equal("panel k=2: Fleiss equals Scott's pi",
                               bool(np.isclose(fleiss.point, scott.point,
                                               rtol=0, atol=1e-12)), True)


def check_abstracts(checks: PageChecks) -> None:
    page = PAGE["abstracts"]
    ratings, report = load_ratings("abstracts_data.csv", "abstracts_semantics.yaml",
                                   cleaner="abstracts_cleaner_config.yaml")
    balance = BalanceProfile.build_from_ratings(ratings)
    overlap = RaterProfile.build_from_ratings(ratings).rater_overlap
    checks.check_equal("abstracts: dropped by the cleaner (rated once)",
                       report.n_subjects_lost, page["rated_once"])
    checks.check_equal("abstracts: rated twice",
                       balance.subject_counts.get(2, 0), page["rated_twice"])
    checks.check_equal("abstracts: rated three times",
                       balance.subject_counts.get(3, 0), page["rated_three"])
    checks.check_equal("abstracts: A and B both saw", overlap.value_for(0, 1), page["A_B"])
    checks.check_equal("abstracts: A and C both saw", overlap.value_for(0, 2), page["A_C"])
    checks.check_equal("abstracts: B and C both saw", overlap.value_for(1, 2), page["B_C"])


SECTIONS = (
    ("Percent agreement and Cohen's kappa", check_professors),
    ("Pooled or individual rates", check_strict_grader),
    ("When the condition is rare", check_pathologists),
    ("When the condition is rare: prevalence sweep", check_sweep),
    ("More than two raters: two readers", check_biopsies),
    ("More than two raters: the panel (simulated, seed 42)", check_panel),
    ("When not everyone rates everything (simulated, seed 42)", check_abstracts),
)


def run(verbose: bool = True) -> PageChecks:
    """Run every section; the returned checks hold any failures."""
    checks = PageChecks(verbose=verbose)
    for title, section in SECTIONS:
        if verbose:
            print(f"\n{title}")
        section(checks)
    return checks


def main() -> int:
    checks = run()
    print(f"\n{checks.n_checks - len(checks.failures)} of {checks.n_checks} checks passed")
    for line in checks.failures:
        print(f"FAIL {line}")
    return 1 if checks.failures else 0


if __name__ == "__main__":
    sys.exit(main())
