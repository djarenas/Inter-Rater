#!/usr/bin/env python3
"""
Reproduces every number in the `statistical_background` docs page.

Standard library only (fractions, itertools, random, statistics) so it runs
anywhere with no install step.  Exact arithmetic via `fractions.Fraction`
wherever the quantity is exact; Monte Carlo only where the page reports a
sampling property.

Run:  python3 statistical_background_numbers.py
Every `check()` is an identity the `interrater` package should also satisfy.
"""

from fractions import Fraction as F
from itertools import combinations, product
import random
import statistics as st

FAILURES = []


def check(label, got, want, tol=0):
    ok = (got == want) if tol == 0 else (abs(float(got) - float(want)) <= tol)
    if not ok:
        FAILURES.append(f"{label}: got {got} want {want}")
    print(f"   {'ok  ' if ok else 'FAIL'} {label}")


# ---------------------------------------------------------------- coefficients

def rates(rows, cats):
    """Per-rater category rates. rows = subjects x raters."""
    n = len(rows)
    k = len(rows[0])
    return [[F(sum(1 for s in rows if s[r] == c), n) for c in cats] for r in range(k)]


def pooled(p, cats):
    k = len(p)
    return [sum(p[r][j] for r in range(k)) / k for j in range(len(cats))]


def p_o_pairs(rows, cats):
    """Fleiss' observed agreement: mean over subjects of the fraction of
    agreeing rater pairs. Equals mean pairwise p_o."""
    n, k = len(rows), len(rows[0])
    npair = F(k * (k - 1), 2)
    tot = F(0)
    for s in rows:
        tot += sum(F(m * (m - 1), 2) for m in (s.count(c) for c in cats)) / npair
    return tot / n


def pe_pooled(p, cats):
    """Scott / Fleiss chance term: pooled marginals."""
    return sum(x * x for x in pooled(p, cats))


def pe_individual(p, cats):
    """Cohen / Conger chance term: mean of the pairwise individual terms."""
    k = len(p)
    pairs = list(combinations(range(k), 2))
    return sum(sum(p[r][j] * p[s][j] for j in range(len(cats))) for r, s in pairs) / len(pairs)


def pe_bennett(cats):
    """Bennett's S / Brennan-Prediger: 1/q, ignores the observed rates."""
    return F(1, len(cats))


def pe_gwet(p, cats):
    """Gwet's AC1."""
    q = len(cats)
    pi = pooled(p, cats)
    return sum(x * (1 - x) for x in pi) / (q - 1)


def coef(p_o, p_e):
    return (p_o - p_e) / (1 - p_e)


def var_across_raters(p, cats):
    """Unbiased across-rater variance of the usage rate, per category."""
    k = len(p)
    pb = pooled(p, cats)
    return [sum((p[r][j] - pb[j]) ** 2 for r in range(k)) / (k - 1) for j in range(len(cats))]


def krippendorff_alpha(rows, cats):
    """Nominal alpha, k coders, complete data, from the coincidence matrix.

    Exact identity (any k, any q, complete data):
        alpha = pi + (1 - pi) / (N * k)
    where pi is the pooled-marginal coefficient (Scott at k=2, Fleiss above).
    The 1/(2N) of the two-rater case is just k=2 here.
    """
    nsub, k = len(rows), len(rows[0])
    n = nsub * k
    o = {(c, d): F(0) for c in cats for d in cats}
    for s in rows:
        cnt = {c: s.count(c) for c in cats}
        for c in cats:
            for d in cats:
                num = cnt[c] * (cnt[c] - 1) if c == d else cnt[c] * cnt[d]
                o[(c, d)] += F(num, k - 1)
    nc = {c: sum(o[(c, d)] for d in cats) for c in cats}
    assert sum(nc.values()) == n
    D_o = sum(o[(c, d)] for c in cats for d in cats if c != d) / n
    D_e = sum(nc[c] * nc[d] for c in cats for d in cats if c != d) / (n * (n - 1))
    return 1 - D_o / D_e


def from_2x2(a, b, c, d):
    """rows from a 2x2 count table: a=both 1, b=A 1 B 0, c=A 0 B 1, d=both 0."""
    return [[1, 1]] * a + [[1, 0]] * b + [[0, 1]] * c + [[0, 0]] * d


def report(name, rows, cats=(1, 0), show=("p_o", "pooled", "indiv")):
    p = rates(rows, cats)
    p_o = p_o_pairs(rows, cats)
    out = {"p_o": p_o, "p": p,
           "pe_pooled": pe_pooled(p, cats), "pe_indiv": pe_individual(p, cats),
           "pe_S": pe_bennett(cats), "pe_AC1": pe_gwet(p, cats)}
    out["pooled"] = coef(p_o, out["pe_pooled"])
    out["indiv"] = coef(p_o, out["pe_indiv"])
    out["S"] = coef(p_o, out["pe_S"])
    out["AC1"] = coef(p_o, out["pe_AC1"])
    print(f"\n{name}")
    print(f"   rater rates      {[[str(x) for x in pr] for pr in p]}")
    print(f"   p_o              {float(p_o):.4f}")
    for key, lab in (("pe_pooled", "p_e pooled    "), ("pe_indiv", "p_e individual"),
                     ("pe_S", "p_e Bennett S "), ("pe_AC1", "p_e Gwet AC1  ")):
        print(f"   {lab}   {float(out[key]):.5f}")
    for key, lab in (("pooled", "Scott/Fleiss  "), ("indiv", "Cohen/Conger  "),
                     ("S", "Bennett S     "), ("AC1", "Gwet AC1      ")):
        print(f"   {lab}   {float(out[key]):+.4f}   ({out[key]})")
    return out


# =============================================================== sections

print("=" * 72)
print("S1-S2  grammar/creativity professors, 100 essays")
print("=" * 72)
profs = from_2x2(82, 8, 8, 2)
r1 = report("82 / 8 / 8 / 2", profs)
check("p_o = 0.84", r1["p_o"], F(21, 25))
check("chance term = 0.82", r1["pe_indiv"], F(41, 50))
check("Cohen kappa = 1/9 = 0.111", r1["indiv"], F(1, 9))
check("equal marginals -> pooled == individual", r1["pe_pooled"], r1["pe_indiv"])
check("S2 counterfactual: 95/100 agreed -> 0.722",
      coef(F(95, 100), F(41, 50)), F(13, 18))
check("S2 counterfactual: 98/100 agreed -> 0.889",
      coef(F(98, 100), F(41, 50)), F(8, 9))
check("S4 callback: professors under Bennett S = 0.680", r1["S"], F(17, 25))
check("S4 callback: professors under Gwet AC1 = 33/41 = 0.805", r1["AC1"], F(33, 41))

print("\n" + "=" * 72)
print("S3  the strict grader, 100 essays")
print("=" * 72)
strict = from_2x2(50, 0, 40, 10)
r3 = report("50 / 0 / 40 / 10", strict)
check("p_o = 0.60", r3["p_o"], F(3, 5))
check("Cohen p_e = 0.50", r3["pe_indiv"], F(1, 2))
check("Cohen kappa = 0.200", r3["indiv"], F(1, 5))
check("Scott p_e = 0.58", r3["pe_pooled"], F(29, 50))
check("Scott pi = 0.0476", r3["pooled"], F(1, 21))
check("coefficient gap = 0.152", r3["indiv"] - r3["pooled"], F(1, 5) - F(1, 21))

pA, pB = r3["p"][0], r3["p"][1]
quarter = sum((pA[j] - pB[j]) ** 2 for j in range(2)) / 4
check("S3 identity  pooled - individual = (1/4) sum (p_A - p_B)^2",
      r3["pe_pooled"] - r3["pe_indiv"], quarter)
check("   ... and it equals 0.08", quarter, F(2, 25))

alpha = krippendorff_alpha(strict, (1, 0))
print(f"\n   Krippendorff alpha = {float(alpha):.6f}   ({alpha})")
pi = r3["pooled"]
check("S3: alpha = pi + (1-pi)/(2N), the k=2 case   [page says 1/(2N-1)]",
      alpha, pi + (1 - pi) / 200)
print(f"   -> correction is (1-pi)/2N = {float((1-pi)/200):.6f}, moving "
      f"{float(pi):.4f} to {float(alpha):.4f}")

print("\n" + "=" * 72)
print("S4  rare cancer, 1000 slides")
print("=" * 72)
path = from_2x2(15, 50, 50, 885)
r4 = report("15 / 50 / 50 / 885", path)
check("p_o = 0.900", r4["p_o"], F(9, 10))
check("chance term = 0.87845", r4["pe_indiv"], F(87845, 100000))
check("equal marginals -> pooled == individual", r4["pe_pooled"], r4["pe_indiv"])
check("Cohen kappa = 431/2431 = 0.177", r4["indiv"], F(431, 2431))
check("Bennett S = 0.800", r4["S"], F(4, 5))
check("Gwet AC1 = 15569/17569 = 0.886", r4["AC1"], F(15569, 17569))
print(f"   headroom 1 - p_e = {float(1 - r4['pe_indiv']):.5f} "
      f"-> {round(float(r4['pe_indiv'])*1000)} slides by chance, "
      f"{round(float(r4['p_o'])*1000)} achieved, "
      f"{round(float(r4['p_o'])*1000) - round(float(r4['pe_indiv'])*1000)} ahead "
      f"of a possible {round(float(1-r4['pe_indiv'])*1000)}")
print(f"   of A's 65 cancer calls, B confirmed 15 ({15/65:.1%})")

print("\n   prevalence sweep: sensitivity 80%, specificity 95% fixed,")
print("   two raters independent given the truth; only prevalence moves")
print(f"   {'prev':>6} {'p_o':>9} {'kappa':>7} {'S':>7} {'AC1':>7}  {'+call rate':>10}")
SENS, SPEC = F(8, 10), F(95, 100)
prev_kappa = []
for prev in (F(1, 2), F(3, 10), F(1, 5), F(1, 10), F(1, 20), F(1, 50)):
    both1 = prev * SENS ** 2 + (1 - prev) * (1 - SPEC) ** 2
    both0 = prev * (1 - SENS) ** 2 + (1 - prev) * SPEC ** 2
    p_o = both1 + both0
    marg = prev * SENS + (1 - prev) * (1 - SPEC)
    p_e = marg ** 2 + (1 - marg) ** 2
    pe_ac1 = 2 * marg * (1 - marg)
    prev_kappa.append((prev, coef(p_o, p_e)))
    print(f"   {float(prev):>6.0%} {float(p_o):>9.4f} {float(coef(p_o, p_e)):>7.3f} "
          f"{float(coef(p_o, F(1,2))):>7.3f} {float(coef(p_o, pe_ac1)):>7.3f}  {float(marg):>10.4f}")
    if prev == F(1, 50):
        check("sweep 2% row reproduces the S4 marginal rate of 6.5%", marg, F(13, 200))
print("   NB p_o at 30% and 10% prevalence is exactly 0.8375 and 0.8825 -- ties at the")
print("      third decimal, so the page's .837 / .883 and a round-half-up .838 / .883")
print("      both 'look wrong'. Show four decimals in the page table.")
print(f"   NB kappa is NOT monotone in prevalence: {float(prev_kappa[1][1]):.3f} at 30% "
      f"exceeds {float(prev_kappa[0][1]):.3f} at 50%; it peaks between.")
check("kappa at 30% prevalence exceeds kappa at 50%", prev_kappa[1][1] > prev_kappa[0][1], True)

print("\n" + "=" * 72)
print("S5  more than two raters: one study read by 2, 4, 6, 16 and 32 readers")
print("=" * 72)

# Propensity model (no notion of truth): subject i has a propensity p_i that a
# reader calls it positive; readers are iid Bernoulli(p_i) given the subject.
# This is exactly rater exchangeability -- the condition under which Fleiss'
# kappa is invariant to panel size -- so the example cannot be rigged toward
# the claim. p_i ~ Beta(A, B) is U-shaped, mean 0.402.
S5_SEED, S5_N, S5_A, S5_B = 256, 200, 0.35, 0.52
S5_KS = (2, 4, 6, 16, 32)


def s5_panel(seed=S5_SEED, k=32, n=S5_N):
    rng = random.Random(seed)
    rows = []
    for _ in range(n):
        p = rng.betavariate(S5_A, S5_B)
        rows.append([1 if rng.random() < p else 0 for _ in range(k)])
    return rows


def s5_fleiss(sub, k):
    """float Fleiss kappa, binary; used inside the bootstrap loop for speed."""
    Po = 0.0
    tot = 0
    for s in sub:
        a = sum(s)
        b = k - a
        Po += (a * (a - 1) + b * (b - 1)) / (k * (k - 1))
        tot += a
    Po /= len(sub)
    p = tot / (len(sub) * k)
    Pe = p * p + (1 - p) * (1 - p)
    return (Po - Pe) / (1 - Pe)


rows32 = s5_panel()
print(f"\n   propensity p_i ~ Beta({S5_A}, {S5_B}), mean {S5_A/(S5_A+S5_B):.3f};"
      f"  {S5_N} biopsies, 32 readers, seed {S5_SEED}")
print(f"   pooled positive rate {sum(sum(s) for s in rows32)/(S5_N*32):.4f}")
print(f"   per-reader positive rate, first eight: "
      f"{[round(sum(s[r] for s in rows32)/S5_N, 3) for r in range(8)]}")

rng = random.Random(9)
s5 = {}
print(f"\n   {'k':>3} {'pairs':>6} {'P_o':>8} {'P_e':>8} {'kappa':>8} "
      f"{'all agree':>10} {'boot SE':>9} {'SE/SE(2)':>9} {'1/sqrt(k/2)':>12}")
for k in S5_KS:
    sub = [s[:k] for s in rows32]
    tot = sum(sum(s) for s in sub)
    pbar = F(tot, S5_N * k)
    Pe = pbar ** 2 + (1 - pbar) ** 2
    Po = sum(F(sum(s) * (sum(s) - 1) + (k - sum(s)) * (k - sum(s) - 1), k * (k - 1))
             for s in sub) / S5_N
    boot = [s5_fleiss([sub[rng.randrange(S5_N)] for _ in range(S5_N)], k)
            for _ in range(4000)]
    s5[k] = dict(Po=Po, Pe=Pe, kap=coef(Po, Pe), se=st.stdev(boot),
                 cons=F(sum(1 for s in sub if len(set(s)) == 1), S5_N))
    d = s5[k]
    print(f"   {k:>3} {k*(k-1)//2:>6} {float(d['Po']):>8.4f} {float(d['Pe']):>8.4f} "
          f"{float(d['kap']):>8.4f} {float(d['cons']):>9.1%} {d['se']:>9.4f} "
          f"{d['se']/s5[2]['se']:>9.3f} {(2/k)**.5:>12.3f}")

# --- claim 1: Fleiss at k = 2 IS Scott's pi -----------------------------------
sub2 = [s[:2] for s in rows32]
po2 = F(sum(1 for s in sub2 if s[0] == s[1]), S5_N)
pb2 = F(sum(sum(s) for s in sub2), 2 * S5_N)
check("S5: Fleiss P_o at k=2 is plain percent agreement", s5[2]['Po'], po2)
check("S5: Fleiss P_e at k=2 is Scott's pooled chance term",
      s5[2]['Pe'], pb2 ** 2 + (1 - pb2) ** 2)
check("S5: so Fleiss kappa at k=2 IS Scott's pi exactly",
      s5[2]['kap'], coef(po2, pb2 ** 2 + (1 - pb2) ** 2))

# --- claim 2: kappa is flat in k ----------------------------------------------
kaps = [float(s5[k]['kap']) for k in S5_KS]
print(f"\n   kappa across k: {[round(x, 4) for x in kaps]}  "
      f"spread {max(kaps)-min(kaps):.4f}")
check("S5 claim 1: kappa flat in k (spread < 0.01 over a 16-fold change)",
      max(kaps) - min(kaps) < 0.01, True)

# --- claim 3: consensus collapses ---------------------------------------------
cons = [float(s5[k]['cons']) for k in S5_KS]
print(f"   unanimity across k: {[f'{x:.1%}' for x in cons]}")
check("S5 claim 2: unanimity strictly decreasing in k",
      all(a > b for a, b in zip(cons, cons[1:])), True)
check("S5 claim 2: unanimity at k=32 is under half its k=2 value",
      cons[-1] < cons[0] / 2, True)

# --- claim 4: SE falls like 1/sqrt(k), then plateaus --------------------------
print(f"   SE across k:        {[round(s5[k]['se'], 4) for k in S5_KS]}")
r6 = s5[6]['se'] / s5[2]['se']
r32 = s5[32]['se'] / s5[2]['se']
check("S5 claim 3: SE tracks 1/sqrt(k) up to k=6 (within 0.03)",
      abs(r6 - (2/6) ** .5) < 0.03, True)
check("S5 claim 3: SE departs sharply from 1/sqrt(k) by k=32 (ratio > 1.5x)",
      r32 / (2/32) ** .5 > 1.5, True)
print(f"      k=6:  observed {r6:.3f} vs 1/sqrt(k/2) {(2/6)**.5:.3f}  -> tracks")
print(f"      k=32: observed {r32:.3f} vs 1/sqrt(k/2) {(2/32)**.5:.3f}  -> plateau")
print("      the floor is the 200 biopsies: more readers measure THESE subjects")
print("      better, and say nothing about subjects that were never collected.")

print("\n   how many of the 32 readers called it positive:")
dist = {}
for s in rows32:
    dist[sum(s)] = dist.get(sum(s), 0) + 1
for lo, hi in [(0, 0), (1, 4), (5, 11), (12, 20), (21, 27), (28, 31), (32, 32)]:
    print(f"      {lo:>2}-{hi:<2}: {sum(v for kk, v in dist.items() if lo <= kk <= hi):>3} biopsies")

print("\n" + "=" * 72)
print("S6  four raters, hijacked")
print("=" * 72)

def independent_panel(rs):
    """Exact 2^k product table for mutually independent raters."""
    rs = [F(x).limit_denominator(1000) for x in rs]
    k = len(rs)
    cells = []
    for pat in product((1, 0), repeat=k):
        pr = F(1)
        for r, v in zip(rs, pat):
            pr *= r if v == 1 else 1 - r
        cells.append((pat, pr))
    p = [[rs[r], 1 - rs[r]] for r in range(k)]
    npair = F(k * (k - 1), 2)
    p_o = sum(pr * (sum(F(m * (m - 1), 2) for m in (pat.count(1), pat.count(0))) / npair)
              for pat, pr in cells)
    return p, p_o, cells

for lbl, rs in (("all four aiming at 90%", [F(9,10)]*4),
                ("the hacker's RNGs, 90/70/50/30", [F(9,10), F(7,10), F(1,2), F(3,10)])):
    p, p_o, cells = independent_panel(rs)
    cats = (0, 1)
    pe_f = sum(x * x for x in pooled(p, cats))
    pe_c = pe_individual(p, cats)
    s2 = var_across_raters(p, cats)
    print(f"\n   {lbl}  (exact, no sampling noise)")
    print(f"      p_o {float(p_o):.4f}   p_e pooled {float(pe_f):.4f}   p_e individual {float(pe_c):.4f}")
    print(f"      Fleiss {float(coef(p_o, pe_f)):+.4f}   Conger {float(coef(p_o, pe_c)):+.4f}")
    check(f"      [{lbl}] independence -> Conger exactly 0", coef(p_o, pe_c), 0)
    check(f"      [{lbl}] Conger p_e = sum pbar^2 - (1/k) sum s2",
          pe_c, pe_f - sum(s2) / len(p))
    print("      distribution of how many of the four passed: " +
          ", ".join(f"{n}:{float(sum(pr for pat,pr in cells if pat.count(1)==n)):.4%}"
                    for n in (4, 3, 2, 1, 0)))

check("S5 four independent raters at 90% land on S1's chance term 0.82",
      independent_panel([F(9,10)]*4)[1], F(41, 50))
check("S5 hacker panel Fleiss = -5/72", 
      coef(independent_panel([F(9,10),F(7,10),F(1,2),F(3,10)])[1],
           sum(x*x for x in pooled(independent_panel([F(9,10),F(7,10),F(1,2),F(3,10)])[0], (0,1)))),
      F(-5, 72))

# ---- the realized panels shown on the page -------------------------------
SEED, N = 434, 500

def draw(rs, n, rng):
    return [[1 if rng.random() < r else 0 for r in rs] for _ in range(n)]

rng = random.Random(SEED)
panelA = draw([.9] * 4, N, rng)
panelB = draw([.9, .7, .5, .3], N, rng)

for lbl, rows in (("A  four professors, each aiming at 90%", panelA),
                  ("B  the hacker's four RNGs at 90/70/50/30", panelB)):
    cats = (0, 1)
    p = rates(rows, cats)
    p_o = p_o_pairs(rows, cats)
    pe_f, pe_c = pe_pooled(p, cats), pe_individual(p, cats)
    print(f"\n   {lbl}   (seed {SEED}, N={N})")
    print(f"      realized pass counts {[int(p[r][1]*N) for r in range(4)]} "
          f"= {[f'{float(p[r][1]):.1%}' for r in range(4)]}")
    dist = {}
    for s in rows:
        dist[s.count(1)] = dist.get(s.count(1), 0) + 1
    print(f"      essays by number passing  " +
          " / ".join(str(dist.get(n, 0)) for n in (4, 3, 2, 1, 0)))
    print(f"      p_o {float(p_o):.4f}   p_e pooled {float(pe_f):.4f}   "
          f"p_e individual {float(pe_c):.4f}")
    print(f"      Fleiss {float(coef(p_o, pe_f)):+.4f}   Conger {float(coef(p_o, pe_c)):+.4f}")
    pw = []
    for r, s in combinations(range(4), 2):
        po = F(sum(1 for x in rows if x[r] == x[s]), N)
        pe = sum(p[r][j] * p[s][j] for j in range(2))
        pw.append(float(coef(po, pe)))
    print(f"      six pairwise Cohen kappas  {[f'{x:+.3f}' for x in pw]}")
    check(f"      [{lbl[0]}] Conger p_e = sum pbar^2 - (1/k) sum s2",
          pe_c, pe_f - sum(var_across_raters(p, cats)) / 4)

# ---- sampling behaviour reported alongside the table ---------------------
print(f"\n   sampling behaviour of the hacker panel, {N} essays, 1500 draws")
rng = random.Random(20251024)
fl, cg = [], []
for _ in range(1500):
    rows = draw([.9, .7, .5, .3], N, rng)
    p = rates(rows, (0, 1))
    p_o = p_o_pairs(rows, (0, 1))
    fl.append(float(coef(p_o, pe_pooled(p, (0, 1)))))
    cg.append(float(coef(p_o, pe_individual(p, (0, 1)))))
print(f"      Fleiss  mean {st.mean(fl):+.4f}  sd {st.stdev(fl):.4f}   "
      f"(exact value {float(F(-5,72)):+.4f})")
print(f"      Conger  mean {st.mean(cg):+.4f}  sd {st.stdev(cg):.4f}   (exact value +0.0000)")

print("\n" + "=" * 72)
print(f"{len(FAILURES)} failure(s)" if FAILURES else "all checks passed")
for f in FAILURES:
    print("  " + f)
print("=" * 72)
