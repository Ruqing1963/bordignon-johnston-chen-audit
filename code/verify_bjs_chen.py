r"""verify_bjs_chen.py -- self-contained numerical audit of

  M. Bordignon, D. R. Johnston, V. Starichkova, "An explicit version of Chen's theorem and the
  linear sieve", Int. J. Number Theory 21 (2025) 2497-2572; arXiv:2207.09452 (v6).

All line numbers ("l.NNNN") refer to the arXiv v6 LaTeX source; results are numbered as in the arXiv v6
PDF (one global counter: main theorem = Theorem 3, linear sieve = Theorem 6, S(A) = Theorem 44, ...).  Every formula used below is
transcribed from that source; no other input is used.

Parts
  A. constants: twin-prime constant, the factor (2 e^{+-gamma} C_2)^{-1}
  B. the explicit linear sieve (Section 2): f_n, c_n (Table 2), tau_n (Lemma 11, Table 3),
     tau'_n (Lemma 12), J(k) and C_1, C_2 (Lemma 14, Table 1)
  C. the final inequality (Section 9): Theorems 44, 50, 51 combined through Lemma 37,
     in units of U_N N / log^2 N; printed computation, corrections, threshold scans and the
     variant with the sharpened Lemma 41.

Outputs (relative to the audit folder): data/verify_bjs_chen.txt, data/margin_scan.csv,
data/sieve_constants.csv.   Run:  python code/verify_bjs_chen.py      (about 2-4 minutes)
"""
import math, os, sys
import numpy as np
from mpmath import mp, mpf, log, exp, sqrt, quad, euler, primezeta, e as E

mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
os.makedirs(DATA, exist_ok=True)
OUT = []


def say(s=""):
    print(s)
    OUT.append(str(s))


def nstr(x, n=8):
    return mp.nstr(mpf(x), n)


# =====================================================================================
# A. constants
# =====================================================================================
EG = exp(euler)


def twin_prime_constant():
    # log C_2 = -sum_{n>=2} (2^n - 2)/n * sum_{p>2} p^{-n}
    s = mpf(0)
    for n in range(2, 400):
        t = (mpf(2) ** n - 2) / n * (primezeta(n) - mpf(2) ** (-n))
        s += t
        if abs(t) < mpf(10) ** (-28):
            break
    return exp(-s)


C2TW = twin_prime_constant()
U_PRINTED = 1 / (2 * EG * C2TW)      # (2 e^{+gamma} C_2)^{-1}: used at l.2103, 2335, 2353, 2357, 2510
U_CORRECT = EG / (2 * C2TW)          # (2 e^{-gamma} C_2)^{-1}: the true bound for 1/U_N, by (l.280)


def part_A():
    say("=" * 90)
    say("A. Constants")
    say("=" * 90)
    say("gamma = %s, e^gamma = %s" % (nstr(euler, 12), nstr(EG, 12)))
    say("twin prime constant C_2 = prod_{p>2}(1-(p-1)^-2) = %s   (paper l.2070: 0.66016...)" % nstr(C2TW, 15))
    say("min_N U_N = 2 e^{-gamma} C_2 = %s  (U_N >= this since prod_{p|N,p>2}(p-1)/(p-2) >= 1)" % nstr(2 * C2TW / EG, 10))
    say("printed factor (2 e^{+gamma} C_2)^{-1} = %s" % nstr(U_PRINTED, 10))
    say("correct factor (2 e^{-gamma} C_2)^{-1} = %s" % nstr(U_CORRECT, 10))
    say("ratio = e^{2 gamma} = %s" % nstr(U_CORRECT / U_PRINTED, 10))


# =====================================================================================
# B. the explicit linear sieve
# =====================================================================================
def h_fun(s):
    s = np.asarray(s, dtype=float)
    return np.where(s <= 2, math.exp(-2), np.where(s <= 3, np.exp(-s), 3 / s * np.exp(-s)))


C_TABLE2 = {}
for rng, val in [((2, 2), .33), ((3, 3), .39), ((4, 4), .45), ((5, 5), .51), ((6, 6), .54), ((7, 7), .57),
                 ((8, 8), .58), ((9, 10), .61), ((11, 12), .63), ((13, 13), .64), ((14, 14), .65),
                 ((15, 18), .66), ((19, 20), .67), ((21, 26), .68), ((27, 34), .69), ((35, 46), .70),
                 ((47, 82), .71), ((83, 345), .72), ((346, 500), .73)]:
    for n in range(rng[0], rng[1] + 1):
        C_TABLE2[n] = val
C_TABLE2[1] = 1.0
KAPPA2 = 0.9607  # Lemma 8 (Nathanson, Lemma 9.7): c_n = 0.9607 for all n >= 2


def c_of(n, ctab):
    return ctab[n] if n in ctab else KAPPA2


def compute_fn(nmax=500, H=2000):
    """f_n on the grid s = 1 + i/H (trapezoid rule), definitions (l.465-481).  Returns the accurate
    c_n = max_s (f_n(s) / (2 e^2 h(s)))^{1/(n-1)} and, as a check, f(4), F(4) via (l.459-461)."""
    smax = nmax + 4
    s = 1 + np.arange(int((smax - 1) * H) + 1) / H
    hs = h_fun(s)
    shift = H  # index offset for s -> s-1
    f = np.where(s <= 3, (3 - s) / s, 0.0)  # f_1, valid for s >= 1
    cn = {1: 1.0}
    i4 = 3 * H  # index of s = 4
    F4 = 1 + f[i4]; f4 = 1.0
    i3 = 2 * H
    for n in range(2, nmax + 1):
        lo = 0 if (n - 1) % 2 == 1 else H             # valid region of f_{n-1}: s>=1 (odd) or s>=2 (even)
        g = f.copy(); g[:lo] = 0.0
        # tail integral T(u) = int_u^inf f_{n-1}, trapezoid
        seg = (g[1:] + g[:-1]) / (2 * H)
        T = np.zeros_like(g); T[:-1] = np.cumsum(seg[::-1])[::-1]
        new = np.zeros_like(f)
        new[shift:] = T[:-shift] / s[shift:]          # s f_n(s) = int_{s-1}^inf f_{n-1}
        if n % 2 == 1:
            new[:i3] = 3 * new[i3] / s[:i3]          # s f_n(s) = 3 f_n(3) for 1 <= s <= 3
            valid = slice(0, None)
            F4 += new[i4]
        else:
            valid = slice(H, None)                    # s >= 2
            f4 -= new[i4]
        f = new
        ratio = f[valid] / (2 * math.e ** 2 * hs[valid])
        cn[n] = float(ratio.max()) ** (1.0 / (n - 1))
    return cn, f4, F4


GAMMA_S0 = {"printed": lambda s0: math.exp(s0 - 2), "correct": lambda s0: math.exp(s0 - 1.8)}
KAP = {2.0: .9607, 2.2: .9557, 2.4: .9457, 2.6: .9261, 2.8: .8914}
KAP3, KTIL = .8349, .9214
GAM3 = 4 * math.e / 3


def tau_seq(eps, gam="correct", xi3=True, nmax=500, ctab=C_TABLE2):
    """tau_n of Lemma 11, (l.756-768).  gam: gamma_{s0} = e^{s0-2} (printed l.608) or e^{s0-1.8};
    xi3=False drops xi_3 from the maxima (only used to test reproducibility of Table 3)."""
    g = GAMMA_S0[gam]
    t = [None, 3.0]
    for n in range(2, nmax + 1):
        cm, cc = c_of(n - 1, ctab) ** (n - 2), c_of(n, ctab) ** (n - 1)
        prev = t[-1]
        cands = []
        if n % 2 == 0:
            for s0, k in KAP.items():
                gm = g(s0)
                cands.append(prev * (k + (gm + k) * eps) + 2 * gm * cm + 2 * cc)
        else:
            cands.append(prev * (KTIL + (5 + 4 * eps) * eps) + 6 * (1 + eps) * cm + 2 * (2 + eps) * cc)
        if xi3:
            cands.append(prev * (KAP3 + (GAM3 + KAP3) * eps) + 2 * GAM3 * cm + 2 * cc)
        t.append(max(cands))
    return t


def tau_prime_seq(eps, gam="correct", nmax=502, ctab=C_TABLE2):
    """tau'_n of Lemma 12, (l.928-932)."""
    g2 = GAMMA_S0[gam](2.0)
    A = KAPPA2 + (g2 + KAPPA2) * eps
    t = [None, 3.0]
    for n in range(2, nmax + 1):
        t.append(t[-1] * A + 8 * math.e / 3 * c_of(n - 1, ctab) ** (n - 2) + 2 * (2 + eps) * c_of(n, ctab) ** (n - 1))
    return t


def J_of(k, eps, gam="correct"):
    """J(k) of Lemma 14 (l.1004)."""
    g2 = GAMMA_S0[gam](2.0)
    A = KAPPA2 + (g2 + KAPPA2) * eps
    return A + 8 * math.e / 6 * (KAPPA2 / A) ** (k - 2) + 2 * (2 + eps) / 3 * KAPPA2 ** (k - 1) / A ** (k - 2)


def C_lemma28(eps, gam="correct", ko=249, ke=250):
    """C_1, C_2 by Lemma 14 with k_o = 249, k_e = 250 (l.1038); inf if the geometric series diverges."""
    t = tau_seq(eps, gam)
    tp = tau_prime_seq(eps, gam)
    Jo, Je = J_of(2 * ko + 1, eps, gam), J_of(2 * ke + 2, eps, gam)
    C1 = sum(t[2 * n - 1] for n in range(1, ko + 1)) + (tp[2 * ko] * Jo / (1 - Jo ** 2) if Jo < 1 else math.inf)
    C2 = sum(t[2 * n] for n in range(1, ke + 1)) + (tp[2 * ke + 1] * Je / (1 - Je ** 2) if Je < 1 else math.inf)
    return C1, C2, Jo, Je


def C_direct(eps, gam="correct", xi3=True, ctab=C_TABLE2, M=40000):
    """Sum of tau_n using Lemma 11 for n <= 500 and, for n > 500, the one-step bound from the proof
    of Lemma 12, tau_n <= A tau_{n-1} + (8e/3) c_{n-1}^{n-2} + 2(2+eps) c_n^{n-1}, with c_n = 0.9607
    (Lemma 8) for n > 500; explicit summation to M, then the closed-form geometric tail."""
    t = tau_seq(eps, gam, xi3, 500, ctab)
    g2 = GAMMA_S0[gam](2.0)
    A = KAPPA2 + (g2 + KAPPA2) * eps
    assert A < 1
    odd = sum(t[n] for n in range(1, 501, 2)); even = sum(t[n] for n in range(2, 501, 2))
    prev = t[500]
    for n in range(501, M + 1):
        prev = A * prev + 8 * math.e / 3 * c_of(n - 1, ctab) ** (n - 2) + 2 * (2 + eps) * KAPPA2 ** (n - 1)
        if n % 2: odd += prev
        else: even += prev
    # tail n > M: tau_n <= A tau_{n-1} + B kappa^{n-2}, B = 8e/3 + 2(2+eps) kappa; sum_{n>M} tau_n <= T
    B = 8 * math.e / 3 + 2 * (2 + eps) * KAPPA2
    T = (A * prev + B * KAPPA2 ** (M - 1) / (1 - KAPPA2)) / (1 - A)
    return odd + T, even + T


TABLE1 = {74: (631, 630), 76: (559, 559), 78: (504, 504), 80: (461, 461), 85: (386, 386), 90: (336, 337),
          95: (302, 302), 100: (276, 277), 120: (218, 219), 140: (189, 190), 160: (172, 173), 180: (161, 162),
          200: (153, 154), 300: (133, 134), 400: (125, 126), 500: (121, 122), 600: (118, 119), 700: (116, 117),
          800: (115, 116), 900: (114, 115), 1000: (113, 114), 2000: (109, 110), 10000: (106, 108),
          100000: (106, 107)}
TABLE3 = [((1, 1), 3), ((2, 2), 8), ((3, 10), 10), ((11, 13), 9), ((14, 16), 8), ((17, 20), 7), ((21, 24), 6),
          ((25, 30), 5), ((31, 37), 4), ((38, 46), 3), ((47, 63), 2), ((64, 118), 1), ((119, 173), 1e-1),
          ((174, 228), 1e-2), ((229, 283), 1e-3), ((284, 338), 1e-4), ((339, 393), 1e-5), ((394, 448), 1e-6),
          ((449, 500), 1e-7)]


def part_B():
    say("\n" + "=" * 90)
    say("B. The explicit linear sieve (Section 2)")
    say("=" * 90)
    # --- Lemma 9 (l.607-621)
    say("Lemma 9 (l.608): sup_{2<=s<=s0+0.2} h(s-1)/h(s) versus gamma_{s0}=e^{s0-2} (printed) and e^{s0-1.8}")
    ss = np.linspace(2, 3, 100001)
    for s0 in [2.0, 2.2, 2.4, 2.6, 2.8]:
        m = ss[ss <= s0 + 0.2 + 1e-12]
        sup = float((h_fun(m - 1) / h_fun(m)).max())
        say("  s0=%.1f  sup=%.6f   e^{s0-2}=%.6f   e^{s0-1.8}=%.6f" % (s0, sup, math.exp(s0 - 2), math.exp(s0 - 1.8)))
    say("  (Remark after Lemma 14, l.1045, itself uses gamma_2 = e^{0.2}.)")
    # --- c_n (Table 2)
    cn, f4, F4 = compute_fn()
    say("\nTable 2 (l.495-559): c_n from f_n on a grid of step 1/2000 (trapezoid rule)")
    say("  checks: f(4) via f_n = %.7f vs 2e^g log3/4 = %.7f;  F(4) via f_n = %.7f" % (f4, float(2 * EG * log(3) / 4), F4))
    bad = [(n, cn[n], C_TABLE2[n]) for n in range(2, 501) if cn[n] > C_TABLE2[n]]
    say("  n with computed c_n > tabulated c_n: " + ", ".join("n=%d: %.5f > %.2f" % b for b in bad))
    say("  sample: " + ", ".join("c_%d=%.4f" % (n, cn[n]) for n in [2, 3, 5, 10, 20, 33, 50, 100, 345, 500]))
    ctab_fix = dict(C_TABLE2)
    for n, v, _ in bad:
        ctab_fix[n] = math.ceil(v * 100) / 100
    # --- Table 3
    say("\nTable 3 (l.860-924): tau_n for eps = 1/200 from (l.756-768) with Table 2 c_n")
    for gam in ["printed", "correct"]:
        for xi3 in [True, False]:
            t = tau_seq(1 / 200, gam, xi3)
            viol = [n for (a, b), v in TABLE3 for n in range(a, b + 1) if t[n] > v]
            d1, d2 = C_direct(1 / 200, gam, xi3)
            say("  gamma %-7s xi_3 %-8s tau_2=%.3f max_{3<=n<=10} tau_n=%.3f ; Table 3 exceeded at %2d values of n%s;"
                " C_1, C_2 = %.1f, %.1f (Table 1: 153, 154)"
                % (gam, "incl." if xi3 else "omitted", t[2], max(t[3:11]), len(viol),
                   (" (first n=%d)" % viol[0]) if viol else "", d1, d2))
    # --- Lemma 12 hypotheses with the correct gamma
    say("\nLemma 12 (l.940-953) hypotheses, gamma_{s0} = e^{s0-1.8}:")
    lin = lambda s0, k: (k, GAMMA_S0['correct'](s0) + k)
    worst = None
    for x in np.linspace(1e-6, 1 / 40, 20000):
        vals = [k + (g + k) * x for (k, g) in [(KAP[s], GAMMA_S0['correct'](s)) for s in KAP]] + [KAP3 + (GAM3 + KAP3) * x]
        if max(vals) > vals[0] + 1e-15:
            worst = x; break
    say("  max_{s0} kappa_{s0}+(gamma_{s0}+kappa_{s0})eps = value at s0=2 holds for eps < %.6f (1/%.2f) [needs 1/74]" % (worst, 1 / worst))
    for gname in ["correct", "printed"]:
        b_ = 5 - (GAMMA_S0[gname](2.0) + KAPPA2); c_ = KAPPA2 - KTIL     # 4 eps^2 + b eps - c <= 0
        xt = (-b_ + math.sqrt(b_ * b_ + 16 * c_)) / 8
        say("  kappa~+(5+4eps)eps <= kappa_2+(gamma_2+kappa_2)eps  [gamma_2 %s] holds iff eps <= %.6f (1/%.2f)"
            % (gname, xt, 1 / xt))
    # --- Table 1
    say("\nTable 1 (l.376-448): C_1, C_2.  Lemma 14 with k_o=249, k_e=250, and direct summation")
    say("  %7s | %-10s | %-34s | %-19s | %-19s" % ("1/eps", "printed", "Lemma 14 (gamma corr.)  J(499) J(502)", "direct, gamma corr.", "direct, gamma printed"))
    rows = []
    for inv in sorted(TABLE1) + [round(1 / 1.452e-7)]:
        eps = 1 / inv
        L1, L2, Jo, Je = C_lemma28(eps, "correct")
        D1, D2 = C_direct(eps, "correct")
        P1, P2 = C_direct(eps, "printed")
        F1, F2 = C_direct(eps, "correct", ctab=ctab_fix)
        pr = TABLE1.get(inv, (106, 107) if inv > 10 ** 6 else (float('nan'),) * 2)
        say("  %7d | %4s %4s | %8s %8s  %6.4f %6.4f | %8.2f %8.2f | %8.2f %8.2f"
            % (inv, pr[0], pr[1], "%.1f" % L1 if L1 < math.inf else "inf", "%.1f" % L2 if L2 < math.inf else "inf",
               Jo, Je, D1, D2, P1, P2))
        rows.append((inv, pr[0], pr[1], L1, L2, Jo, Je, D1, D2, P1, P2, F1, F2))
    inv = round(1 / 1.452e-7)
    r = rows[-1]
    say("  eps = 1.452e-7 (used in Sections 6-8 with C_1=106, C_2=107, l.2052): direct C_1=%.3f C_2=%.3f;"
        " with c_n corrected where Table 2 is exceeded: C_1=%.3f C_2=%.3f" % (r[7], r[8], r[11], r[12]))
    say("  first 1/eps in Table 1 with J(499) >= 1 (Lemma 14 vacuous): %s"
        % next(row[0] for row in rows if row[5] >= 1))
    # exact boundary
    lo, hi = 200.0, 300.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if J_of(499, 1 / mid) < 1: lo = mid
        else: hi = mid
    say("  J(499) = 1 at 1/eps = %.2f ; J(502) = 1 at 1/eps = %.2f"
        % (lo, next(x for x in np.arange(200, 400, 0.01) if J_of(502, 1 / x) >= 1)))
    say("  ratios direct/printed: " + ", ".join("1/eps=%d: %.3f/%.3f" % (row[0], row[7] / row[1], row[8] / row[2])
                                                for row in rows if row[0] in (100, 200, 1000, 100000)))
    with open(os.path.join(DATA, "sieve_constants.csv"), "w") as fh:
        fh.write("inv_eps,C1_printed,C2_printed,C1_lemma28,C2_lemma28,J499,J502,C1_direct,C2_direct,"
                 "C1_direct_printedgamma,C2_direct_printedgamma,C1_direct_cfix,C2_direct_cfix\n")
        for row in rows:
            fh.write(",".join(str(x) for x in row) + "\n")
    # Lemma 14 J(k) as eps -> 0
    say("  J(k) -> kappa_2 + 8e/6 + (4/3)kappa_2 = %.4f as eps -> 0 (any fixed k)" % (KAPPA2 + 8 * math.e / 6 + 4 / 3 * KAPPA2))
    say("  (2.5) check (l.454): 196.1 * e^2 h(s)/200 <= %.4f" % (196.1 / 200))
    return r[7], r[8]


# =====================================================================================
# C. The final inequality (Section 9)
# =====================================================================================
R1 = mpf('2.0452')
PRIMES = [p for p in range(3, 400) if all(p % q for q in range(2, int(p ** .5) + 1))]
LOGPRIM = []  # (p, log prod_{2<q<=p} q)
_acc = mpf(0)
for _p in PRIMES:
    _acc += log(_p)
    LOGPRIM.append((_p, _acc))


def eps_lemma63(logK, sharp=False):
    """Lemma 41 (l.1897-1922).  paper: pbar = largest prime with prod_{2<p<=pbar} p <= K, eps = 1/(pbar-2).
    sharp: p+ = smallest prime with prod_{2<p<=p+} p >= K; then q_1 >= p+ and eps = 1/(p+ - 2)."""
    if not sharp:
        pb = max(p for p, lp in LOGPRIM if lp <= logK)
    else:
        pb = min(p for p, lp in LOGPRIM if lp >= logK)
    return 1 / mpf(pb - 2), pb


def lx1(LX): return LX - 5 * log(LX)        # log x_1(X), x_1(X) = X/log^5 X   (l.1361)
def lx2(LX): return LX - 15 * log(LX)       # log x_2(Y), x_2(Y) = Y/log^15 Y


def nu(Lx, delta):
    """(l.1402), (l.1545): nu(x) = min{100/(sqrt(K) log^2 K), 1/(2 R_1 log Q_1(x))}, K = log^delta x."""
    logK = delta * log(Lx)
    return min(100 / (exp(logK / 2) * logK ** 2), 1 / (2 * R1 * 10 * log(Lx)))


GRID = [mpf(1)] + [1 + mpf(10) ** k for k in range(-9, 7)]


def sup_from(fun, L0):
    """sup of fun(L) over L >= L0 (L = log of the variable), on a geometric grid; also reports whether the
    maximum sits at the left end point (all functions below are decreasing there)."""
    vals = [fun(L0 * r) for r in GRID]
    m = max(vals)
    return m, vals.index(m) == 0


ENDPOINT_OK = [True]


def smax(fun, L0):
    m, ok = sup_from(fun, L0)
    if not ok: ENDPOINT_OK[0] = False
    return m


def p_const(LX2, delta, star=False):
    """p(X_2), p*(X_2): Lemmas 31 and 34 (l.1478-1483, l.1674-1678)."""
    L1 = lx1(LX2)
    b0 = 1 - 1 / (2 * R1 * 10 * log(L1)) if star else 1 - nu(L1, delta)

    def inner(ly):
        siegel = mpf('1.1') * 10 * log(ly) * (mpf('3.2e-8') / ly ** 6 + ly ** 2 * exp((b0 - 1) * ly) / b0)
        Ey = 27 * (4 * ly ** mpf('6.5') / lx1(ly) ** 10 + 4 / ly ** mpf('3.5') + 18 * exp(-ly / 12) * ly ** mpf('1.5')
                   + mpf('2.5') * exp(-ly / 6) * ly ** mpf('7.5'))
        rest = exp(-ly / 2) / (2 * log(2) * ly ** 6) + mpf('0.4') * ly ** 5 * exp(-ly)
        return siegel + Ey + rest
    p2 = smax(inner, L1)
    p1 = p2 + (mpf('0.67') + 2 * exp(-L1 / 6)) / L1 ** 8
    return p1 * (1 + 1 / (LX2 ** 2 * L1 ** 3) + 1 / ((1 - 4 / L1) * LX2)) + mpf('2.2') / LX2 ** 2


def c_const(LX2):
    """c(X_2): Lemma 33 (l.1612-1617)."""
    L1 = lx1(LX2)

    def inner(ly):
        nu_ = 1 / (2 * R1 * 10 * log(ly))
        return (mpf('3.2e-8') / ly ** 6 + ly ** 2 * (exp(-nu_ * ly) / (1 - nu_)
                + ly ** 10 * (mpf('1.02') * exp(-ly / 2) + 3 * exp(-2 * ly / 3))
                + mpf('9.4') * ly ** mpf('1.515') * exp(-mpf('0.8274') * sqrt(ly))))
    c1 = smax(inner, L1)
    return c1 * (1 + 1 / (LX2 ** 2 * L1 ** 3) + 1 / ((1 - 4 / L1) * LX2)) + 1 / LX2 ** 2


def a_const(LX2, delta):
    """c_2, c_3 (Lemma 42, l.1948-1962), a_1, a (l.2039-2041)."""
    c = c_const(LX2)
    c2 = c + mpf('1.3841') * LX2 ** 4 * exp(-LX2) / log(LX2)

    def c3_in(LN):
        L1 = lx1(LN); lQ = 10 * log(L1); lll = log(log(LN))
        t = (mpf('1.5') / LN + (LN + lQ) / (LN - lQ)) * (LN / L1) ** delta
        u = log(delta * log(L1))
        return t * (EG * u + mpf('2.5') / u) / lll + mpf('1.3841') * LN ** (2 + delta) * exp(-LN) / (log(LN) * lll)
    c3 = smax(c3_in, LX2)

    def a1_in(LN):
        lQ = 10 * log(lx1(LN))
        return c2 / (LN ** (2 - delta) * log(log(LN))) * mpf('1.3841') * lQ / log(lQ)
    a1 = smax(a1_in, LX2) + c3

    def a_in(LN):
        lQ = 10 * log(lx1(LN))
        return log(log(LN)) / LN ** delta / C2TW * (EG * log(lQ) + mpf('2.5') / log(lQ))
    return a1 * smax(a_in, LX2), c3


def m_const(LX3, delta, star=False):
    """m(X_3): Lemma 32 (l.1540-1546) with v_0 from Lemma 29 (l.1311-1319), beta_0(x_2) at x_2 = x_2(X_3)."""
    L2 = lx2(LX3)
    b0 = 1 - 1 / (2 * R1 * 10 * log(L2)) if star else 1 - nu(L2, delta)

    def vp_in(ly):
        return mpf('3.2e-8') / ly ** 4 + ly ** 4 * (exp((b0 - 1) * ly) / b0 + mpf('1.02') * ly ** 10 * exp(-ly / 2)
                                                    + 3 * ly ** 10 * exp(-2 * ly / 3))
    vp = smax(vp_in, L2)
    v0 = vp * (1 + 1 / (LX3 ** 10 * L2 ** 5) + 1 / ((1 - 6 / L2) * LX3)) + 3 / LX3
    m = (39 * v0 + 108 * LX3 ** 16 * exp(-LX3) / log(LX3) + 26 * LX3 ** 5 / L2 ** 10
         + 88 * LX3 ** 5 * (exp(-8 * LX3 / 3) + exp(-LX3 / 2)) + 106 / LX3 ** 6)
    return m, vp, 1 - b0


CBAR = quad(lambda b: log(2 - 3 * b) / (b * (1 - b)), [mpf(1) / 8, mpf(1) / 3])


def f_lin(s): return 2 * EG * log(s - 1) / s                                              # (l.1186)
def F_lin(s): return 2 * EG / s * (1 + quad(lambda t: log(t - 1) / t, [2, s - 1]))      # (l.1184)
def h_mp(s):
    s = mpf(s)
    return exp(-2) if s <= 2 else (exp(-s) if s <= 3 else 3 / s * exp(-s))


def chen(LLX2=mpf('32.7'), LX3=None, delta=mpf('1.478'), al=mpf('1e-5'), eps0=mpf('1e-4'), eps=mpf('1.452e-7'),
         C1=mpf(106), C2=mpf(107), ufix=False, eps2_arg="X3^8", sharp=False):
    """Lower bound for pi_2(N) / (U_N N / log^2 N) at N = X_2 (all error terms decrease in N).
    ufix: replace (2e^{+g}C_2)^{-1} by (2e^{-g}C_2)^{-1} wherever the former is used.
    eps2_arg: 'X3^8' (as printed, l.2508, 2692, 2701) or 'X3' (= the bound actually available, since
    k_2 >= K_delta(x_2(Y_j)) with Y_j >= N^{1/8} >= X_3).   sharp: sharpened Lemma 41."""
    LX2 = exp(LLX2)
    if LX3 is None: LX3 = LX2 / 8
    L = LX2
    Uf = U_CORRECT if ufix else U_PRINTED
    Cb = max(C1, C2)
    e1, pb1 = eps_lemma63(delta * log(lx1(LX2)), sharp)
    LXe2 = 8 * LX3 if eps2_arg == "X3^8" else LX3
    e2, pb2 = eps_lemma63(delta * log(lx2(LXe2)), sharp)
    p = p_const(LX2, delta); ps = p_const(LX2, delta, True)
    a, c3 = a_const(LX2, delta)
    L1 = lx1(LX2)
    tail = mpf('0.9') * exp(L1 / 2) * LX2 ** 4 / (exp(LX2) * L1 ** 10 * log(LX2))
    c4, c4s = p + tail, ps + tail                                                       # (l.2016-2017)
    lQ = 10 * log(L1)
    eta = mpf('32.02') / L ** 2
    # ---- Theorem 44 (l.2044-2067); the 1/U_N factor in the statement is (2e^{-g}C_2)^{-1}
    sa = 4 - 8 * al
    SA_a = 8 * (1 - eta) * (f_lin(sa) - C2 * eps * E ** 2 * h_mp(sa) - (1 - eta) ** -1 * U_CORRECT * c4 / (8 * L))
    cal = 4 - 8 * al - 160 * log(LX2) / LX2
    mbar = max(1 - f_lin(cal), F_lin(cal) - 1)
    SA_b = 8 * (1 + eta) * (f_lin(cal) - e1 * (1 - f_lin(cal)) - (1 + e1) * eps * C2 * E ** 2 * h_mp(cal)
                            - (3 * e1 + a) * (mbar + eps * Cb * E ** 2 * h_mp(cal)) - a - mpf('64.04') / L ** 2
                            - (1 + eta) ** -1 * U_CORRECT * c4s / (8 * L) * mpf('1.3841') * lQ / log(lQ))
    # ---- Theorem 50 (l.2345-2370)
    k = 8 * (mpf(1) / 6 - al)
    r8 = 1 / (1 - exp(-LX2 / 8))
    l2 = r8 * (EG / 4 * ((log(6) + log((3 - 8 * al) / (3 - 18 * al))) / (mpf(1) / 2 - al) + 512 / (k * LX2 ** 2))
               + (log(mpf(8) / 3) + 64 / LX2 ** 2) * eps * C1 * E ** 2 * h_mp(k))

    def l1f(pp):
        return (LX2 ** 2 * exp(-LX2)) * (2 * EG / k + eps * C1 * E ** 2 * h_mp(k)) * (
            r8 * (log(mpf(8) / 3) + 64 / LX2 ** 2) * mpf('1.3841') / log(LX2) + pp * exp(LX2) / LX2 ** 4)

    def lf(pp):
        return pp * (1 / LX2 + mpf('0.55')) + mpf('1.3841') * (log(mpf(8) / 3) + 64 / LX2 ** 2) * exp(L1 / 2) * LX2 ** 3 / (
            exp(LX2) * log(LX2) * L1 ** 10)
    SQ_a = 8 * mpf('1.00005') * (1 + eta) * (l1f(p) + l2) + lf(p) * Uf / L
    SQ_b = 8 * mpf('1.00005') * (1 + eta) * (1 + e1) * (l1f(ps) + l2) + lf(ps) * Uf / L
    # ---- Theorem 51 (l.2499-2512)
    m, vp, onemb = m_const(LX3, delta)
    SB_main = (mpf('1.00005') * (1 + e2) * (1 + mpf('4.51') / L ** 2)
               * (2 * EG / (mpf(1) / 2 - al) + 3 * eps * C1 * E ** 2 * h_mp(mpf('1.5') - 3 * al))
               * (1 + eps0 + 9 / L) * (CBAR + 36 / L ** 2 + (log(mpf(8) / 3) + 64 / L ** 2) * (10 * log(1 + eps0) / L + 27 / L ** 2)))
    SB_err = Uf * (320 * m * (1 + eps0) / (3 * log(1 + eps0))
                   + mpf('0.13') * (1 + eps0) * L ** 5 * exp(-L / 8) / log(1 + eps0))
    SB = SB_main + SB_err
    # ---- Lemma 37 (l.1782-1784); 2N^{7/8}+N^{1/3} in units U_N N/log^2 N is < 3 L^2 e^{-L/8} / min U_N
    extra = 3 * L ** 2 * exp(-L / 8) * U_CORRECT
    res_a = SA_a - SQ_a / 2 - SB / 2 - extra
    res_b = SA_b - SQ_b / 2 - SB / 2 - extra
    return dict(e1=e1, pb1=pb1, e2=e2, pb2=pb2, p=p, ps=ps, a=a, c3=c3, m=m, vp=vp, onemb=onemb, SA_a=SA_a,
                SA_b=SA_b, SQ_a=SQ_a, SQ_b=SQ_b, SB_main=SB_main, SB_err=SB_err, res_a=res_a, res_b=res_b)


def conditions_ok(LLX2, LX3, delta, al=mpf('1e-5')):
    """Hypotheses of Theorems 44, 50, 51 at N = X_2 (all monotone in N)."""
    LX2 = exp(LLX2); L = LX2; u0 = mpf(10) ** 9
    rhs = u0 * (1 + mpf('9e-7') / log(u0))
    L1 = lx1(LX2)
    ok = [al * L - 10 * log(L1) - mpf('2.5') * log(L) >= rhs,           # (l.2049), (l.2349)
          (mpf('0.5') - al) * L - 20 * log(L) >= L / 4,                    # (l.2049)
          8 * al + 160 * log(L) / L < 1,                                   # (l.2050)
          delta * log(L1) >= log(3022),                                    # K_delta(x_1(X_2)) >= 3022
          10 * log(3) + al * L - 10 * log(L) >= rhs,                       # (l.2503)
          log(lx1(LX2)) >= mpf('10.4'), log(lx2(LX3)) >= mpf('10.4'),
          8 * LX3 <= LX2 + mpf('1e-20') * LX2]                             # N >= X_2 >= X_3^8
    return all(ok)


def delta_candidates(LLX2, LX3, eps2_arg, sharp, lo=mpf('1.30'), hi=mpf('1.70')):
    """The margin is piecewise decreasing in delta, with upward jumps where K_delta crosses an odd
    primorial; its supremum over delta is attained just to the right of a jump.  Return those points,
    rounded up to 4 decimals."""
    LX2 = exp(LLX2)
    LXe2 = 8 * LX3 if eps2_arg == "X3^8" else LX3
    base = [log(lx1(LX2)), log(lx2(LXe2))]
    out = set()
    for b in base:
        for p, lp in LOGPRIM:
            d = lp / b
            if lo <= d <= hi:
                out.add(mpf(math.ceil(float(d) * 10 ** 4 + 1e-9)) / 10 ** 4)
    return sorted(out)


def best_delta(LLX2, LX3=None, **kw):
    LX2 = exp(LLX2)
    if LX3 is None: LX3 = LX2 / 8
    best = None
    for d in delta_candidates(LLX2, LX3, kw.get("eps2_arg", "X3^8"), kw.get("sharp", False)):
        if not conditions_ok(LLX2, LX3, d): continue
        r = chen(LLX2=LLX2, LX3=LX3, delta=d, **kw)
        if best is None or r['res_b'] > best[1]['res_b']:
            best = (d, r)
    return best


def show(label, r, delta):
    say("  [%s]  delta=%s" % (label, nstr(delta, 5)))
    say("    eps_1 = 1/%d (prime %d), eps_2 = 1/%d (prime %d); m(X_3) = %s (v' = %s, 1-beta_0(x_2) = %s)"
        % (r['pb1'] - 2, r['pb1'], r['pb2'] - 2, r['pb2'], nstr(r['m'], 6), nstr(r['vp'], 6), nstr(r['onemb'], 4)))
    say("    p = %s, p* = %s, a = %s, c_3 = %s" % (nstr(r['p'], 3), nstr(r['ps'], 3), nstr(r['a'], 3), nstr(r['c3'], 5)))
    say("    S(A) >= %s (a), %s (b);  sum S(A_q) <= %s (a), %s (b);  S(B) <= %s + %s"
        % tuple(nstr(r[k], 10) for k in ['SA_a', 'SA_b', 'SQ_a', 'SQ_b', 'SB_main', 'SB_err']))
    say("    ==> pi_2(N) / (U_N N/log^2 N) >  case (a): %s   case (b): %s" % (nstr(r['res_a'], 8), nstr(r['res_b'], 8)))


def part_C(C1c, C2c):
    say("\n" + "=" * 90)
    say("C. The final inequality (Section 9, l.2725-2735), units U_N N / log^2 N, evaluated at N = X_2")
    say("=" * 90)
    say("cbar = %s (l.2553: < 0.363084);  f(4) = %s, F(4) = %s" % (nstr(CBAR, 12), nstr(f_lin(4), 10), nstr(F_lin(4), 10)))
    LX3p = exp(mpf('30.62'))
    Cc = dict(C1=mpf(math.ceil(C1c * 10) / 10), C2=mpf(math.ceil(C2c * 10) / 10))
    say("recomputed sieve constants used below at eps = 1.452e-7: C_1 = %s, C_2 = %s" % (nstr(Cc['C1'], 5), nstr(Cc['C2'], 5)))
    ENDPOINT_OK[0] = True
    say("\n-- C.1 the paper's parameters: X_2 = exp(exp(32.7)), X_3 = exp(exp(30.62)), delta = 1.478,"
        " alpha_i = 1e-5, eps_0 = 1e-4 (l.2727)")
    d0 = mpf('1.478')
    runs = [("as printed: C=106/107, (2e^{+g}C_2)^{-1}, eps_2((X_3)^8)", dict()),
            ("only the factor corrected to (2e^{-g}C_2)^{-1}", dict(ufix=True)),
            ("only eps_2 corrected to eps_2(X_3)", dict(eps2_arg="X3")),
            ("only C_1, C_2 recomputed", dict(**Cc)),
            ("all three corrected", dict(ufix=True, eps2_arg="X3", **Cc)),
            ("all three corrected + sharpened Lemma 41", dict(ufix=True, eps2_arg="X3", sharp=True, **Cc))]
    summary = {}
    for lab, kw in runs:
        r = chen(LX3=LX3p, delta=d0, **kw)
        show(lab, r, d0)
        summary[lab] = r['res_b']
    r0 = chen(LX3=LX3p, delta=d0)
    term = 320 * r0['m'] * (1 + mpf('1e-4')) / (3 * log(1 + mpf('1e-4')))
    say("\n  bracket 320 m(X_3)(1+eps_0)/(3 log(1+eps_0)) = %s; x printed factor = %s; x correct factor = %s"
        % (nstr(term, 6), nstr(term * U_PRINTED, 6), nstr(term * U_CORRECT, 6)))
    say("  change of the case (b) margin from the factor alone: -(1/2)(%s - %s) = %s"
        % (nstr(term * U_CORRECT, 6), nstr(term * U_PRINTED, 6), nstr(-(term * (U_CORRECT - U_PRINTED)) / 2, 6)))
    say("  suprema in C.1 attained at the left end point: %s" % ENDPOINT_OK[0])

    # ---- threshold scans
    say("\n-- C.2 scans in loglog X_2 (X_3 = X_2^{1/8}, best delta at each point, alpha_i = 1e-5, eps_0 = 1e-4)")
    variants = [("printed", dict()),
                ("ufix_only", dict(ufix=True)),
                ("corrected", dict(ufix=True, eps2_arg="X3", **Cc)),
                ("corrected_sharp", dict(ufix=True, eps2_arg="X3", sharp=True, **Cc))]
    grid = [mpf(32400 + 2 * i) / 1000 for i in range(401)]        # 32.400 .. 33.200, step 0.002
    table = {}
    for name, kw in variants:
        for LL in grid:
            table[(name, LL)] = best_delta(LL, **kw)
    with open(os.path.join(DATA, "margin_scan.csv"), "w") as fh:
        fh.write("loglogX2," + ",".join("%s_margin,%s_delta" % (n, n) for n, _ in variants) + "\n")
        for LL in grid:
            fh.write(nstr(LL, 6) + "," + ",".join("%s,%s" % (nstr(table[(n, LL)][1]['res_b'], 10), nstr(table[(n, LL)][0], 5))
                                                   for n, _ in variants) + "\n")
    show_pts = [mpf(x) / 1000 for x in (32500, 32600, 32700, 32740, 32750, 32800, 32900, 33000)]
    for name, _ in variants:
        say("  %s:" % name)
        for LL in show_pts:
            d, r = table[(name, LL)]
            say("    loglog X2 = %s  best delta = %s  eps1=1/%d eps2=1/%d  m=%s  case(b) = %s  case(a) = %s"
                % (nstr(LL, 5), nstr(d, 5), r['pb1'] - 2, r['pb2'] - 2, nstr(r['m'], 3), nstr(r['res_b'], 6), nstr(r['res_a'], 5)))

    # ---- thresholds
    say("\n-- C.3 thresholds (grid step 0.002 on [32.4, 33.2], refined to 0.0005; X_3 = X_2^{1/8}; best delta)")
    for name, kw in variants:
        for target in [mpf(0), mpf('2e-4')]:
            above = [LL for LL in grid if table[(name, LL)][1]['res_b'] > target]
            if not above:
                say("  %-16s target %s: not reached by 33.2" % (name, nstr(target, 2))); continue
            # last grid point at which the margin is <= target
            below = [LL for LL in grid if table[(name, LL)][1]['res_b'] <= target]
            if not below:
                say("  %-16s target %s: already above at 32.4" % (name, nstr(target, 2))); continue
            lastbad = max(below)
            LL = lastbad
            while best_delta(LL, **kw)[1]['res_b'] <= target:
                LL += mpf('0.0005')
            d, r = best_delta(LL, **kw)
            say("  %-16s margin > %s for every grid point >= %s; first such point at resolution 0.0005: loglog X_2 = %s"
                " (delta = %s, margin %s)" % (name, nstr(target, 2), nstr(min(x for x in above if x > lastbad), 6),
                                              nstr(LL, 6), nstr(d, 5), nstr(r['res_b'], 6)))

    # ---- statements actually certified
    say("\n-- C.4 certified points")
    ENDPOINT_OK[0] = True
    LL = mpf('32.72')
    b = best_delta(LL, ufix=True, eps2_arg="X3", **Cc)
    show("corrected, existing lemmas, loglog X2 = 32.72, X3 = X2^{1/8}, best delta", b[1], b[0])
    LL = mpf('32.75')
    for d in [mpf('1.476')]:
        r = chen(LLX2=LL, delta=d, ufix=True, eps2_arg="X3", **Cc)
        show("corrected, existing lemmas, loglog X2 = 32.75, X3 = X2^{1/8}", r, d)
        r = chen(LLX2=LL, LX3=exp(mpf('30.67')), delta=d, ufix=True, eps2_arg="X3", **Cc)
        show("same with X3 = exp(exp(30.67))", r, d)
    LL = mpf('32.7')
    b = best_delta(LL, ufix=True, eps2_arg="X3", sharp=True, **Cc)
    show("corrected + sharpened Lemma 41, loglog X2 = 32.7, X3 = X2^{1/8}, best delta", b[1], b[0])
    r = chen(LLX2=LL, LX3=exp(mpf('30.62')), delta=b[0], ufix=True, eps2_arg="X3", sharp=True, **Cc)
    show("same with X3 = exp(exp(30.62)) as in the paper", r, b[0])
    for d in [mpf('1.462')]:
        r = chen(LLX2=LL, delta=d, ufix=True, eps2_arg="X3", sharp=True, **Cc)
        show("corrected + sharpened Lemma 41 at delta = 1.462", r, d)
        r = chen(LLX2=LL, delta=d, ufix=True, eps2_arg="X3", sharp=True)
        show("as above but with printed C = 106/107", r, d)
    say("  all suprema over y >= x (or N >= X_2) in C.4 attained at the left end point: %s" % ENDPOINT_OK[0])
    say("  sharpened Lemma 41, lower range (alpha_i = 1e-5 kept; None = hypotheses of Thm 44 fail):")
    for x in [3220, 3224, 3230, 3240]:
        LL = mpf(x) / 100
        bb = best_delta(LL, ufix=True, eps2_arg="X3", sharp=True, **Cc)
        say("    loglog X2 = %s: %s" % (nstr(LL, 5), "None" if bb is None else "delta %s, case (b) %s, case (a) %s"
                                        % (nstr(bb[0], 5), nstr(bb[1]['res_b'], 6), nstr(bb[1]['res_a'], 6))))
    say("  hypotheses of Theorems 44, 50, 51 at (32.72, 1.4767), (32.75, 1.476) and (32.7, %s): %s, %s, %s"
        % (nstr(b[0], 5), conditions_ok(mpf('32.72'), exp(mpf('32.72')) / 8, mpf('1.4767')),
           conditions_ok(mpf('32.75'), exp(mpf('32.75')) / 8, mpf('1.476')),
           conditions_ok(mpf('32.7'), exp(mpf('32.7')) / 8, b[0])))
    # Theorem 5 with threshold 32.75: largest k with prod_{i<=k} p_i <= e^{e^{32.75}}: theta(p_k) ~ p_k
    for LL in [mpf('32.7'), mpf('32.75')]:
        x = exp(LL)   # theta(p_k) <= e^{LL};  k = pi(x) ~ li(x)
        say("  Theorem 5 check: loglog bound %s -> log pi(e^{%s}) ~ log li = %s (< 29.3)" % (nstr(LL, 4), nstr(LL, 4), nstr(log(mp.li(x)), 6)))


if __name__ == "__main__":
    part_A()
    C1c, C2c = part_B()
    part_C(C1c, C2c)
    with open(os.path.join(DATA, "verify_bjs_chen.txt"), "w") as fh:
        fh.write("\n".join(OUT) + "\n")
