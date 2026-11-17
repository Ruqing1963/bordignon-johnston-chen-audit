# A note on an explicit version of Chen's theorem

Ruqing Chen (GUT Geoservice Inc., Montreal). Repository: <https://github.com/Ruqing1963/bordignon-johnston-chen-audit>.
DOI: [10.5281/zenodo.24323781](https://doi.org/10.5281/zenodo.24323781).

This is a 5-page audit note on M. Bordignon, D. R. Johnston and V. Starichkova,
*An explicit version of Chen's theorem and the linear sieve*, Int. J. Number Theory 21 (2025), no. 10,
2497–2572, doi:10.1142/S1793042125501192, arXiv:2207.09452.

Line numbers (l.NNNN) refer to the arXiv v6 LaTeX source of the target paper (download it from
<https://arxiv.org/abs/2207.09452v6>; it is not redistributed here). Result numbers follow the arXiv v6 PDF, which
uses one global counter for all results. The journal version numbers results by section; section numbers and
Tables 1–3 are the same in both.

| short name used in the audit | arXiv v6 number | journal number |
|---|---|---|
| main theorem (exp(exp(32.7))) | Theorem 3 | Theorem 1.4 |
| square-free corollary | Corollary 4 | Corollary 1.5 |
| at most e^29.3 prime factors | Theorem 5 | Theorem 1.6 |
| linear sieve | Theorem 6 | Theorem 2.1 |
| c_n table, c_n <= 0.9607 | Lemmas 7, 8 | Lemmas 2.2, 2.3 |
| Lemmas on gamma, kappa, tau, tau', J | Lemmas 9, 10, 11, 12, 14 | Lemmas 2.4, 2.5, 2.7, 2.8, 2.10 |
| primes in progressions | Lemma 27 | Lemma 4.2 |
| Lemmas not re-derived | Lemmas 29–36, 38–43 | Lemmas 4.4–4.11, 5.2–5.7 |
| m(X_3) | Lemma 32 | Lemma 4.7 |
| pi_2 >= S(A) - ... | Lemma 37 | Lemma 5.1 |
| eps_i, q_1 >= pbar | Lemma 41 | Lemma 5.5 |
| S(A) | Theorem 44 | Theorem 6.1 |
| case (b) lemma using 1/(q_1-2) <= eps_i | Lemma 49 | Lemma 6.7 |
| sum S(A_q) | Theorem 50 | Theorem 7.1 |
| S(B) | Theorem 51 | Theorem 8.1 |

The journal numbers were read off a copy of the publisher's web page. In that copy, the passages discussed in
Findings 1, 2 and 4 read the same as in arXiv v6, apart from copy-editing; the entries of Tables 1–3 were not
included in the copy and were not compared.

## Files
- `bordignon_johnston_chen_audit.tex` and `.pdf`: the note. pdflatex is run twice and gives no warnings and no overfull or underfull boxes.
- `code/verify_bjs_chen.py`: a self-contained reproduction using mpmath and numpy. It takes about 4 minutes and writes `data/verify_bjs_chen.txt`,
  `data/margin_scan.csv` and `data/sieve_constants.csv`.
- `code/make_figures.py`: writes `figures/fig_bjs_chen_audit.pdf` from the CSV files. Run it after the verify script.

## Findings
1. **The factor e^{2 gamma} (material).** Theorem 51 (l.2510) bounds U_N^{-1} by (2e^{+gamma}C_2)^{-1} = 0.42524.
   The definition (l.280) only gives (2e^{-gamma}C_2)^{-1} = 1.34897. The term affected is 320 m(X_3)(1+eps_0)/(3 log(1+eps_0)) = 0.0131041,
   and it enters Lemma 37 with weight 1/2. With the correct factor the case (b) margin falls from +4.4753e-4 to -5.6047e-3.
   The same slip appears at l.2103, 2335, 2353 and 2357, but there it is multiplied by 1/log N, so it has no numerical effect.
   Nothing elsewhere compensates for it.
2. **eps_2 is evaluated at (X_3)^8 (l.2508, 2692, 2701).** The hypothesis of Lemma 41 only holds at X_3, because Y_j >= N^{1/8}.
   With the correct argument, eps_2 changes from 1/57 to 1/51, which costs another 2.67e-3. The paper itself writes eps_2(X_2^{1/8}) at l.2742.
3. **Linear-sieve constants.** gamma_{s0} at l.608 should be e^{s0-1.8}, not e^{s0-2}. No variant of the recursion that we tried reproduces Table 3.
   Lemma 14 as stated gives C(1/200) = 192.9/196.1 (printed: 153/154). J(499) >= 1 for 1/eps >= 219.4, so for those
   eps Lemma 14 gives no bound. Direct summation gives C = 132.7/136.1 at 1/eps = 1e5 (printed: 106/107), about 25% more.
   The effect on Chen's theorem is only -2.1e-5. Numerically, c_33 = 0.69008, slightly above the tabulated 0.69.
4. **Restoration.**
   (A) Keep all of the paper's lemmas and apply the three corrections, with X_3 = X_2^{1/8} and delta = 1.4767. Then the case (b) bound is 3.1e-4 and the case (a) bound is 0.124.
   So Theorem 3 holds for N > exp(exp(32.72)); on the scan grid, the target 2e-4 is first reached at loglog X_2 = 32.716.
   (B) Lemma 41 can be sharpened by a one-line argument: q_1 >= p+, the smallest prime whose odd primorial is >= K.
   With this, the three corrections and delta = 1.4776, the case (b) bound at exp(exp(32.7)) is 4.10e-3, so Theorem 3 holds as stated.
   The bound stays above 3.2e-3 down to loglog X_2 = 32.24.
5. **Not checked:** the external inputs. These are Bordignon's Siegel-zero and PNT-in-AP bounds, Kadiri, Akbary–Hambrook,
   Broadbent et al., Büthe, Rosser–Schoenfeld, Montgomery–Vaughan, Johnston–Yang, and Nathanson's Lemma 9.7 and Theorems 9.3 and 10.2.
   The proofs of Lemmas 29–36 and 38–43 were also not re-derived; their constants were used as printed.

## Status of the earlier exploratory claims
- (i) Confirmed. With delta = 1.476 the corrected case (b) bound at loglog N = 32.75 is +5.2e-4, which is above 2e-4. The actual threshold is lower, at 32.716.
- (ii) Confirmed. The sharpened lemma is proved in the note. At 32.7 the bound is +3.21e-3 with delta = 1.462, and the best value is +4.10e-3 at delta = 1.4776.
- (iii) Partly confirmed. The C values, the gamma misprint, the divergence of Lemma 14 and an effect of about 1e-5 all hold. Not confirmed: the claim that Table 3 is reproduced with the printed gamma and xi_3 omitted. That variant still gives tau_n up to 10.5 for 3 <= n <= 10, above the tabulated 10.
