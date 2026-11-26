#!/usr/bin/env python3
"""IRMA-like simulation (the 5 genes and the 7 signed regulations of IRMA, mRNA measured) with three switches that separate what makes the signs of GREMA wrong on real data.

usage: sim_irma_like.py OUT_DIR SEED [--protein 0|1] [--input 0|1] [--posttrans 0|1] [--series 3] [--points 20] [--noise 0.05]
Genes G1 CBF1, G2 GAL4, G3 SWI5, G4 GAL80, G5 ASH1. Signed regulations (as input/IRMA/IRMA_signs.tsv): CBF1->GAL4 +, GAL4->SWI5 +, SWI5->ASH1 +, SWI5->CBF1 +, SWI5->GAL80 +, ASH1->CBF1 -, GAL80->GAL4 -.
Rules (fixed before any run). mRNA m_i: dm_i/dt = prod_i - degrade_i m_i (Euler, step 1), prod_i = TransMax_i prod_j 1/(1 + (k_ij/r_j)^(s_ij n_ij)) + b_i over the regulators j of gene i,
r_j = the regulator signal: the mRNA m_j (switch off), the protein p_j (--protein 1: dp_j/dt = d_j (m_j - p_j), d_j in U(0.1, 0.3)). n in U(2, 4), k in U(0.2, 0.6), TransMax in U(0.6, 1.0), b in U(0, 0.05), degrade in U(0.2, 0.5).
  --protein 1   the regulators act through the protein that follows the mRNA with a delay (the mRNA is measured)
  --input 1     an external input u(t) = 0 for t < 5 and 1 from t = 5 (the shift to galactose): the regulation of SWI5 by GAL4 acts with the factor 0.2 + 0.8 u(t)
  --posttrans 1 GAL80 does not repress the transcription of GAL4; it lowers the activity of the Gal4 protein: the signal of Gal4 in the regulation of SWI5 is r_Gal4 (1 - 1/(1 + (k/r_Gal80)^(-n))) (the more Gal80, the less Gal4 acts),
                so GAL80->GAL4 has no direct mRNA-level counterpart
The series start from random m_j (= p_j) in U(0.05, 0.9). The measured data are the mRNAs, with multiplicative Gaussian noise and one global scaling to [0.01, 1]. The draw is repeated until every gene varies by more than 0.15.
Writes OUT_DIR/expression.txt, knowledge.txt (every pair '?', diagonal '0'), gold.tsv and signs.tsv (the 7 regulations).
"""
import sys, os
import numpy as np

out, seed = sys.argv[1], int(sys.argv[2])
opt = {"--protein": 0, "--input": 0, "--posttrans": 0, "--series": 3, "--points": 20, "--noise": 0.05}
for i in range(3, len(sys.argv), 2):
    opt[sys.argv[i]] = type(opt[sys.argv[i]])(sys.argv[i + 1])
R, T, noise = opt["--series"], opt["--points"], opt["--noise"]
N = 5
# (regulator, target): sign
EDGES = {(0, 1): 1, (1, 2): 1, (2, 4): 1, (2, 0): 1, (2, 3): 1, (4, 0): -1, (3, 1): -1}
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(seed)
for _ in range(1000):
    tm = rng.uniform(.6, 1.0, N); b = rng.uniform(0, .05, N); deg = rng.uniform(.2, .5, N)
    k = rng.uniform(.2, .6, (N, N)); n = rng.uniform(2, 4, (N, N)); dp = rng.uniform(.1, .3, N)
    def hill(x, kk, e):
        return 1 / (1 + (kk / max(x, 0.01)) ** e)
    def sim(m0):
        m = m0.copy(); p = m0.copy(); M = np.zeros((N, T))
        for t in range(T):
            M[:, t] = m
            r = p if opt["--protein"] else m
            u = (1.0 if t >= 5 else 0.0)
            new = np.zeros(N)
            for i in range(N):
                pr = 1.0
                for (j, ii), s in EDGES.items():
                    if ii != i:
                        continue
                    if opt["--posttrans"] and (j, ii) == (3, 1):
                        continue                                   # GAL80 does not repress the transcription of GAL4
                    sig = r[j]
                    if opt["--posttrans"] and (j, ii) == (1, 2):
                        sig = r[1] * (1 - hill(r[3], k[3, 1], n[3, 1]))   # Gal4 activity lowered by Gal80
                    h = hill(sig, k[i, j], s * n[i, j])
                    if opt["--input"] and (j, ii) == (1, 2):
                        h = h * (0.2 + 0.8 * u)
                    pr *= h
                new[i] = m[i] + tm[i] * pr + b[i] - deg[i] * m[i]
            if opt["--protein"]:
                p = p + dp * (m - p)
            m = np.clip(new, 0.01, 1.5)
        return M
    runs = [sim(rng.uniform(.05, .9, N)) for _ in range(R)]
    allx = np.concatenate(runs, axis=1)
    if (allx.max(1) - allx.min(1)).min() > 0.15:
        break
lo, hi = allx.min(), allx.max()
with open(out + "/expression.txt", "w") as f:
    f.write("repeat_number=%d\ntimepoint_number=%d\n" % (R, T))
    for r_, A in enumerate(runs):
        A = A * (1 + noise * rng.standard_normal(A.shape))
        A = np.clip(0.01 + (A - lo) / (hi - lo) * 0.99, 0.01, 1)
        for i in range(N):
            f.write("rep%d\tG%d\t" % (r_ + 1, i + 1) + "\t".join("%.6f" % v for v in A[i]) + "\n")
names = ["G%d" % (i + 1) for i in range(N)]
with open(out + "/knowledge.txt", "w") as f:
    f.write("TF\tGENE\tREGULATORY\n")
    for a in range(N):
        for b2 in range(N):
            f.write("%s\t%s\t%s\n" % (names[a], names[b2], "0" if a == b2 else "?"))
with open(out + "/gold.tsv", "w") as f:
    for a in range(N):
        for b2 in range(N):
            f.write("%s\t%s\t%d\n" % (names[a], names[b2], 1 if (a, b2) in EDGES else 0))
with open(out + "/signs.tsv", "w") as f:
    for (a, b2), s in sorted(EDGES.items()):
        f.write("%s\t%s\t%+d\n" % (names[a], names[b2], s))
print("seed %d done" % seed)
