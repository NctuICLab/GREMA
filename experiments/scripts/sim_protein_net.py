#!/usr/bin/env python3
"""Simulated signed network in which the regulator is a protein, but the measured variable is the promoter activity (the rate of production).

usage: sim_protein_net.py OUT_DIR SEED [--genes 10] [--series 4] [--points 20] [--noise 0.03] [--driver none|transcription|cleavage]
Rules (fixed before any run):
  Network as in sim_signed_net.py: every gene has 1 or 2 regulators from the other genes, the sign of every edge is + or - (probability 1/2).
  Every gene j has a protein P_j that follows its promoter activity:  dP_j/dt = d_j (a_j - P_j),  d_j in U(0.1, 0.3)  (Euler, step 1, time constant 3-10 steps);
  the promoter activity of gene i is a function of the PROTEINS of its regulators:
      a_i = TransMax_i * prod_j 1 / (1 + (k_ij / P_j)^(s_ij n_ij)) + b_i,   n in U(2, 4), k in U(0.2, 0.6), TransMax in U(0.6, 1.0), b in U(0, 0.05).
  Several series start from random initial proteins P_j(0) in U(0.05, 0.9). The measured data are the promoter activities a_i(t) (the proteins are not
  measured), with multiplicative Gaussian noise and one global scaling to [0.01, 1]. The draw is repeated until every gene varies by more than 0.15.
  --driver (a shared external pulse u_s(t) = s exp(-t/6), s = 0.5, 1, 1.5, 2 for the 4 series; all series start at the equilibrium of the proteins):
    transcription: every promoter activity is multiplied by (1 + u_s(t)); the proteins still follow the promoter activity;
    cleavage: every protein is also cleaved, dP_j/dt = d_j (a_j - P_j) - 0.4 u_s(t) P_j: the protein falls, the promoters it represses rise (like LexA after UV);
    none (default): no pulse, the series differ by the random initial proteins (as before).
Writes OUT_DIR/expression.txt, knowledge.txt (every pair '?', diagonal '0'), gold.tsv, signs.tsv and decay.tsv (the d_j).
"""
import sys, os
import numpy as np

out, seed = sys.argv[1], int(sys.argv[2])
opt = {"--genes": 10, "--series": 4, "--points": 20, "--noise": 0.03, "--driver": "none"}
for i in range(3, len(sys.argv), 2):
    opt[sys.argv[i]] = type(opt[sys.argv[i]])(sys.argv[i + 1])
N, R, T, noise = opt["--genes"], opt["--series"], opt["--points"], opt["--noise"]
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(seed)
edges = {}
for i in range(N):
    others = [j for j in range(N) if j != i]
    for j in rng.choice(others, size=int(rng.integers(1, 3)), replace=False):
        edges[(int(j), i)] = int(rng.choice([-1, 1]))
for _ in range(1000):
    tm = rng.uniform(.6, 1.0, N); b = rng.uniform(0, .05, N); d = rng.uniform(.1, .3, N)
    k = rng.uniform(.2, .6, (N, N)); n = rng.uniform(2, 4, (N, N))
    drv = opt["--driver"]
    def step(P, amp, t):
        u = amp * np.exp(-t / 6.0) if drv != "none" else 0.0
        a = np.zeros(N)
        for i in range(N):
            pr = 1.0
            for (j, ii), s in edges.items():
                if ii == i:
                    pr *= 1 / (1 + (k[i, j] / max(P[j], 0.01)) ** (s * n[i, j]))
            a[i] = tm[i] * pr * ((1 + u) if drv == "transcription" else 1.0) + b[i]
        Pn = P + d * (a - P) - (0.4 * u * P if drv == "cleavage" else 0.0)
        return a, np.maximum(Pn, 0.0)
    def sim(p0, amp):
        P = p0.copy(); A = np.zeros((N, T))
        for t in range(T):
            A[:, t], P = step(P, amp, t)
        return A
    if drv == "none":
        runs = [sim(rng.uniform(.05, .9, N), 0.0) for _ in range(R)]
    else:
        P0 = np.full(N, 0.5)
        for _ in range(100):                       # equilibrium of the proteins without the pulse (amplitude 0)
            a0, P0 = step(P0, 0.0, 0)
        runs = [sim(P0, [0.5, 1.0, 1.5, 2.0, 2.5, 3.0][r]) for r in range(R)]
    allx = np.concatenate(runs, axis=1)
    if (allx.max(1) - allx.min(1)).min() > 0.15:
        break
lo, hi = allx.min(), allx.max()
with open(out + "/expression.txt", "w") as f:
    f.write("repeat_number=%d\ntimepoint_number=%d\n" % (R, T))
    for r, A in enumerate(runs):
        A = A * (1 + noise * rng.standard_normal(A.shape))
        A = np.clip(0.01 + (A - lo) / (hi - lo) * 0.99, 0.01, 1)
        for i in range(N):
            f.write("rep%d\tG%d\t" % (r + 1, i + 1) + "\t".join("%.6f" % v for v in A[i]) + "\n")
names = ["G%d" % (i + 1) for i in range(N)]
with open(out + "/knowledge.txt", "w") as f:
    f.write("TF\tGENE\tREGULATORY\n")
    for a in range(N):
        for b2 in range(N):
            f.write("%s\t%s\t%s\n" % (names[a], names[b2], "0" if a == b2 else "?"))
with open(out + "/gold.tsv", "w") as f:
    for a in range(N):
        for b2 in range(N):
            f.write("%s\t%s\t%d\n" % (names[a], names[b2], 1 if (a, b2) in edges else 0))
with open(out + "/signs.tsv", "w") as f:
    for (a, b2), s in sorted(edges.items()):
        f.write("%s\t%s\t%+d\n" % (names[a], names[b2], s))
with open(out + "/decay.tsv", "w") as f:
    for i in range(N):
        f.write("%s\t%.4f\n" % (names[i], d[i]))
print("seed %d: %d edges (%d activating, %d repressing)" % (seed, len(edges), sum(1 for v in edges.values() if v > 0), sum(1 for v in edges.values() if v < 0)))
