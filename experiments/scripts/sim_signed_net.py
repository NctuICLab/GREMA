#!/usr/bin/env python3
"""Simulated time series of a random signed network (HFODE form), with the true regulations and their signs, to test the signs found by GREMA.

usage: sim_signed_net.py OUT_DIR SEED [--genes 10] [--series 4] [--points 20] [--noise 0.03] [--model product|additive] [--hidden 0]
Rules (fixed before any run):
  Network: every target gene has 1 or 2 regulators drawn from the other genes; the sign of every edge is + (activation) or - (repression) with probability 1/2.
  Dynamics (Euler, step 1, as the solver of GREMA):  dx_i/dt = TransMax_i * prod_j 1/(1 + (k_ij / x_j)^(s_ij n_ij)) + b_i - degrade_i x_i,
    s_ij = +1 / -1 is the sign of the edge, n_ij in U(2, 4), k_ij in U(0.2, 0.6), TransMax_i in U(0.6, 1.0), b_i in U(0, 0.05), degrade_i in U(0.2, 0.5);
    x is kept in [0.01, 1.5]. Several series start from random initial values in U(0.05, 0.9) (the perturbation of the experiments).
  The draw is repeated until every gene varies by more than 0.25 over all series. All values are scaled to [0.01, 1] (one global map), multiplicative
  Gaussian noise (--noise) is added before the scaling of every series.
Variants (misspecification): --model additive: the regulators of a gene act as the MEAN of their Hill terms instead of the product;
  --hidden H: H more genes take part in the network and the regulation but are not written (not measured); gold and signs only list the edges between measured genes.
Writes OUT_DIR/expression.txt (GREMA format), knowledge.txt (every pair '?', diagonal '0'), gold.tsv (TF GENE 0/1) and signs.tsv (TF GENE +1/-1 of the true edges).
"""
import sys, os
import numpy as np

out, seed = sys.argv[1], int(sys.argv[2])
opt = {"--genes": 10, "--series": 4, "--points": 20, "--noise": 0.03, "--model": "product", "--hidden": 0}
for i in range(3, len(sys.argv), 2):
    opt[sys.argv[i]] = type(opt[sys.argv[i]])(sys.argv[i + 1])
N, R, T, noise = opt["--genes"], opt["--series"], opt["--points"], opt["--noise"]
H = opt["--hidden"]
NT = N + H   # all genes of the simulation, the first N are measured
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(seed)
edges = {}
for i in range(NT):
    others = [j for j in range(NT) if j != i]
    for j in rng.choice(others, size=int(rng.integers(1, 3)), replace=False):
        edges[(int(j), i)] = int(rng.choice([-1, 1]))
for _ in range(1000):
    tm = rng.uniform(.6, 1.0, NT); b = rng.uniform(0, .05, NT); deg = rng.uniform(.2, .5, NT)
    k = rng.uniform(.2, .6, (NT, NT)); n = rng.uniform(2, 4, (NT, NT))
    def sim(x0):
        X = np.zeros((NT, T)); X[:, 0] = x0
        for t in range(T - 1):
            for i in range(NT):
                terms = [1 / (1 + (k[i, j] / max(X[j, t], 0.01)) ** (s * n[i, j])) for (j, ii), s in edges.items() if ii == i]
                p = float(np.mean(terms)) if opt["--model"] == "additive" else float(np.prod(terms))
                X[i, t + 1] = min(max(X[i, t] + tm[i] * p + b[i] - deg[i] * X[i, t], 0.01), 1.5)
        return X
    runs = [sim(rng.uniform(.05, .9, NT)) for _ in range(R)]
    allx = np.concatenate(runs, axis=1)
    if (allx.max(1) - allx.min(1)).min() > 0.25:
        break
runs = [X[:N] for X in runs]   # the hidden genes are not measured
allx = np.concatenate(runs, axis=1)
lo, hi = allx.min(), allx.max()
with open(out + "/expression.txt", "w") as f:
    f.write("repeat_number=%d\ntimepoint_number=%d\n" % (R, T))
    for r, X in enumerate(runs):
        X = X * (1 + noise * rng.standard_normal(X.shape))
        X = np.clip(0.01 + (X - lo) / (hi - lo) * 0.99, 0.01, 1)
        for i in range(N):
            f.write("rep%d\tG%d\t" % (r + 1, i + 1) + "\t".join("%.6f" % v for v in X[i]) + "\n")
names = ["G%d" % (i + 1) for i in range(N)]
with open(out + "/knowledge.txt", "w") as f:
    f.write("TF\tGENE\tREGULATORY\n")
    for a in range(N):
        for b2 in range(N):
            f.write("%s\t%s\t%s\n" % (names[a], names[b2], "0" if a == b2 else "?"))
edges = {e: v for e, v in edges.items() if e[0] < N and e[1] < N}
with open(out + "/gold.tsv", "w") as f:
    for a in range(N):
        for b2 in range(N):
            f.write("%s\t%s\t%d\n" % (names[a], names[b2], 1 if (a, b2) in edges else 0))
with open(out + "/signs.tsv", "w") as f:
    for (a, b2), s in sorted(edges.items()):
        f.write("%s\t%s\t%+d\n" % (names[a], names[b2], s))
print("seed %d: %d edges (%d activating, %d repressing)" % (seed, len(edges), sum(1 for v in edges.values() if v > 0), sum(1 for v in edges.values() if v < 0)))
