#!/usr/bin/env python3
"""Simulated data with an external input and a known truth (HFODE, the IRMA topology of sim_gen.py plus a stimulus U1).

usage: sim_input_net.py OUT_DIR SEED [--noise 0.05] [--points 16]
Network: the 8 gene-gene regulations of sim_gen.py (6 activating, 2 repressing) and the input U1 -> G1 (activation), U1 -> G4 (repression).
dx_i/dt = TransMax_i prod_j 1/(1+(k_ij/x_j)^e_ij) prod_m 1/(1+(k_im/u_m(t))^e_im) + b_i - degrade_i x_i, e = n * sign, Euler step 1 like GREMA; random initial values U(0.05, 0.9).
Conditions (series): 3 training series with the stimulus (start, level) = (3, 1.0), (6, 0.6), (4, 0.8) and 1 held-out condition (5, 0.7); u = 0.01 before the start.
Noise: multiplicative, --noise; the genes of all conditions are scaled with one linear map to [0.01, 1]; the input is already on that scale.
Writes (every genes file has the rows of the genes and, where said, U1):
  train_input.txt      3 series, genes + U1          train_noinput.txt   the same without U1
  train_const.txt      U1 = 0.5 in every series (a constant input)       train_shuf.txt   U1 of series 1,2,3 given to series 2,3,1 (the stimulus shuffled between the series)
  heldout.txt          the held-out condition (genes observed + its U1)
  know_input.txt / know_noinput.txt   no prior knowledge (self = 0, every other pair '?', U1 -> every gene '?' in know_input)
  gold.tsv  TF GENE gold (all gene pairs and the U1 pairs)    signs.tsv   the true edges with their signs
"""
import sys, os
import numpy as np

out, seed = sys.argv[1], int(sys.argv[2])
noise = float(sys.argv[sys.argv.index("--noise") + 1]) if "--noise" in sys.argv else 0.05
T = int(sys.argv[sys.argv.index("--points") + 1]) if "--points" in sys.argv else 16
os.makedirs(out, exist_ok=True)
N = 5
edges = {(1, 2): 1, (2, 4): 1, (2, 3): 1, (3, 5): 1, (3, 1): 1, (3, 4): 1, (4, 2): -1, (5, 1): -1}
in_edges = {1: 1, 4: -1}  # U1 -> G1 (+), U1 -> G4 (-)
conds = [(3, 1.0), (6, 0.6), (4, 0.8), (5, 0.7)]  # the last one is held out
rng = np.random.default_rng(seed)


def u_series(start, level):
    u = np.full(T, 0.01); u[start:] = level
    return u


while True:
    tm = rng.uniform(.4, 1, N); b = rng.uniform(0, .05, N); deg = rng.uniform(.2, .5, N)
    k = rng.uniform(.2, .6, (N, N)); n = rng.uniform(2, 4, (N, N))
    ku = rng.uniform(.2, .6, N); nu = rng.uniform(2, 4, N)

    def sim(x0, u):
        X = np.zeros((N, T)); X[:, 0] = x0
        for t in range(T - 1):
            for i in range(N):
                p = 1.0
                for j in range(N):
                    if (j + 1, i + 1) in edges:
                        p *= 1 / (1 + (k[i, j] / X[j, t]) ** (n[i, j] * edges[(j + 1, i + 1)]))
                if (i + 1) in in_edges:
                    p *= 1 / (1 + (ku[i] / u[t]) ** (nu[i] * in_edges[i + 1]))
                X[i, t + 1] = min(max(X[i, t] + tm[i] * p + b[i] - deg[i] * X[i, t], 0.01), 1.5)
        return X
    runs = [sim(rng.uniform(.05, .9, N), u_series(*c)) for c in conds]
    allx = np.concatenate(runs, axis=1)
    if (allx.max(1) - allx.min(1)).min() > 0.25:  # every gene has visible dynamics
        break
lo, hi = allx.min(), allx.max()
obs = []
for X in runs:
    Y = X * (1 + noise * rng.standard_normal(X.shape))
    obs.append(np.clip(0.01 + (Y - lo) / (hi - lo) * 0.99, 0.01, 1))
U = [u_series(*c) for c in conds]


def write(path, series, inputs=None):
    with open(path, "w") as f:
        f.write(f"repeat_number={len(series)}\ntimepoint_number={T}\n")
        for r, X in enumerate(series):
            for i in range(N):
                f.write(f"rep{r + 1}\tG{i + 1}\t" + "\t".join("%.6f" % v for v in X[i]) + "\n")
            if inputs is not None:
                f.write(f"rep{r + 1}\tU1\t" + "\t".join("%.6f" % v for v in inputs[r]) + "\n")


write(f"{out}/train_input.txt", obs[:3], U[:3])
write(f"{out}/train_noinput.txt", obs[:3])
write(f"{out}/train_const.txt", obs[:3], [np.full(T, 0.5)] * 3)
write(f"{out}/train_shuf.txt", obs[:3], [U[1], U[2], U[0]])
write(f"{out}/heldout.txt", obs[3:], U[3:])
for name, withinput in (("know_input", True), ("know_noinput", False)):
    with open(f"{out}/{name}.txt", "w") as f:
        f.write("TF\tGENE\tREGULATORY\n")
        for g in range(1, N + 1):
            for t in range(1, N + 1):
                f.write(f"G{t}\tG{g}\t{'0' if t == g else '?'}\n")
            if withinput:
                f.write(f"U1\tG{g}\t?\n")
with open(f"{out}/gold.tsv", "w") as f:
    for g in range(1, N + 1):
        for t in range(1, N + 1):
            if t != g:
                f.write(f"G{t}\tG{g}\t{int((t, g) in edges)}\n")
        f.write(f"U1\tG{g}\t{int(g in in_edges)}\n")
with open(f"{out}/signs.tsv", "w") as f:
    for (t, g), s in edges.items():
        f.write(f"G{t}\tG{g}\t{s:+d}\n")
    for g, s in in_edges.items():
        f.write(f"U1\tG{g}\t{s:+d}\n")
