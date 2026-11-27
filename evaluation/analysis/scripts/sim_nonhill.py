#!/usr/bin/env python3
"""Simulated time series of the IRMA topology with regulation functions that are NOT the Hill product that GREMA fits (a test of the form of the model).

usage: sim_nonhill.py OUT_DIR FORM SEED --series S [--points 16] [--noise 0.05]
FORM  hill    the product of Hill factors (the model of GREMA: the reference)
      or      the factors of the regulators are combined by an OR gate, 1 - prod(1 - f_j), instead of the product (an AND gate); the same Hill factors
      bell    the edges G2 -> G3 and G3 -> G5 are non-monotone (a bell: 4 r / (1 + r)^2 with r = (x/k)^n, the largest effect at x = k, so the same regulator activates at low and represses at high levels);
              the other edges are Hill factors
      delay   the regulators act through a protein that follows the mRNA, dp/dt = a (x - p), a = 0.1-0.3; only the mRNA is measured
      linear  production = clip(c + sum of +- w x_j, 0, 2): no saturation, no threshold
Network: 8 regulations, G1 -> G2 +, G2 -> G4 +, G2 -> G3 +, G3 -> G5 +, G3 -> G1 +, G3 -> G4 +, G4 -> G2 -, G5 -> G1 - (as in sim_gen.py).
dx_i/dt = TransMax_i * q_i + b_i - degrade_i x_i, Euler step 1 like GREMA, x clipped to [0.01, 1.5]; TransMax 0.4-1, b 0-0.05, degrade 0.2-0.5, k 0.2-0.6, n 2-4;
random initial values U(0.05, 0.9) in every series; parameter sets are drawn again until every gene moves by more than 0.25 over all the series.
S training series and one held-out series (new initial values) of --points points; multiplicative noise --noise; one linear map of all values to [0.01, 1].
Writes train.txt, heldout.txt, know.txt (no prior knowledge: the diagonal 0, every other pair ?), gold.tsv (TF GENE 0/1 for the 20 pairs), signs.tsv (TF GENE sign type), summary.txt.
"""
import sys, os
import numpy as np

out, form, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
a = sys.argv
S = int(a[a.index("--series") + 1])
T = int(a[a.index("--points") + 1]) if "--points" in a else 16
noise = float(a[a.index("--noise") + 1]) if "--noise" in a else 0.05
assert form in ("hill", "or", "bell", "delay", "linear")
os.makedirs(out, exist_ok=True)
N = 5
edges = {(1, 2): 1, (2, 4): 1, (2, 3): 1, (3, 5): 1, (3, 1): 1, (3, 4): 1, (4, 2): -1, (5, 1): -1}
bell_edges = {(2, 3), (3, 5)} if form == "bell" else set()
regs = {i: [(j, s) for (j, t), s in edges.items() if t == i] for i in range(1, N + 1)}
rng = np.random.default_rng(seed)
tries = 0
while True:
    tries += 1
    tm = rng.uniform(.4, 1, N); b = rng.uniform(0, .05, N); deg = rng.uniform(.2, .5, N)
    k = rng.uniform(.2, .6, (N, N)); n = rng.uniform(2, 4, (N, N))
    ap = rng.uniform(.1, .3, N)                       # delay: rate of the protein
    w = rng.uniform(.3, .9, (N, N)); c = rng.uniform(.2, .5, N)  # linear

    def factor(i, j, s, y):
        if (j, i) in bell_edges:
            r = (y / k[i - 1, j - 1]) ** n[i - 1, j - 1]
            return 4 * r / (1 + r) ** 2
        return 1 / (1 + (k[i - 1, j - 1] / y) ** (n[i - 1, j - 1] * s))

    def sim(x0):
        X = np.zeros((N, T)); X[:, 0] = x0
        P = x0.copy()                                # protein levels (delay)
        for t in range(T - 1):
            src = P if form == "delay" else X[:, t]
            nxt = np.zeros(N)
            for i in range(1, N + 1):
                if form == "linear":
                    q = min(max(c[i - 1] + sum(s * w[i - 1, j - 1] * src[j - 1] for j, s in regs[i]), 0), 2)
                else:
                    fs = [factor(i, j, s, src[j - 1]) for j, s in regs[i]]
                    q = 1 - np.prod([1 - f for f in fs]) if form == "or" else np.prod(fs)
                nxt[i - 1] = min(max(X[i - 1, t] + tm[i - 1] * q + b[i - 1] - deg[i - 1] * X[i - 1, t], 0.01), 1.5)
            if form == "delay":
                P = P + ap * (X[:, t] - P)
            X[:, t + 1] = nxt
        return X
    runs = [sim(rng.uniform(.05, .9, N)) for _ in range(S + 1)]
    allx = np.concatenate(runs, axis=1)
    if (allx.max(1) - allx.min(1)).min() > 0.25:
        break
lo, hi = allx.min(), allx.max()
obs = []
for X in runs:
    Y = X * (1 + noise * rng.standard_normal(X.shape))
    obs.append(np.clip(0.01 + (Y - lo) / (hi - lo) * 0.99, 0.01, 1))


def write(path, series):
    with open(path, "w") as f:
        f.write(f"repeat_number={len(series)}\ntimepoint_number={T}\n")
        for r, X in enumerate(series):
            for i in range(N):
                f.write(f"rep{r + 1}\tG{i + 1}\t" + "\t".join("%.6f" % v for v in X[i]) + "\n")


write(f"{out}/train.txt", obs[:S]); write(f"{out}/heldout.txt", obs[S:])
with open(f"{out}/know.txt", "w") as f:
    f.write("TF\tGENE\tREGULATORY\n")
    for g in range(1, N + 1):
        for t in range(1, N + 1):
            f.write(f"G{t}\tG{g}\t{'0' if t == g else '?'}\n")
with open(f"{out}/gold.tsv", "w") as f:
    for g in range(1, N + 1):
        for t in range(1, N + 1):
            if t != g:
                f.write(f"G{t}\tG{g}\t{int((t, g) in edges)}\n")
with open(f"{out}/signs.tsv", "w") as f:
    for (t, g), s in edges.items():
        f.write(f"G{t}\tG{g}\t{s:+d}\t{'bell' if (t, g) in bell_edges else 'monotone'}\n")
with open(f"{out}/summary.txt", "w") as f:
    f.write(f"form {form}, seed {seed}, series {S}, points {T}, noise {noise}, parameter sets drawn {tries}\n")
