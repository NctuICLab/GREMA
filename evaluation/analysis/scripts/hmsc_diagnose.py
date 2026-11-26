#!/usr/bin/env python3
"""Diagnostics for the hMSC design: is GREMA's score only the correlation of the profiles, and where do the gold pairs rank?

usage: hmsc_diagnose.py DESIGN_DIR RUNS_GLOB
Over all TF -> gene pairs of DESIGN_DIR/gold.tsv: Spearman between the mean VOTE_SCORE of the runs and |Pearson| of the two 5-point profiles; between runs (reproducibility);
for every gold pair: the rank of its vote score and of its |Pearson| among all pairs (1 = best), the TRRUST sign and the sign called by GREMA (majority of the runs).
"""
import sys, glob, collections
import numpy as np
from scipy.stats import spearmanr, rankdata

D, g = sys.argv[1], sys.argv[2]
names = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{D}/names.tsv")}
expr = {}
for l in open(f"{D}/expression.txt"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "rep1":
        expr[f[1]] = np.array([float(x) for x in f[2:]])
rows = [l.rstrip("\n").split("\t") for l in open(f"{D}/gold.tsv")][1:]
keys = [(r[0], r[1]) for r in rows]
runs = []
for p in sorted(glob.glob(g)):
    d = {}
    for l in open(p):
        f = l.rstrip("\n").split("\t")
        if f[0] == "TF":
            continue
        d[(f[0], f[1])] = (float(f[4]), f[2])
    runs.append(d)
V = np.array([[r.get(k, (0.0, "0"))[0] for k in keys] for r in runs])
vote = V.mean(0)
pear = np.array([abs(np.corrcoef(expr[a], expr[b])[0, 1]) for a, b in keys])
print(f"runs {len(runs)}, pairs {len(keys)}")
rr = [spearmanr(V[i], V[j])[0] for i in range(len(runs)) for j in range(i + 1, len(runs))]
print(f"Spearman between runs (reproducibility): {np.mean(rr):.3f}")
print(f"Spearman of the mean vote with |Pearson|: {spearmanr(vote, pear)[0]:.3f}")
rv = len(keys) + 1 - rankdata(vote, method="average"); rp = len(keys) + 1 - rankdata(pear, method="average")
print("gold pairs (TF -> gene, TRRUST sign, GREMA sign by majority, rank of vote, rank of |Pearson|, |Pearson|):")
for i, r in enumerate(rows):
    if r[2] == "1":
        sg = collections.Counter(run.get(keys[i], (0, "0"))[1] for run in runs).most_common(1)[0][0]
        print(f"  {names[r[0]]:8s} -> {names[r[1]]:10s} {r[4] or '?':11s} {sg:2s} vote rank {rv[i]:6.0f}  pearson rank {rp[i]:6.0f}  {pear[i]:.2f}")
