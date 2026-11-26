#!/usr/bin/env python3
"""Does GREMA use the time order of a single series? Compare the VOTE_SCORE of the candidate pairs between runs (no gold standard needed).

usage: hmsc_time_control.py KNOWLEDGE.txt REAL_GLOB SHUF_GLOB
Prints, over the candidate pairs ('?' in KNOWLEDGE): the Spearman correlation of the scores between two runs on the same data (reproducibility: mean over all pairs of runs),
between a run on the real data and a run on the time-shuffled data, and between two runs on two different shufflings; and the number of pairs called (REGULATORY + or -).
If real-vs-shuffled is as high as real-vs-real, the scores do not depend on the time order (they come from the levels / the prior); if it is lower, the time order is used.
"""
import sys, glob, itertools
import numpy as np
from scipy.stats import spearmanr


def load(path, cands):
    s = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if (f[0], f[1]) in cands:
            s[(f[0], f[1])] = (float(f[4]), f[2])
    return s


cands = {tuple(l.split("\t")[:2]) for l in open(sys.argv[1]) if l.rstrip("\n").split("\t")[2] == "?"}
real = [load(p, cands) for p in sorted(glob.glob(sys.argv[2]))]
shuf = [load(p, cands) for p in sorted(glob.glob(sys.argv[3]))]
keys = sorted(cands)


def vec(r):
    return np.array([r.get(k, (0.0, "0"))[0] for k in keys])


def mean_rho(A, B, same):
    pairs = list(itertools.combinations(range(len(A)), 2)) if same else list(itertools.product(range(len(A)), range(len(B))))
    return np.mean([spearmanr(vec(A[i]), vec(B[j]))[0] for i, j in pairs])


print(f"candidate pairs {len(keys)}; runs: real {len(real)}, shuffled {len(shuf)}")
print(f"called pairs per run: real {np.mean([sum(1 for v in r.values() if v[1] in '+-' and v[1] != '0') for r in real]):.0f}, shuffled {np.mean([sum(1 for v in r.values() if v[1] in '+-' and v[1] != '0') for r in shuf]):.0f}")
print(f"Spearman of VOTE_SCORE: real vs real {mean_rho(real, real, True):.3f}, shuffled vs shuffled (different shufflings) {mean_rho(shuf, shuf, True):.3f}, real vs shuffled {mean_rho(real, shuf, False):.3f}")
mr = np.mean([vec(r) for r in real], axis=0); ms = np.mean([vec(r) for r in shuf], axis=0)
print(f"Spearman of the mean scores: real vs shuffled {spearmanr(mr, ms)[0]:.3f}")
