#!/usr/bin/env python3
"""Candidate regulators from HSIC-Lasso (nonlinear) with stability selection; the same rules as lasso_prefilter.py with the linear LASSO replaced.

usage: hsic_prefilter.py EXPRESSION.txt KNOW_IN.txt OUT_DIR [--keep-frac 0.5] [--boots 100] [--seed 0]
Needs pyHSICLasso (pixi: `pixi add --pypi pyHSICLasso`).
Rules (fixed before any evaluation, no gold standard is used):
  For every target gene i: y = x_i(t+1) - x_i(t) (inside every series), features = x_j(t) of the candidates j != i (self pairs are not candidates).
  The own x_i(t) (degradation) is removed from y by least squares first, as in lasso_prefilter.py (a linear step; the candidates act non-linearly).
  In each of --boots subsamples (50% of the rows) HSIC-Lasso (Gaussian kernels, B = 0, M = 1) selects K = ceil(--keep-frac * candidates) features;
  the score of a candidate is the fraction of subsamples in which it was selected (ties: |Spearman correlation| with y, scaled to 1e-6).
  The K candidates with the highest score of every target are kept.
Outputs as lasso_prefilter.py: OUT_DIR/prior_scores.txt, OUT_DIR/knowledge_filtered.txt, OUT_DIR/summary.txt.
"""
import sys, os, math
import numpy as np
from scipy.stats import spearmanr
from pyHSICLasso import HSICLasso

args = sys.argv[1:]
expr, know_in, outdir = args[:3]
opt = {"--keep-frac": "0.5", "--boots": "100", "--seed": "0"}
for i in range(3, len(args), 2):
    opt[args[i]] = args[i + 1]
os.makedirs(outdir, exist_ok=True)
rng = np.random.default_rng(int(opt["--seed"]))

genes, series = [], {}
for line in open(expr):
    p = line.rstrip("\n").split("\t")
    if len(p) > 3 and p[0].startswith("rep"):
        series.setdefault(p[0], {})[p[1]] = [float(x) for x in p[2:]]
        if p[1] not in genes:
            genes.append(p[1])
N = len(genes)
idx = {g: i for i, g in enumerate(genes)}
X, DX = [], []
for r in sorted(series, key=lambda s: int(s[3:])):
    m = np.array([series[r][g] for g in genes]).T
    X.append(m[:-1])
    DX.append(m[1:] - m[:-1])
X, DX = np.vstack(X), np.vstack(DX)
Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1)

know = {}
for line in open(know_in):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 3 and p[0] != "TF":
        know[(p[0], p[1])] = p[2].strip()
cand = {g: [t for t in genes if t != g and know.get((t, g), "?") == "?"] for g in genes}

B = int(opt["--boots"])
score, keep = {}, {}
for g in genes:
    cs = cand[g]
    if not cs:
        keep[g] = []
        continue
    y = DX[:, idx[g]]
    xi = Xs[:, idx[g]]
    y = y - xi * (xi @ y) / (xi @ xi) if xi @ xi > 0 else y
    F = Xs[:, [idx[t] for t in cs]]
    K = max(1, math.ceil(float(opt["--keep-frac"]) * len(cs)))
    n = len(y)
    freq = np.zeros(len(cs))
    for b in range(B):
        rows = rng.choice(n, size=max(4, n // 2), replace=False)
        h = HSICLasso()
        h.input(F[rows], y[rows], featname=[str(k) for k in range(len(cs))])
        try:
            h.regression(K, B=0, M=1, n_jobs=1)
        except Exception:
            continue
        for name in h.get_features():
            freq[int(name)] += 1
    freq /= B
    rho = np.abs([spearmanr(F[:, k], y)[0] if F[:, k].std() > 0 and y.std() > 0 else 0 for k in range(len(cs))])
    tot = freq + 1e-6 * np.nan_to_num(rho)
    for k, t in enumerate(cs):
        score[(t, g)] = tot[k]
    keep[g] = [cs[k] for k in np.argsort(-tot)[:K]]

with open(outdir + "/prior_scores.txt", "w") as f:
    f.write("TF\tGENE\tREGULATORY\tCONFIDENCE_LEVEL\tVOTE_SCORE\n")
    for t in genes:
        for g in genes:
            s = score.get((t, g), 0.0)
            f.write("%s\t%s\t%s\t%.6f\t%.6f\n" % (t, g, "+" if t != g else "0", s, s))
with open(outdir + "/knowledge_filtered.txt", "w") as f:
    f.write("TF\tGENE\tREGULATORY\n")
    for t in genes:
        for g in genes:
            v = know.get((t, g), "?")
            if v == "?":
                v = "?" if t in keep.get(g, []) else "0"
            f.write("%s\t%s\t%s\n" % (t, g, v))
with open(outdir + "/summary.txt", "w") as f:
    f.write("genes %d, series %d, rows %d, boots %d, keep-frac %s\n" % (N, len(series), len(y), B, opt["--keep-frac"]))
    f.write("kept candidates per target: " + " ".join("%s:%d/%d" % (g, len(keep[g]), len(cand[g])) for g in genes if cand[g]) + "\n")
print(open(outdir + "/summary.txt").read())
