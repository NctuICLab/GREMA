#!/usr/bin/env python3
"""Candidate regulators from a sparse linear ODE with stability selection (Inferelator / TIGRESS style), to restrict the knowledge file of GREMA.

usage: lasso_prefilter.py EXPRESSION.txt KNOW_IN.txt OUT_DIR [--keep-frac 0.5] [--boots 200] [--seed 0]
  EXPRESSION.txt  GREMA time-series file (every replicate is a series)
  KNOW_IN.txt     knowledge file; the pairs with '?' are the candidate regulators (self pairs are not candidates)
Rules (fixed before any evaluation, no gold standard is used):
  For every target gene i: y = x_i(t+1) - x_i(t) (inside every series), features = x_j(t) of the candidates j != i, all columns standardised.
  The own x_i(t) (degradation) is removed from y by least squares first. In each of --boots subsamples (50% of the rows) every feature is multiplied
  by a random weight in {0.5, 1} (randomised lasso) and the LARS path is followed for L = K steps, K = ceil(--keep-frac * candidates);
  the score of a candidate is the fraction of subsamples in which it entered the path within the first L steps (ties: |correlation| with y, scaled to 1e-6).
  The K candidates with the highest score of every target are kept.
Writes OUT_DIR/prior_scores.txt (final_results format: TF GENE REGULATORY CONFIDENCE_LEVEL VOTE_SCORE, the score in both columns),
OUT_DIR/knowledge_filtered.txt ('?' for the kept pairs, '0' for all other pairs; fixed '+', '-', '0' of the input are kept) and OUT_DIR/summary.txt.
"""
import sys, os, math
import numpy as np
from sklearn.linear_model import lars_path

args = sys.argv[1:]
expr, know_in, outdir = args[:3]
opt = {"--keep-frac": "0.5", "--boots": "200", "--seed": "0"}
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
    m = np.array([series[r][g] for g in genes]).T          # time x genes
    X.append(m[:-1])
    DX.append(m[1:] - m[:-1])
X, DX = np.vstack(X), np.vstack(DX)
Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1)

know = {}
fixed = {}
for line in open(know_in):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 3 and p[0] != "TF":
        know[(p[0], p[1])] = p[2].strip()
cand = {g: [t for t in genes if t != g and know.get((t, g), "?") == "?"] for g in genes}

B = int(opt["--boots"])
score = {}
keep = {}
for g in genes:
    cs = cand[g]
    if not cs:
        keep[g] = []
        continue
    y = DX[:, idx[g]]
    xi = Xs[:, idx[g]]
    y = y - xi * (xi @ y) / (xi @ xi) if xi @ xi > 0 else y            # remove the degradation term
    F = Xs[:, [idx[t] for t in cs]]
    K = max(1, math.ceil(float(opt["--keep-frac"]) * len(cs)))
    n = len(y)
    freq = np.zeros(len(cs))
    for b in range(B):
        rows = rng.choice(n, size=max(2, n // 2), replace=False)
        w = rng.choice([0.5, 1.0], size=len(cs))
        try:
            _, active, _ = lars_path(F[rows] * w, y[rows], method="lasso", max_iter=K)
        except Exception:
            continue
        for a in active[:K]:
            freq[a] += 1
    freq /= B
    corr = np.abs([np.corrcoef(F[:, k], y)[0, 1] if F[:, k].std() > 0 and y.std() > 0 else 0 for k in range(len(cs))])
    tot = freq + 1e-6 * corr
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
