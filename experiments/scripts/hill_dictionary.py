#!/usr/bin/env python3
"""Hill-function dictionary with sparse regression and stability selection: a convex alternative to the genetic search of GREMA (prototype).

usage: hill_dictionary.py EXPRESSION.txt KNOW_IN.txt OUT_DIR [--keep-frac 0.5] [--boots 200] [--seed 0]
Model for every target gene i (additive; HFODE multiplies the regulators):
  dx_i/dt  =  c - d x_i + sum_j sum_(k,n) beta_(j,k,n) h(x_j; k, n),      h(x; k, n) = 1 / (1 + (k/x)^n)
  candidates j: the pairs with '?' of the knowledge file (self pairs are not candidates); n in {-4,-2,-1,1,2,4} (n < 0: repression type),
  k in {25, 50, 75 percentile of x_j over all series}  -> 18 columns per candidate (one group).
Rules (fixed before any evaluation, no gold standard is used):
  y = x_i(t+1) - x_i(t) inside every series. The constant and x_i(t) (degradation) are projected out of y and of every dictionary column by least squares;
  the columns are standardised. In each of --boots subsamples (50% of the rows, every column multiplied by a random weight in {0.5, 1}) the LARS path
  (lasso) is followed until K = ceil(--keep-frac * candidates) different candidates have entered; the score of a candidate is the fraction of subsamples
  in which it was among them (ties: largest |correlation| of one of its columns with y, scaled to 1e-6). The K best candidates of every target are kept.
Outputs as lasso_prefilter.py: OUT_DIR/prior_scores.txt, OUT_DIR/knowledge_filtered.txt, OUT_DIR/summary.txt.
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
NS = [-4, -2, -1, 1, 2, 4]

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
all_x = np.vstack([np.array([series[r][g] for g in genes]).T for r in series])

know = {}
for line in open(know_in):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 3 and p[0] != "TF":
        know[(p[0], p[1])] = p[2].strip()
cand = {g: [t for t in genes if t != g and know.get((t, g), "?") == "?"] for g in genes}


def hill(x, k, n):
    return 1.0 / (1.0 + (k / np.maximum(x, 1e-6)) ** n)


# dictionary columns of every gene j (independent of the target)
dict_cols = {}
for j, t in enumerate(genes):
    ks = np.percentile(all_x[:, j], [25, 50, 75])
    dict_cols[t] = np.column_stack([hill(X[:, j], k, n) for k in ks for n in NS])

B = int(opt["--boots"])
score, keep = {}, {}
for g in genes:
    cs = cand[g]
    if not cs:
        keep[g] = []
        continue
    Z = np.column_stack([np.ones(len(X)), X[:, idx[g]]])                 # constant and degradation
    proj = lambda M: M - Z @ np.linalg.lstsq(Z, M, rcond=None)[0]
    y = proj(DX[:, idx[g]][:, None])[:, 0]
    F = np.column_stack([dict_cols[t] for t in cs])
    F = proj(F)
    grp = np.repeat(np.arange(len(cs)), dict_cols[cs[0]].shape[1])
    sd = F.std(0)
    ok = sd > 1e-9
    F = F[:, ok] / sd[ok]
    grp = grp[ok]
    K = max(1, math.ceil(float(opt["--keep-frac"]) * len(cs)))
    n = len(y)
    freq = np.zeros(len(cs))
    for b in range(B):
        rows = rng.choice(n, size=max(3, n // 2), replace=False)
        w = rng.choice([0.5, 1.0], size=F.shape[1])
        try:
            _, active, _ = lars_path(F[rows] * w, y[rows], method="lasso", max_iter=min(F.shape[1], 6 * K + 6))
        except Exception:
            continue
        chosen = []
        for a in active:
            if grp[a] not in chosen:
                chosen.append(grp[a])
            if len(chosen) == K:
                break
        for c in chosen:
            freq[c] += 1
    freq /= B
    corr = np.zeros(len(cs))
    for c in range(len(cs)):
        cols = F[:, grp == c]
        if cols.shape[1] and y.std() > 0:
            corr[c] = np.nanmax(np.abs([np.corrcoef(cols[:, q], y)[0, 1] for q in range(cols.shape[1])]))
    tot = freq + 1e-6 * np.nan_to_num(corr)
    for c, t in enumerate(cs):
        score[(t, g)] = tot[c]
    keep[g] = [cs[c] for c in np.argsort(-tot)[:K]]

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
    f.write("genes %d, series %d, rows %d, boots %d, keep-frac %s, dictionary columns per candidate %d\n" % (N, len(series), len(y), B, opt["--keep-frac"], 3 * len(NS)))
print(open(outdir + "/summary.txt").read())
