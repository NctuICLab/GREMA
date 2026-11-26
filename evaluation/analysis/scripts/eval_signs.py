#!/usr/bin/env python3
"""Detection and signs of GREMA on the simulated signed networks (sim_signed_net.py).

usage: eval_signs.py DATA_DIR RUN_GLOB [RUN_GLOB ...]
  DATA_DIR  holds gold.tsv and signs.tsv of one simulated network
  RUN_GLOB  final_results.txt of the GREMA runs of that network
Per run: AUROC / AUPRC of VOTE_SCORE; among the true edges, the share called as a regulation (REGULATORY + or -) and, among those, the share with the
true sign (random: 0.5); and the same for the K best edges by vote (K = number of true edges): precision and sign accuracy of the true edges among them.
"""
import sys, glob, os
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

d = sys.argv[1]
gold = {}
for line in open(os.path.join(d, "gold.tsv")):
    f = line.split()
    if f[0] != f[1]:
        gold[(f[0], f[1])] = int(f[2])
sign = {(f[0], f[1]): int(f[2]) for f in (l.split() for l in open(os.path.join(d, "signs.tsv")))}
keys = sorted(gold)
y = np.array([gold[k] for k in keys])
K = int(y.sum())
rows = []
for pat in sys.argv[2:]:
    for fn in sorted(glob.glob(pat)):
        res = {}
        for line in open(fn):
            f = line.rstrip("\n").split("\t")
            if f[0] != "TF":
                res[(f[0], f[1])] = (f[2], float(f[4]))
        vote = np.array([res[k][1] for k in keys])
        called = [k for k in sign if res[k][0] in "+-" and res[k][0] != "0"]
        correct = [k for k in called if (res[k][0] == "+") == (sign[k] > 0)]
        top = [keys[i] for i in np.argsort(-vote, kind="stable")[:K]]
        top_true = [k for k in top if gold[k] == 1]
        top_correct = [k for k in top_true if (res[k][0] == "+") == (sign[k] > 0) and res[k][0] in "+-"]
        rows.append((roc_auc_score(y, vote), average_precision_score(y, vote), len(called) / len(sign),
                     len(correct) / max(len(called), 1), len(top_true) / K, len(top_correct) / max(len(top_true), 1), len(called), len(correct), len(top_true), len(top_correct)))
r = np.array(rows)
print("%s: %d true edges (%d activating), random AUPRC %.3f, %d runs" % (os.path.basename(d.rstrip("/")), K, sum(1 for v in sign.values() if v > 0), y.mean(), len(r)))
print("  vote AUROC %.3f / AUPRC %.3f;  true edges called as a regulation %.2f, of those with the true sign %.2f (%d of %d);  top-%d by vote: %.2f true, sign right %.2f (%d of %d)" % (
    r[:, 0].mean(), r[:, 1].mean(), r[:, 2].mean(), r[:, 3].mean(), r[:, 7].sum(), r[:, 6].sum(), K, r[:, 4].mean(), r[:, 5].mean(), r[:, 9].sum(), r[:, 8].sum()))
