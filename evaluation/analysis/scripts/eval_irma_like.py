#!/usr/bin/env python3
"""Evaluate GREMA on the IRMA-like simulations (sim_irma_like.py): detection and signs for the 8 combinations of the three switches, and the effect of every switch.

usage: eval_irma_like.py SIM_DIR     (SIM_DIR holds p<P>i<I>t<T>/seed<k>/{gold,signs}.tsv and runs/p<P>i<I>t<T>_s<k>/**/final_results.txt)
Per condition: AUROC of VOTE_SCORE over the 20 pairs (7 true regulations); the 6 regulations that are transcriptional in every condition (all but GAL80 -> GAL4)
called as a regulation (`+` or `-`) and the share of those with the true sign; GAL80 -> GAL4 separately (a protein-level inhibition when --posttrans 1).
Then the effect of every switch: the mean of the sign accuracy and of the AUROC over the conditions with the switch on, minus with the switch off.
"""
import glob, itertools, os, sys
import numpy as np
from sklearn.metrics import roc_auc_score

O = sys.argv[1].rstrip("/")
GAL80_GAL4 = ("G4", "G2")
rows = {}
for p, i, t in itertools.product([0, 1], repeat=3):
    cond = "p%di%dt%d" % (p, i, t)
    auc, called, right, g_called, g_right = [], 0, 0, 0, 0
    for s in range(1, 11):
        d = "%s/%s/seed%d" % (O, cond, s)
        gold = {}
        for line in open(d + "/gold.tsv"):
            f = line.split()
            if f[0] != f[1]:
                gold[(f[0], f[1])] = int(f[2])
        sign = {(f[0], f[1]): int(f[2]) for f in (l.split() for l in open(d + "/signs.tsv"))}
        fn = sorted(glob.glob("%s/runs/%s_s%d/**/final_results.txt" % (O, cond, s), recursive=True))[0]
        res = {}
        for line in open(fn):
            f = line.rstrip("\n").split("\t")
            if f[0] != "TF":
                res[(f[0], f[1])] = (f[2], float(f[4]))
        keys = sorted(gold)
        auc.append(roc_auc_score([gold[k] for k in keys], [res[k][1] for k in keys]))
        for e, sg in sign.items():
            r = res[e][0]
            if r in ("+", "-"):
                ok = (r == "+") == (sg > 0)
                if e == GAL80_GAL4:
                    g_called += 1; g_right += ok
                else:
                    called += 1; right += ok
    rows[(p, i, t)] = (np.mean(auc), called, right, g_called, g_right)
print("%-10s %-30s %-10s %-34s %-22s" % ("condition", "(protein, input, post-translational)", "AUROC", "6 transcriptional edges: called, sign right", "GAL80->GAL4: called, right"))
for (p, i, t), (a, c, r, gc, gr) in rows.items():
    print("%-10s %-30s %.3f      %3d / 60 called, sign right %3d (%3.0f%%)       %2d / 10, right %2d" % (
        "p%di%dt%d" % (p, i, t), "(%s, %s, %s)" % ("on" if p else "off", "on" if i else "off", "on" if t else "off"), a, c, r, 100.0 * r / max(c, 1), gc, gr))
print()
names = ["protein stage (delay)", "external input", "post-translational inhibition"]
for k, name in enumerate(names):
    on = [v for key, v in rows.items() if key[k] == 1]
    off = [v for key, v in rows.items() if key[k] == 0]
    acc = lambda vs: 100.0 * sum(v[2] for v in vs) / max(sum(v[1] for v in vs), 1)
    print("%-30s sign accuracy: on %3.0f%%, off %3.0f%%;   AUROC: on %.3f, off %.3f" % (name, acc(on), acc(off), np.mean([v[0] for v in on]), np.mean([v[0] for v in off])))
