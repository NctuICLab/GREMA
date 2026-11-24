#!/usr/bin/env python3
"""Does combining CONFIDENCE_LEVEL and VOTE_SCORE rank the edges better than either alone?

usage: combine_scores.py --gold GOLD.tsv [--self-pairs] LABEL=GLOB [LABEL=GLOB ...]

The two rules were fixed before looking at any result (no tuning on the gold standard):
  mean  = (conf + vote) / 2          (both are in [0, 1]; conf is 0 for a pair without regulation)
  rank  = mean of the two average ranks (ties get the average rank)
Printed for each group: per-run mean AUROC / AUPRC (sd of the AUROC) for conf, vote, mean, rank.
"""
import argparse, glob, os, sys
import numpy as np
from scipy.stats import rankdata
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_runs import load_gold, load_run, metrics, fmt

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--gold", required=True)
ap.add_argument("--self-pairs", action="store_true")
ap.add_argument("groups", nargs="+", metavar="LABEL=GLOB")
a = ap.parse_args()
keys, y = load_gold(a.gold, a.self_pairs)
print("pairs=%d true=%d random AUPRC=%.3f" % (len(keys), y.sum(), y.mean()))
for g in a.groups:
    label, pat = g.split("=", 1)
    files = sorted(glob.glob(pat))
    if not files:
        print(label, ": no files for", pat); continue
    res = {n: [] for n in ("conf", "vote", "mean", "rank")}
    for f in files:
        r, _ = load_run(f, keys)
        conf, vote = r[1], r[2]
        sc = {"conf": conf, "vote": vote, "mean": (conf + vote) / 2,
              "rank": (rankdata(conf) + rankdata(vote)) / 2}
        for n, s in sc.items():
            res[n].append(metrics(y, s))
    print("\n%s  (%d runs)" % (label, len(files)))
    for n, v in res.items():
        v = np.array(v)
        print("  %-5s %s (sd AUROC %.3f)" % (n, fmt(v.mean(0)), v[:, 0].std()))
