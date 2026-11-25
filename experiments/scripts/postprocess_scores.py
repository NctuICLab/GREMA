#!/usr/bin/env python3
"""Post-processing of the edge scores of GREMA (VOTE_SCORE) and combination with dynGENIE3; rules fixed before the evaluation.

usage: postprocess_scores.py --gold GOLD.tsv --grema GLOB --dyn-rf GLOB --dyn-et GLOB [--self-pairs]
  --grema    final_results.txt of the GREMA runs (sorted by name; run i is paired with dynGENIE3 seed i)
  --dyn-rf / --dyn-et   dynGENIE3 results written by dyngenie3_irma.py (final_results.txt of seed0, seed1, ...)

Scores of every run, AUROC / AUPRC averaged over the runs (the gold standard only defines the evaluated pairs and the metrics):
  vote       VOTE_SCORE of GREMA
  A1 z       vote standardised within every target gene (over its candidate regulators; 0 if the sd is 0)
  A2 sum     vote divided by the sum of the votes of the target gene (the GENIE3 convention)
  dynRF/ET   dynGENIE3 score (Random Forest / Extra-Trees)
  B  vote+dynRF, vote+dynET   mean of the two average ranks (over all evaluated pairs)
"""
import argparse, glob, os, sys
import collections
import numpy as np
from scipy.stats import rankdata
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_runs import load_gold, load_run, metrics, fmt

ap = argparse.ArgumentParser()
ap.add_argument("--gold", required=True)
ap.add_argument("--grema", required=True)
ap.add_argument("--dyn-rf", required=True)
ap.add_argument("--dyn-et", required=True)
ap.add_argument("--self-pairs", action="store_true")
a = ap.parse_args()
keys, y = load_gold(a.gold, a.self_pairs)
targets = collections.defaultdict(list)
for i, k in enumerate(keys):
    targets[k[1]].append(i)


def files(g):
    return sorted(glob.glob(g))


def vote(f):
    return load_run(f, keys)[0][2]


def per_target(v, how):
    out = np.zeros_like(v)
    for t, idx in targets.items():
        x = v[idx]
        if how == "z":
            sd = x.std()
            out[idx] = (x - x.mean()) / sd if sd > 0 else 0.0
        else:
            s = x.sum()
            out[idx] = x / s if s > 0 else 0.0
    return out


def mean_rank(u, v):
    return (rankdata(u) + rankdata(v)) / 2


grema, rf, et = files(a.grema), files(a.dyn_rf), files(a.dyn_et)
n = min(len(grema), len(rf), len(et))
res = collections.defaultdict(list)
for i in range(n):
    v = vote(grema[i])
    d_rf, d_et = vote(rf[i]), vote(et[i])
    sc = {"vote": v, "A1 z": per_target(v, "z"), "A2 sum": per_target(v, "sum"),
          "dynRF": d_rf, "dynET": d_et,
          "B vote+dynRF": mean_rank(v, d_rf), "B vote+dynET": mean_rank(v, d_et)}
    for k, s in sc.items():
        res[k].append(metrics(y, s))
print("pairs=%d true=%d random AUPRC=%.3f   runs paired: %d" % (len(keys), y.sum(), y.mean(), n))
for k in ["vote", "A1 z", "A2 sum", "dynRF", "dynET", "B vote+dynRF", "B vote+dynET"]:
    v = np.array(res[k])
    print("  %-14s %s (sd AUROC %.3f)" % (k, fmt(v.mean(0)), v[:, 0].std()))
