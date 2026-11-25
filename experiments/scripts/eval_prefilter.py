#!/usr/bin/env python3
"""Evaluate GREMA runs on a prefiltered knowledge file (lasso_prefilter.py).

usage: eval_prefilter.py --gold GOLD.tsv --prior prior_scores.txt --runs GLOB [--self-pairs]
Rule fixed before the evaluation: the score of a pair is the VOTE_SCORE of the GREMA run (0 for the pairs that were filtered out) plus
1e-4 * the prior score (the tie-breaker among the pairs with the same vote).
Printed (AUROC / AUPRC, mean over the runs): the vote alone, the prior alone and the combined score; and how many true edges the filter kept.
"""
import argparse, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_runs import load_gold, load_run, metrics, fmt

ap = argparse.ArgumentParser()
ap.add_argument("--gold", required=True)
ap.add_argument("--prior", required=True)
ap.add_argument("--runs", required=True)
ap.add_argument("--self-pairs", action="store_true")
a = ap.parse_args()
keys, y = load_gold(a.gold, a.self_pairs)
prior = load_run(a.prior, keys)[0][2]
files = sorted(glob.glob(a.runs))
res = {"vote alone": [], "prior alone": [metrics(y, prior)], "vote + 1e-4 prior": []}
kept = []
for f in files:
    v = load_run(f, keys)[0][2]
    res["vote alone"].append(metrics(y, v))
    res["vote + 1e-4 prior"].append(metrics(y, v + 1e-4 * prior))
    kept.append(int(((v > 0) & (y == 1)).sum()))
print("pairs=%d true=%d random AUPRC=%.3f   runs: %d   true edges with a vote > 0 (mean): %.1f" % (len(keys), y.sum(), y.mean(), len(files), np.mean(kept)))
for k, v in res.items():
    v = np.array(v)
    print("  %-18s %s (sd AUROC %.3f)" % (k, fmt(v.mean(0)), v[:, 0].std()))
