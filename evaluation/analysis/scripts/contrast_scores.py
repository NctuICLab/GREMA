#!/usr/bin/env python3
"""Edge scores from the GA runs of the first iteration: do the runs that contain a candidate regulator fit better than those without it?

usage: contrast_scores.py RUN_OUT_DIR OUT_FILE
  RUN_OUT_DIR  the output directory of a GREMA run made with GREMA_KEEP_ITER=1 (it holds iga_results/<gene>_step1.txt and final_results.txt)
  OUT_FILE     written in the format of final_results.txt (TF GENE REGULATORY CONFIDENCE_LEVEL VOTE_SCORE), so that eval_runs.py can rank it;
               the score is written in both score columns and REGULATORY is '+' for every pair (the score has no sign)

HFODE with a random subset of the candidates in every GA run: a candidate that is not in the subset has the Hill coefficient n exactly 0.
For every target gene and candidate j the 30 runs of the first iteration are split into the runs with n_j != 0 and the runs with n_j == 0,
and the score is the Welch t statistic of log(fitness) (without j) - (with j): positive when the runs that contain j fit better. The subsets are
drawn at random, so the comparison is not biased by the choice of the candidates. A candidate that is in no run or in all runs gets 0.
(VOTE_SCORE counts how often a candidate is in the subset of the better runs; its variance comes from the random membership.)
Result file columns of one run: b, TransMax, n(1..N), k(1..N), degrade, ..., fitness, one trailing field (separated by spaces).
"""
import sys, os, math
import numpy as np

outdir, outfile = sys.argv[1], sys.argv[2]
genes = []
for line in open(os.path.join(outdir, "final_results.txt")):
    f = line.rstrip("\n").split("\t")
    if f[0] != "TF" and f[0] not in genes:
        genes.append(f[0])
N = len(genes)
rows = []
for t, tgt in enumerate(genes):
    a = []
    for line in open(os.path.join(outdir, "iga_results", "%s_step1.txt" % tgt)):
        e = line.split()
        if len(e) < 2 + N + 2:
            continue
        v = [float(x) for x in e]
        if all(math.isfinite(x) for x in v[2:2 + N]) and math.isfinite(v[-2]) and v[-2] > 0:
            a.append(v)
    a = np.array(a)
    logf = np.log(a[:, -2])
    for j, tf in enumerate(genes):
        present = a[:, 2 + j] != 0
        x_with, x_without = logf[present], logf[~present]
        if len(x_with) < 2 or len(x_without) < 2:
            score = 0.0
        else:
            se = math.sqrt(x_with.var(ddof=1) / len(x_with) + x_without.var(ddof=1) / len(x_without))
            score = (x_without.mean() - x_with.mean()) / se if se > 0 else 0.0
        rows.append((tf, tgt, score))
with open(outfile, "w") as f:
    f.write("TF\tGENE\tREGULATORY\tCONFIDENCE_LEVEL\tVOTE_SCORE\n")
    for tf, tgt, s in rows:
        f.write("%s\t%s\t%s\t%.6f\t%.6f\n" % (tf, tgt, "+" if tf != tgt else "0", s, s))
