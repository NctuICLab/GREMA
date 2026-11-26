#!/usr/bin/env python3
"""Run dynGENIE3 (Huynh-Thu and Geurts 2018) on a GREMA time-series file and write a GREMA-like final_results.txt.

usage: python3 -I dyngenie3_runner.py DYNGENIE3_DIR IN.txt OUT_DIR --seeds 10 [--tree RF|ET] [--ntrees 1000] [--know KNOW.txt] [--nthreads 1]
  --know  GREMA knowledge file: only the genes with at least one '?' as regulator are candidate regulators (default: every gene)
  IN.txt  GREMA format (repeat_number=, timepoint_number=, then 'rep<k> <gene> v1 v2 ...')
  Every replicate is one time series (time points 0, 1, 2, ...); the genes are G1..Gn as in the file.
  One file OUT_DIR/seed<s>/final_results.txt per seed with the columns TF GENE REGULATORY CONFIDENCE_LEVEL VOTE_SCORE:
  the dynGENIE3 score of the edge TF -> GENE divided by the largest score, written in both score columns; REGULATORY is '+' for every pair
  (dynGENIE3 gives no sign), so that eval_runs.py can rank the edges (binary is not meaningful).
"""
import sys, os
import numpy as np

args = sys.argv[1:]
src, infile, outdir = args[:3]
opt = {"--seeds": "10", "--tree": "RF", "--ntrees": "1000", "--know": "", "--nthreads": "1"}
for i in range(3, len(args), 2):
    opt[args[i]] = args[i + 1]
sys.path.insert(0, src)
import dynGENIE3 as dg

genes, series = [], {}
for line in open(infile):
    p = line.rstrip("\n").split("\t")
    if len(p) > 3 and p[0].startswith("rep"):
        series.setdefault(p[0], {})[p[1]] = [float(x) for x in p[2:]]
        if p[1] not in genes:
            genes.append(p[1])
ts = [np.array([series[r][g] for g in genes]).T for r in sorted(series, key=lambda s: int(s[3:]))]
times = [np.arange(t.shape[0], dtype=float) for t in ts]
regs = "all"
if opt["--know"]:
    cand = []
    for line in open(opt["--know"]):
        q = line.rstrip("\n").split("\t")
        if len(q) >= 3 and q[2].strip() == "?" and q[0] not in cand:
            cand.append(q[0])
    regs = [g for g in genes if g in cand]
print("genes", len(genes), "series", len(ts), "points", ts[0].shape[0], "candidate regulators", regs if regs == "all" else len(regs))
for s in range(int(opt["--seeds"])):
    np.random.seed(s)
    VIM, alphas, _, _, _ = dg.dynGENIE3(ts, times, gene_names=genes, regulators=regs, tree_method=opt["--tree"], ntrees=int(opt["--ntrees"]), nthreads=int(opt["--nthreads"]))
    m = VIM.max() if VIM.max() > 0 else 1.0
    d = os.path.join(outdir, "seed%d" % s)
    os.makedirs(d, exist_ok=True)
    with open(d + "/final_results.txt", "w") as f:
        f.write("TF\tGENE\tREGULATORY\tCONFIDENCE_LEVEL\tVOTE_SCORE\n")
        for i, tf in enumerate(genes):
            for j, g in enumerate(genes):
                sc = VIM[i, j] / m
                f.write("%s\t%s\t%s\t%.6f\t%.6f\n" % (tf, g, "+" if i != j else "0", sc, sc))
    print("seed", s, "done")
