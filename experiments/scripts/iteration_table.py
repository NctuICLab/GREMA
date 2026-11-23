#!/usr/bin/env python3
"""Per-iteration table of GREMA runs, like Table 5 of the GREMA paper (SOS data).

usage: iteration_table.py --gold GOLD.tsv [--self-pairs] FINAL_RESULTS.txt [FINAL_RESULTS.txt ...]

GREMA fixes regulations iteration by iteration; CONFIDENCE_LEVEL = 1 / iteration in which a pair was determined
(a pair determined as "no regulation" counts too). For each confidence level the mean over the files of:
inferred pairs, TP (found +/- and true), TN (found 0 and false), FP, FN, accuracy and the ratio of determined pairs.
Paper (SOS, 64 pairs): conf 1.00 / 0.50 / 0.33 / 0.25: inferred 49 / 9 / 3 / 3, TP 5 / 1 / 0 / 1, TN 33 / 6 / 2 / 1,
FP 10 / 2 / 1 / 1, FN 1 / 0 / 0 / 0, accuracy 0.78 / 0.78 / 0.66 / 0.66.
"""
import argparse
import collections
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--gold", required=True)
ap.add_argument("--self-pairs", action="store_true")
ap.add_argument("files", nargs="+")
a = ap.parse_args()
gold = {}
for line in open(a.gold):
    f = line.split()
    if len(f) >= 3 and (a.self_pairs or f[0] != f[1]):
        gold[(f[0], f[1])] = int(f[2])
levels = collections.defaultdict(lambda: np.zeros((len(a.files), 5)))  # inferred, TP, TN, FP, FN
for i, path in enumerate(a.files):
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "TF" or (f[0], f[1]) not in gold:
            continue
        pred, true = f[2] != "0", gold[(f[0], f[1])] == 1
        row = levels[round(float(f[3]), 2)][i]
        row[0] += 1
        row[1 + (0 if pred and true else 1 if not pred and not true else 2 if pred else 3)] += 1
n = len(gold)
print("%d pairs, %d files (mean per file)" % (n, len(a.files)))
print("%-6s %8s %5s %5s %5s %5s %8s %6s" % ("conf", "inferred", "TP", "TN", "FP", "FN", "accuracy", "ratio"))
for lv in sorted(levels, reverse=True):
    m = levels[lv].mean(0)
    print("%-6.2f %8.1f %5.1f %5.1f %5.1f %5.1f %8.2f %6.2f" % (lv, m[0], m[1], m[2], m[3], m[4], (m[1] + m[2]) / m[0], m[0] / n))
t = np.array([sum(levels[l][i] for l in levels) for i in range(len(a.files))]).mean(0)
print("all    %8.1f %5.1f %5.1f %5.1f %5.1f   SEN %.2f SPE %.2f ACC %.2f" % (
    t[0], t[1], t[2], t[3], t[4], t[1] / max(t[1] + t[4], 1e-9), t[2] / max(t[2] + t[3], 1e-9), (t[1] + t[2]) / t[0]))
