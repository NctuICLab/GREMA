#!/usr/bin/env python3
"""Signs found by GREMA on a real data set with known signs of the regulations.

usage: eval_signs_real.py SIGNS.tsv LABEL=GLOB [LABEL=GLOB ...]
  SIGNS.tsv  TF GENE +1/-1 of the true regulations whose sign is known (input/IRMA/IRMA_signs.tsv)
  GLOB       final_results.txt of the GREMA runs of one setting
Per setting: how often the signed true edges are called a regulation (`+` or `-`), and among those how often the sign is right, per edge and in all;
also the baseline "every edge `+`" (the share of activating edges among the signed edges) and random (0.5).
"""
import sys, glob, collections

signs = {}
for line in open(sys.argv[1]):
    f = line.split()
    signs[(f[0], f[1])] = int(f[2])
n_act = sum(1 for v in signs.values() if v > 0)
print("%d signed true edges, %d activating: the baseline 'every edge +' is right in %.0f%%, random 50%%" % (len(signs), n_act, 100.0 * n_act / len(signs)))
for arg in sys.argv[2:]:
    label, pat = arg.split("=", 1)
    files = sorted(glob.glob(pat))
    called = collections.Counter(); right = collections.Counter()
    for fn in files:
        res = {}
        for line in open(fn):
            f = line.rstrip("\n").split("\t")
            if f[0] != "TF":
                res[(f[0], f[1])] = f[2]
        for e, s in signs.items():
            r = res.get(e, "0")
            if r in ("+", "-"):
                called[e] += 1
                if (r == "+") == (s > 0):
                    right[e] += 1
    tc, tr = sum(called.values()), sum(right.values())
    print("\n%s (%d runs): signed edges called %d / %d (%.0f%%), sign right %d / %d (%.0f%%)" % (
        label, len(files), tc, len(signs) * len(files), 100.0 * tc / (len(signs) * len(files)), tr, tc, 100.0 * tr / max(tc, 1)))
    print("   per edge (called / sign right, of %d runs): " % len(files) + ";  ".join(
        "%s->%s %+d: %d/%d" % (e[0], e[1], signs[e], called[e], right[e]) for e in signs))
