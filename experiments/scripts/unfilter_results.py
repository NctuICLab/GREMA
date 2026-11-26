#!/usr/bin/env python3
"""Map the results of a run on the augmented data (augment_filtered.py) back to the original genes: the pair F_j -> G_i becomes G_j -> G_i.

usage: unfilter_results.py FINAL_RESULTS_IN FINAL_RESULTS_OUT
All other pairs (they were fixed '0') are dropped; the diagonal is written as '0'.
"""
import sys

fin, fout = sys.argv[1], sys.argv[2]
res = {}
genes = []
for line in open(fin):
    f = line.rstrip("\n").split("\t")
    if f[0] == "TF":
        continue
    if f[0].startswith("F") and f[1].startswith("G"):
        res[("G" + f[0][1:], f[1])] = f
        if f[1] not in genes:
            genes.append(f[1])
with open(fout, "w") as o:
    o.write("TF\tGENE\tREGULATORY\tCONFIDENCE_LEVEL\tVOTE_SCORE\n")
    for a in genes:
        for b in genes:
            if a == b:
                o.write("%s\t%s\t0\t1.000\t0.000\n" % (a, b))
            else:
                f = res[(a, b)]
                o.write("%s\t%s\t%s\t%s\t%s\n" % (a, b, f[2], f[3], f[4]))
