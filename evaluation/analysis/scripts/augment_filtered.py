#!/usr/bin/env python3
"""Add a filtered (accumulated) copy of every gene as an extra variable that is the only possible regulator (no change of the solver).

usage: augment_filtered.py IN_DIR OUT_DIR D
  IN_DIR  has expression.txt and knowledge.txt (GREMA format, genes G1..Gn)
  D       decay rate of the exponential filter  F(t+1) = F(t) + D (x(t) - F(t)),  F(0) = x(0)   (inside every series, time step 1)
Writes OUT_DIR/expression.txt with the genes G1..Gn (targets) and F1..Fn (filtered copies), and OUT_DIR/knowledge.txt in which
  only F_j may regulate G_i (j != i: '?'), the measured genes G_j never regulate ('0'), the F variables have no regulator ('0').
The regulation of a pair F_j -> G_i is the regulation G_j -> G_i of the original genes (`unfilter_results.py`).
"""
import sys, os

indir, outdir, D = sys.argv[1], sys.argv[2], float(sys.argv[3])
os.makedirs(outdir, exist_ok=True)
header, rows, genes = [], [], []
for line in open(os.path.join(indir, "expression.txt")):
    p = line.rstrip("\n").split("\t")
    if len(p) > 3 and p[0].startswith("rep"):
        rows.append((p[0], p[1], [float(x) for x in p[2:]]))
        if p[1] not in genes:
            genes.append(p[1])
    else:
        header.append(line.rstrip("\n"))
with open(os.path.join(outdir, "expression.txt"), "w") as f:
    for h in header:
        f.write(h + "\n")
    for rep, g, v in rows:
        f.write("%s\t%s\t%s\n" % (rep, g, "\t".join("%.6f" % x for x in v)))
    for rep, g, v in rows:
        fv = [v[0]]
        for t in range(len(v) - 1):
            fv.append(fv[-1] + D * (v[t] - fv[-1]))
        f.write("%s\tF%s\t%s\n" % (rep, g[1:], "\t".join("%.6f" % x for x in fv)))
allg = genes + ["F" + g[1:] for g in genes]
with open(os.path.join(outdir, "knowledge.txt"), "w") as f:
    f.write("TF\tGENE\tREGULATORY\n")
    for a in allg:
        for b in allg:
            ok = a.startswith("F") and b.startswith("G") and a[1:] != b[1:]
            f.write("%s\t%s\t%s\n" % (a, b, "?" if ok else "0"))
