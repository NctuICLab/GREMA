#!/usr/bin/env python3
"""Control: permute the experiments (replicates) of every gene independently.

usage: permute_replicates.py IN.txt OUT.txt SEED
In the output, gene g of replicate r holds the series that gene g has in replicate sigma_g(r), where sigma_g is a random permutation of the
replicates drawn independently for every gene. The time points stay aligned, so what all experiments share (the common time course, the "phase")
is kept, while the coupling of the genes within one experiment is destroyed. A method whose score stays high on such data does not use
the regulation between the genes of an experiment.
"""
import sys
import numpy as np

infile, outfile, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
rng = np.random.default_rng(seed)
header, rows = [], {}
genes = []
for line in open(infile):
    p = line.rstrip("\n").split("\t")
    if len(p) > 3 and p[0].startswith("rep"):
        rows[(int(p[0][3:]), p[1])] = p[2:]
        if p[1] not in genes:
            genes.append(p[1])
    else:
        header.append(line.rstrip("\n"))
reps = sorted(set(r for r, g in rows))
perm = {g: rng.permutation(reps) for g in genes}
with open(outfile, "w") as f:
    for h in header:
        f.write(h + "\n")
    for r_i, r in enumerate(reps):
        for g in genes:
            src = int(perm[g][r_i])
            f.write("rep%d\t%s\t%s\n" % (r, g, "\t".join(rows[(src, g)])))
