#!/usr/bin/env python3
"""Convert the IRMA yeast synthetic network data (Cantone et al. 2009) to the GREMA format.

Usage: python3 convert_irma.py <dir with data_on_original.mat and data_off_original.mat> <output dir>

Needs numpy and scipy. See README.md in this directory for the source and the preprocessing.
"""
import os
import sys

import numpy as np
from scipy.io import loadmat

GENES = ["CBF1", "GAL4", "SWI5", "GAL80", "ASH1"]  # row order of the .mat files -> G1..G5
# true network (Cantone et al. 2009), 1-based gene indexes: regulator -> target
EDGES = [(1, 2), (2, 4), (2, 3), (3, 5), (3, 1), (3, 4), (4, 2), (5, 1)]
LOW, HIGH = 0.01, 1.0  # the expression is scaled to [LOW, HIGH]; GREMA cannot use 0


def main(src, out):
    on = loadmat(os.path.join(src, "data_on_original.mat"))["data_on_original"]
    off = loadmat(os.path.join(src, "data_off_original.mat"))["data_off_original"]
    lo = min(on.min(), off.min())
    hi = max(on.max(), off.max())

    def scale(x):  # the same linear scaling for both data sets
        return LOW + (x - lo) / (hi - lo) * (HIGH - LOW)

    os.makedirs(out, exist_ok=True)
    for name, data in (("switch_on", on), ("switch_off", off)):
        data = scale(data)
        with open(os.path.join(out, "IRMA_%s_timeseries_expression.txt" % name), "w") as f:
            f.write("repeat_number=1\ntimepoint_number=%d\n" % data.shape[1])
            for i in range(len(GENES)):
                f.write("rep1\tG%d\t%s\n" % (i + 1, "\t".join("%.6f" % v for v in data[i])))

    with open(os.path.join(out, "IRMA_no_prior_knowledge.txt"), "w") as f:
        f.write("TF\tGENE\tREGULATORY\n")
        for gene in range(1, 6):
            for tf in range(1, 6):
                f.write("G%d\tG%d\t%s\n" % (tf, gene, "0" if tf == gene else "?"))

    true = set(EDGES)
    with open(os.path.join(out, "IRMA_GoldStandard.tsv"), "w") as f:  # directed
        for a in range(1, 6):
            for b in range(1, 6):
                if a != b:
                    f.write("G%d\tG%d\t%d\n" % (a, b, int((a, b) in true)))
    undirected = {tuple(sorted(e)) for e in EDGES}
    with open(os.path.join(out, "IRMA_GoldStandardUndirected.tsv"), "w") as f:
        for a in range(1, 6):
            for b in range(a + 1, 6):
                f.write("G%d\tG%d\t%d\n" % (a, b, int((a, b) in undirected)))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
