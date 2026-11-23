#!/usr/bin/env python3
"""Convert the SOS DNA repair data (SPEM package, from Ronen et al. 2002) to the GREMA format.

Usage: python3 convert_sos.py <sos csv from export_sos.R> <output dir>
Only the standard library is needed.
"""
import csv
import os
import sys

GENES = ["uvrD", "lexA", "umuD", "recA", "uvrA", "uvrY", "ruvA", "polB"]  # G1..G8
LOW, HIGH = 0.01, 1.0  # the expression is scaled to [LOW, HIGH]; GREMA cannot use 0


def main(src, out):
    with open(src) as f:
        rows = list(csv.reader(f))
    header, body = rows[0][1:], rows[1:]
    names = [r[0] for r in body]
    assert names == GENES, names
    data = [[float(v) for v in r[1:]] for r in body]
    lo = min(min(r) for r in data)
    hi = max(max(r) for r in data)
    n = len(header)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "SOS_timeseries_expression.txt"), "w") as f:
        f.write("repeat_number=1\ntimepoint_number=%d\n" % n)
        for i, row in enumerate(data):
            f.write("rep1\tG%d\t%s\n" % (i + 1, "\t".join("%.6f" % (LOW + (v - lo) / (hi - lo) * (HIGH - LOW)) for v in row)))
    with open(os.path.join(out, "SOS_no_prior_knowledge.txt"), "w") as f:
        f.write("TF\tGENE\tREGULATORY\n")
        for gene in range(1, len(GENES) + 1):
            for tf in range(1, len(GENES) + 1):
                f.write("G%d\tG%d\t%s\n" % (tf, gene, "0" if tf == gene else "?"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
