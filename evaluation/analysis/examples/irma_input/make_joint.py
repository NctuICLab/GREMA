#!/usr/bin/env python3
"""IRMA switch-on and switch-off as the two series of one profile, cut to the first N time points (default 16, the length of switch-on).

usage: make_joint.py IRMA_DIR OUT.txt [N]
"""
import sys

d, out = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 16
series = []
for c in ("on", "off"):
    rows = {}
    for l in open(f"{d}/IRMA_switch_{c}_timeseries_expression.txt"):
        f = l.rstrip("\n").split("\t")
        if len(f) > 3 and f[0] == "rep1":
            rows[f[1]] = f[2:2 + n]
    series.append(rows)
with open(out, "w") as f:
    f.write(f"repeat_number=2\ntimepoint_number={n}\n")
    for r, rows in enumerate(series, 1):
        for g in sorted(rows, key=lambda x: int(x[1:])):
            f.write(f"rep{r}\t{g}\t" + "\t".join(rows[g]) + "\n")
