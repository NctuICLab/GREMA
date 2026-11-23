#!/usr/bin/env python3
"""Convert the original SOS DNA repair data of Ronen et al. (2002) to the GREMA format, 4 replicates.

The data (Uri Alon's lab, "Data on SOS system", sosdata.zip) has 4 experiments, 8 genes, 50 time points
(6 min step); the values are promoter activities (GFP/min/OD). Experiments 1 and 2: UV = 5 J/m^2,
experiments 3 and 4: UV = 20 J/m^2.

As in convert_sos.py (SPEM, which holds Exp3 only), the values are transformed with log2(x + 1) and then
scaled linearly to [0.01, 1] with one map for all experiments and genes (GREMA cannot use 0).

Usage: python3 convert_sos_raw.py <dir with Exp1.txt..Exp4.txt> <output dir>
Only the standard library is needed.
"""
import math
import os
import sys

GENES = ["uvrD", "lexA", "umuD", "recA", "uvrA", "uvrY", "ruvA", "polB"]  # G1..G8 (umuD is called umuDC in the files)
LOW, HIGH = 0.01, 1.0


def read_exp(path):
    with open(path) as f:
        lines = [l.rstrip("\r\n").split("\t") for l in f if l.strip()]
    assert lines[0][0].startswith("time"), lines[0][:2]
    names = [l[0] for l in lines[1:]]
    assert [n.replace("umuDC", "umuD") for n in names] == GENES, names
    data = [[math.log2(float(v) + 1) for v in l[1:]] for l in lines[1:]]
    assert all(len(r) == 50 for r in data)
    return data


def main(src, out):
    exps = [read_exp(os.path.join(src, "Exp%d.txt" % i)) for i in range(1, 5)]
    lo = min(min(min(r) for r in e) for e in exps)
    hi = max(max(max(r) for r in e) for e in exps)
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "SOS_4exp_timeseries_expression.txt"), "w") as f:
        f.write("repeat_number=4\ntimepoint_number=50\n")
        for k, e in enumerate(exps, 1):
            for i, row in enumerate(e):
                f.write("rep%d\tG%d\t%s\n" % (k, i + 1, "\t".join("%.6f" % (LOW + (v - lo) / (hi - lo) * (HIGH - LOW)) for v in row)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
