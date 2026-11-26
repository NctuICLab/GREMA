#!/usr/bin/env python3
"""SOS (4 experiments, UV 5 J/m^2 in experiments 1 and 2, 20 J/m^2 in 3 and 4) with the UV dose as an external input: the data of the held-out dose test.

usage: prep_sos_input.py SOS_4exp_timeseries_expression.txt OUT_DIR
For every fold (loeo1..loeo4: leave experiment k out, fit on the other three, predict k = the main design, the input has both doses in the training data;
to20: fit on the two experiments of 5, predict the two of 20; to5: the other way: the input is then the same in both training series, check_data.py calls it not identifiable)
and every input shape it writes, in OUT_DIR/<fold>/:
  train_<shape>.txt, test_<shape>.txt   the experiments as rep1, rep2, ..., genes + the input row U1
  train_noinput.txt, test_noinput.txt   without the input;  know_input.txt, know_noinput.txt (no prior knowledge: the diagonal 0, every other pair '?'; U1 -> every gene '?')
Shapes (the input of an experiment; raw dose 5 or 20, scaled with u = 0.01 + 0.99 * raw / 20, so 5 -> 0.2575 and 20 -> 1; the exposure is at the first time point):
  decay10 (the main one: the dose falls to 0.01 with a half-time of 10 time points = 60 min), decay5 (30 min), decay20 (120 min), step (the dose, constant),
  rev10 (decay10 in the reverse order of time: a control, the stimulus grows instead of decaying).
"""
import sys, os

src, out = sys.argv[1], sys.argv[2]
hdr, rows = [], {}
for l in open(src):
    l = l.rstrip("\n")
    if l.startswith("repeat_number") or l.startswith("timepoint_number"):
        hdr.append(l)
    elif l.startswith("rep"):
        f = l.split("\t")
        rows.setdefault(int(f[0][3:]), {})[f[1]] = f[2:]
genes = sorted(rows[1], key=lambda g: int(g[1:]))
T = len(rows[1][genes[0]])
dose = {1: 5, 2: 5, 3: 20, 4: 20}


def u_of(raw, shape):
    base = 0.0
    if shape == "step":
        v = [raw] * T
    else:
        h = {"decay10": 10, "decay5": 5, "decay20": 20, "rev10": 10}[shape]
        v = [base + (raw - base) * 2.0 ** (-t / h) for t in range(T)]
        if shape == "rev10":
            v = v[::-1]
    return [0.01 + 0.99 * x / 20.0 for x in v]


def write(path, exps, shape):
    with open(path, "w") as f:
        f.write(f"repeat_number={len(exps)}\ntimepoint_number={T}\n")
        for r, e in enumerate(exps, 1):
            for g in genes:
                f.write(f"rep{r}\t{g}\t" + "\t".join(rows[e][g]) + "\n")
            if shape:
                f.write(f"rep{r}\tU1\t" + "\t".join("%.6f" % x for x in u_of(dose[e], shape)) + "\n")


# to20 / to5: fit on the two experiments of one dose, predict the two of the other (the input is then the same in both training series: not identifiable, check_data says so);
# loeo<k>: leave experiment k out (fit on the other three, which have both doses), predict experiment k (the identifiable design, the main one)
folds = {"to20": ([1, 2], [3, 4]), "to5": ([3, 4], [1, 2])}
for k in (1, 2, 3, 4):
    folds["loeo%d" % k] = ([e for e in (1, 2, 3, 4) if e != k], [k])
for d, (tr, te) in folds.items():
    os.makedirs(f"{out}/{d}", exist_ok=True)
    for shape in ("decay10", "decay5", "decay20", "step", "rev10", None):
        name = shape or "noinput"
        write(f"{out}/{d}/train_{name}.txt", tr, shape)
        write(f"{out}/{d}/test_{name}.txt", te, shape)
    for kname, withinput in (("know_input", True), ("know_noinput", False)):
        with open(f"{out}/{d}/{kname}.txt", "w") as f:
            f.write("TF\tGENE\tREGULATORY\n")
            for g in genes:
                for t in genes:
                    f.write(f"{t}\t{g}\t{'0' if t == g else '?'}\n")
                if withinput:
                    f.write(f"U1\t{g}\t?\n")
print("written", out)
