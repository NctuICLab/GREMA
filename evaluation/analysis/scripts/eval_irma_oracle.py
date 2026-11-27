#!/usr/bin/env python3
"""Signs found by GREMA on IRMA when the true structure is told (D2 of the sign diagnosis; irma_oracle_knowledge.py + GREMA runs).

usage: eval_irma_oracle.py RUNS_DIR [IRMA_DIR]
RUNS_DIR has <cond>_<setting>_r<k>/out/final_results.txt with cond = on, off and setting = unknown_signs (the 8 true edges are '?', all other pairs 0) or
loo_<TF>_<GENE> (that signed edge is '?', the other signed edges carry their true sign).
For every signed edge and condition: in how many of the runs the edge was called with the right sign, with the wrong sign, or not called (0); for loo the edge that was left open,
for unknown_signs every edge. Also the mean vote score. The reference of the earlier report (structure unknown, 6 settings of 10 runs): 25-50% of the called signs right.
"""
import sys, os, glob, collections

D = sys.argv[1]
IR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "input", "IRMA")
signs = {}
for l in open(f"{IR}/IRMA_signs.tsv"):
    f = l.split(); signs[(f[0], f[1])] = int(f[2])


def res(path):
    r = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if f[0] != "TF":
            r[(f[0], f[1])] = (f[2], float(f[4]))
    return r


def tally(files, edge):
    c = collections.Counter(); votes = []
    for fn in files:
        sg, vt = res(fn)[edge]
        votes.append(vt)
        c["not called" if sg == "0" else ("right" if (sg == "+") == (signs[edge] > 0) else "wrong")] += 1
    return c, sum(votes) / max(len(votes), 1)


def line(c, v):
    n = sum(c.values())
    called = c["right"] + c["wrong"]
    return "right %2d  wrong %2d  not called %2d   (right among called: %s)  vote %.2f" % (c["right"], c["wrong"], c["not called"], "%d%%" % round(100.0 * c["right"] / called) if called else "n/a", v)


print("edge     known  | structure known, this sign open (the others told)           | structure known, no sign told")
tot = {"loo": collections.Counter(), "unk": collections.Counter()}
for cond in ("on", "off"):
    print(f"\n-- switch-{cond}")
    for edge, s in signs.items():
        lo = sorted(glob.glob(f"{D}/{cond}_loo_{edge[0]}_{edge[1]}_r*/out/final_results.txt"))
        un = sorted(glob.glob(f"{D}/{cond}_unknown_signs_r*/out/final_results.txt"))
        if not lo or not un:
            continue
        cl, vl = tally(lo, edge); cu, vu = tally(un, edge)
        tot["loo"].update(cl); tot["unk"].update(cu)
        print("%-8s %-5s  | %-72s | %s" % (f"{edge[0]}->{edge[1]}", "+" if s > 0 else "-", line(cl, vl), line(cu, vu)))
for k, lab in (("loo", "one sign open, the structure and the other signs told"), ("unk", "structure told, no sign told")):
    c = tot[k]; called = c["right"] + c["wrong"]
    print(f"\nall edges, both conditions, {lab}: right {c['right']}, wrong {c['wrong']}, not called {c['not called']}; right among the called: {100.0 * c['right'] / max(called, 1):.0f}%")
