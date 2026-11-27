#!/usr/bin/env python3
"""Do the real IRMA data carry the expected sign of each known regulation when the structure is known? (D1, no GREMA)

usage: irma_sign_oracle_linear.py [IRMA_DIR] [--grema RUNS_DIR]   (IRMA_DIR default: input/IRMA of the repository)
--grema  also compare the sign of the regression at lag 0 with the sign that GREMA finds when the structure is told but no sign (RUNS_DIR of irma_oracle_knowledge.py + GREMA, runs <cond>_unknown_signs_r<k>)
For every gene i with known regulators R_i (the gold standard of IRMA; the signs are those of IRMA_signs.tsv, the edge G2 -> G4 has no known sign and is only a covariate):
OLS of the change dx_i(t) = x_i(t+1) - x_i(t) on [1, x_i(t), x_j(t - L) for j in R_i] with the lag L = 0, 1, 2, on switch-on, on switch-off and on both series together
(with a constant for the condition). The sign of the coefficient of a regulator, with its t value, against the known sign of the regulation.
This is the best a linear model of the one-step change can do with the TRUE structure; with 14-20 points per series the t values are small, so read the signs as a question
"does the data point the right way at all", not as a test. Only numpy is used.
"""
import sys, os
import numpy as np

args = [a for a in sys.argv[1:]]
GR = None
if "--grema" in args:
    i = args.index("--grema"); GR = args[i + 1]; del args[i:i + 2]
IR = args[0] if args else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "input", "IRMA")


def read(path):
    rows = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) > 3 and f[0] == "rep1":
            rows[f[1]] = np.array([float(x) for x in f[2:]])
    return np.array([rows["G%d" % i] for i in range(1, 6)])


data = {c: read(f"{IR}/IRMA_switch_{c}_timeseries_expression.txt") for c in ("on", "off")}
signs = {}
for l in open(f"{IR}/IRMA_signs.tsv"):
    f = l.split()
    signs[(f[0], f[1])] = int(f[2])
gold = [tuple(l.split("\t")[:2]) for l in open(f"{IR}/IRMA_GoldStandard.tsv") if l.rstrip("\n").endswith("\t1")]


def design(conds, gene, L):
    """rows and columns for gene (1-based): returns y, X, names of the regulators"""
    regs = [t for (t, g) in gold if g == f"G{gene}"]
    ys, Xs = [], []
    for ci, c in enumerate(conds):
        x = data[c]
        T = x.shape[1]
        for t in range(L, T - 1):
            row = [1.0, x[gene - 1, t]] + [x[int(r[1:]) - 1, t - L] for r in regs]
            if len(conds) > 1:
                row.append(float(ci))
            Xs.append(row); ys.append(x[gene - 1, t + 1] - x[gene - 1, t])
    return np.array(ys), np.array(Xs), regs


def fit(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ b
    dof = max(len(y) - X.shape[1], 1)
    cov = (res @ res / dof) * np.linalg.pinv(X.T @ X)
    return b, b / np.sqrt(np.maximum(np.diag(cov), 1e-30))


print("known edges with a sign:", ", ".join(f"{t}->{g} {'+' if s > 0 else '-'}" for (t, g), s in signs.items()))
for label, conds in (("switch-on", ["on"]), ("switch-off", ["off"]), ("both series, with a constant for the condition", ["on", "off"])):
    print(f"\n== {label}")
    print("%-8s %-5s %s" % ("edge", "known", "   ".join("lag %d: coef (t)" % L for L in (0, 1, 2))))
    right = {L: 0 for L in (0, 1, 2)}
    for (t, g), s in signs.items():
        cells = []
        for L in (0, 1, 2):
            y, X, regs = design(conds, int(g[1:]), L)
            b, tv = fit(y, X)
            k = 2 + regs.index(t)
            ok = np.sign(b[k]) == s
            right[L] += int(ok)
            cells.append("%+.2f (%+.1f) %s" % (b[k], tv[k], "ok" if ok else "WRONG"))
        print("%-8s %-5s %s" % (f"{t}->{g}", "+" if s > 0 else "-", "   ".join(cells)))
    print("right signs of the 7 edges: " + ", ".join(f"lag {L}: {right[L]}" for L in (0, 1, 2)))

if GR:
    import glob, collections
    print("\n== the sign of the regression at lag 0 against the sign GREMA finds when the structure is told and no sign (majority of the runs)")
    agree = tot = 0
    for c in ("on", "off"):
        for (t, g), s in signs.items():
            y, X, regs = design([c], int(g[1:]), 0)
            b, tv = fit(y, X)
            lin = int(np.sign(b[2 + regs.index(t)]))
            votes = collections.Counter()
            for fn in glob.glob(f"{GR}/{c}_unknown_signs_r*/out/final_results.txt"):
                for l in open(fn):
                    f = l.rstrip("\n").split("\t")
                    if (f[0], f[1]) == (t, g) and f[2] in "+-" and f[2] != "0":
                        votes[1 if f[2] == "+" else -1] += 1
            gr = votes.most_common(1)[0][0] if votes else 0
            agree += int(gr == lin); tot += 1
            print("  switch-%-3s %s->%s  known %s  regression %s  GREMA %s (%s of %d runs)  %s" % (c, t, g, "+" if s > 0 else "-", "+" if lin > 0 else "-", "+" if gr > 0 else "-", votes[gr], sum(votes.values()), "same" if gr == lin else "DIFFERENT"))
    print("the regression and GREMA give the same sign for %d of %d edges and conditions" % (agree, tot))
