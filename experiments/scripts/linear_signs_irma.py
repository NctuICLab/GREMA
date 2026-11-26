"""Signs of a linear model on the real IRMA data (a check that does not use the Hill function).

usage: linear_signs_irma.py      (run in the repository root)
For every target gene i: the change x_i(t+1) - x_i(t) is regressed on the standardised expression x_j(t) of the other genes (the constant and x_i(t) are projected out first);
the sign of the coefficient (ridge, lambda 0.1, all candidates together) or of the correlation (one candidate at a time) is the predicted sign of the regulation j -> i.
Prints how many of the 7 signed true regulations of input/IRMA/IRMA_signs.tsv get the right sign (switch-on and switch-off separately).
"""
import numpy as np
W = "."
signs = {}
for line in open(W + "/input/IRMA/IRMA_signs.tsv"):
    f = line.split(); signs[(f[0], f[1])] = int(f[2])
def load(fn):
    series, genes = {}, []
    for line in open(fn):
        p = line.rstrip("\n").split("\t")
        if len(p) > 3 and p[0].startswith("rep"):
            series.setdefault(p[0], {})[p[1]] = np.array(p[2:], float)
            if p[1] not in genes: genes.append(p[1])
    return np.array([series["rep1"][g] for g in genes]), genes
for name in ["on", "off"]:
    X, genes = load("%s/input/IRMA/IRMA_switch_%s_timeseries_expression.txt" % (W, name))
    N, T = X.shape
    dX = X[:, 1:] - X[:, :-1]; Xc = X[:, :-1]
    Z = (Xc - Xc.mean(1, keepdims=True)) / Xc.std(1, keepdims=True)
    res = {"ridge": {}, "correlation": {}}
    for i in range(N):
        y = dX[i]
        # remove the degradation term x_i (and the constant) by least squares
        D = np.column_stack([np.ones(T - 1), Xc[i]])
        y_res = y - D @ np.linalg.lstsq(D, y, rcond=None)[0]
        others = [j for j in range(N) if j != i]
        F = np.column_stack([Z[j] - D @ np.linalg.lstsq(D, Z[j], rcond=None)[0] for j in others])
        beta = np.linalg.solve(F.T @ F + 0.1 * np.eye(len(others)), F.T @ y_res)
        for k, j in enumerate(others):
            res["ridge"][("G%d" % (j + 1), "G%d" % (i + 1))] = np.sign(beta[k])
            res["correlation"][("G%d" % (j + 1), "G%d" % (i + 1))] = np.sign(np.corrcoef(Z[j], y)[0, 1])
    for method, r in res.items():
        right = [e for e, s in signs.items() if r[e] == s]
        print("switch-%-3s %-12s sign right for %d of 7 signed true edges (%s)" % (name, method, len(right), ", ".join("%s>%s %s" % (e[0], e[1], "ok" if e in right else "WRONG") for e in signs)))
print("\nbaselines: every edge + : 5 of 7; random 3.5 of 7")
