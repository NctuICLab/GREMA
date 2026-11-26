#!/usr/bin/env python3
"""Simulate the HFODE model that GREMA fitted (with or without external inputs) for a new condition.

usage: grema_simulate.py --runs RUN_DIR [RUN_DIR ...] --data DATA.txt [--inputs U1,U2] [--out OUT.tsv] [--score]
  RUN_DIR   a GREMA output directory (it has out/G*.txt or G*.txt; the genes are G1..Gn in the order of run.log)
  DATA.txt  GREMA profile of the new condition (repeat_number=, timepoint_number=, rows 'rep<k> <name> v1 v2 ...'); the first value of every gene
            is the initial value, the rows named in --inputs are the input u(t) of the new condition (same scale as the genes, 0.01 to 1)
Model (as in EMA_HFODE/src/Model.cpp, equal steps of 1):
  x_i(t+1) = clip(x_i + TransMax_i * prod_j 1/(1+(k_ij/x_j)^n_ij) + b_i - degrade_i * x_i, 0.0001, 100000), j over the genes (simulated) and the inputs (given);
  a factor with n_ij = 0 is left out; the parameters of a gene are the run with the lowest fitness of its runs (columns: b, TransMax, n(1..V), k(1..V), degrade, fitness, one more).
  V = genes + inputs; the order of the variables is: the genes in natural order of their names, then the inputs.
--score  the observed genes of DATA.txt are compared with the simulation: RMSE and skill (1 - MSE / MSE of persistence) per run, mean and sd over the runs.
--out    the simulated time courses of every run, TSV: run, gene, v1 ... vT.
Importable: load_profile, load_model, free_run, metrics.
"""
import sys, os, re, math, argparse
import numpy as np

XMIN, XMAX = 0.0001, 100000.0


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def load_profile(path, inputs=()):
    """-> genes (names), inputs (names), array [series][variable][time] with the genes first (natural order) then the inputs (natural order)."""
    rows, nrep, npt = {}, 0, 0
    for l in open(path):
        l = l.rstrip("\n").rstrip("\r")
        if l.startswith("repeat_number="):
            nrep = int(l.split("=")[1])
        elif l.startswith("timepoint_number="):
            npt = int(l.split("=")[1])
        elif re.match(r"rep\d+\t", l):
            f = l.split("\t")
            rows.setdefault(int(f[0][3:]), {})[f[1]] = [float(x) for x in f[2:]]
    names = sorted(rows[1], key=natural_key)
    inputs = [i for i in names if i in set(inputs)]
    genes = [g for g in names if g not in set(inputs)]
    order = genes + inputs
    X = np.array([[rows[r][n] for n in order] for r in sorted(rows)])
    return genes, inputs, X


def load_model(rundir, n_genes, n_vars):
    d = rundir + "/out" if os.path.isdir(rundir + "/out") else rundir
    P = []
    for i in range(1, n_genes + 1):
        best = None
        for l in open(f"{d}/G{i}.txt"):
            f = l.split()
            if len(f) < 2 * n_vars + 5:
                continue
            v = [float(x) for x in f]
            fit = v[2 * n_vars + 3]
            if math.isfinite(fit) and (best is None or fit < best[0]):
                best = (fit, v)
        v = best[1]
        P.append(dict(b=v[0], tm=v[1], n=np.array(v[2:2 + n_vars]), k=np.array(v[2 + n_vars:2 + 2 * n_vars]), deg=v[2 + 2 * n_vars]))
    return P


def step_gene(p, x, xi):
    temp = 1.0
    for j in range(len(x)):
        if p["n"][j] != 0:
            temp *= 1.0 / (1.0 + (p["k"][j] / x[j]) ** p["n"][j])
    act = max(p["tm"] * temp + p["b"], 0.0)
    return min(max(xi + act - p["deg"] * xi, XMIN), XMAX)


def free_run(P, x0, U, T=None):
    """x0: initial values of the genes; U: inputs [input][time] (empty for a model without inputs, then T is needed). Returns [gene][time]."""
    U = np.asarray(U, float)
    if U.size == 0:
        U = np.zeros((0, T))
    T = U.shape[1]
    G = len(P)
    X = np.zeros((G, T)); X[:, 0] = x0
    for t in range(T - 1):
        v = np.concatenate([X[:, t], U[:, t]])
        for i in range(G):
            X[i, t + 1] = step_gene(P[i], v, X[i, t])
    return X


def metrics(pred, obs):
    pers = np.repeat(obs[:, :1], obs.shape[1], axis=1)
    mse = ((pred - obs) ** 2).mean(1); mse0 = ((pers - obs) ** 2).mean(1)
    return float(np.sqrt(mse).mean()), float((1 - mse / np.maximum(mse0, 1e-12)).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--inputs", default="")
    ap.add_argument("--out")
    ap.add_argument("--score", action="store_true")
    a = ap.parse_args()
    genes, inputs, X = load_profile(a.data, [s for s in a.inputs.split(",") if s])
    G, V = len(genes), len(genes) + len(inputs)
    series = X[0]
    T = series.shape[1]
    out = open(a.out, "w") if a.out else None
    res = []
    for rd in a.runs:
        P = load_model(rd, G, V)
        S = free_run(P, series[:G, 0], series[G:], T)
        if out:
            for g, name in enumerate(genes):
                out.write(rd + "\t" + name + "\t" + "\t".join("%.6f" % v for v in S[g]) + "\n")
        if a.score:
            res.append(metrics(S, series[:G]))
    if a.score:
        r = np.array(res)
        print(f"{len(res)} runs, {G} genes, {len(inputs)} inputs, {T} time points: RMSE {r[:, 0].mean():.3f} +- {r[:, 0].std():.3f}, skill {r[:, 1].mean():.3f} +- {r[:, 1].std():.3f}")


if __name__ == "__main__":
    main()
