#!/usr/bin/env python3
"""Does the HFODE model fitted by GREMA on one IRMA condition predict the time course of the other condition? (a fixed-network model, a new external input)

usage: predict_irma.py RUNS_DIR [--min-fitness-rows 1]
RUNS_DIR holds <cond>_<setting>_r<k>/out/G<i>.txt for cond in on, off (GREMA runs on IRMA_switch_<cond>_timeseries_expression.txt, input/IRMA).
Rules (fixed before looking at any prediction):
  model     for every gene the run with the lowest fitness of its 30 runs (columns: b, TransMax, n(1..N), k(1..N), degrade, fitness, trailing field);
            dx_i/dt = TransMax_i * prod_j 1/(1+(k_ij/x_j)^n_ij) + b_i - degrade_i * x_i (a term with n_ij = 0 is left out), explicit Euler with step 1 per time point like GREMA,
            x clipped to [0.0001, 100000]; the same functions as Model.cpp (EMA_HFODE/src)
  predict   free-running simulation of all genes together from the observed first time point of the TARGET series, for the length of the target series;
            also the one-step prediction from the observed values (what GREMA fits)
  baselines persistence (x(0) constant); training mean (the mean level of every gene in the training series, constant);
            no-regulation ODE (dx = c - d x fitted by least squares on the training series, free-running)
  metrics   RMSE per gene averaged over the genes; skill = 1 - MSE / MSE of persistence (mean over the genes; 0 = as good as persistence, < 0 worse);
            Pearson correlation of the predicted and the observed time course (mean over the genes; undefined for a constant prediction)
Prints the means and standard deviations over the GA runs of each setting, for train on -> predict off, train off -> predict on, and the training condition itself.
"""
import sys, glob, os, math
import numpy as np

XMIN, XMAX = 0.0001, 100000.0


def read_series(path):
    rows = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if f[0].startswith("rep") and len(f) > 2:
            rows[f[1]] = np.array([float(x) for x in f[2:]])
    names = sorted(rows, key=lambda g: int(g[1:]))
    return names, np.array([rows[g] for g in names])


def load_model(rundir, N):
    P = []
    for i in range(1, N + 1):
        best = None
        for l in open(f"{rundir}/out/G{i}.txt"):
            f = l.split()
            if len(f) < 2 * N + 5:
                continue
            v = [float(x) for x in f]
            fit = v[2 * N + 3]
            if math.isfinite(fit) and (best is None or fit < best[0]):
                best = (fit, v)
        v = best[1]
        P.append(dict(b=v[0], tm=v[1], n=np.array(v[2:2 + N]), k=np.array(v[2 + N:2 + 2 * N]), deg=v[2 + 2 * N]))
    return P


def step_gene(p, x, xi):
    temp = 1.0
    for j in range(len(x)):
        if p["n"][j] != 0:
            temp *= 1.0 / (1.0 + (p["k"][j] / x[j]) ** p["n"][j])
    act = max(p["tm"] * temp + p["b"], 0.0)
    return min(max(xi + act - p["deg"] * xi, XMIN), XMAX)


def free_run(P, x0, T):
    N = len(P)
    X = np.zeros((N, T)); X[:, 0] = x0
    for t in range(T - 1):
        for i in range(N):
            X[i, t + 1] = step_gene(P[i], X[:, t], X[i, t])
    return X


def one_step(P, Xobs):
    N, T = Xobs.shape
    X = np.zeros((N, T)); X[:, 0] = Xobs[:, 0]
    for t in range(T - 1):
        for i in range(N):
            X[i, t + 1] = step_gene(P[i], Xobs[:, t], Xobs[i, t])
    return X


def no_reg_model(Xtr):
    P = []
    for i in range(Xtr.shape[0]):
        x = Xtr[i, :-1]; d = np.diff(Xtr[i])
        A = np.vstack([np.ones_like(x), -x]).T
        c, dg = np.linalg.lstsq(A, d, rcond=None)[0]
        P.append(dict(b=c, tm=0.0, n=np.zeros(Xtr.shape[0]), k=np.zeros(Xtr.shape[0]), deg=dg))
    return P


def metrics(pred, obs):
    pers = np.repeat(obs[:, :1], obs.shape[1], axis=1)
    mse = ((pred - obs) ** 2).mean(1); mse0 = ((pers - obs) ** 2).mean(1)
    rm = float(np.sqrt(mse).mean())
    skill = float((1 - mse / np.maximum(mse0, 1e-12)).mean())
    cors = []
    for i in range(obs.shape[0]):
        if pred[i].std() > 1e-9 and obs[i].std() > 1e-9:
            cors.append(np.corrcoef(pred[i], obs[i])[0, 1])
    return rm, skill, (float(np.mean(cors)) if cors else float("nan"))


def fmt(vals):
    a = np.array(vals, dtype=float)
    return f"{np.nanmean(a):6.3f} +- {np.nanstd(a):5.3f}"


D = sys.argv[1]
IR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "input", "IRMA")
data = {c: read_series(f"{IR}/IRMA_switch_{c}_timeseries_expression.txt")[1] for c in ("on", "off")}
N = data["on"].shape[0]
print(f"genes {N}; time points on {data['on'].shape[1]}, off {data['off'].shape[1]}")
for train, target in [("on", "off"), ("off", "on"), ("on", "on"), ("off", "off")]:
    Xtr, Xtg = data[train], data[target]
    T = Xtg.shape[1]
    kind = "in-sample" if train == target else "PREDICT"
    print(f"\n== train {train} -> {kind} {target}  (RMSE, skill vs persistence, Pearson r; mean +- sd over the GA runs)")
    base = {
        "persistence": np.repeat(Xtg[:, :1], T, axis=1),
        "training mean": np.repeat(Xtr.mean(1, keepdims=True), T, axis=1),
        "no-regulation ODE": free_run(no_reg_model(Xtr), Xtg[:, 0], T),
    }
    for name, pr in base.items():
        rm, sk, co = metrics(pr, Xtg)
        print(f"  {name:22s} RMSE {rm:6.3f}   skill {sk:6.3f}   r {co:6.3f}")
    for setting in ("def", "div2"):
        dirs = sorted(glob.glob(f"{D}/{train}_{setting}_r*"))
        if not dirs:
            continue
        fr, os_ = [], []
        for d in dirs:
            P = load_model(d, N)
            fr.append(metrics(free_run(P, Xtg[:, 0], T), Xtg)); os_.append(metrics(one_step(P, Xtg), Xtg))
        for lab, res in (("free-running", fr), ("one-step", os_)):
            r = np.array(res)
            print(f"  GREMA {setting:4s} {lab:12s} RMSE {fmt(r[:, 0])}   skill {fmt(r[:, 1])}   r {fmt(r[:, 2])}   ({len(dirs)} runs)")
