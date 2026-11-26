#!/usr/bin/env python3
"""Forecast of the last time points of both IRMA series from a model fitted on the first points, with and without the input.

usage: forecast.py JOINT_INPUT_16.txt RUNS_DIR [N_TRAIN=12]
JOINT_INPUT_16.txt  the 16-point profile of both series with the input row U1 (made by run_example.sh)
RUNS_DIR            has in12_r<k>/ (GREMA with -input U1 on the first N_TRAIN points) and no12_r<k>/ (without the input, both series as replicates)
The model simulates every series free-running from its observed point N_TRAIN-1 to the end (the input of the series is given), and is compared with the observed points.
Prints RMSE and skill (1 - MSE / MSE of persistence = the last observed value kept) mean +- sd over the GA runs, with two baselines: the mean level of the training points
and a no-regulation ODE (dx = c - d x per gene fitted on the training points of both series).
"""
import sys, os, glob
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "scripts"))
import grema_simulate as gs

data, runs = sys.argv[1], sys.argv[2]
n0 = int(sys.argv[3]) if len(sys.argv) > 3 else 12
genes, inputs, X = gs.load_profile(data, ["U1"])          # series x (variables x time)
G = len(genes); S, _, T = X.shape
obs = X[:, :G, n0 - 1:]                                      # the observed points n0-1 .. T-1
u = X[:, G:, n0 - 1:]
L = obs.shape[2]


def score(sim):
    r, sk = [], []
    for s in range(S):
        a, b = gs.metrics(sim[s][:, 1:], obs[s][:, 1:])        # points n0 .. T-1 (the first one is the start)
        r.append(a)
        p = np.repeat(obs[s][:, :1], L - 1, axis=1)
        sk.append(1 - ((sim[s][:, 1:] - obs[s][:, 1:]) ** 2).mean() / max(((p - obs[s][:, 1:]) ** 2).mean(), 1e-12))
    return float(np.mean(r)), float(np.mean(sk))


def ms(v):
    return "%.3f +- %.3f" % (np.mean(v), np.std(v))


print(f"{S} series, {G} genes, fitted on the first {n0} points, forecast of the points {n0 + 1}..{T}")
for lab, pre, V, use_u in (("with the input U1", "in", G + 1, True), ("without the input (both series as replicates)", "no", G, False)):
    res = []
    for d in sorted(glob.glob(f"{runs}/{pre}{n0}_r*")):
        P = gs.load_model(d, G, V)
        sim = [gs.free_run(P, obs[s][:, 0], u[s] if use_u else np.zeros((0, L)), L) for s in range(S)]
        res.append(score(sim))
    r = np.array(res)
    print(f"  GREMA {lab:46s} RMSE {ms(r[:, 0])}  skill {ms(r[:, 1])}   ({len(res)} runs)")
tr = X[:, :G, :n0]
mean_level = tr.mean(axis=(0, 2))
sim = [np.repeat(mean_level[:, None], L, axis=1) for _ in range(S)]
print(f"  training mean level                                        RMSE {score(sim)[0]:.3f}  skill {score(sim)[1]:.3f}")
c = np.zeros(G); dg = np.zeros(G)
for i in range(G):
    xs = tr[:, i, :-1].ravel(); ds = np.diff(tr[:, i, :], axis=1).ravel()
    c[i], dg[i] = np.linalg.lstsq(np.vstack([np.ones_like(xs), -xs]).T, ds, rcond=None)[0]
sim = []
for s in range(S):
    Y = np.zeros((G, L)); Y[:, 0] = obs[s][:, 0]
    for t in range(L - 1):
        Y[:, t + 1] = np.clip(Y[:, t] + c - dg * Y[:, t], gs.XMIN, gs.XMAX)
    sim.append(Y)
print(f"  no-regulation ODE                                          RMSE {score(sim)[0]:.3f}  skill {score(sim)[1]:.3f}")
