#!/usr/bin/env python3
"""Can augmented data take the place of real series? Evaluation of sim_aug_prep.py + GREMA runs (INPUT: SIM_DIR with seed<k>/ and runs/seed<k>_<cond>/out/).

usage: eval_augment.py SIM_DIR [--seeds 5]
Conditions per network (seed): s1, s2, s3, s4 = 1-4 real series; p3, p9 = s1 with -perturb 3 / 9 -noise 0.03 (the built-in augmentation: 4 and 10 series, 3 and 9 of them noisy copies);
gp = 10 GP samples of s1 (gp_v2.py), GREMA run on each one alone, an edge scored by the fraction of the 10 runs that call it (also the mean vote score).
Measures: AUROC of the vote score (gp: the fraction of runs) against the 14-18 true edges of the 90 pairs; the free-running skill (1 - MSE / MSE of persistence) on the held-out series 5
(gp: the mean over the 10 models) and whether the prediction is stable (skill > 0).
Rules fixed before the runs: an augmentation takes the place of real series if it closes at least half of the gap in AUROC between s1 and s4 (paired over the seeds) AND its held-out
prediction is stable (skill > 0) in at least as many networks as s4 minus one.
"""
import sys, os, glob, math
import numpy as np
from sklearn.metrics import roc_auc_score

D = sys.argv[1]
nseed = int(sys.argv[sys.argv.index("--seeds") + 1]) if "--seeds" in sys.argv else 5
XMIN, XMAX = 0.0001, 100000.0
CONDS = ["s1", "p3", "p9", "gp", "s2", "s3", "s4"]


def read_final(path):
    r = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if f[0] != "TF":
            r[(f[0], f[1])] = (f[2], float(f[4]))
    return r


def read_series(path):
    rows = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) > 3 and f[0].startswith("rep"):
            rows.setdefault(f[0], {})[f[1]] = [float(x) for x in f[2:]]
    return rows


def load_model(rundir, N):
    P = []
    for i in range(1, N + 1):
        best = None
        for l in open(f"{rundir}/out/G{i}.txt"):
            f = l.split()
            if len(f) < 2 * N + 5:
                continue
            v = [float(x) for x in f]
            if math.isfinite(v[2 * N + 3]) and (best is None or v[2 * N + 3] < best[0]):
                best = (v[2 * N + 3], v)
        v = best[1]
        P.append(dict(b=v[0], tm=v[1], n=np.array(v[2:2 + N]), k=np.array(v[2 + N:2 + 2 * N]), deg=v[2 + 2 * N]))
    return P


def free_run(P, x0, T):
    N = len(P)
    X = np.zeros((N, T)); X[:, 0] = x0
    for t in range(T - 1):
        for i in range(N):
            p = P[i]; temp = 1.0
            for j in range(N):
                if p["n"][j] != 0:
                    temp *= 1.0 / (1.0 + (p["k"][j] / X[j, t]) ** p["n"][j])
            X[i, t + 1] = min(max(X[i, t] + max(p["tm"] * temp + p["b"], 0.0) - p["deg"] * X[i, t], XMIN), XMAX)
    return X


def skill(pred, obs):
    pers = np.repeat(obs[:, :1], obs.shape[1], axis=1)
    return float((1 - ((pred - obs) ** 2).mean(1) / np.maximum(((pers - obs) ** 2).mean(1), 1e-12)).mean())


res = {c: {"auc": {}, "auc_vote": {}, "sk": {}} for c in CONDS}
for s in range(1, nseed + 1):
    sd = f"{D}/seed{s}"
    gold = {(l.split("\t")[0], l.split("\t")[1]): int(l.split("\t")[2]) for l in open(f"{sd}/gold.tsv")}
    keys = [k for k in gold if k[0] != k[1]]; y = [gold[k] for k in keys]   # the 90 non-self pairs
    ho = read_series(f"{sd}/heldout.txt")["rep1"]
    genes = sorted(ho, key=lambda g: int(g[1:])); N = len(genes)
    H = np.array([ho[g] for g in genes])
    for c in CONDS:
        if c != "gp":
            rd = f"{D}/runs/seed{s}_{c}"
            if not os.path.exists(f"{rd}/out/final_results.txt"):
                continue
            fin = read_final(f"{rd}/out/final_results.txt")
            res[c]["auc"][s] = roc_auc_score(y, [fin[k][1] for k in keys])
            res[c]["sk"][s] = skill(free_run(load_model(rd, N), H[:, 0], H.shape[1]), H)
        else:
            rds = sorted(glob.glob(f"{D}/runs/seed{s}_gp*"))
            rds = [r for r in rds if os.path.exists(f"{r}/out/final_results.txt")]
            if len(rds) < 10:
                continue
            fr = np.zeros(len(keys)); vt = np.zeros(len(keys)); sk = []
            for rd in rds:
                fin = read_final(f"{rd}/out/final_results.txt")
                fr += np.array([int(fin[k][0] in "+-" and fin[k][0] != "0") for k in keys]) / len(rds)
                vt += np.array([fin[k][1] for k in keys]) / len(rds)
                sk.append(skill(free_run(load_model(rd, N), H[:, 0], H.shape[1]), H))
            res[c]["auc"][s] = roc_auc_score(y, fr); res[c]["auc_vote"][s] = roc_auc_score(y, vt)
            res[c]["sk"][s] = float(np.mean(sk))

LAB = {"s1": "1 real series", "p3": "1 series + -perturb 3 (4 series)", "p9": "1 series + -perturb 9 (10 series)", "gp": "1 series + 10 GP samples", "s2": "2 real series", "s3": "3 real series", "s4": "4 real series"}
print("%-34s %-17s %-22s %-10s" % ("condition (%d networks)" % nseed, "AUROC", "held-out skill (median)", "stable"))
for c in CONDS:
    a = list(res[c]["auc"].values()); k = list(res[c]["sk"].values())
    if a:
        print("%-34s %.3f +- %.3f     %-22s %d of %d" % (LAB[c], np.mean(a), np.std(a), "%.2f" % np.median(k), sum(x > 0 for x in k), len(k)))
if res["gp"]["auc_vote"]:
    print("  (GP, mean vote score instead of the fraction of runs: AUROC %.3f)" % np.mean(list(res["gp"]["auc_vote"].values())))
print("\npaired by network: AUROC minus that of 1 real series, mean +- standard error")
ref = res["s1"]["auc"]
gap = None
for c in CONDS[1:]:
    d = [res[c]["auc"][s] - ref[s] for s in res[c]["auc"] if s in ref]
    if len(d) > 1:
        print("  %-34s %+.3f +- %.3f   (higher in %d of %d)" % (LAB[c], np.mean(d), np.std(d, ddof=1) / np.sqrt(len(d)), sum(x > 0 for x in d), len(d)))
        if c == "s4":
            gap = np.mean(d)
print("\nrules")
if gap:
    s4_stable = sum(x > 0 for x in res["s4"]["sk"].values())
    for c in ("p3", "p9", "gp"):
        d = [res[c]["auc"][s] - ref[s] for s in res[c]["auc"] if s in ref]
        if not d:
            continue
        st = sum(x > 0 for x in res[c]["sk"].values())
        closes = np.mean(d) >= 0.5 * gap
        stable = st >= s4_stable - 1
        print("  %-34s closes >= half of the gap (%.3f of %.3f): %s; stable held-out prediction %d vs s4 %d: %s -> %s" % (LAB[c], np.mean(d), gap, "yes" if closes else "NO", st, s4_stable, "yes" if stable else "NO", "TAKES THE PLACE" if (closes and stable) else "does not take the place"))
