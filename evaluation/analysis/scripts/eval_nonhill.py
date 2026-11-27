#!/usr/bin/env python3
"""Evaluate the GREMA runs on data with non-Hill regulation functions (sim_nonhill.py): is the limit of the Hill model the form or the data?

usage: eval_nonhill.py SIM_DIR [--seeds 6] [--dyn DYN_DIR]
SIM_DIR has <form>_s<series>/seed<k>/ (sim_nonhill.py) and runs/<form>_s<series>_seed<k>/out/ (GREMA, no prior knowledge, -f 1 -g 500); DYN_DIR (optional) has
<form>_s<series>/seed<k>/seed<m>/final_results.txt from dyngenie3_runner.py (5 seeds of the random forest, averaged).
Per form and number of series (mean +- sd over the seeds):
  AUROC     vote score of GREMA against the 8 true regulations of 20 pairs; dynGENIE3 on the same data; for bell the AUROC of all 20 pairs too
  sign      fraction of the true monotone edges called by GREMA (+ or -) that have the right sign
  held-out  free-running skill (1 - MSE / MSE of persistence) on a held-out series of the same system, GREMA against a no-regulation ODE fitted on the training series
Rules fixed before the runs (see the report): the AUROC drop of a form = AUROC(hill) - AUROC(form) with the same number of series;
< 0.05 and sign accuracy >= 85%: the Hill fit is robust to this form; drop >= 0.10 or sign accuracy < 70%: the form limits the Hill fit; between: moderate.
The effect of the data = AUROC(4 series) - AUROC(1 series) for the Hill form is compared with the largest drop of a form.
"""
import sys, os, glob, math
import numpy as np
from sklearn.metrics import roc_auc_score

D = sys.argv[1]
nseed = int(sys.argv[sys.argv.index("--seeds") + 1]) if "--seeds" in sys.argv else 6
DYN = sys.argv[sys.argv.index("--dyn") + 1] if "--dyn" in sys.argv else None
FORMS = ["hill", "or", "bell", "delay", "linear"]
XMIN, XMAX = 0.0001, 100000.0
N = 5


def read_final(path):
    r = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if f[0] != "TF":
            r[(f[0], f[1])] = (f[2], float(f[3]), float(f[4]))
    return r


def read_series(path):
    rows = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) > 3 and f[0].startswith("rep"):
            rows.setdefault(f[0], {})[f[1]] = [float(x) for x in f[2:]]
    return np.array([[rows[r]["G%d" % (i + 1)] for i in range(N)] for r in sorted(rows, key=lambda x: int(x[3:]))])


def load_model(rundir):
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
    X = np.zeros((N, T)); X[:, 0] = x0
    for t in range(T - 1):
        for i in range(N):
            p = P[i]; temp = 1.0
            for j in range(N):
                if p["n"][j] != 0:
                    temp *= 1.0 / (1.0 + (p["k"][j] / X[j, t]) ** p["n"][j])
            act = max(p["tm"] * temp + p["b"], 0.0)
            X[i, t + 1] = min(max(X[i, t] + act - p["deg"] * X[i, t], XMIN), XMAX)
    return X


def skill(pred, obs):
    pers = np.repeat(obs[:, :1], obs.shape[1], axis=1)
    return float((1 - ((pred - obs) ** 2).mean(1) / np.maximum(((pers - obs) ** 2).mean(1), 1e-12)).mean())


def noreg_skill(tr, obs):
    c = np.zeros(N); dg = np.zeros(N)
    for i in range(N):
        xs = tr[:, i, :-1].ravel(); ds = np.diff(tr[:, i, :], axis=1).ravel()
        c[i], dg[i] = np.linalg.lstsq(np.vstack([np.ones_like(xs), -xs]).T, ds, rcond=None)[0]
    Y = np.zeros((N, obs.shape[1])); Y[:, 0] = obs[:, 0]
    for t in range(obs.shape[1] - 1):
        Y[:, t + 1] = np.clip(Y[:, t] + c - dg * Y[:, t], XMIN, XMAX)
    return skill(Y, obs)


def ms(v):
    return "%.3f +- %.3f" % (np.mean(v), np.std(v)) if len(v) else "n/a"


R = {}
for size in (4, 1):
    for form in FORMS:
        a = dict(auc=[], auc_c=[], dyn=[], sign=[], sk=[], sk0=[], bell=[], hit=[], by_seed={})
        for s in range(1, nseed + 1):
            sd = f"{D}/{form}_s{size}/seed{s}"; rd = f"{D}/runs/{form}_s{size}_seed{s}"
            if not os.path.exists(f"{rd}/out/final_results.txt"):
                continue
            gold = {(l.split("\t")[0], l.split("\t")[1]): int(l.split("\t")[2]) for l in open(f"{sd}/gold.tsv")}
            typ = {(l.split("\t")[0], l.split("\t")[1]): (l.split("\t")[2], l.split("\t")[3].strip()) for l in open(f"{sd}/signs.tsv")}
            fin = read_final(f"{rd}/out/final_results.txt")
            keys = list(gold)
            a["auc"].append(roc_auc_score([gold[k] for k in keys], [fin[k][2] for k in keys]))
            a["by_seed"][s] = a["auc"][-1]
            a["auc_c"].append(roc_auc_score([gold[k] for k in keys], [fin[k][1] if fin[k][0] != "0" else 0 for k in keys]))
            for k in keys:
                if gold[k] and typ[k][1] == "monotone" and fin[k][0] in "+-" and fin[k][0] != "0":
                    a["sign"].append(int(fin[k][0] == ("+" if typ[k][0] == "+1" else "-")))
            for k in keys:
                if gold[k] and typ[k][1] == "bell":
                    a["bell"].append(fin[k][2])
            tr = read_series(f"{sd}/train.txt"); ho = read_series(f"{sd}/heldout.txt")[0]
            a["sk"].append(skill(free_run(load_model(rd), ho[:, 0], ho.shape[1]), ho))
            a["sk0"].append(noreg_skill(tr, ho))
            if DYN:
                ds = []
                for f in glob.glob(f"{DYN}/{form}_s{size}/seed{s}/seed*/final_results.txt"):
                    fr = read_final(f)
                    ds.append([fr[k][2] for k in keys])
                if ds:
                    a["dyn"].append(roc_auc_score([gold[k] for k in keys], np.mean(ds, axis=0)))
        R[(form, size)] = a

for size in (4, 1):
    print(f"\n== {size} training series x 16 points")
    print("%-7s %-17s %-17s %-17s %-8s %-17s %-17s" % ("form", "AUROC GREMA", "AUROC dynGENIE3", "held-out skill", "(no-reg)", "sign right", "bell-edge vote"))
    for form in FORMS:
        a = R[(form, size)]
        sg = "%d/%d" % (sum(a["sign"]), len(a["sign"])) if a["sign"] else "n/a"
        print("%-7s %-17s %-17s %-17s %-8s %-17s %-17s  (%d seeds)" % (form, ms(a["auc"]), ms(a["dyn"]), ms(a["sk"]), "%.3f" % np.mean(a["sk0"]) if a["sk0"] else "n/a", sg, ms(a["bell"]) if a["bell"] else "-", len(a["auc"])))
print("\nrules")
for size in (4, 1):
    h = np.mean(R[("hill", size)]["auc"])
    for form in FORMS[1:]:
        a = R[(form, size)]
        if not a["auc"]:
            continue
        drop = h - np.mean(a["auc"])
        sg = np.mean(a["sign"]) if a["sign"] else float("nan")
        verdict = "robust" if (drop < 0.05 and sg >= 0.85) else ("LIMITS the Hill fit" if (drop >= 0.10 or sg < 0.70) else "moderate")
        print(f"  {size} series, {form:6s}: AUROC drop {drop:+.3f}, sign accuracy {sg:.2f}  -> {verdict}")
print("\npaired by seed (the parameters and the initial values of a seed are the same in every form): AUROC of the form minus AUROC of hill, mean +- standard error")
for size in (4, 1):
    for form in FORMS[1:]:
        d = [R[(form, size)]["by_seed"][s] - R[("hill", size)]["by_seed"][s] for s in R[(form, size)]["by_seed"] if s in R[("hill", size)]["by_seed"]]
        if len(d) > 1:
            print(f"  {size} series, {form:6s}: {np.mean(d):+.3f} +- {np.std(d, ddof=1) / np.sqrt(len(d)):.3f}   ({sum(x < 0 for x in d)} of {len(d)} seeds lower than hill)")
d = [R[("hill", 4)]["by_seed"][s] - R[("hill", 1)]["by_seed"][s] for s in R[("hill", 4)]["by_seed"] if s in R[("hill", 1)]["by_seed"]]
print(f"  hill, 4 against 1 series (different data of the same seed): {np.mean(d):+.3f} +- {np.std(d, ddof=1) / np.sqrt(len(d)):.3f}")
g4, g1 = np.mean(R[("hill", 4)]["auc"]), np.mean(R[("hill", 1)]["auc"])
drops = [np.mean(R[("hill", 4)]["auc"]) - np.mean(R[(f, 4)]["auc"]) for f in FORMS[1:] if R[(f, 4)]["auc"]]
print(f"  effect of the data (hill, 4 against 1 series): {g4 - g1:+.3f}; the largest drop of a form with 4 series: {max(drops):+.3f}")
