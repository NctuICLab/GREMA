#!/usr/bin/env python3
"""Evaluate the simulation with a known input (sim_input_net.py + GREMA runs): the criteria of evaluation/analysis/INPUT_TERM_DESIGN.md, section 5, test 1.

usage: eval_input_sim.py SIM_DIR [--seeds 10]
SIM_DIR has seed<k>/ (sim_input_net.py) and runs/seed<k>_<mode>/out/ (GREMA, mode = input, noinput, const, shuf).
Prints, over the seeds:
  gene -> gene   AUROC of the vote and of the confidence (20 pairs per seed, 8 true), per mode; sign accuracy of the called true edges
  input -> gene  pooled AUROC of the vote over the seeds (5 pairs per seed, 2 true), true edges called, sign accuracy; for input, const and shuf
  held-out       free-running skill (1 - MSE / MSE of persistence) on the held-out condition: input model, no-input model, training mean, no-regulation ODE
Criteria (fixed before the runs): input-edge AUROC >= 0.8 and sign accuracy >= 90%; gene-gene AUROC (input) not more than 0.05 below the no-input fit;
held-out skill of the input model above the no-regulation ODE by at least 0.1; the controls (const, shuf) give an input-edge AUROC near 0.5.
"""
import sys, os
import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grema_simulate as gs

D = sys.argv[1]
nseed = int(sys.argv[sys.argv.index("--seeds") + 1]) if "--seeds" in sys.argv else 10
MODES = ["input", "noinput", "const", "shuf"]


def read_final(path):
    r = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if f[0] != "TF":
            r[(f[0], f[1])] = (f[2], float(f[3]), float(f[4]))
    return r


def read_gold(path):
    return {(l.split("\t")[0], l.split("\t")[1]): int(l.split("\t")[2]) for l in open(path)}


def read_signs(path):
    return {(l.split("\t")[0], l.split("\t")[1]): l.split("\t")[2].strip() for l in open(path)}


res = {m: dict(gauc_v=[], gauc_c=[], gsign=[], iv=[], il=[], icalled=[], isign=[]) for m in MODES}
for s in range(1, nseed + 1):
    gold = read_gold(f"{D}/seed{s}/gold.tsv"); signs = read_signs(f"{D}/seed{s}/signs.tsv")
    for m in MODES:
        f = f"{D}/runs/seed{s}_{m}/out/final_results.txt"
        if not os.path.exists(f):
            continue
        fin = read_final(f)
        gp = [k for k in gold if k[0].startswith("G")]
        y = [gold[k] for k in gp]
        res[m]["gauc_v"].append(roc_auc_score(y, [fin[k][2] for k in gp]))
        res[m]["gauc_c"].append(roc_auc_score(y, [fin[k][1] if fin[k][0] != "0" else 0 for k in gp]))
        called = [k for k in gp if gold[k] and fin[k][0] in "+-" and fin[k][0] != "0"]
        res[m]["gsign"] += [int(fin[k][0] == ("+" if signs[k] == "+1" else "-")) for k in called]
        if m != "noinput":
            ip = [k for k in gold if k[0] == "U1"]
            res[m]["iv"] += [(gold[k], fin[k][2]) for k in ip]
            res[m]["il"] += [(gold[k], fin[k][1] if fin[k][0] != "0" else 0) for k in ip]
            tc = [k for k in ip if gold[k] and fin[k][0] in "+-"]
            res[m]["icalled"].append((len(tc), sum(gold[k] for k in ip)))
            res[m]["isign"] += [int(fin[k][0] == ("+" if signs[k] == "+1" else "-")) for k in tc]


def ms(v):
    return f"{np.mean(v):.3f} +- {np.std(v):.3f}" if len(v) else "n/a"


print(f"seeds {nseed}\n")
print("gene -> gene (20 pairs, 8 true per seed; mean +- sd over the seeds)")
for m in MODES:
    r = res[m]
    if r["gauc_v"]:
        print(f"  {m:8s} AUROC vote {ms(r['gauc_v'])}  conf {ms(r['gauc_c'])}   sign right {sum(r['gsign'])}/{len(r['gsign'])} of the called true edges  ({len(r['gauc_v'])} seeds)")
print("\ninput -> gene (5 pairs, 2 true per seed; pooled over the seeds)")
for m in ("input", "const", "shuf"):
    r = res[m]
    if r["iv"]:
        y = [a for a, _ in r["iv"]]
        tc = sum(a for a, _ in r["icalled"]); tt = sum(b for _, b in r["icalled"])
        print(f"  {m:8s} AUROC vote {roc_auc_score(y, [b for _, b in r['iv']]):.3f}  conf {roc_auc_score(y, [b for _, b in r['il']]):.3f}   true edges called {tc}/{tt}   sign right {sum(r['isign'])}/{len(r['isign'])}")

# held-out prediction
rows = []
for s in range(1, nseed + 1):
    hf = f"{D}/seed{s}/heldout.txt"
    genes, inputs, H = gs.load_profile(hf, ["U1"])
    Hs = H[0]; G = len(genes); T = Hs.shape[1]
    obs, u = Hs[:G], Hs[G:]
    _, _, Tr = gs.load_profile(f"{D}/seed{s}/train_input.txt", ["U1"])
    tr = Tr[:, :G, :]
    out = {}
    for m in ("input", "noinput"):
        rd = f"{D}/runs/seed{s}_{m}"
        if os.path.isdir(rd):
            V = G + (1 if m == "input" else 0)
            P = gs.load_model(rd, G, V)
            out[m] = gs.metrics(gs.free_run(P, obs[:, 0], u if m == "input" else np.zeros((0, T)), T), obs)[1]
    out["training mean"] = gs.metrics(np.repeat(tr.mean(axis=(0, 2))[:, None], T, axis=1), obs)[1]
    c = np.zeros(G); dg = np.zeros(G)
    for i in range(G):
        xs = tr[:, i, :-1].ravel(); ds = np.diff(tr[:, i, :], axis=1).ravel()
        c[i], dg[i] = np.linalg.lstsq(np.vstack([np.ones_like(xs), -xs]).T, ds, rcond=None)[0]
    X = np.zeros((G, T)); X[:, 0] = obs[:, 0]
    for t in range(T - 1):
        X[:, t + 1] = np.clip(X[:, t] + c - dg * X[:, t], gs.XMIN, gs.XMAX)
    out["no-regulation ODE"] = gs.metrics(X, obs)[1]
    rows.append(out)
print("\nheld-out condition, free-running skill against persistence (mean +- sd over the seeds)")
for k in ("input", "noinput", "training mean", "no-regulation ODE"):
    v = [r[k] for r in rows if k in r]
    print(f"  {k:18s} {ms(v)}   ({len(v)} seeds)")
d = [r["input"] - r["no-regulation ODE"] for r in rows if "input" in r]
if d:
    print(f"  input - no-regulation ODE: mean {np.mean(d):.3f}, input better in {sum(x > 0 for x in d)}/{len(d)} seeds")

print("\ncriteria")
r = res["input"]
if r["iv"]:
    y = [a for a, _ in r["iv"]]
    auc = roc_auc_score(y, [b for _, b in r["iv"]])
    sg = np.mean(r["isign"]) if r["isign"] else float("nan")
    print(f"  input-edge AUROC >= 0.8:                    {auc:.3f}   {'PASS' if auc >= 0.8 else 'FAIL'}")
    print(f"  input-edge sign accuracy >= 90%:            {sg:.3f}   {'PASS' if sg >= 0.9 else 'FAIL'}")
if res["input"]["gauc_v"] and res["noinput"]["gauc_v"]:
    dd = np.mean(res["input"]["gauc_v"]) - np.mean(res["noinput"]["gauc_v"])
    print(f"  gene-gene AUROC (input) - (no input) >= -0.05: {dd:+.3f}   {'PASS' if dd >= -0.05 else 'FAIL'}")
if d:
    print(f"  held-out skill above the no-regulation ODE by >= 0.1: {np.mean(d):+.3f}   {'PASS' if np.mean(d) >= 0.1 else 'FAIL'}")
for m in ("const", "shuf"):
    if res[m]["iv"]:
        y = [a for a, _ in res[m]["iv"]]
        a = roc_auc_score(y, [b for _, b in res[m]["iv"]])
        print(f"  control {m}: input-edge AUROC near 0.5 (0.4-0.6): {a:.3f}   {'PASS' if 0.4 <= a <= 0.6 else 'FAIL'}")
