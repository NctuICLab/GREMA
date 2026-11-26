#!/usr/bin/env python3
"""Evaluate the SOS runs with the UV dose as input (prep_sos_input.py + GREMA runs): held-out experiment prediction.

usage: eval_sos_input.py INPUT_DOSE_DIR
INPUT_DOSE_DIR has <fold>/ (prep_sos_input.py) and runs/<fold>_<mode>_r<k>/out/ (GREMA; mode = noinput, decay10, rev10, decay5, decay20, step).
For every fold and mode: the model of every GA run (the run of lowest fitness for each gene) is simulated free-running from the observed first point of every held-out experiment,
with the input of that experiment, and compared with the observation: skill = 1 - MSE / MSE of persistence (mean over the genes and the held-out experiments).
Baselines: persistence (0), the training mean, a no-regulation ODE (dx = c - d x per gene, fitted on the training experiments).
Main design: loeo1..loeo4 (mean over the four folds). Also printed: the single-dose folds to20 and to5 (not identifiable), the input edges (U1 -> gene, vote score, sign) and the LexA (G2) edges.
"""
import sys, os, glob, collections
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grema_simulate as gs

D = sys.argv[1]
MODES = ["noinput", "decay10", "rev10", "decay5", "decay20", "step"]


def baselines(tr, te_obs):
    G, T = te_obs.shape[1], te_obs.shape[2]
    out = {}
    out["training mean"] = np.mean([gs.metrics(np.repeat(tr.mean(axis=(0, 2))[:, None], T, axis=1), te_obs[s])[1] for s in range(te_obs.shape[0])])
    c = np.zeros(G); dg = np.zeros(G)
    for i in range(G):
        xs = tr[:, i, :-1].ravel(); ds = np.diff(tr[:, i, :], axis=1).ravel()
        c[i], dg[i] = np.linalg.lstsq(np.vstack([np.ones_like(xs), -xs]).T, ds, rcond=None)[0]
    sk = []
    for s in range(te_obs.shape[0]):
        Y = np.zeros((G, T)); Y[:, 0] = te_obs[s][:, 0]
        for t in range(T - 1):
            Y[:, t + 1] = np.clip(Y[:, t] + c - dg * Y[:, t], gs.XMIN, gs.XMAX)
        sk.append(gs.metrics(Y, te_obs[s])[1])
    out["no-regulation ODE"] = float(np.mean(sk))
    return out


table = collections.defaultdict(dict)       # fold -> mode -> [skill per run]
base = {}
edges = collections.defaultdict(lambda: collections.defaultdict(list))
for fold in sorted(os.listdir(D)):
    if not os.path.isdir(f"{D}/{fold}") or not os.path.exists(f"{D}/{fold}/train_noinput.txt"):
        continue
    g0, _, Tr = gs.load_profile(f"{D}/{fold}/train_noinput.txt")
    G = len(g0)
    _, _, Te = gs.load_profile(f"{D}/{fold}/test_noinput.txt")
    base[fold] = baselines(Tr, Te)
    for mode in MODES:
        dirs = sorted(glob.glob(f"{D}/runs/{fold}_{mode}_r*"))
        if not dirs:
            continue
        has_in = mode != "noinput"
        _, _, Tt = gs.load_profile(f"{D}/{fold}/test_{mode}.txt", ["U1"] if has_in else [])
        sk = []
        for d in dirs:
            if not os.path.exists(f"{d}/out/final_results.txt"):
                continue
            P = gs.load_model(d, G, G + (1 if has_in else 0))
            T = Tt.shape[2]
            sk.append(np.mean([gs.metrics(gs.free_run(P, Tt[s][:G, 0], Tt[s][G:] if has_in else np.zeros((0, T)), T), Tt[s][:G])[1] for s in range(Tt.shape[0])]))
            for l in open(f"{d}/out/final_results.txt"):
                p = l.rstrip("\n").split("\t")
                if p[0] == "U1" or p[0] == "G2":
                    edges[(fold.startswith("loeo"), mode, p[0])][p[1]].append((p[2], float(p[4])))
        table[fold][mode] = sk


def ms(v):
    return "%7.3f +- %5.3f" % (np.mean(v), np.std(v)) if len(v) else "      n/a      "


print("held-out free-running skill against persistence (mean +- sd over the GA runs)")
print("%-7s %-17s %-17s %-17s %-17s %-17s %-9s %-9s" % ("fold", "no input", "decay10", "rev10 (control)", "decay5", "decay20", "train mean", "no-reg ODE"))
for fold in sorted(table):
    t = table[fold]
    print("%-7s %-17s %-17s %-17s %-17s %-17s %9.3f %9.3f" % (fold, ms(t.get("noinput", [])), ms(t.get("decay10", [])), ms(t.get("rev10", [])), ms(t.get("decay5", [])), ms(t.get("decay20", [])), base[fold]["training mean"], base[fold]["no-regulation ODE"]))
lo = [f for f in table if f.startswith("loeo")]
if lo:
    print("\nmain design (leave one experiment out), mean over the folds of the mean over the GA runs:")
    for mode in MODES:
        v = [np.mean(table[f][mode]) for f in lo if mode in table[f] and table[f][mode]]
        if v:
            print(f"  {mode:8s} {np.mean(v):7.3f}  ({len(v)} folds)")
    print(f"  training mean {np.mean([base[f]['training mean'] for f in lo]):7.3f}")
    print(f"  no-regulation ODE {np.mean([base[f]['no-regulation ODE'] for f in lo]):7.3f}")
    a = [np.mean(table[f]["decay10"]) for f in lo if table[f].get("decay10")]
    b = [np.mean(table[f]["noinput"]) for f in lo if table[f].get("noinput")]
    c = [np.mean(table[f]["rev10"]) for f in lo if table[f].get("rev10")]
    r = [base[f]["no-regulation ODE"] for f in lo]
    if len(a) == len(lo) and len(b) == len(lo):
        print("criteria: decay10 above no input: %s (%+.3f); above the no-regulation ODE: %s (%+.3f)" % (
            "PASS" if np.mean(a) > np.mean(b) else "FAIL", np.mean(a) - np.mean(b), "PASS" if np.mean(a) > np.mean(r) else "FAIL", np.mean(a) - np.mean(r)))
    if len(a) == len(lo) and len(c) == len(lo):
        print("          rev10 worse than decay10 by >= 0.1: %s (%+.3f)" % ("PASS" if np.mean(a) - np.mean(c) >= 0.1 else "FAIL", np.mean(a) - np.mean(c)))
print("\nedges (main design, vote score averaged over the folds and runs; the signs of all runs):")
for (isloeo, mode, tf), d in sorted(edges.items()):
    if not isloeo or mode not in ("decay10", "noinput"):
        continue
    if tf == "U1" and mode == "noinput":
        continue
    items = ", ".join("%s %.2f %s" % (g, np.mean([v for _, v in d[g]]), "".join(sorted(s for s, _ in d[g]))) for g in sorted(d, key=lambda x: int(x[1:])))
    print(f"  {mode:8s} {tf} -> : {items}")
