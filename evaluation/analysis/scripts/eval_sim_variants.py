"""Summary of GREMA on the simulated signed networks and their misspecified variants (sim_signed_net.py).

usage: eval_sim_variants.py SIM_DIR      (SIM_DIR holds net1-5, runs/, var_<name>/net1-5 and var_<name>/runs/; see the archive of the report)
Per setting: AUROC / AUPRC of VOTE_SCORE, the true edges called as a regulation, the share of those with the true sign and the same for the K best edges by vote.
"""
import glob, os, sys
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
D = sys.argv[1].rstrip("/")

def one(netdir, fn):
    gold, sign = {}, {}
    for line in open(os.path.join(netdir, "gold.tsv")):
        f = line.split()
        if f[0] != f[1]:
            gold[(f[0], f[1])] = int(f[2])
    for line in open(os.path.join(netdir, "signs.tsv")):
        f = line.split(); sign[(f[0], f[1])] = int(f[2])
    keys = sorted(gold); y = np.array([gold[k] for k in keys]); K = int(y.sum())
    res = {}
    for line in open(fn):
        f = line.rstrip("\n").split("\t")
        if f[0] != "TF":
            res[(f[0], f[1])] = (f[2], float(f[4]))
    vote = np.array([res[k][1] for k in keys])
    called = [k for k in sign if res[k][0] in ("+", "-")]
    correct = [k for k in called if (res[k][0] == "+") == (sign[k] > 0)]
    top = [keys[i] for i in np.argsort(-vote, kind="stable")[:K]]
    top_true = [k for k in top if gold[k] == 1]
    top_correct = [k for k in top_true if res[k][0] in ("+", "-") and (res[k][0] == "+") == (sign[k] > 0)]
    return (roc_auc_score(y, vote), average_precision_score(y, vote), y.mean(), len(sign), len(called), len(correct), len(top_true), len(top_correct), K)

def summarize(name, netdirs_files):
    rows = np.array([one(a, b) for a, b in netdirs_files])
    edges, called, correct = rows[:, 3].sum(), rows[:, 4].sum(), rows[:, 5].sum()
    tt, tc, kk = rows[:, 6].sum(), rows[:, 7].sum(), rows[:, 8].sum()
    print("%-9s runs %2d  AUROC %.3f  AUPRC %.3f (random %.2f)  true edges called %3d/%3d (%.0f%%)  sign right %3d/%3d (%.0f%%)  top-K true %3d/%3d (%.0f%%), sign right %3d/%3d" % (
        name, len(rows), rows[:, 0].mean(), rows[:, 1].mean(), rows[:, 2].mean(), called, edges, 100 * called / edges, correct, called, 100 * correct / max(called, 1), tt, kk, 100 * tt / kk, tc, tt))

summarize("baseline", [("%s/net%d" % (D, n), fn) for n in range(1, 6) for fn in sorted(glob.glob("%s/runs/net%d_r*/**/final_results.txt" % (D, n), recursive=True))])
for v in ["noise10", "noise20", "points8", "points5", "series1", "additive", "hidden2"]:
    summarize(v, [("%s/var_%s/net%d" % (D, v, n), fn) for n in range(1, 6) for fn in sorted(glob.glob("%s/var_%s/runs/net%d_r*/**/final_results.txt" % (D, v, n), recursive=True))])
