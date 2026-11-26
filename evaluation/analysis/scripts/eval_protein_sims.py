#!/usr/bin/env python3
"""Detection and signs of GREMA (and dynGENIE3) on the simulated networks with protein regulators (sim_protein_net.py, augment_filtered.py).

usage: eval_protein_sims.py SIM_DIR     (the unpacked sim_protein_results.tgz)
Layout of SIM_DIR: net1-5 (no driver), drv_transcription/net1-5 and drv_cleavage/net1-5 (data and gold.tsv / signs.tsv), runs/base_net<k>, runs/d<D>_net<k> (GREMA on the
data with the filtered copies of the genes as the only regulators), drv_*/runs/net<k>, dyn/<cond>_<RF|ET>_net<k>/seed<s> (dynGENIE3), each holding final_results.txt.
Prints for GREMA: AUROC / AUPRC of VOTE_SCORE, the true edges called as a regulation, among them the share with the true sign, and the same for the K best pairs by vote
(K = number of true edges); for dynGENIE3 AUROC / AUPRC.
"""
import glob, os, sys
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

O = sys.argv[1].rstrip("/")


def load_truth(netdir):
    gold, sign = {}, {}
    for line in open(os.path.join(netdir, "gold.tsv")):
        f = line.split()
        if f[0] != f[1]:
            gold[(f[0], f[1])] = int(f[2])
    for line in open(os.path.join(netdir, "signs.tsv")):
        f = line.split()
        sign[(f[0], f[1])] = int(f[2])
    return gold, sign


def load_results(fn, filtered=False):
    res = {}
    for line in open(fn):
        f = line.rstrip("\n").split("\t")
        if f[0] == "TF":
            continue
        if filtered:           # pair F_j -> G_i is the regulation G_j -> G_i
            if f[0].startswith("F") and f[1].startswith("G"):
                res[("G" + f[0][1:], f[1])] = (f[2], float(f[4]))
        else:
            res[(f[0], f[1])] = (f[2], float(f[4]))
    return res


def one(netdir, fn, filtered=False):
    gold, sign = load_truth(netdir)
    keys = sorted(gold)
    y = np.array([gold[k] for k in keys])
    K = int(y.sum())
    res = load_results(fn, filtered)
    vote = np.array([res[k][1] for k in keys])
    called = [k for k in sign if res[k][0] in ("+", "-")]
    correct = [k for k in called if (res[k][0] == "+") == (sign[k] > 0)]
    top = [keys[i] for i in np.argsort(-vote, kind="stable")[:K]]
    tt = [k for k in top if gold[k] == 1]
    tc = [k for k in tt if res[k][0] in ("+", "-") and (res[k][0] == "+") == (sign[k] > 0)]
    return roc_auc_score(y, vote), average_precision_score(y, vote), y.mean(), len(sign), len(called), len(correct), len(tt), len(tc), K


def show(name, items, filtered=False):
    r = np.array([one(a, b, filtered) for a, b in items])
    print("%-36s runs %2d  AUROC %.3f / AUPRC %.3f (random %.2f)  true edges called %3d/%3d  sign right %3d/%3d (%3.0f%%)  top-K: true %3d/%3d, sign right %3d/%3d" % (
        name, len(r), r[:, 0].mean(), r[:, 1].mean(), r[:, 2].mean(), r[:, 4].sum(), r[:, 3].sum(), r[:, 5].sum(), r[:, 4].sum(),
        100 * r[:, 5].sum() / max(r[:, 4].sum(), 1), r[:, 6].sum(), r[:, 8].sum(), r[:, 7].sum(), r[:, 6].sum()))


def fr(pat):
    return sorted(glob.glob(pat, recursive=True))


print("== GREMA, regulator = the measured signal (default settings)")
show("no driver (random initial proteins)", [("%s/net%d" % (O, s), fr("%s/runs/base_net%d/**/final_results.txt" % (O, s))[0]) for s in range(1, 6)])
show("transcription driver", [("%s/drv_transcription/net%d" % (O, s), fr("%s/drv_transcription/runs/net%d/**/final_results.txt" % (O, s))[0]) for s in range(1, 6)])
show("cleavage driver", [("%s/drv_cleavage/net%d" % (O, s), fr("%s/drv_cleavage/runs/net%d/**/final_results.txt" % (O, s))[0]) for s in range(1, 6)])
print("== GREMA, no driver, the filtered copy of every gene is the only regulator")
for d in ["0.1", "0.2", "0.3"]:
    show("filtered, d = %s" % d, [("%s/net%d" % (O, s), fr("%s/runs/d%s_net%d/**/final_results.txt" % (O, d, s))[0]) for s in range(1, 6)], filtered=True)
print("== dynGENIE3 (no signs): AUROC / AUPRC")
for cond, base in [("none", "net"), ("transcription", "drv_transcription/net"), ("cleavage", "drv_cleavage/net")]:
    out = []
    for tree in ["RF", "ET"]:
        per = []
        for s in range(1, 6):
            gold, _ = load_truth("%s/%s%d" % (O, base, s))
            keys = sorted(gold)
            y = np.array([gold[k] for k in keys])
            for fn in fr("%s/dyn/%s_%s_net%d/seed*/final_results.txt" % (O, cond, tree, s)):
                res = load_results(fn)
                sc = np.array([res[k][1] for k in keys])
                per.append((roc_auc_score(y, sc), average_precision_score(y, sc)))
        a = np.array(per)
        out.append("%s %.3f / %.3f (%d runs)" % (tree, a[:, 0].mean(), a[:, 1].mean(), len(a)))
    print("%-14s %s" % (cond, "   ".join(out)))
