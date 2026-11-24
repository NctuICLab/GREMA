#!/usr/bin/env python3
"""Evaluate GREMA runs against a gold standard with AUROC / AUPRC.

usage: eval_runs.py --gold GOLD.tsv [--self-pairs] [--by REGEX] LABEL=GLOB [LABEL=GLOB ...]

GLOB selects final_results.txt files (quote it so the shell does not expand it). For every group three
scores are used for ranking the edges: binary (regulation found or not), conf (CONFIDENCE_LEVEL, 0 for no
regulation) and vote (VOTE_SCORE). Printed for each score:
  per-run  mean AUROC / AUPRC over the files (and the standard deviation of the AUROC)
  pooled   AUROC / AUPRC of the mean score over all files (the "vote fraction" of the GP experiments)
  --by REGEX: the files are split into subgroups by the first capture group of REGEX applied to the path
              (e.g. 'seed(\\d+)'); pooled scores are computed per subgroup and their mean and sd are printed
              ("seed-mean" in the GP experiments).
The gold file has one line per pair: TF<TAB>GENE<TAB>0/1. Self pairs are skipped unless --self-pairs.
Random AUROC is 0.5; random AUPRC is the fraction of true edges.
"""
import argparse, glob, re, sys
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

NAMES = ["binary", "conf", "vote"]


def load_gold(path, self_pairs):
    gold = {}
    for line in open(path):
        f = line.split()
        if len(f) >= 3 and (self_pairs or f[0] != f[1]):
            gold[(f[0], f[1])] = int(f[2])
    keys = sorted(gold)
    return keys, np.array([gold[k] for k in keys])


def load_run(path, keys):
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "TF" or len(f) < 5:
            continue
        reg = f[2] != "0"
        d[(f[0], f[1])] = (float(reg), float(f[3]) if reg else 0.0, float(f[4]), f[2])
    missing = [k for k in keys if k not in d]
    if missing:
        sys.exit("%s: pairs missing: %s" % (path, missing[:3]))
    return np.array([[d[k][i] for k in keys] for i in range(3)]), [d[k][3] for k in keys]


def metrics(y, s):
    return roc_auc_score(y, s), average_precision_score(y, s)


def fmt(a):
    return "%.3f / %.3f" % (a[0], a[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gold", required=True)
    ap.add_argument("--self-pairs", action="store_true")
    ap.add_argument("--by")
    ap.add_argument("groups", nargs="+", metavar="LABEL=GLOB")
    a = ap.parse_args()
    keys, y = load_gold(a.gold, a.self_pairs)
    print("pairs=%d true=%d random AUPRC=%.3f" % (len(keys), y.sum(), y.mean()))
    for g in a.groups:
        label, pat = g.split("=", 1)
        files = sorted(glob.glob(pat))
        if not files:
            print(label, ": no files for", pat); continue
        runs, signs = zip(*[load_run(f, keys) for f in files])
        runs = np.array(runs)  # (files, 3, pairs)
        tp = [int(((r[0] == 1) & (y == 1)).sum()) for r in runs]
        fp = [int(((r[0] == 1) & (y == 0)).sum()) for r in runs]
        tps = [s for sg in signs for s, yy in zip(sg, y) if yy == 1 and s != "0"]
        print("\n%s  (%d runs)  TP %.1f  FP %.1f  signs of TP edges %s" % (
            label, len(files), np.mean(tp), np.mean(fp), {s: tps.count(s) for s in sorted(set(tps))}))
        for i, name in enumerate(NAMES):
            per = np.array([metrics(y, r[i]) for r in runs])
            line = "  %-7s per-run %s (sd AUROC %.3f)  pooled %s" % (
                name, fmt(per.mean(0)), per[:, 0].std(), fmt(metrics(y, runs[:, i].mean(0))))
            if a.by:
                sub = {}
                for f, r in zip(files, runs):
                    m = re.search(a.by, f)
                    sub.setdefault(m.group(1) if m else "?", []).append(r[i])
                ps = np.array([metrics(y, np.mean(v, axis=0)) for v in sub.values()])
                line += "  | %d subgroups: mean %s (sd AUROC %.3f)" % (len(sub), fmt(ps.mean(0)), ps[:, 0].std())
            print(line)


if __name__ == "__main__":
    main()
