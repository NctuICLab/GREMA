#!/usr/bin/env python3
"""Baselines for the hMSC design, scored on the same pairs and the same gold standard as GREMA (eval_runs.py).

usage: hmsc_baselines.py DESIGN_DIR
Scores for every TF -> gene pair of DESIGN_DIR/gold.tsv (columns TF GENE gold atac_candidate sign):
  footprint   1 if the TF has a footprint near the gene (the ATAC prior alone), else 0
  |pearson|   absolute Pearson correlation of the two 5-point profiles; |pearson lag|  of the TF at t and the gene change from t to t+1 (4 points)
  |spearman|  absolute Spearman correlation
  random      mean AUROC of 200 random scores
Also: AUROC among the candidates with a footprint only (does the ranking inside the prior help?), and the number of gold pairs.
"""
import sys, random
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import spearmanr

D = sys.argv[1]
expr = {}
for l in open(f"{D}/expression.txt"):
    f = l.rstrip("\n").split("\t")
    if f[0] == "rep1":
        expr[f[1]] = np.array([float(x) for x in f[2:]])
rows = [l.rstrip("\n").split("\t") for l in open(f"{D}/gold.tsv")][1:]
y = np.array([int(r[2]) for r in rows]); fp = np.array([int(r[3]) for r in rows])
pearson = np.array([abs(np.corrcoef(expr[r[0]], expr[r[1]])[0, 1]) for r in rows])
spear = np.array([abs(spearmanr(expr[r[0]], expr[r[1]])[0]) for r in rows])
lag = np.array([abs(np.corrcoef(expr[r[0]][:-1], np.diff(expr[r[1]]))[0, 1]) for r in rows])
lag = np.nan_to_num(lag)
print(f"pairs {len(y)}, gold {y.sum()}, with a footprint {fp.sum()}, gold with a footprint {int((y * fp).sum())}")
for name, s in [("footprint", fp), ("|pearson|", pearson), ("|spearman|", spear), ("|pearson lag|", lag)]:
    print(f"{name:14s} AUROC {roc_auc_score(y, s):.3f}  AUPRC {average_precision_score(y, s):.3f}")
random.seed(1)
r = [roc_auc_score(y, np.random.RandomState(k).rand(len(y))) for k in range(200)]
print(f"random         AUROC {np.mean(r):.3f} (sd {np.std(r):.3f}), AUPRC {y.mean():.3f}")
m = fp == 1
if y[m].sum() and (1 - y[m]).sum():
    print(f"inside the footprint candidates ({m.sum()} pairs, {int(y[m].sum())} gold): |pearson| AUROC {roc_auc_score(y[m], pearson[m]):.3f}")
