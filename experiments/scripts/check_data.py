#!/usr/bin/env python3
"""Check a GREMA input file before the run: is it a data set from which a method can infer regulation?

usage: check_data.py EXPRESSION.txt
Prints the size of the data, three diagnostics and warnings. The thresholds are rules of thumb from the data sets of this repository (IRMA, SOS, the
Hackett TF induction data, simulations; see experiments/METHODS_REPORT.md and experiments/PROTEIN_REGULATOR_REPORT.md), not tests with a guarantee.
  series, points, samples   samples = series x (points - 1) are the rows from which the rate of change is learnt
  mean |correlation|        mean absolute correlation between the genes (all series pooled): high = the genes move together (a shared driver)
  series similarity         mean correlation of the same gene between two different series: high = the series have the same time course (all
                            experiments are one stimulus, no different perturbations); scores can then come from the shared time course and not from the
                            regulation (SOS: 0.66, the experiment-permutation control left the scores unchanged; simulated shared drivers 0.80-0.96;
                            data with different perturbations: Hackett 0.03-0.05, DREAM4 0.16, simulations with random initial values 0.10-0.23)
Only numpy is used.
"""
import sys
import numpy as np

series = {}
genes = []
for line in open(sys.argv[1]):
    p = line.rstrip("\n").split("\t")
    if len(p) > 3 and p[0].startswith("rep"):
        series.setdefault(p[0], {})[p[1]] = np.array([float(x) for x in p[2:]])
        if p[1] not in genes:
            genes.append(p[1])
reps = sorted(series, key=lambda s: int(s[3:]))
X = [np.array([series[r][g] for g in genes]) for r in reps]      # series x (genes x points)
R, G, T = len(X), len(genes), X[0].shape[1]
print("genes %d, series %d, points per series %d, samples (series x (points - 1)) %d" % (G, R, T, R * (T - 1)))
pooled = np.concatenate(X, axis=1)
C = np.corrcoef(pooled)
mc = float(np.abs(C[np.triu_indices(G, 1)]).mean())
print("mean |correlation| between the genes: %.2f" % mc)
sim = None
if R >= 2:
    cs = []
    for g in range(G):
        for i in range(R):
            for j in range(i + 1, R):
                a, b = X[i][g], X[j][g]
                if a.std() > 1e-9 and b.std() > 1e-9:
                    cs.append(np.corrcoef(a, b)[0, 1])
    sim = float(np.mean(cs)) if cs else None
    print("series similarity (mean correlation of the same gene between series): %s" % ("%.2f" % sim if sim is not None else "n/a"))
else:
    print("series similarity: only one series")
warn = []
if R < 3:
    warn.append("only %d series: the detection of the regulations is much weaker with one or two series (simulations: AUROC 0.91 with 4 series, 0.82 with 1)." % R)
if T < 8:
    warn.append("only %d time points per series." % T)
if R * (T - 1) < 30:
    warn.append("only %d samples for %d genes: models with many parameters are not determined; use few slots and few candidate regulators." % (R * (T - 1), G))
if mc >= 0.55:
    warn.append("the genes are strongly correlated (mean |r| %.2f): a shared driver; the detection is weak, any gene can stand in for another." % mc)
if sim is not None and sim >= 0.5:
    warn.append("the series have the same time course (similarity %.2f): there are no different perturbations; a high score can come from the shared time course and not from the regulation; run the experiment-permutation control (permute_replicates.py)." % sim)
print()
if warn:
    for w in warn:
        print("WARNING:", w)
else:
    print("no warning from these rules of thumb.")
print("\nAlways: the measured variable should be the activity of the regulator (mRNA of a TF that acts as mRNA; not a promoter activity or a level that is changed after the translation);")
print("the signs found were right in simulations (92-100% when the model is right) but only 25-50% on IRMA (known signs): do not trust them without a check.")
