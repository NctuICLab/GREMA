#!/usr/bin/env python3
"""Check a GREMA input file before the run: is it a data set from which a method can infer regulation?

usage: check_data.py EXPRESSION.txt [--input U1,U2]
--input  the rows that are external inputs (GREMA option -input): they are left out of the gene statistics and checked as inputs (see below)
Prints the size of the data, three diagnostics and warnings. The thresholds are rules of thumb from the data sets of this repository (IRMA, SOS, the
Hackett TF induction data, simulations; see evaluation/analysis/METHODS_REPORT.md and evaluation/analysis/PROTEIN_REGULATOR_REPORT.md), not tests with a guarantee.
  series, points, samples   samples = series x (points - 1) are the rows from which the rate of change is learnt
  mean |correlation|        mean absolute correlation between the genes (all series pooled): high = the genes move together (a shared driver)
  series similarity         mean correlation of the same gene between two different series: high = the series have the same time course (all
                            experiments are one stimulus, no different perturbations); scores can then come from the shared time course and not from the
                            regulation (SOS: 0.66, the experiment-permutation control left the scores unchanged; simulated shared drivers 0.80-0.96;
                            data with different perturbations: Hackett 0.03-0.05, DREAM4 0.16, simulations with random initial values 0.10-0.23)
Inputs (with --input): warnings when an input is not on the scale of the genes (0.01 to 1), is constant, is the same in every series (it then has the shape of the
time index and cannot be told from the programme of the cell: a single stimulus in one series, as in a differentiation time course), when the series have
different numbers of points, and a hint on the design (at least two different levels or start times and at least 8 points).
Only numpy is used.
"""
import sys
import numpy as np

args = sys.argv[1:]
input_names = []
if "--input" in args:
    i = args.index("--input")
    input_names = [x for x in args[i + 1].split(",") if x]
    del args[i:i + 2]
series = {}
genes = []
inputs = []
lengths = set()
for line in open(args[0]):
    p = line.rstrip("\n").split("\t")
    if len(p) > 3 and p[0].startswith("rep"):
        lengths.add(len(p) - 2)
        series.setdefault(p[0], {})[p[1]] = np.array([float(x) for x in p[2:]])
        if p[1] in input_names:
            if p[1] not in inputs:
                inputs.append(p[1])
        elif p[1] not in genes:
            genes.append(p[1])
for nm in input_names:
    if nm not in inputs:
        print("WARNING: -input %s is not a row of the file" % nm)
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
note = []
if R < 3:
    warn.append("only %d series: the detection of the regulations is much weaker with one or two series (simulations: AUROC 0.91 with 4 series, 0.82 with 1)." % R)
if T < 8:
    warn.append("only %d time points per series." % T)
if R * (T - 1) < 30:
    warn.append("only %d samples for %d genes: models with many parameters are not determined; use few slots and few candidate regulators." % (R * (T - 1), G))
if mc >= 0.55:
    warn.append("the genes are strongly correlated (mean |r| %.2f): a shared driver; the detection is weak, any gene can stand in for another." % mc)
varying_input = False
for nm in inputs:
    U = np.array([series[r][nm] for r in reps])
    if len({tuple(np.round(u, 4)) for u in U}) >= 2 and not all(u.std() < 1e-9 for u in U):
        varying_input = True
if sim is not None and sim >= 0.5 and varying_input:
    note.append("the series have a similar time course (similarity %.2f), as expected when one stimulus drives them; the input differs between the series, so the input term can explain this common part (check that the input edges are the same when the input is shuffled between the series)." % sim)
elif sim is not None and sim >= 0.5:
    warn.append("the series have the same time course (similarity %.2f): there are no different perturbations; a high score can come from the shared time course and not from the regulation; run the experiment-permutation control (permute_replicates.py)." % sim)
if inputs:
    print("\ninputs: %s" % ", ".join(inputs))
    if len(lengths) > 1:
        warn.append("the rows have different numbers of points (%s): GREMA needs the same number of time points in every series." % sorted(lengths))
    for nm in inputs:
        U = np.array([series[r][nm] for r in reps])                    # series x points
        lo, hi = float(U.min()), float(U.max())
        print("  %s: range %.3f to %.3f, %d distinct profiles in %d series" % (nm, lo, hi, len({tuple(np.round(u, 4)) for u in U}), R))
        if lo < 0.01 - 1e-9 or hi > 1 + 1e-9:
            warn.append("%s is not on the scale of the genes (0.01 to 1): scale it with the same linear map as the genes (GREMA cannot use 0)." % nm)
        if all(u.std() < 1e-9 for u in U) and len({round(float(u[0]), 4) for u in U}) == 1:
            warn.append("%s is constant and the same in every series: it carries no information; leave it out." % nm)
        elif len({tuple(np.round(u, 4)) for u in U}) == 1:
            warn.append("%s is the same in every series: it has the shape of the time index and cannot be told from the programme of the cell (a differentiation time course, a single stimulus); the edges from it are not identifiable. Use at least two series with a different level or start time of the stimulus." % nm)
        elif R == 1:
            warn.append("%s: only one series, the input has the shape of the time index." % nm)
        if T < 8:
            warn.append("%s: only %d time points; the effect of the input on the genes is hard to tell from the change in time." % (nm, T))
print()
for n_ in note:
    print("NOTE:", n_)
if warn:
    for w in warn:
        print("WARNING:", w)
else:
    print("no warning from these rules of thumb.")
print("\nAlways: the measured variable should be the activity of the regulator (mRNA of a TF that acts as mRNA; not a promoter activity or a level that is changed after the translation);")
print("the signs found were right in simulations (92-100% when the model is right) but only 25-50% on IRMA (known signs): do not trust them without a check.")
