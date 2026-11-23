#!/usr/bin/env python3
"""How far has the genetic algorithm converged? Uses the progress files written by the solver.

usage: convergence.py 'GLOB' [GLOB ...]        e.g. 'runs/irma_on_f1_*/out/tmp/progress*_step1.txt'

Each progress file has one line per run and checkpoint (every -g/10 generations):
  run  generations_since_last_improvement  generation  evaluations  best_fitness  cc
The last checkpoint is generation 0.9*g (450 for -g 500). Prints, over all runs found:
  - the share of runs that still improved in the last 100 / 50 generations before the last checkpoint
  - the relative fitness gain from generation 0.5*g and from 0.2*g to the last checkpoint (median, 90th pct, max)
Use step1 files (the first GREMA iteration, all genes) to compare settings. g = 500 is assumed.
"""
import glob, sys
import numpy as np

G_LAST, G_MID, G_EARLY = 450, 250, 100
rows = []
for pat in sys.argv[1:]:
    for f in glob.glob(pat):
        d = {}
        for line in open(f):
            p = line.split()
            if len(p) >= 6:
                d.setdefault(int(p[0]), {})[int(p[2])] = (int(p[1]), float(p[4]))
        for g in d.values():
            if G_LAST in g and G_MID in g and G_EARLY in g:
                rows.append((g[G_LAST][0], g[G_EARLY][1], g[G_MID][1], g[G_LAST][1]))
if not rows:
    sys.exit("no progress data found (needs checkpoints at generations %d, %d and %d, i.e. -g 500)" % (G_EARLY, G_MID, G_LAST))
a = np.array(rows)
same, f100, f250, f450 = a.T
imp = (f250 - f450) / np.maximum(f250, 1e-12)
imp100 = (f100 - f450) / np.maximum(f100, 1e-12)
print("runs: %d" % len(a))
print("still improving in the last 100 generations: %.1f%%, in the last 50: %.1f%%" % (np.mean(same < 100) * 100, np.mean(same < 50) * 100))
print("relative fitness gain gen 250 -> 450: median %.2f%%, 90th pct %.1f%%, max %.1f%%" % (np.median(imp) * 100, np.percentile(imp, 90) * 100, imp.max() * 100))
print("relative fitness gain gen 100 -> 450: median %.2f%%, 90th pct %.1f%%" % (np.median(imp100) * 100, np.percentile(imp100, 90) * 100))
