#!/usr/bin/env python3
"""How far has the genetic algorithm converged? Uses the progress files written by the solver.

usage: convergence.py 'GLOB' [GLOB ...]        e.g. 'runs/irma_on_f1_*/out/tmp/progress*_step1.txt'

Each progress file has one line per run and checkpoint (every -g/10 generations):
  run  generations_since_last_improvement  generation  evaluations  best_fitness  cc
The last checkpoint is generation 0.9*g (450 for -g 500). Prints, over all runs found:
  - the share of runs that still improved in the last 100 / 50 generations before the last checkpoint
  - the relative fitness gain from generation 0.5*g and from 0.2*g to the last checkpoint (median, 90th pct, max)
Use step1 files (the first GREMA iteration, all genes) to compare settings. g is taken from the files (the last
checkpoint is 0.9*g); the early / mid / last checkpoints are 0.2*g / 0.5*g / 0.9*g, and "still improving" is
measured over the last 0.2*g and 0.1*g generations (100 / 50 for -g 500).
"""
import glob, sys
import numpy as np

rows = []
G = 0
for pat in sys.argv[1:]:
    for f in glob.glob(pat):
        d = {}
        for line in open(f):
            p = line.split()
            if len(p) >= 6:
                d.setdefault(int(p[0]), {})[int(p[2])] = (int(p[1]), float(p[4]))
        for g in d.values():
            last = max(g)
            G = max(G, round(last / 0.9))
            mid, early = round(last * 5 / 9), round(last * 2 / 9)
            if mid in g and early in g:
                rows.append((g[last][0], g[early][1], g[mid][1], g[last][1]))
if not rows:
    sys.exit("no progress data found")
a = np.array(rows)
same, f100, f250, f450 = a.T
imp = (f250 - f450) / np.maximum(f250, 1e-12)
imp100 = (f100 - f450) / np.maximum(f100, 1e-12)
print("runs: %d, g = %d" % (len(a), G))
print("still improving in the last %d generations: %.1f%%, in the last %d: %.1f%%" % (0.2 * G, np.mean(same < 0.2 * G) * 100, 0.1 * G, np.mean(same < 0.1 * G) * 100))
print("relative fitness gain 0.5g -> 0.9g: median %.2f%%, 90th pct %.1f%%, max %.1f%%" % (np.median(imp) * 100, np.percentile(imp, 90) * 100, imp.max() * 100))
print("relative fitness gain 0.2g -> 0.9g: median %.2f%%, 90th pct %.1f%%" % (np.median(imp100) * 100, np.percentile(imp100, 90) * 100))
