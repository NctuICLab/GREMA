#!/usr/bin/env python3
"""Make a GREMA input from the Hackett et al. (2020, GSE142864) TF-induction time courses and the Ma et al. (2014) gold standard.

usage: python3 -I hackett_prep.py DATA_DIR GOLD_XLSX OUT_DIR [--tfs 10] [--genes 30] [--network 1] [--per-tf N] [--step 5] [--skip-tfs N]
Rules (fixed before any GREMA run, no gold-standard edge is used to choose anything):
  TFs    the --tfs induced TFs that are regulators in the gold standard (network N) and have the 8 standard time points
         (0,5,10,15,20,30,45,90 min), with the most responding genes (max |value| over time > 1)
  genes  the --genes other genes with the largest mean (over the chosen experiments) of the max |value| over time
         (genes with a missing value in a chosen experiment are skipped)
  data   each chosen experiment = one replicate (repeat_number); linear interpolation to every --step min (default 5: 19 points; 15: 7 points);
         one global linear scaling of all values to [0.01, 1]
  know   only the chosen TFs can regulate (all other rows are 0), every TF -> gene pair is '?'
  gold   TF -> gene pairs of the chosen TFs (self pairs skipped) with 1 if the edge is in the gold network
Writes expression.txt, knowledge.txt, gold.tsv, names.tsv (G<i> -> ORF, symbol, TF or not) and summary.txt in OUT_DIR.
"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hackett_explore import read_gold, read_platform, read_matrix, experiments, TIMES

args = sys.argv[1:]
opt = {"--tfs": 10, "--genes": 30, "--network": 1, "--per-tf": 0, "--step": 5, "--skip-tfs": 0}
pos = []
i = 0
while i < len(args):
    if args[i] in opt:
        opt[args[i]] = int(args[i + 1]); i += 2
    else:
        pos.append(args[i]); i += 1
data_dir, gold_xlsx, out_dir = pos
os.makedirs(out_dir, exist_ok=True)
CUT = 1.0
GRID = list(range(0, 91, opt["--step"]))

gold = read_gold(gold_xlsx, opt["--network"])
gold_set = set(gold)
gold_regs = set(t for t, g in gold)
probe2orf, sym2orf = read_platform(data_dir + "/GPL16244_annot.txt")
titles, rows = read_matrix(data_dir + "/GPL16244.txt.gz")
exps = experiments(titles)
orf2sym = {}
for s, o in sym2orf.items():
    orf2sym.setdefault(o, s)
orf_rows = collections.defaultdict(list)
for pid, vals in rows.items():
    if pid in probe2orf:
        orf_rows[probe2orf[pid]].append(vals)


def profile(orf, idx):
    """mean over the probes of the gene at the 8 time points, None if any value is missing"""
    out = []
    for t in idx:
        xs = [v[t] for v in orf_rows[orf] if v[t] == v[t]]
        if not xs:
            return None
        out.append(sum(xs) / len(xs))
    return out


# 1. TFs: induced TFs that are gold regulators, ranked by the number of responding genes
cand = []
for (strain, tf), cols in sorted(exps.items()):
    orf = sym2orf.get(tf.upper())
    if orf is None or orf not in gold_regs:
        continue
    idx = [cols[t] for t in TIMES]
    n = 0
    for g in orf_rows:
        p = profile(g, idx)
        if p is not None and max(abs(x) for x in p) > CUT:
            n += 1
    cand.append((n, tf, orf, strain, idx))
cand.sort(key=lambda c: (-c[0], c[1]))
chosen = []
seen = set()
skipped = 0
for c in cand:
    if c[2] in seen:
        continue
    seen.add(c[2])
    if skipped < opt["--skip-tfs"]:   # --skip-tfs N: leave out the N best-ranked TFs (an independent second subset)
        skipped += 1
        continue
    chosen.append(c)
    if len(chosen) == opt["--tfs"]:
        break
tf_orfs = [c[2] for c in chosen]

# 2. genes: strongest mean response over the chosen experiments
score = {}
for g in orf_rows:
    if g in tf_orfs:
        continue
    ps = [profile(g, c[4]) for c in chosen]
    if any(p is None for p in ps):
        continue
    score[g] = sum(max(abs(x) for x in p) for p in ps) / len(ps)
if opt["--per-tf"]:
    # each chosen TF contributes the genes that respond most strongly in its own experiment
    genes = []
    for c in chosen:
        own = {g: max(abs(x) for x in profile(g, c[4])) for g in score}
        for g, _ in sorted(own.items(), key=lambda kv: (-kv[1], kv[0]))[: opt["--per-tf"]]:
            if g not in genes:
                genes.append(g)
else:
    genes = [g for g, s in sorted(score.items(), key=lambda kv: (-kv[1], kv[0]))[: opt["--genes"]]]
variables = tf_orfs + genes   # G1..G<n_tf> are the TFs


def interp(p):
    out = []
    for t in GRID:
        for k in range(len(TIMES) - 1):
            if TIMES[k] <= t <= TIMES[k + 1]:
                w = (t - TIMES[k]) / (TIMES[k + 1] - TIMES[k])
                out.append(p[k] * (1 - w) + p[k + 1] * w)
                break
    return out


series = []   # [replicate][variable] -> 19 values (TF orf may be missing in its own experiment -> skipped below)
for c in chosen:
    series.append([interp(profile(v, c[4])) if profile(v, c[4]) is not None else None for v in variables])
bad = sorted(set(v for s in series for v, p in zip(variables, s) if p is None))
if bad:
    sys.exit("missing profiles for %s" % bad)
allv = [x for s in series for p in s for x in p]
lo, hi = min(allv), max(allv)

with open(out_dir + "/expression.txt", "w") as f:
    f.write("repeat_number=%d\ntimepoint_number=%d\n" % (len(chosen), len(GRID)))
    for r, s in enumerate(series):
        for j, p in enumerate(s):
            f.write("rep%d\tG%d\t" % (r + 1, j + 1) + "\t".join("%.6f" % (0.01 + (x - lo) / (hi - lo) * 0.99) for x in p) + "\n")
names = ["G%d" % (j + 1) for j in range(len(variables))]
with open(out_dir + "/knowledge.txt", "w") as f:
    f.write("TF\tGENE\tREGULATORY\n")
    for a in range(len(variables)):
        for b in range(len(variables)):
            f.write("%s\t%s\t%s\n" % (names[a], names[b], "?" if a < len(tf_orfs) else "0"))
n_true = 0
with open(out_dir + "/gold.tsv", "w") as f:
    for a in range(len(tf_orfs)):
        for b in range(len(variables)):
            y = 1 if (variables[a], variables[b]) in gold_set else 0
            n_true += y if a != b else 0
            f.write("%s\t%s\t%d\n" % (names[a], names[b], y))
with open(out_dir + "/names.tsv", "w") as f:
    for j, v in enumerate(variables):
        f.write("%s\t%s\t%s\t%s\n" % (names[j], v, orf2sym.get(v, ""), "TF" if j < len(tf_orfs) else "gene"))
with open(out_dir + "/summary.txt", "w") as f:
    f.write("gold network %d: %d edges\n" % (opt["--network"], len(gold)))
    f.write("TFs (n responding genes): " + ", ".join("%s/%s (%d)" % (c[1], c[2], c[0]) for c in chosen) + "\n")
    f.write("variables: %d (%d TFs + %d genes), replicates %d, time points %d\n" % (len(variables), len(tf_orfs), len(genes), len(chosen), len(GRID)))
    f.write("evaluated pairs (TF -> variable, no self): %d, gold edges among them: %d\n" % (len(tf_orfs) * len(variables) - len(tf_orfs), n_true))
print(open(out_dir + "/summary.txt").read())
