#!/usr/bin/env python3
"""Does an ATAC footprint near a gene's TSS tell which TF -> gene pairs are in TRRUST? (no GREMA run)

usage: python3 hmsc_prior_enrichment.py [--window 10000] [--min-tpm 5]
Pairs: every TRRUST TF that is expressed (TPM >= min-tpm on some day) and has a JASPAR motif in the footprint files x every expressed protein-coding gene (self pairs skipped).
Footprint = a HINT-ATAC footprint of the motif within +-window of the TSS on any of the 5 days (same rule as hmsc_prep.py).
Prints: pairs, TRRUST pairs, the fraction of pairs and of TRRUST pairs with a footprint, the odds ratio, the AUROC of the footprint (0/1) per TF (mean over TFs with
>= 5 TRRUST targets) and for all pairs, and the same when the gene labels are shuffled within each TF (control, 20 shuffles).
"""
import sys, random, bisect, collections, tempfile
sys.argv = [sys.argv[0], tempfile.mkdtemp()] + sys.argv[1:]  # reuse the loaders of hmsc_prep (its part before the design) with a scratch output directory
src = open(__file__.replace("hmsc_prior_enrichment.py", "hmsc_prep.py")).read()
head = src[: src.index('if opt["--design"] == "runx2"')]
exec(compile(head, "hmsc_prep_head", "exec"))
from itertools import islice
allsym = {v[2]: g for g, v in expr.items()}
tfs = [s for s in trrust if s in allsym]
idx = load_footprints(set(tfs))
tfs = [s for s in tfs if s in idx]
genes = [(g, v) for g, v in expr.items()]
W = opt["--window"]
rows = []
for s in tfs:
    for g, v in genes:
        if v[2] == s:
            continue
        rows.append((s, v[2], has_fp(idx, s, v[0], v[1], W), v[2] in trrust[s]))
n = len(rows); pos = sum(r[3] for r in rows); fp = sum(r[2] for r in rows)
pf = sum(1 for r in rows if r[3] and r[2])
print(f"TFs {len(tfs)}, pairs {n}, TRRUST pairs {pos}, pairs with a footprint {fp} ({100.0*fp/n:.1f}%), TRRUST pairs with a footprint {pf} ({100.0*pf/pos:.1f}%)")
a, b, c, d = pf, pos - pf, fp - pf, n - pos - (fp - pf)
print(f"odds ratio {(a * d) / (b * c):.2f}")


def auroc(items):
    p = [x for x in items if x[1]]; q = [x for x in items if not x[1]]
    if not p or not q:
        return None
    # binary score: AUROC = P(score_pos > score_neg) + 0.5 P(equal)
    pp = sum(x[0] for x in p) / len(p); qq = sum(x[0] for x in q) / len(q)
    return pp * (1 - qq) + 0.5 * (pp * qq + (1 - pp) * (1 - qq))


by = collections.defaultdict(list)
for s, g, f, t in rows:
    by[s].append((int(f), t))
per = [auroc(v) for v in by.values() if sum(1 for x in v if x[1]) >= 5]
print(f"per-TF AUROC of the footprint (TFs with >= 5 TRRUST targets, {len(per)} TFs): mean {sum(per)/len(per):.3f}")
print("all pairs together AUROC", round(auroc([(int(r[2]), r[3]) for r in rows]), 3))
random.seed(1)
perm = []
for _ in range(20):
    ps = []
    for s, v in by.items():
        if sum(1 for x in v if x[1]) >= 5:
            lab = [x[1] for x in v]; random.shuffle(lab)
            ps.append(auroc([(x[0], l) for x, l in zip(v, lab)]))
    perm.append(sum(ps) / len(ps))
print(f"shuffled gene labels within each TF: mean per-TF AUROC {sum(perm)/20:.3f} (range {min(perm):.3f}-{max(perm):.3f})")
