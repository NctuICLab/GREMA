#!/usr/bin/env python3
"""Make GREMA inputs from the hMSC RNA-seq (RSEM, day 1/4/7/10/15) and ATAC-seq footprints (HINT-ATAC, JASPAR motifs), and a gold standard from TRRUST.

usage: python3 hmsc_prep.py OUT_DIR [--design runx2|multi] [--tfs 15] [--genes 60] [--window 10000] [--min-tpm 5] [--min-fold 2]
Rules (fixed before any GREMA run; no gold edge is used to choose a gene):
  genes     protein-coding genes of GENCODE v35 (gene location file of MAGMA) with TPM >= --min-tpm on some day; TSS = start (+) or end (-)
  data      log2(TPM + 1), one sample per day (1 replicate, 5 time points), one global linear scaling of all values to [0.01, 1]; the real times are 1,4,7,10,15 days (-time)
  footprint a TF has a footprint near a gene when a HINT-ATAC footprint of its JASPAR motif (name after the second '.', upper case, '::' splits dimers) lies within +-window of the TSS on some day
  design runx2   the only regulator is RUNX2; candidates = the genes with a RUNX2 footprint near the TSS and fold >= --min-fold, the --genes with the largest fold
  design multi   regulators = the --tfs TRRUST TFs (TFs of the gold standard) with the largest fold (>= --min-fold) that are expressed and have a motif in the footprint files;
                 genes = the --genes other genes with the largest fold (>= --min-fold) (fold = (max TPM + 1) / (min TPM + 1) over the 5 days)
  know      only the regulators can regulate; TF -> gene = '?' ('atac' knowledge: only when the TF has a footprint near the gene; 'all' knowledge: every pair)
  gold      TRRUST (v2, human, Han et al. 2018) TF -> gene pairs of the candidate pairs (self pairs skipped), any sign; the sign (Activation / Repression) is kept in signs.tsv
Writes expression.txt, knowledge_atac.txt, knowledge_all.txt, gold.tsv (pairs with 1 / 0 for the atac candidates and for the others), signs.tsv, names.tsv, summary.txt.
"""
import sys, os, csv, bisect, collections, math

A = "/mnt/10T-backup-mingju/ifar"
LOC = f"{A}/tools/magma_v109/geneloc_38/gencode.v35.annotation.gtf_protein_coding_gene.txt_magma_gene_loc.txt"
TRRUST = "/mnt/10T-backup-mingju/GREMA_analysis_runs/external/trrust/trrust_rawdata.human.tsv"
DAYS = [1, 4, 7, 10, 15]
FOOT = {1: "hMSC_Day1/hint-atac/jaspar/footprints_mpbs.bed", 4: "hMSC_Day4/hint-atac/footprints_mpbs.bed", 7: "hMSC_Day7/hint-atac/footprints_mpbs.bed",
        10: "hMSC_Day10/hint-atac/footprints_mpbs.bed", 15: "hMSC_Day15/hint-atac/jaspar/footprints_mpbs.bed"}

args = sys.argv[1:]
opt = {"--design": "multi", "--tfs": 15, "--genes": 60, "--window": 10000, "--min-tpm": 5.0, "--min-fold": 2.0}
pos = []
i = 0
while i < len(args):
    if args[i] in opt:
        opt[args[i]] = type(opt[args[i]])(args[i + 1]); i += 2
    else:
        pos.append(args[i]); i += 1
out = pos[0]
os.makedirs(out, exist_ok=True)

loc = {}
for l in open(LOC):
    f = l.split()
    if f[0] != "ENSGID":
        loc[f[0].split(".")[0]] = ("chr" + f[1], int(f[2]) if f[4] == "+" else int(f[3]), f[5])
tpm = {}
for d in DAYS:
    for r in csv.DictReader(open(f"{A}/RNA_seq/hMSC/day{d}/nf_results/star_rsem/hMSC_day{d}.genes.results"), delimiter="\t"):
        tpm.setdefault(r["gene_id"], []).append(float(r["TPM"]))
expr = {}
for g, (c, tss, sym) in loc.items():
    t = tpm.get(g)
    if t and max(t) >= opt["--min-tpm"] and sym not in {v[2] for v in expr.values()}:
        expr[g] = (c, tss, sym, t, (max(t) + 1) / (min(t) + 1))
trrust = collections.defaultdict(dict)
for l in open(TRRUST):
    f = l.rstrip("\n").split("\t")
    if f[2] in ("Activation", "Repression") or f[1] not in trrust[f[0]]:
        trrust[f[0]][f[1]] = f[2] if f[2] != "Unknown" else trrust[f[0]].get(f[1], "Unknown")

# footprint index: motif TF symbol -> chrom -> sorted starts
wanted_tfs = None


def load_footprints(tf_set):
    idx = collections.defaultdict(lambda: collections.defaultdict(list))
    for d in DAYS:
        for l in open(f"{A}/ATAC_seq/hMSC/{FOOT[d]}"):
            f = l.split("\t")
            if len(f) < 4:  # Day 7 has one truncated line
                continue
            name = f[3].split(".", 2)[-1].upper()
            for tf in name.split("::"):
                if tf in tf_set:
                    idx[tf][f[0]].append(int(f[1]))
    for tf in idx:
        for c in idx[tf]:
            idx[tf][c].sort()
    return idx


def has_fp(idx, tf, chrom, tss, w):
    s = idx.get(tf, {}).get(chrom)
    if not s:
        return False
    j = bisect.bisect_left(s, tss - w)
    return j < len(s) and s[j] <= tss + w


if opt["--design"] == "runx2":
    tfs = ["RUNX2"]
    idx = load_footprints(set(tfs))
    cand = [(v[4], g) for g, v in expr.items() if v[2] != "RUNX2" and has_fp(idx, "RUNX2", v[0], v[1], opt["--window"]) and v[4] >= opt["--min-fold"]]
else:
    allsym = {v[2]: g for g, v in expr.items()}
    idx0 = load_footprints({s for s in trrust if s in allsym})
    ranked = sorted([(allsym[s] and expr[allsym[s]][4], s) for s in trrust if s in allsym and s in idx0 and expr[allsym[s]][4] >= opt["--min-fold"]], reverse=True)
    tfs = [s for _, s in ranked[:opt["--tfs"]]]
    idx = {t: idx0[t] for t in tfs}
    tf_ids = {allsym[s] for s in tfs}
    cand = [(v[4], g) for g, v in expr.items() if g not in tf_ids and v[4] >= opt["--min-fold"]]
cand.sort(reverse=True)
sym2id = {v[2]: g for g, v in expr.items()}
tf_ids = [sym2id[s] for s in tfs]
gene_ids = [g for _, g in cand[:opt["--genes"]]]
vars_ = tf_ids + gene_ids
names = [expr[g][2] for g in vars_]

# expression: log2(TPM+1), global scaling to [0.01, 1]
vals = [[math.log2(x + 1) for x in expr[g][3]] for g in vars_]
lo = min(min(v) for v in vals); hi = max(max(v) for v in vals)
with open(f"{out}/expression.txt", "w") as f:
    f.write("repeat_number=1\ntimepoint_number=5\n")
    for k, v in enumerate(vals):
        f.write("rep1\tG%d\t%s\n" % (k + 1, "\t".join("%.6f" % (0.01 + 0.99 * (x - lo) / (hi - lo)) for x in v)))
n = len(vars_)
fp_pairs = set()
gold, signs = [], []
with open(f"{out}/knowledge_atac.txt", "w") as ka, open(f"{out}/knowledge_all.txt", "w") as kall, open(f"{out}/gold.tsv", "w") as fg:
    ka.write("TF\tGENE\tREGULATORY\n"); kall.write("TF\tGENE\tREGULATORY\n"); fg.write("TF\tGENE\tgold\tatac_candidate\tsign\n")
    for a in range(n):
        for b in range(n):
            tfa = a < len(tf_ids)
            fp = tfa and a != b and has_fp(idx, names[a], expr[vars_[b]][0], expr[vars_[b]][1], opt["--window"])
            kall.write("G%d\tG%d\t%s\n" % (a + 1, b + 1, "?" if tfa and a != b else "0"))
            ka.write("G%d\tG%d\t%s\n" % (a + 1, b + 1, "?" if fp else "0"))
            if tfa and a != b:
                g1 = int(names[b] in trrust.get(names[a], {}))
                fg.write("G%d\tG%d\t%d\t%d\t%s\n" % (a + 1, b + 1, g1, int(fp), trrust.get(names[a], {}).get(names[b], "") if g1 else ""))
with open(f"{out}/names.tsv", "w") as f:
    for k, s in enumerate(names):
        f.write("G%d\t%s\t%s\t%s\n" % (k + 1, s, vars_[k], "TF" if k < len(tf_ids) else "gene"))
rows = [l.rstrip("\n").split("\t") for l in open(f"{out}/gold.tsv")][1:]
pairs = len(rows); pos_ = sum(int(r[2]) for r in rows); atac = sum(int(r[3]) for r in rows)
pos_atac = sum(1 for r in rows if r[2] == "1" and r[3] == "1")
with open(f"{out}/summary.txt", "w") as f:
    f.write(f"design {opt['--design']}, window +-{opt['--window']}, regulators {len(tf_ids)}: {' '.join(tfs)}\nvariables {n} (genes {len(gene_ids)})\n")
    f.write(f"candidate pairs {pairs}, with a footprint {atac} ({100.0 * atac / max(pairs, 1):.0f}%)\ngold (TRRUST) pairs among the candidates {pos_}, of them with a footprint {pos_atac}\n")
    f.write(f"signed gold pairs: {sum(1 for r in rows if r[4] in ('Activation', 'Repression'))}\n")
print(open(f"{out}/summary.txt").read())
