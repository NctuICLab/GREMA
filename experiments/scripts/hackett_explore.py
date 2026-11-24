#!/usr/bin/env python3
"""Explore the Hackett et al. (2020, GSE142864) TF-induction time courses against the Ma et al. (2014) gold-standard networks.

usage: python3 -I hackett_explore.py DATA_DIR GOLD_XLSX [network 1|2|3]
  DATA_DIR  holds GPL16244.txt.gz (series matrix) and GPL16244_annot.txt (platform table)
  GOLD_XLSX is Ma et al. 2014, PLoS ONE, Table S? "Gold_standard_network1..3" (pone.0106479.s020.xlsx)
Prints, for every induced TF that is a regulator in the gold standard: the number of responding genes (max |value| over time > cutoff)
and of gold targets, and their overlap. Only the standard library is used (the data are untrusted downloads).
"""
import gzip, re, sys, zipfile, collections
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
TIMES = (0, 5, 10, 15, 20, 30, 45, 90)


def read_gold(xlsx, network):
    z = zipfile.ZipFile(xlsx)
    ss = []
    for ev, el in ET.iterparse(z.open("xl/sharedStrings.xml")):
        if el.tag.endswith("}si"):
            ss.append("".join(t.text or "" for t in el.iter() if t.tag.endswith("}t")))
            el.clear()
    edges = []
    for ev, el in ET.iterparse(z.open("xl/worksheets/sheet%d.xml" % network)):
        if el.tag.endswith("}row"):
            vals = []
            for c in el:
                v = c.find("m:v", NS)
                t = v.text if v is not None else ""
                vals.append(ss[int(t)] if c.get("t") == "s" and t else t)
            if len(vals) >= 2 and vals[0] != "Transcription factor":
                edges.append((vals[0], vals[1]))
            el.clear()
    return edges


def read_platform(path):
    probe2orf, sym2orf = {}, {}
    inside = False
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("!platform_table_begin"):
            inside = True
            next(open(path)) if False else None
            continue
        if line.startswith("!platform_table_end"):
            break
        if inside:
            f = line.rstrip("\n").split("\t")
            if f[0] == "ID" or len(f) < 5:
                continue
            if f[2] and f[3].upper() == "FALSE":   # real probes (the others are controls)
                probe2orf[f[0]] = f[2]
                if f[4]:
                    sym2orf.setdefault(f[4].upper(), f[2])
    return probe2orf, sym2orf


def read_matrix(path):
    titles, cols = [], None
    rows = {}
    with gzip.open(path, "rt") as f:
        inmat = False
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [t.strip('"') for t in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                inmat = True
                next(f)
            elif line.startswith("!series_matrix_table_end"):
                break
            elif inmat:
                p = line.rstrip("\n").split("\t")
                vals = []
                for x in p[1:]:
                    try:
                        vals.append(float(x))
                    except ValueError:
                        vals.append(float("nan"))
                rows[p[0].strip('"')] = vals
    return titles, rows


def experiments(titles):
    """{(strain, TF): {time: column index}} for titles like SMY257_GCN4_t15 (replicate t0 columns are skipped)."""
    exp = collections.defaultdict(dict)
    for i, t in enumerate(titles):
        m = re.match(r"^([A-Za-z0-9]+)_([^_]+)_t(\d+)$", t)
        if m and int(m.group(3)) not in exp[(m.group(1), m.group(2))]:
            exp[(m.group(1), m.group(2))][int(m.group(3))] = i
    return {k: v for k, v in exp.items() if tuple(sorted(v)) == TIMES}


if __name__ == "__main__":
    data_dir, gold_xlsx = sys.argv[1], sys.argv[2]
    network = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    cutoff = 1.0
    gold = read_gold(gold_xlsx, network)
    gold_tfs = collections.defaultdict(set)
    for tf, g in gold:
        gold_tfs[tf].add(g)
    print("gold network %d: %d edges, %d regulators" % (network, len(gold), len(gold_tfs)))
    probe2orf, sym2orf = read_platform(data_dir + "/GPL16244_annot.txt")
    titles, rows = read_matrix(data_dir + "/GPL16244.txt.gz")
    exps = experiments(titles)
    print("platform: %d probes with ORF, %d symbols; experiments with the 8 standard time points: %d" % (len(probe2orf), len(sym2orf), len(exps)))
    orf_rows = collections.defaultdict(list)
    for pid, vals in rows.items():
        if pid in probe2orf:
            orf_rows[probe2orf[pid]].append(vals)
    out = []
    for (strain, tf), cols in sorted(exps.items()):
        orf = sym2orf.get(tf.upper())
        if orf is None or orf not in gold_tfs:
            continue
        idx = [cols[t] for t in TIMES]
        resp = set()
        for g, vl in orf_rows.items():
            best = 0.0
            for v in vl:
                xs = [v[i] for i in idx]
                if all(x == x for x in xs):
                    best = max(best, max(abs(x) for x in xs))
            if best > cutoff:
                resp.add(g)
        tg = gold_tfs[orf]
        out.append((len(resp), tf, orf, len(tg), len(resp & tg)))
    print("\ninduced TFs that are gold-standard regulators: %d (cutoff |value| > %.1f)" % (len(out), cutoff))
    print("%-8s %-9s %9s %8s %8s" % ("TF", "ORF", "responding", "targets", "overlap"))
    for n, tf, orf, nt, ov in sorted(out, reverse=True)[:25]:
        print("%-8s %-9s %9d %8d %8d" % (tf, orf, n, nt, ov))
