#!/usr/bin/env python3
"""Knowledge files of IRMA with the TRUE structure, to ask whether GREMA gets the sign right when the structure is known (D2 of the sign diagnosis).

usage: irma_oracle_knowledge.py OUT_DIR [IRMA_DIR]
The structure is the gold standard of IRMA (8 regulations; input/IRMA/IRMA_GoldStandard.tsv), the signs of 7 of them are known (input/IRMA/IRMA_signs.tsv); G2 -> G4 has no known sign.
Writes, in OUT_DIR:
  know_unknown_signs.txt   the 8 true regulations are '?', every other pair is 0 (structure known, no sign known)
  know_loo_<TF>_<GENE>.txt for each of the 7 signed edges: that edge is '?', the other 6 signed edges carry their true sign (+ or -), G2 -> G4 is '?', every other pair is 0
                            (everything else is told to GREMA; does it find the sign of the one edge that is left open?)
This uses the gold standard on purpose: it is a diagnosis of why the signs are wrong, not a method.
"""
import sys, os

out = sys.argv[1]
IR = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "input", "IRMA")
os.makedirs(out, exist_ok=True)
signs = {}
for l in open(f"{IR}/IRMA_signs.tsv"):
    f = l.split(); signs[(f[0], f[1])] = int(f[2])
gold = [tuple(l.split("\t")[:2]) for l in open(f"{IR}/IRMA_GoldStandard.tsv") if l.rstrip("\n").endswith("\t1")]
genes = ["G%d" % i for i in range(1, 6)]


def write(path, value):
    with open(path, "w") as f:
        f.write("TF\tGENE\tREGULATORY\n")
        for g in genes:
            for t in genes:
                f.write(f"{t}\t{g}\t{value(t, g)}\n")


write(f"{out}/know_unknown_signs.txt", lambda t, g: "?" if (t, g) in gold else "0")
for (tt, gg) in signs:
    def v(t, g, tt=tt, gg=gg):
        if (t, g) not in gold:
            return "0"
        if (t, g) == (tt, gg) or (t, g) not in signs:
            return "?"
        return "+" if signs[(t, g)] > 0 else "-"
    write(f"{out}/know_loo_{tt}_{gg}.txt", v)
print("written", out)
