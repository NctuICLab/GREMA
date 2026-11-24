import sys
# usage: make_sos_p_inputs.py OUTDIR   -> know9.txt (9x9 all '?') and gold9.tsv (G9 = P -> G1..G8 is true)
d = sys.argv[1]
g = ["G%d" % i for i in range(1, 10)]
open(d + "/know9.txt", "w").write("TF\tGENE\tREGULATORY\n" + "".join("%s\t%s\t?\n" % (a, b) for a in g for b in g))
open(d + "/gold9.tsv", "w").write("".join("%s\t%s\t%d\n" % (a, b, 1 if a == "G9" and b != "G9" else 0) for a in g for b in g))
