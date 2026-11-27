#!/usr/bin/env python3
"""Can augmented data take the place of real series? Data sets for the comparison (simulation with a known truth).

usage: sim_aug_prep.py OUT_DIR SEED
Simulates a random signed network of 10 genes with 5 series of 20 points (sim_signed_net.py, 3% noise, the model of GREMA is the generating model) and writes in OUT_DIR:
  s1.txt, s2.txt, s3.txt, s4.txt   the first 1, 2, 3, 4 series (the "real" data); the series 5 is never used for the fit
  heldout.txt                      series 5 (a new initial state of the same system), for the prediction test
  know.txt, gold.tsv, signs.tsv    knowledge (every pair '?', diagonal 0), the 10 x 9 pairs with 0/1, the signs of the true edges
Conditions that are made from s1 (see sim_aug_launch in the report): s1 with `-perturb K -noise 0.03` (built-in augmentation, GREMA_main.pl) and GP samples of s1 (gp_v2.py).
"""
import sys, os, subprocess, shutil

out, seed = sys.argv[1], sys.argv[2]
here = os.path.dirname(os.path.abspath(__file__))
os.makedirs(out, exist_ok=True)
subprocess.check_call([sys.executable, f"{here}/sim_signed_net.py", out, seed, "--genes", "10", "--series", "5", "--points", "20", "--noise", "0.03"], stdout=subprocess.DEVNULL)
lines = open(f"{out}/expression.txt").read().rstrip("\n").split("\n")
T = lines[1].split("=")[1]
rows = lines[2:]
byrep = {}
for l in rows:
    byrep.setdefault(int(l.split("\t")[0][3:]), []).append(l)


def write(path, reps):
    with open(path, "w") as f:
        f.write(f"repeat_number={len(reps)}\ntimepoint_number={T}\n")
        for new, r in enumerate(reps, 1):
            for l in byrep[r]:
                p = l.split("\t")
                f.write("\t".join([f"rep{new}"] + p[1:]) + "\n")


for k in (1, 2, 3, 4):
    write(f"{out}/s{k}.txt", list(range(1, k + 1)))
write(f"{out}/heldout.txt", [5])
os.rename(f"{out}/knowledge.txt", f"{out}/know.txt")
os.remove(f"{out}/expression.txt")
print("written", out)
