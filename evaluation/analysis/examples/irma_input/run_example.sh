#!/bin/bash
# IRMA switch-on and switch-off as two series with the medium as external input: the data, the check, GREMA with and without the input, a forecast.
# usage (from anywhere): bash evaluation/analysis/examples/irma_input/run_example.sh OUT_DIR [THREADS]
# Needs the solver (cd EMA_HFODE; make). Takes about 20-30 minutes with 4 threads.
set -eu
R=$(cd "$(dirname "$0")/../../../.." && pwd)
E=$R/evaluation/analysis
O=${1:?usage: run_example.sh OUT_DIR [THREADS]}; T=${2:-4}
mkdir -p "$O"; O=$(cd "$O" && pwd)
X=$E/examples/irma_input

echo "== 1. the profile: switch-on and switch-off (the first 16 points) as two series, and the input U1"
python3 $X/make_joint.py $R/input/IRMA $O/joint16.txt 16
python3 $E/scripts/make_input.py $O/joint16.txt $X/table.tsv U1 --out $O/joint16_input.txt
python3 $X/make_joint.py $R/input/IRMA $O/joint12.txt 12
python3 $E/scripts/make_input.py $O/joint12.txt $X/table.tsv U1 --out $O/joint12_input.txt

echo "== 2. the data check"
python3 $E/scripts/check_data.py $O/joint16_input.txt --input U1 | sed -n 1,12p

echo "== 3. the knowledge files (no prior knowledge; with the input: U1 may regulate every gene)"
cp $R/input/IRMA/IRMA_no_prior_knowledge.txt $O/know_no.txt
cp $O/know_no.txt $O/know_in.txt
for g in 1 2 3 4 5; do printf 'U1\tG%s\t?\n' $g >> $O/know_in.txt; done

echo "== 4. GREMA: 5 runs with the input and 5 without, on 12 points (for the forecast) and on 16 points"
export GREMA_ROOT=$R DONELOG=$O/done.log
{
for k in 1 2 3 4 5; do
  echo "in12_r$k $O/joint12_input.txt $O/know_in.txt -input U1"
  echo "no12_r$k $O/joint12.txt $O/know_no.txt"
  echo "in16_r$k $O/joint16_input.txt $O/know_in.txt -input U1"
  echo "no16_r$k $O/joint16.txt $O/know_no.txt"
done
} > $O/jobs.txt
mkdir -p $O/runs
cat $O/jobs.txt | xargs -P 4 -L 1 bash -c "$E/scripts/run_one.sh $O/runs/\$0 \$1 \$2 -g 500 -f 1 -t $T \$3 \$4" > /dev/null 2>&1
echo "runs finished: $(grep -c 'exit=0' $O/done.log) of 20 ok"

echo "== 5. the effect of the input on the genes (16 points, vote score = fraction of the runs that use it; mean over the 5 runs)"
python3 - "$O" <<'E'
import sys, glob, collections
O = sys.argv[1]
v = collections.defaultdict(list); s = collections.defaultdict(collections.Counter)
for f in glob.glob(f"{O}/runs/in16_r*/out/final_results.txt"):
    for l in open(f):
        p = l.rstrip("\n").split("\t")
        if p[0] == "U1":
            v[p[1]].append(float(p[4])); s[p[1]][p[2]] += 1
for g in sorted(v):
    print(f"  U1 -> {g}: vote {sum(v[g]) / len(v[g]):.2f}   sign in the 5 runs: {dict(s[g])}")
E

echo "== 6. forecast of points 13-16 from a model fitted on points 1-12"
python3 $X/forecast.py $O/joint16_input.txt $O/runs 12
