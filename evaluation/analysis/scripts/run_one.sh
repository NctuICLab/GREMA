#!/bin/bash
# Run GREMA once in an isolated directory.
# usage: run_one.sh OUTDIR INPUT.txt KNOWLEDGE.txt [GREMA_main.pl options...]
#   e.g. run_one.sh runs/irma_on_f3_r1 input/IRMA/IRMA_switch_on_timeseries_expression.txt \
#            input/IRMA/IRMA_no_prior_knowledge.txt -g 500 -f 3 -t 5
# OUTDIR gets work/ (input copy), out/ (GREMA output, final_results.txt is in out/) and log.txt.
# A line "<name> exit=<code>" is appended to $DONELOG (default: OUTDIR/../done.log).
# GREMA_ROOT defaults to the repository that holds this script.
ROOT=${GREMA_ROOT:-$(cd "$(dirname "$0")/../../.." && pwd)}
out=$1; in=$2; know=$3; shift 3
rm -rf "$out"; mkdir -p "$out/work" "$out/out" || exit 1
out=$(cd "$out" && pwd)
cp "$in" "$out/work/in.txt" && cp "$know" "$out/work/know.txt" || exit 1
cd "$out/work" || exit 1
perl "$ROOT/GREMA_main.pl" -i in.txt -k know.txt -o "$out/out" "$@" > "$out/log.txt" 2>&1
rc=$?
echo "$(basename "$out") exit=$rc" >> "${DONELOG:-$out/../done.log}"
exit $rc
