#!/bin/bash
# Without -input the files that GREMA_main.pl generates must be identical in two versions of GREMA (the runs themselves are random, so only the files are compared).
# usage: check_input_regression.sh OLD_GREMA_DIR NEW_GREMA_DIR [PROFILE.txt KNOWLEDGE.txt]   (default: IRMA switch-on, no prior knowledge, of NEW_GREMA_DIR)
# Compared: tmp/<profile>_config, tmp/<knowledge>_ForStep1.txt and the mapping of run.log (first 7 lines). Needs both solvers built (make in EMA_HFODE).
set -u
OLD=${1:?usage: check_input_regression.sh OLD_GREMA_DIR NEW_GREMA_DIR [PROFILE KNOWLEDGE]}; NEW=${2:?}
P=${3:-$NEW/input/IRMA/IRMA_switch_on_timeseries_expression.txt}; K=${4:-$NEW/input/IRMA/IRMA_no_prior_knowledge.txt}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
for v in old new; do
  if [ $v = old ]; then R=$OLD; else R=$NEW; fi
  mkdir -p "$T/$v"
  ( cd "$T/$v" && perl "$R/GREMA_main.pl" -i "$P" -k "$K" -o "$T/$v" -g 100 -f 1 -t 4 > /dev/null 2>&1 ) || { echo "the run of $v failed"; exit 1; }
done
bad=0
for f in "$(basename "$P")_config" "$(basename "${K%.*}")_ForStep1.txt"; do
  if cmp -s "$T/old/tmp/$f" "$T/new/tmp/$f"; then echo "identical: $f"; else echo "DIFFERENT: $f"; bad=1; fi
done
if cmp -s <(head -7 "$T/old/run.log") <(head -7 "$T/new/run.log"); then echo "identical: run.log mapping"; else echo "DIFFERENT: run.log mapping"; bad=1; fi
exit $bad
