# Hand-off: GREMA analyses to continue on a server

This file is written for an agent (or person) that continues the work on a faster machine. Everything needed is in the
repository: code, data, scripts and the results of the finished runs. The user's language is Chinese; report to them in Chinese.

## 1. Goal and context

GREMA infers a gene regulatory network from time-series data with an ODE model (HFODE = Hill-function ODE, or S-system)
fitted by an evolutionary algorithm (30 GA runs per gene, voting over the best 5, iterations that fix the regulations that
agree). Output `final_results.txt`: `TF GENE REGULATORY CONFIDENCE_LEVEL VOTE_SCORE`. See `README.md`.

The user's question is how accurate GREMA is on **real data with few time series** (IRMA, SOS) and how the accuracy can be
improved: better fitness, more data (augmentation), more replicates. Accuracy = AUROC / AUPRC of the edge ranking against a gold standard.

Repository state: `master` contains PR #18 (a code review fixed ~30 bugs in the Perl driver and both C++ solvers; results
produced **before** it, e.g. `-f 3`, are not comparable). Branch `handoff/analysis` adds: the fix that makes the solver use the
perturbation sets (`-perturb`, `-noise`), the SOS data (`input/SOS/`), and `experiments/` (this file, scripts, archived results).

## 2. Setup

```shell
git clone https://github.com/NctuICLab/GREMA.git && cd GREMA && git checkout handoff/analysis
cd EMA_HFODE && make && cd ..           # HFODE solver (the one that was evaluated)
cd EMA_SSYSTEM && make && cd ..         # S-system solver (only smoke-tested, see section 6)
python3 -m pip install numpy scipy scikit-learn
perl -V:useithreads                     # must print useithreads='define'
tar xzf experiments/results/final_results.tgz -C /some/dir    # archived results, see experiments/results/README.md
```

Scripts (`experiments/scripts/`): `run_one.sh` (one GREMA run in an isolated directory), `eval_runs.py` (AUROC/AUPRC of
groups of runs; `--by` for seed-wise means), `convergence.py` (GA convergence from the progress files), `gp_auto.py`,
`gp_v2.py` (GP augmentation: sample trajectories from the data), `sim_gen.py` (simulated data with the IRMA topology).

Practical notes
- One run uses `-t` processes (one per gene). Keep `parallel runs x t <= cores`, otherwise everything slows down (the laptop was
  at load 100-170 on 8 cores).
- `final_results.txt` is written at the very end, so its existence means the run finished. GREMA runs inside `<out>/tmp`
  (config, progress files `progress<gene>_step<k>.txt`, `calprofile*`).
- The output directory must exist and must not contain a `final_results.txt`; `run_one.sh` handles this.
- Timing on the (loaded) laptop: IRMA single series `-g 500 -t 5` about 20 s; SOS 4 experiments 30-60 min per run.

## 3. Rules the user set (do not break them)

- Design choices must **not use the gold standard** ("不能偷看答案的做法"). Fix settings first, then evaluate.
- The target is **real data with few replicates** (IRMA, SOS), not DREAM4 simulations.
- Keep IRMA switch-on and switch-off **separate** (16 and 21 time points, one replicate each); no merging of time points.
- Do **not** generate augmented data from the HFODE model (it made results worse).
- Keep the first time point of IRMA (removing it made the results worse).
- Differences below about 0.05 AUROC are noise (8 true edges; sd over runs 0.01-0.04, over GP seeds 0.04-0.09). Say so when reporting.
- Do not commit or open PRs without being asked; do not merge.

## 4. Results so far (HFODE, `-g 500`, after PR #18, mean of runs; AUROC / AUPRC)

IRMA: 20 non-self directed pairs, 8 true (random 0.5 / 0.4). `conf` = CONFIDENCE_LEVEL, `vote` = VOTE_SCORE.

| Setting (10 runs) | switch-on conf | switch-on vote | switch-off conf | switch-off vote |
|---|---|---|---|---|
| `-f 1` (default) | 0.705 / 0.596 | 0.781 / 0.754 | 0.683 / 0.572 | 0.605 / 0.544 |
| `-f 3` | 0.779 / 0.732 | 0.762 / 0.667 | 0.708 / 0.660 | 0.609 / 0.574 |

- Before the fix `-f 3` was at random level (AUROC 0.50) because the derivative term read an uninitialised value at the last time point.
- `-c 1` (5 runs) did not help: conf switch-on 0.646 (`-f 1`) / 0.675 (`-f 3`), switch-off 0.663 / 0.708.
- GP augmentation (5 seeds x 10 trajectories, score = mean over the 50 trajectories): helped switch-on only for the v2 latent curve with
  `-f 1` and vote (0.823 / 0.834); first GP version and `-f 3` not better than the original data. Switch-off: every GP setting
  was worse than the original data (AUROC 0.50-0.60 vs 0.68-0.71).
- GA convergence on IRMA (progress files): 85% of the GA runs still improved in the last 100 generations, but the median relative
  gain from generation 250 to 450 is only 0.4-1.5% (90th percentile 3-13%). Earlier (5 runs) `-g 500` and `-g 1000` gave the same
  TP/FP on IRMA. Not yet checked: `-g 2000`, other data sets.

SOS (see `input/SOS/README.md`): gold standard = LexA represses all 8 genes incl. itself (8 true of 64 pairs; 7 of 56 without self).
Paper (Table 5): TP 7, TN 42, FP 14, FN 1 (SEN 0.88, SPE 0.75, ACC 0.77); 4 GREMA iterations.

| Setting | pairs | runs | TP / FP | binary | conf | vote |
|---|---|---|---|---|---|---|
| 1 experiment, diagonal `0`, `-f 1` | 56 | 5 | 2.2 / 14.6 | 0.508 / 0.128 | 0.509 / 0.128 | 0.607 / 0.182 |
| 1 experiment, diagonal `0`, `-f 3` | 56 | 5 | 2.8 / 15.0 | 0.547 / 0.139 | 0.519 / 0.128 | 0.557 / 0.195 |
| 1 experiment, diagonal `?`, `-f 1` | 64 | 5 | 4.2 / 21.4 | 0.571 / 0.147 | 0.578 / 0.154 | 0.542 / 0.155 |
| 4 experiments as replicates, diagonal `?`, `-f 1` | 64 | 3 (partial) | 5.7 / 21.7 | 0.661 / 0.184 | 0.680 / 0.205 | 0.754 / 0.290 |

- With 1 experiment the result is at random level; the true edges found all have sign `+` although LexA is a repressor.
- With 4 experiments it is clearly better (preliminary, 3 runs); LexA -> LexA comes out as `-` (correct); still far from the paper (TP 7 / FP 14).
- The paper's data are "4 experiments"; how it fed them to GREMA, and its preprocessing, are not stated (we use `log2(x+1)`, one global linear scaling, 4 replicates).

## 5. Tasks to run (in this order of priority)

Define once: `R=$PWD; OUT=$HOME/grema_runs; mkdir -p $OUT` and use `bash`. A job list is piped to `xargs`; each line is `outdir fitness [more]`.

**T1. SOS with 4 experiments (paper protocol), `-f 1` and `-f 3`, 5 runs each** (3 `-f 1` runs finished on the laptop; rerun all for a clean comparison).
```shell
(for f in 1 3; do for r in 1 2 3 4 5; do echo "$OUT/sos4/f${f}_r$r $f"; done; done) | xargs -P 3 -L1 sh -c \
  "$R/experiments/scripts/run_one.sh \"\$0\" $R/input/SOS/SOS_4exp_timeseries_expression.txt $R/input/SOS/SOS_all_unknown_knowledge.txt -g 500 -f \"\$1\" -t 4"
python3 experiments/scripts/eval_runs.py --gold input/SOS/SOS_GoldStandard.tsv --self-pairs \
   "f1=$OUT/sos4/f1_r*/out/final_results.txt" "f3=$OUT/sos4/f3_r*/out/final_results.txt"      # 64 pairs, 8 true
python3 experiments/scripts/eval_runs.py --gold input/SOS/SOS_GoldStandard.tsv \
   "f1=$OUT/sos4/f1_r*/out/final_results.txt" "f3=$OUT/sos4/f3_r*/out/final_results.txt"      # 56 pairs, 7 true
```
Report: AUROC/AUPRC, TP/FP vs the paper (TP 7, FP 14), signs of the true edges, number of GREMA iterations (`ls <run>/out/tmp | grep ForStep`).
Optional: `-t 8` when the server has the cores; 2 distinct UV doses are in the data (5 and 20 J/m^2) - could also be tried separately.

**T2. IRMA with the built-in perturbation sets** (`-perturb 9 -noise 0.05` and `0.1`, `-f 1` and `-f 3`, switch-on/off, 5 runs each = 40 runs, 6 finished).
The solver adds noise (standard deviation = noise x value, independent per time point and gene, also on the initial value) to
copies of the data. Before the fix on this branch the sets were generated but ignored by the fitness.
```shell
(for n in 0.05 0.1; do for f in 1 3; do for sw in on off; do for r in 1 2 3 4 5; do echo "$OUT/pt/n${n}_f${f}_${sw}_r$r $n $f $sw"; done; done; done; done) | xargs -P 4 -L1 sh -c \
  "$R/experiments/scripts/run_one.sh \"\$0\" $R/input/IRMA/IRMA_switch_\$3_timeseries_expression.txt $R/input/IRMA/IRMA_no_prior_knowledge.txt -g 500 -f \"\$2\" -t 5 -perturb 9 -noise \"\$1\""
python3 experiments/scripts/eval_runs.py --gold input/IRMA/IRMA_GoldStandard.tsv "pt0.05-f1-on=$OUT/pt/n0.05_f1_on_r*/out/final_results.txt" \
   "pt0.1-f3-off=$OUT/pt/n0.1_f3_off_r*/out/final_results.txt"      # one LABEL=GLOB per setting
```
Compare with the baseline (`irma_fixed`: `f1_*`, `f3_*`) and with GP (`irma_gp`). Hypothesis to test, not a claim: independent noise per time point adds no
information and is less plausible than GP trajectories, so a gain is not expected.

**T3. GA convergence, `-g 500` vs `-g 1000` vs `-g 2000`** on IRMA (switch-on/off, `-f 1` and `-f 3`, 5-10 runs). Compare AUROC/AUPRC and the
TP/FP of the voted network (`eval_runs.py`), and the fitness convergence with `convergence.py` (assumes `-g 500`; for other `-g` the checkpoints are every `g/10` generations, adapt the constants).
Question: is the result stable between 500 and 2000 generations? Do the same check on SOS 4 experiments if T1 shows that runs are still improving.

**T4 (optional).** More runs (20) of `-f 1` vs `-f 3` on IRMA to decide whether `-f 3` is really better (the switch-off difference of 0.025 AUROC is within noise);
S-system accuracy on IRMA (`-m s-system`, never evaluated, see below); GP augmentation on SOS 4 experiments (one GP per experiment).

## 6. Known limits and open points

- S-system: fixed in PR #18 but run end-to-end only once (IRMA, `-g 100`); accuracy not evaluated.
- Gold standards: IRMA (8 regulations, from the literature/paper), SOS (from Figure 5 of the GREMA paper; LexA star incl. self-regulation). Signs are not used in the evaluation.
- The GP results used 5 seeds x 10 trajectories; the seed-to-seed sd is large. The "best" GP setting was chosen after seeing the results (selection bias).
- Perturbation and GP data contain no new information about the network (only the uncertainty of the existing series).
- `manuscript/` in the user's working copy (paper PDF, supplement) is not in the repository on purpose.
