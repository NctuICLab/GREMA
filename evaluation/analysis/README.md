# Analyses of GREMA: reports, scripts and archived results

Everything about the evaluation of GREMA beyond the DREAM4 scoring script (`evaluation/dream4_evaluate.py`) is kept here: the reports, the scripts that make and score the data, and the archived results.
All paths in the reports are relative to the root of the repository.

```
evaluation/
  dream4_evaluate.py            scoring of a result against a DREAM4-format gold standard
  gold_standards/, gold_standards_undirected/
  analysis/
    README.md                   this file
    *_REPORT.md, HANDOFF.md     the reports (table below)
    scripts/                    scripts of the analyses (table below)
    results/                    archived results (*.tgz) and results/README.md (what is in each archive)
```

## Reports

| Report | Question and answer in one line | Archive in `results/` |
|---|---|---|
| [HANDOFF.md](HANDOFF.md) | The hand-off file the analyses started from (tasks T1-T4, rules, first results). Historical. | `final_results.tgz`, `irma_gp_inputs.tgz` |
| [SOS_REPORT.md](SOS_REPORT.md) | SOS (4 experiments): `-f 1` is better than `-f 3`, 500 generations are enough, the signs of the LexA edges are wrong, a simulation with a hidden LexA protein reproduces this, a restricted number of regulator slots trades TP against FP. | `sos_results.tgz` |
| [IRMA_REPORT.md](IRMA_REPORT.md) | IRMA: the built-in perturbation sets (`-perturb`) do not help; 500 / 1000 / 2000 generations give the same result. | `irma_t2_t3_results.tgz` |
| [HACKETT_REPORT.md](HACKETT_REPORT.md) | TF-induction time courses of Hackett et al. (a microarray) with a gold standard of direct regulations: design, 5 min vs 15 min time grid, 10 / 3 / 2 regulator slots per gene. | `hackett_results.tgz` |
| [SSYSTEM_SLOTS_REPORT.md](SSYSTEM_SLOTS_REPORT.md) | Default number of regulator slots (HFODE and S-system), 1 / 2 / 4 slots on IRMA, and the S-system output is not reversed. | `ssystem_slots_results.tgz` |
| [METHODS_REPORT.md](METHODS_REPORT.md) | GREMA against dynGENIE3; why SOS cannot rank methods (permutation control); edge scores, filters before GREMA, a second TF subset, real time points (`-time`). | `methods_comparison_results.tgz`, `hackett_set2_results.tgz`, `real_time_points_results.tgz` |
| [SIGN_SIMULATION_REPORT.md](SIGN_SIMULATION_REPORT.md) | Signs of GREMA on simulated signed networks: right (92-100%) with the model right and with seven mild misspecifications. | `sign_simulation_results.tgz` |
| [PROTEIN_REGULATOR_REPORT.md](PROTEIN_REGULATOR_REPORT.md) | Promoter activity measured, protein regulators, a shared driver: collinearity lowers the detection; an opposite direction (cleavage) makes the signs wrong. | `sim_protein_results.tgz` |
| [SIGNS_REAL_DATA_REPORT.md](SIGNS_REAL_DATA_REPORT.md) | On the real IRMA data only 25-50% of the signs are right; not caused by the amount of data, a delay, an input or a protein-level inhibition, and a linear model is as wrong; a check of the data before the run. | `irma_sign_factors_results.tgz` |

## Main findings (details and limits in the reports)

- A high score on data with **one shared stimulus** (SOS) does not come from the regulation (the experiment-permutation control leaves it unchanged, for GREMA and for dynGENIE3); use data with **different perturbations** (Hackett) and run the control (`permute_replicates.py`).
- Fewer regulator slots per gene (`GREMA_SLOT_DIV`) give fewer false regulations; larger Hill-coefficient bounds (`-nmax`), thresholds (`-min`), a filtered regulator and real time points (`-time`) did not improve the ranking.
- On the data with different perturbations the fast methods (dynGENIE3, a Hill-function dictionary with sparse regression) are as good as or better than GREMA and need minutes instead of an hour; no setting was found in which GREMA is better.
- The signs of GREMA are right in simulations, but not on real IRMA data (25-50%); the cause is not known (not specific to the Hill function). `scripts/check_data.py` checks an input file before a run.
- Everything rests on few data sets (IRMA, SOS, Hackett, simulations), few runs per setting and few true edges; differences below about 0.05 AUROC are noise. The rules of every comparison were fixed before the results were looked at.

## Scripts (`scripts/`; each has a usage line in its docstring)

| Group | Scripts |
|---|---|
| Run and score GREMA | `run_one.sh` (one run in its own directory), `eval_runs.py` (AUROC / AUPRC of groups of runs), `iteration_table.py`, `convergence.py`, `eval_prefilter.py`, `combine_scores.py`, `postprocess_scores.py`, `contrast_scores.py` |
| Make or change the input data | `hackett_explore.py`, `hackett_prep.py` (Hackett data and the gold standard of Ma et al.), `make_sos_p_inputs.py`, `permute_replicates.py` (the permutation control), `augment_filtered.py`, `unfilter_results.py` |
| Simulations | `sim_gen.py` (IRMA topology), `sim_sos.py`, `sim_sos_p.py` (SOS with a hidden / measured LexA protein), `sim_signed_net.py` (random signed networks and variants), `sim_protein_net.py` (protein regulators, shared driver), `sim_irma_like.py` (IRMA-like with three switches) |
| Other methods and filters | `dyngenie3_runner.py` (needs the dynGENIE3 code, see the report), `lasso_prefilter.py`, `hsic_prefilter.py`, `hill_dictionary.py`, `gp_auto.py`, `gp_v2.py` (Gaussian-process augmentation) |
| Signs and checks of the data | `eval_signs.py`, `eval_sim_variants.py`, `eval_protein_sims.py`, `eval_irma_like.py`, `eval_signs_real.py`, `linear_signs_irma.py`, `check_data.py` |

Python scripts need numpy, scipy and scikit-learn; `hsic_prefilter.py` also pyHSICLasso; `dyngenie3_runner.py` needs the dynGENIE3 code and was run in a pixi environment (see `METHODS_REPORT.md`).

## Archived results (`results/`)

Each `*.tgz` holds the `final_results.txt` of the runs (columns `TF GENE REGULATORY CONFIDENCE_LEVEL VOTE_SCORE`) and the input files that are needed to score them, not the full run directories (10-30 MB each).
`results/README.md` lists the content of every archive. To score one:

```shell
mkdir -p /tmp/res && tar xzf evaluation/analysis/results/final_results.tgz -C /tmp/res
python3 evaluation/analysis/scripts/eval_runs.py --gold input/IRMA/IRMA_GoldStandard.tsv \
    "f1-on=/tmp/res/irma_fixed/f1_on_r*/final_results.txt"
```

The archives of the Hackett data hold the inputs made by `hackett_prep.py`; the raw data (GEO GSE142864, the gold standard of Ma et al. 2014) are public and not in the repository, the download commands are in `HACKETT_REPORT.md`.
