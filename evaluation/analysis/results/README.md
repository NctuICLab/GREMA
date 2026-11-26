# Archived results

`final_results.tgz` holds the `final_results.txt` of every finished run (columns `TF GENE REGULATORY
CONFIDENCE_LEVEL VOTE_SCORE`), nothing else (the full run directories are 10-30 MB each). Evaluate with
`evaluation/analysis/scripts/eval_runs.py`, e.g.

```shell
mkdir -p /tmp/res && tar xzf evaluation/analysis/results/final_results.tgz -C /tmp/res
python3 evaluation/analysis/scripts/eval_runs.py --gold input/IRMA/IRMA_GoldStandard.tsv \
    "f3-on=/tmp/res/irma_fixed/f3_on_r*/final_results.txt" "f1-on=/tmp/res/irma_fixed/f1_on_r*/final_results.txt"
```

| Directory | Content | Runs |
|-----------|---------|------|
| `irma_prefix_f1`, `irma_prefix_f3` | IRMA, original data, **before** the code-review fixes (PR #18), `-f 1` / `-f 3`, 5 runs per switch | 10 + 10 |
| `irma_fixed` | IRMA, original data, after the fixes, `-g 500`: `f1_*`, `f3_*` (`on`/`off`, r1-r10) and `ss1_on_r1` (S-system smoke test, `-g 100`) | 41 |
| `irma_fixed_c1` | as above with `-c 1` (`c1f1_*`, `c1f3_*`, r1-r5) | 20 |
| `irma_gp` | IRMA, GP augmentation after the fixes: `{auto,lat}/f{1,3}/seed{0-4}/{on,off}/m{0-9}` (one run per sampled trajectory; `auto` = first automatic GP, `lat` = v2 latent curve) | 400 |
| `irma_perturb` | IRMA, `-perturb 9 -noise N` (names `n<noise>_f<fitness>_<on|off>_r<rep>`), **unfinished** | 6 |
| `sos_1exp_diag0` | SOS, 1 experiment (SPEM, = Ronen Exp3), knowledge file with the diagonal = `0` (`input/SOS/SOS_no_prior_knowledge.txt`) | 10 |
| `sos_1exp_diagU` | SOS, 1 experiment, diagonal = `?` (`SOS_all_unknown_knowledge.txt`), **unfinished** (some `-f 1` runs only) | 6 |
| `sos_4exp_diagU` | SOS, 4 experiments as 4 replicates, diagonal = `?`, **unfinished** | 3 |

`irma_gp_inputs.tgz` holds the sampled GP trajectories that were used as inputs (`auto/seed<S>/<on|off>/in_<M>.txt`,
`lat/...`); they were made with `evaluation/analysis/scripts/gp_auto.py` and `gp_v2.py` (seeds 0-4, 10 trajectories each).
Use the same input files if you want to re-run with different options.

## `sos_results.tgz` (SOS follow-up analyses, see `evaluation/analysis/SOS_REPORT.md`)

`final_results.txt` of every run, one directory per setting, `<setting>/<run>/final_results.txt`; all runs finished.

| Directory | Content | Runs |
|-----------|---------|------|
| `sos4` | SOS, 4 experiments as replicates, diagonal `?`, `-g 500`: `f1_r1-5`, `f3_r1-5` (T1; gold `input/SOS/SOS_GoldStandard.tsv`, `--self-pairs`) | 10 |
| `sos4_g1000` | as `sos4` with `-g 1000` | 10 |
| `sos4_selfrep` | as `sos4` `-f 1`, knowledge file `know_G2G2_minus.txt` (LexA -> LexA fixed to `-`; diagnostic only) | 5 |
| `sos4_slot1`, `sos4_slot4` | as `sos4` `-f 1` with `GREMA_SLOT_DIV=1` (8 regulator slots) and `=4` (2 slots) | 5 + 5 |
| `sim_sos` | simulated SOS data, hidden LexA protein (`sim<seed>.txt` made by `sim_sos.py <seed>`), `-f 1`, runs `f1_s1-5`, evaluate with the SOS gold standard | 5 |
| `sim_sos_p` | same simulation with the protein measured as G9 (`sim_sos_p.py`), 9 x 9 knowledge `know9.txt`, gold `gold9.tsv` (G9 -> G1..G8), runs `f1_s1-5`, `f3_s1-5` | 10 |

## `irma_t2_t3_results.tgz` (IRMA perturbation sets and number of generations, see `evaluation/analysis/IRMA_REPORT.md`)

`<setting>/<run>/final_results.txt`; for `gen` and `gen500` also `tmp/progress*_step1.txt` (input of `convergence.py`). All runs finished.

| Directory | Content | Runs |
|-----------|---------|------|
| `pt` | `-perturb 9 -noise N`, `-g 500`, names `n<noise>_f<fitness>_<on|off>_r<rep>`: noise 0.05 and 0.1, `-f 1` / `-f 3`, on / off, r1-r5; noise 0.1 switch-off also r6-r20 | 70 |
| `gen` | `g<G>_f<fitness>_<on|off>_r<rep>`, `-g 1000` and `-g 2000`, `-f 1` / `-f 3`, on / off, r1-r10 | 80 |
| `gen500` | the same with `-g 500` (rerun, to have the progress files) | 40 |

## `hackett_results.tgz` (TF-induction time courses of Hackett et al. 2020, see `evaluation/analysis/HACKETT_REPORT.md`)

Inputs of the three designs (`expression.txt`, `knowledge.txt`, `gold.tsv`, `names.tsv` G<i> -> ORF / symbol / TF, `summary.txt`) and `final_results.txt` of every run; all runs finished. Evaluate with `--gold <design>/gold.tsv`.

| Directory | Content | Runs |
|-----------|---------|------|
| `small` | 4 TFs, 23 variables, 5 min grid (19 points), gold network 1 (4 edges), `-f 1`: `runs` = 4 slots per gene (all TFs), `runs_slot2` = 2 slots (`GREMA_SLOT_DIV=2`) | 3 + 3 |
| `small15` | the same on the 15 min grid (7 points) | 3 + 3 |
| `net3_p10_s15` | 10 TFs, 56 variables, 15 min grid, gold network 3 (18 edges): `runs` = 10 regulator slots per gene, `runs_slot3` = 3 slots (`GREMA_SLOT_DIV=3`), `runs_slot2` = 2 slots (`GREMA_SLOT_DIV=5`) | 3 + 3 + 3 |

## `ssystem_slots_results.tgz` (regulator slots on IRMA and the direction of the S-system output, see `evaluation/analysis/SSYSTEM_SLOTS_REPORT.md`)

`<setting>/<run>/final_results.txt`; all runs finished. Evaluate with `input/IRMA/IRMA_GoldStandard.tsv`.

| Directory | Content | Runs |
|-----------|---------|------|
| `hfode_slots` | HFODE, IRMA, `-f 1`, `-g 500`: `div<d>_f1_<on|off>_r<rep>`, `GREMA_SLOT_DIV=d` (`d` = 1: 4 slots, 2: 2 slots; the default 1 slot is `irma_fixed` of `final_results.tgz`) | 40 |
| `ssystem_slots` | S-system with the slot rule of HFODE, IRMA, `-f 1`: `div0` = default (1 slot), `div2` (2 slots), `div1` (4 slots) | 60 |
| `ssystem_old` | S-system with the old rule `round(N/3)` (2 slots on IRMA), `-f 1`, `ss_f1_<on|off>_r<rep>` | 20 |
| `orient` | simulated data with the IRMA topology (`sim<seed>.txt` from `sim_gen.py <seed> 3 0.02`), `hfode_s<seed>`, `ssystem_s<seed>` (old rule), `-f 1`; evaluate with the IRMA gold standard | 10 |

## `methods_comparison_results.tgz` (dynGENIE3, edge scores, filters before GREMA; see `evaluation/analysis/METHODS_REPORT.md`)

Result files of the comparison, inputs of the controls; GREMA runs as `<setting>/<run>/final_results.txt` (the per-iteration `iga_results/` are not archived), dynGENIE3 as `<data>/seed<k>/final_results.txt`
(same columns: the dynGENIE3 score is written in both score columns, REGULATORY is `+` for every pair). Evaluate with `eval_runs.py` (gold standards: `input/IRMA`, `input/SOS`, Hackett `gold.tsv` in `hackett_results.tgz`).

| Directory | Content |
|-----------|---------|
| `dyngenie` | dynGENIE3 RF / ET, 10 seeds: `RF_on`, `ET_on`, `RF_off`, `ET_off` (IRMA), `sos_RF`, `sos_ET`, `hackett_RF`, `hackett_ET` |
| `dyngenie_perm` | dynGENIE3 on the 10 experiment-permuted data sets of SOS and Hackett (`dyn_<data>_<RF|ET>_p<k>`) |
| `permuted_data` | the permuted data sets (`permute_replicates.py <input> <output> <seed>`, seeds 1-10) |
| `perm_grema` | GREMA on permuted data: `grema_sos_p1-5`, `grema_hackett_p1-3` (2 slots) |
| `contrast_runs`, `contrast_scores` | GREMA runs that kept every iteration (SOS, IRMA on / off, Hackett with 3 and 2 slots) and the contrast scores of them |
| `nmin_runs` | Hackett, 2 slots: `-nmax 8`, `-min 0.5`, `-min 1` (3 runs each) |
| `prior_lasso`, `prior_hsic`, `prior_hill`, `prior_hill_perm` | prior scores, filtered knowledge files and summary of the three filters (`prior_hill_perm`: Hill dictionary on the permuted data) |
| `lasso_grema_runs`, `hsic_grema_runs`, `hill_grema_runs` | GREMA on the filtered knowledge files |

## `hackett_set2_results.tgz` (second TF subset of the Hackett data, see `evaluation/analysis/METHODS_REPORT.md` section 7)

The TFs ranked 11-20 (`hackett_prep.py ... --skip-tfs 10`): `inputs/` (`expression.txt`, `knowledge.txt`, `gold.tsv`, `names.tsv`, `summary.txt`), `grema/` (`grema_2slots_r1-3`: 2 slots on the full candidate list; `grema_hill_r1-3`: on the Hill dictionary filtered knowledge file `hill/s0/knowledge_filtered.txt`),
`dyngenie/RF`, `dyngenie/ET` (10 seeds), `hill/s0-2` (prior scores of the Hill dictionary), `perm/` (the permuted data sets `data_p<k>.txt` and the Hill dictionary / dynGENIE3 scores on them).

## `real_time_points_results.tgz` (option `-time`, see `evaluation/analysis/METHODS_REPORT.md` section 8)

`inputs/set1`, `inputs/set2`: the two Hackett designs with the 8 measured time points (`hackett_prep.py ... --step 0`, `time_list.txt` = the argument of `-time`, `-tunit 5` was used); `rawtime_runs/set1_r1-3`, `set2_r1-3`: GREMA, 2 slots (`GREMA_SLOT_DIV=5`), `-f 1`, `-g 500`, `-time <time_list.txt> -tunit 5`;
`equal_step_test/`: IRMA switch-on, 5 runs without and 5 runs with `-time 0,1,...,15` (the equal-step check).

## `sign_simulation_results.tgz` (signs of GREMA on simulated signed networks, see `evaluation/analysis/SIGN_SIMULATION_REPORT.md`)

`net1-5/` (`expression.txt`, `knowledge.txt`, `gold.tsv`, `signs.tsv`; made with `sim_signed_net.py <dir> <seed>`), `runs/net<k>_r<1,2>/final_results.txt` (GREMA, default settings, `-f 1`, `-g 500`),
`var_<noise10|noise20|points8|points5|series1|additive|hidden2>/net1-5` and `var_*/runs/net<k>_r1` (the misspecified variants, 1 run per network), `dyn/<RF|ET>_net<k>/seed0-4/final_results.txt` (dynGENIE3 on the baseline networks).
`python3 evaluation/analysis/scripts/eval_sim_variants.py <the unpacked archive>` prints the table of the report.

## `sim_protein_results.tgz` (protein regulators, promoter activity measured, shared driver, filtered regulator; see `evaluation/analysis/PROTEIN_REGULATOR_REPORT.md`)

`net1-5/` (no driver), `drv_transcription/net1-5`, `drv_cleavage/net1-5` (`expression.txt`, `knowledge.txt`, `gold.tsv`, `signs.tsv`, `decay.tsv`; `sim_protein_net.py <dir> <seed> [--driver transcription|cleavage]`),
`aug_d0.1|0.2|0.3/net1-5` (data with the filtered copies of the genes and the knowledge file in which only they may regulate; `augment_filtered.py`), `runs/base_net<k>` and `runs/d<D>_net<k>` (GREMA without a driver; the latter on the `aug_d<D>` data),
`drv_*/runs/net<k>` (GREMA with a driver), `dyn/<none|transcription|cleavage>_<RF|ET>_net<k>/seed0-4` (dynGENIE3). `python3 evaluation/analysis/scripts/eval_protein_sims.py <the unpacked archive>` prints the tables of the report.

## `irma_sign_factors_results.tgz` (why the signs on real IRMA are wrong, see `evaluation/analysis/SIGNS_REAL_DATA_REPORT.md` section 2)

`irma_topology_hfode/`: data simulated with HFODE on the IRMA topology (`data_r<1|5>_s<k>.txt`, `sim_gen.py <seed> <replicates> 0.05`), `signs.tsv` (8 signed edges), `runs/r<replicates>_s<k>/final_results.txt` (GREMA default settings);
`irma_like/p<P>i<I>t<T>/seed<k>/` (data, knowledge, gold and signs of the IRMA-like simulation, `sim_irma_like.py <dir> <seed> --protein P --input I --posttrans T`) and `irma_like/runs/p<P>i<I>t<T>_s<k>/final_results.txt`.
`python3 evaluation/analysis/scripts/eval_irma_like.py <the unpacked archive>/irma_like` prints the table of the report; `eval_signs_real.py <signs.tsv> label=<glob>` gives the signs of the first experiment; `linear_signs_irma.py` the linear model on the real data.

## `hmsc_feasibility_results.tgz` (hMSC time series, see `evaluation/analysis/HMSC_FEASIBILITY_REPORT.md`)

`hmsc/runx2/`, `hmsc/multi150/`: `expression.txt` (log2 TPM scaled to [0.01, 1], 5 days), `knowledge_atac.txt` (footprint prior), `knowledge_all.txt` (multi150), `gold.tsv` (TF GENE gold atac_candidate sign), `gold3.tsv` (3 columns for `eval_runs.py`), `names.tsv`, `summary.txt`;
`hmsc/multi150/runs/<atac_r1-5|all_r1-3|shuf1-3_atac>/out/final_results.txt`, `shuf/expression_shuf*.txt` (time-shuffled inputs), `dyn_rf/`, `dyn_rf_atac/` (dynGENIE3, 5 seeds).
`python3 evaluation/analysis/scripts/eval_runs.py --gold hmsc/multi150/gold3.tsv "all=hmsc/multi150/runs/all_r*/out/final_results.txt" ...` gives the table of the report; `hmsc_baselines.py hmsc/multi150`, `hmsc_time_control.py` and `hmsc_diagnose.py` the other numbers.
The raw RNA-seq and ATAC-seq data are internal and not in the repository.

