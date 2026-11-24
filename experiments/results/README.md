# Archived results

`final_results.tgz` holds the `final_results.txt` of every finished run (columns `TF GENE REGULATORY
CONFIDENCE_LEVEL VOTE_SCORE`), nothing else (the full run directories are 10-30 MB each). Evaluate with
`experiments/scripts/eval_runs.py`, e.g.

```shell
mkdir -p /tmp/res && tar xzf experiments/results/final_results.tgz -C /tmp/res
python3 experiments/scripts/eval_runs.py --gold input/IRMA/IRMA_GoldStandard.tsv \
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
`lat/...`); they were made with `experiments/scripts/gp_auto.py` and `gp_v2.py` (seeds 0-4, 10 trajectories each).
Use the same input files if you want to re-run with different options.

## `sos_results.tgz` (SOS follow-up analyses, see `experiments/SOS_REPORT.md`)

`final_results.txt` of every run, one directory per setting, `<setting>/<run>/final_results.txt`; all runs finished.

| Directory | Content | Runs |
|-----------|---------|------|
| `sos4` | SOS, 4 experiments as replicates, diagonal `?`, `-g 500`: `f1_r1-5`, `f3_r1-5` (T1; gold `input/SOS/SOS_GoldStandard.tsv`, `--self-pairs`) | 10 |
| `sos4_g1000` | as `sos4` with `-g 1000` | 10 |
| `sos4_selfrep` | as `sos4` `-f 1`, knowledge file `know_G2G2_minus.txt` (LexA -> LexA fixed to `-`; diagnostic only) | 5 |
| `sos4_slot1`, `sos4_slot4` | as `sos4` `-f 1` with `GREMA_SLOT_DIV=1` (8 regulator slots) and `=4` (2 slots) | 5 + 5 |
| `sim_sos` | simulated SOS data, hidden LexA protein (`sim<seed>.txt` made by `sim_sos.py <seed>`), `-f 1`, runs `f1_s1-5`, evaluate with the SOS gold standard | 5 |
| `sim_sos_p` | same simulation with the protein measured as G9 (`sim_sos_p.py`), 9 x 9 knowledge `know9.txt`, gold `gold9.tsv` (G9 -> G1..G8), runs `f1_s1-5`, `f3_s1-5` | 10 |
