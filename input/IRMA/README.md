# IRMA (yeast synthetic network) in the GREMA format

**Source.** Cantone et al. (2009), *A Yeast Synthetic Network for In Vivo Assessment of Reverse-Engineering and Modeling Approaches*, Cell 137:172-181. The data are the copy published by Grzegorczyk and Husmeier on Figshare (<https://figshare.com/articles/dataset/Data_for_paper_submitted_to_BMC_Bioinformatics/5363545>, file `BMC_data.zip`, folder `YEAST/`), licensed CC BY 4.0.

**Genes.** G1 CBF1, G2 GAL4, G3 SWI5, G4 GAL80, G5 ASH1.

**True network (8 regulations).** CBF1->GAL4, GAL4->GAL80, GAL4->SWI5, SWI5->ASH1, SWI5->CBF1, SWI5->GAL80, GAL80->GAL4, ASH1->CBF1.

| File | Content |
|------|---------|
| `IRMA_switch_on_timeseries_expression.txt` | glucose -> galactose, 1 replicate, 16 time points |
| `IRMA_switch_off_timeseries_expression.txt` | galactose -> glucose, 1 replicate, 21 time points |
| `IRMA_no_prior_knowledge.txt` | knowledge file (all `?`, self-regulation `0`) |
| `IRMA_GoldStandard.tsv`, `IRMA_GoldStandardUndirected.tsv` | gold standards for `evaluation/dream4_evaluate.py` |

**Preprocessing** (`convert_irma.py`): both data sets are scaled with the same linear map to [0.01, 1] (GREMA cannot use 0). All time points are kept; the two experiments are separate inputs because they have different numbers of time points. GREMA treats the time points as unit time steps. According to the literature the sampling interval is constant within each experiment (about 20 min for switch-on and 10 min for switch-off; not verified against the original paper), and each experiment is analysed separately, so this only changes the time scale of the fitted rates.

The first time point (before the shift, the washing period) is kept. Removing it, as done in some other studies, made the results worse (5 runs, HFODE, `-g 500`):

| Data | First point | TP | FP | AUROC | AUPRC |
|------|-------------|----|----|-------|-------|
| switch-on | kept | 4.0 | 1.0 | 0.708 | 0.600 |
| switch-on | removed | 1.6 | 3.4 | 0.460 | 0.393 |
| switch-off | kept | 4.0 | 2.0 | 0.677 | 0.565 |
| switch-off | removed | 3.4 | 2.0 | 0.640 | 0.537 |

AUROC and AUPRC are computed over the 25 possible directed edges (8 true), scoring an edge with its `CONFIDENCE_LEVEL` and a non-edge with 0 (random: AUROC 0.5, AUPRC 0.4).

```shell
mkdir -p output/IRMA_switch_on
perl GREMA_main.pl -i input/IRMA/IRMA_switch_on_timeseries_expression.txt \
    -k input/IRMA/IRMA_no_prior_knowledge.txt -o output/IRMA_switch_on -g 500 -t 5
python3 evaluation/dream4_evaluate.py -d input/IRMA/IRMA_GoldStandard.tsv \
    output/IRMA_switch_on/final_results.txt
```

**Reference results** (HFODE, 5 repeated runs per setting; 8 true edges, 25 possible directed edges). Results are essentially identical for `-g 500` and `-g 1000`, while `-g 100` is noisier, so `-g 500` is a good choice for this data set.

| Data | `-g` | TP | FP | Time |
|------|------|----|----|------|
| switch-off | 100 | 3.4 +- 0.8 | 2.4 +- 1.0 | 6 s |
| switch-off | 500 | 4.0 +- 0.0 | 2.0 +- 0.0 | 24 s |
| switch-off | 1000 | 4.0 +- 0.0 | 1.6 +- 0.5 | 36 s |
| switch-on | 100 | 4.4 +- 0.5 | 1.0 +- 0.0 | 4 s |
| switch-on | 500 | 4.0 +- 0.0 | 1.0 +- 0.0 | 13 s |
| switch-on | 1000 | 4.0 +- 0.0 | 1.0 +- 0.0 | 24 s |

Each run already performs the 30-run voting internally, so repeating the whole run is not required. About half of the true edges are missed, and the missed and false edges are the same in every run, so the limit is the information in the data (1 replicate, 16-21 time points), not randomness.

## Experiments on IRMA

All numbers below are from HFODE with `-g 500` unless stated otherwise, evaluated against `IRMA_GoldStandard.tsv` (8 true edges among 25 possible directed edges; random AUROC 0.5, AUPRC 0.4). The scripts used for the experiments below (simulation, augmentation, evaluation) are not included in this repository. There are only 2 data sets and 8 true edges, so differences of about 0.05 in AUROC should not be over-interpreted.

### HFODE vs S-system

S-system (`-m s-system`, 3 runs) was clearly worse and 5-8 times slower (about 100-120 s):

| Data | Model | TP | FP | AUROC | AUPRC |
|------|-------|----|----|-------|-------|
| switch-on | HFODE | 4 | 1 | 0.708 | 0.600 |
| switch-on | S-system | 1-2 | 6-7 | 0.307 | 0.365 |
| switch-off | HFODE | 4 | 2 | 0.677 | 0.565 |
| switch-off | S-system | 4 | 4-5 | 0.568 | 0.451 |

For `-g 100/500/1000` HFODE had AUROC 0.74/0.71/0.71 and AUPRC 0.66/0.60/0.60 (switch-on), and AUROC 0.62/0.68/0.69 and AUPRC 0.52/0.57/0.58 (switch-off). Scores are `CONFIDENCE_LEVEL` for an edge and 0 for a non-edge; using the binary prediction instead changes AUROC by at most 0.04.

### Is the limited accuracy caused by the model or by the amount of data?

Data were simulated with HFODE itself on the IRMA topology (8 edges and signs, 5 genes, 16 time points, 5% noise, random initial values for each replicate, 3 random parameter sets), so the model is exactly right:

| Replicates | TP | FP | AUROC | AUPRC |
|------------|----|----|-------|-------|
| 1 | 3.0 | 2.7 | 0.58 | 0.51 |
| 5 | 4.7 | 0.3 | 0.78 | 0.71 |
| 10 | 5.0 | 0.3 | 0.81 | 0.77 |

Even with a correct model and no noise, one replicate gives about the same accuracy as on the real IRMA data. Independent replicates (different initial states) help the most; going from 5 to 10 gives little extra, and about 3 of the 8 true edges were still missed. So the main limit is the amount of data, not the model. The noise level had almost no effect (1 replicate without noise: TP 3.0, FP 2.7, AUROC 0.57).

### Can the data be augmented?

All experiments below keep switch-on and switch-off separate. An augmented data set is the original series plus extra series, or extra series only, and GREMA is run on each.

- **Simulating extra replicates from a fitted HFODE** (parameters taken from a GREMA run on the same data) made the results **worse** (switch-on AUROC 0.47, switch-off 0.48-0.54; 4 or 9 extra replicates). The simulated data only repeat what the fitted model already says, including its false edges.
- **Gaussian-process (GP) sampling and voting.** Each gene is fitted with a GP, 20 trajectories are sampled (same number of time points, 1 replicate each), GREMA is run on each trajectory, and every edge is scored by the fraction of the 20 runs in which it was predicted. Five random seeds per setting (mean over seeds; the scores below use the vote fraction, original = binary prediction of a single run):

| Method | switch-on AUROC | switch-on AUPRC | switch-off AUROC | switch-off AUPRC |
|--------|-----------------|-----------------|------------------|------------------|
| Original data, single run | 0.708 | 0.600 | 0.667 | 0.533 |
| GP v1: kernel, length scale and noise chosen by leave-one-out (LOO) prediction | 0.752 | 0.703 | 0.545 | 0.503 |
| GP v2, latent curve only (hyperparameters sampled from a grid weighted by marginal likelihood) | 0.782 | 0.774 | 0.565 | 0.540 |
| GP v2, latent curve + observation noise | 0.786 | 0.740 | 0.575 | 0.541 |

GP voting improved switch-on (AUROC +0.04 to +0.08; the true edges G5->G1, G1->G2, G3->G5 got the most votes) but not switch-off, where the most voted edge was the false positive G3->G2. The seed-to-seed range was small for v1 (switch-on 0.745-0.771) and larger for v2 with noise (switch-on 0.69-0.92). Without using the gold standard, the agreement between the switch-on and switch-off vote scores (Spearman 0.22-0.40) was not clearly higher than for the original data (0.38), so we cannot say that augmentation makes the inference more reliable. The GP settings were fixed before running and the gold standard was not used to choose them, but the first GP settings (a fixed length-scale/noise range, not data-driven) gave switch-on AUROC 0.82 in a single run and switch-off 0.48; we did not repeat that setting with several seeds, and relaxing its range afterwards (switch-off 0.64) was done after seeing the scores, so that number is not a fair result.

GP samples contain no new information about the network, only the uncertainty of the existing series, so they cannot replace additional experiments. The effect on real data is inconsistent across the two experiments.

### After the code-review fixes

A review of the code found and fixed several bugs (see the pull request): the derivative term of fitness 2/3/4 read an uninitialised value at the last time point, regulations given as `+`/`-` by the knowledge file were not always kept in the model, runs with equal fitness overwrote each other in the voting, and others. The numbers below were obtained after the fixes (HFODE, `-g 500`, `-c 0`, 10 repeated runs of the original data; mean over runs, AUROC / AUPRC over the possible directed edges, random AUROC 0.5 and AUPRC 0.4). The earlier tables in this README were produced before the fixes; the default fitness (`-f 1`) gives about the same results as before.

| Data | Fitness | `CONFIDENCE_LEVEL` | `VOTE_SCORE` |
|------|---------|--------------------|--------------|
| switch-on | `-f 1` (default) | 0.705 / 0.596 | 0.781 / 0.754 |
| switch-on | `-f 3` | 0.779 / 0.732 | 0.762 / 0.667 |
| switch-off | `-f 1` (default) | 0.683 / 0.572 | 0.605 / 0.544 |
| switch-off | `-f 3` | 0.708 / 0.660 | 0.609 / 0.574 |

- `-f 3` (squared error of the derivative + squared error) was better than `-f 1` for switch-on (the standard deviation of the AUROC over runs is 0.04 for `-f 3`) and slightly better for switch-off (the difference is as large as the standard deviation). Before the fix `-f 3` was at the random level.
- Adding the correlation term (`-c 1`, 5 runs) did not help: switch-on AUROC 0.646 (`-f 1`) and 0.675 (`-f 3`), switch-off 0.663 and 0.708, measured with `CONFIDENCE_LEVEL`.
- GP augmentation (5 seeds x 10 trajectories, score = mean over the 50 trajectories) helped switch-on only for the v2 latent curve with `-f 1` and `VOTE_SCORE` (0.823 / 0.834); the first GP version and `-f 3` were not better than the original data. For switch-off every GP setting was worse than the original data (AUROC 0.50-0.60 against 0.68-0.71). The seed-to-seed standard deviation of the AUROC is 0.04-0.09 and there are only 8 true edges, so these differences are uncertain.
