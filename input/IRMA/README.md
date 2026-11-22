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
