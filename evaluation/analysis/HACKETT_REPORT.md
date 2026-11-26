# TF-induction time courses of Hackett et al. (2020) with the gold standard of Ma et al. (2014)

A test of GREMA (HFODE, `-f 1`, `-g 500`) on time courses of a larger system with a gold standard of direct regulations.
The final_results.txt of all runs and the input files are in `evaluation/analysis/results/hackett_results.tgz` (see `evaluation/analysis/results/README.md`).

## Data

- **Expression**: GEO [GSE142864](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142864) (Hackett et al., Mol Syst Biol 2020): one TF of *S. cerevisiae* is induced with
  beta-estradiol and the transcriptome is followed. **It is a microarray (Agilent 2-colour) data set, not RNA-seq.** Platform GPL16244: 1548 samples, 182 induction experiments,
  163 of them with the same 8 time points (0, 5, 10, 15, 20, 30, 45, 90 min), 6256 probes; values are log ratios. The processed series matrices were used
  (`GSE142864-GPL16244_series_matrix.txt.gz`, annotation of GPL16244 for the probe -> ORF map); the 3.2 GB raw archive was not downloaded.
- **Gold standard**: Ma, Kemmeren, Gresham, Statnikov, PLoS ONE 2014 ([doi:10.1371/journal.pone.0106479](https://doi.org/10.1371/journal.pone.0106479)), supplementary file
  `pone.0106479.s020.xlsx` (sheets `Gold_standard_network1/2/3`: direct regulatory targets, 1083 / 1785 / 2403 edges, 72 TFs in network 1; network 1 is the most stringent, 3 the most inclusive).
  The files were taken from the journal site (the PMC copy is behind a reCAPTCHA page).
- No sign is stored in the gold standard; signs are not evaluated.

## Design (rules fixed before any GREMA run; `evaluation/analysis/scripts/hackett_prep.py`)

- TFs: the induced TFs that are regulators in the gold standard, with the most responding genes (max |value| over time > 1): SFP1, GLN3, STE12, YOX1, PHO4, MOT3, YAP1, SUM1, ACE2, RLM1
  (`hackett_explore.py` lists the 49 candidates).
- Genes: for every chosen TF the genes that respond most strongly in its own experiment (`--per-tf N`).
- Every chosen experiment is one replicate (10 replicates); linear interpolation to an equal step (`--step`), one global linear scaling to [0.01, 1].
- Knowledge file: only the chosen TFs may regulate (all other rows `0`), every TF -> gene pair `?`; only those pairs are evaluated (self pairs skipped).
- Gold edges among the evaluated pairs, by design (this count was used to choose the design, no GREMA result was):

| Design | variables | pairs | gold edges (network 1) | (network 3) |
|---|---|---|---|---|
| 30 genes with the strongest mean response | 40 | 390 | 2 | |
| 6 genes per TF | 34 | 330 | 6 | |
| 10 genes per TF | 56 | 550 | 9 | **18** (network 2: 12) |
| 15 genes per TF | 86 | 850 | | 23 |
| 20 genes per TF | 120 | 1190 | 16 | 32 |

The first design (strongest shared responders) has too few gold edges to evaluate; network 3 and 10 genes per TF (56 variables, 18 edges) were chosen for the main test because of the number of edges.
Network 3 is inclusive, so many true regulations are missing from the answer and FP are over-estimated.

## Results

Mean AUROC / AUPRC over the runs; random AUROC 0.5.

**Small network** (4 TFs: SFP1, GLN3, STE12, YOX1; 8 genes per TF, 23 variables, 88 pairs, 4 gold edges of network 1; 3 runs; only to check the pipeline and the time step):

| Time grid | binary | conf | vote | time per run |
|---|---|---|---|---|
| 5 min (19 points, 11 interpolated) | 0.562 / 0.052 | 0.631 / 0.061 | 0.758 / 0.091 | about 21 min |
| 15 min (7 points, 2 interpolated) | 0.575 / 0.053 | 0.619 / 0.059 | 0.734 / 0.082 | about 12 min |

All 4 TFs are candidates for every gene (4 slots, no limit). The same network with 2 slots per gene (`GREMA_SLOT_DIV=2`, `round(4/2)`; 3 runs each, `runs_slot2`):

| Time grid, 2 slots | TP / FP | binary | conf | vote |
|---|---|---|---|---|
| 5 min | 4.0 / 42.0 | 0.750 / 0.087 | 0.766 / 0.092 | 0.812 / 0.125 |
| 15 min | 4.0 / 42.0 | 0.750 / 0.087 | 0.760 / 0.090 | 0.859 / 0.145 |
| (4 slots, 5 min / 15 min) | 4.0 / 73.7 and 71.3 | | | |

- With 4 slots the two grids do not differ beyond noise; with 2 slots the 15 min grid is ahead in vote (0.859 vs 0.812) and equal in conf. The difference (0.05) is at the border of noise, and a single true edge of the 4 moves the AUROC by about 0.08 (its rank was 39 on the 5 min grid and 12 on the 15 min grid for one edge), so this is no evidence for a grid. The 15 min grid has fewer made-up points and is faster, so it was used for the main test.
- 2 slots instead of 4: all 4 true edges are still found (TP 4.0), FP fall from about 72 to 42, and binary / conf AUROC rise by 0.13-0.19, vote by 0.05-0.13. The same direction as in the main test below.
  The 2-slot run was made after the 56-variable results were seen, so it is not an independent confirmation.
- 4 true edges cannot show accuracy.

**Main test** (10 TFs, 56 variables, 10 replicates, 15 min grid, network 3: 550 pairs, 18 gold edges, random AUPRC 0.033; 3 runs, about 1.5 h per run with 10 slots):

| regulator slots per gene | TP / FP | binary | conf | vote | SEN / SPE |
|---|---|---|---|---|---|
| 10 (every TF has a slot) | 13.7 / 267.3 | 0.628 / 0.045 | 0.631 / 0.046 | 0.671 / 0.058 | 0.76 / 0.50 |
| 3 (`GREMA_SLOT_DIV=3`) | 10.3 / 153.7 | 0.643 / 0.050 | 0.644 / 0.051 | 0.701 / 0.080 | 0.57 / 0.71 |
| 2 (`GREMA_SLOT_DIV=5`) | 9.0 / 98.3 | 0.658 / 0.058 | 0.661 / 0.062 | 0.720 / 0.101 | |

- With the candidate list restricted to the 10 TFs, GREMA gave every candidate its own slot (10 per gene) and called 86% of the pairs a regulation (475 of 550; the true share is 3%).
  Fewer slots per gene remove FP faster than TP (FP -63%, TP -34% from 10 to 2 slots) and the ranking improves steadily: vote AUROC 0.671 -> 0.701 -> 0.720, AUPRC 0.058 -> 0.080 -> 0.101 (random 0.033).
  The runs of one setting agree closely (sd of the AUROC 0.000-0.02), but only 3 runs and 18 edges were made, and the 10 -> 2 slot difference (0.05) is at the border of noise.
- Even the best setting is weak (vote AUROC 0.72), below the best IRMA and SOS results (about 0.75-0.78; IRMA switch-off is 0.6-0.7).
- The direction is the same as on SOS (fewer slots, fewer FP), but on SOS 2 slots lost many TP and did not improve the AUROC; here 2 slots were better than 3. The reason was not tested
  (the density of true edges is low here, about 0.3 TFs per gene among the 10 candidates).
- Slot 1 was not tried; the signs of the TF -> gene edges are not evaluated.

## Code change: `GREMA_SLOT_DIV` also for restricted candidate lists (`EMA_HFODE/src/Getknowledge.cpp`)

Before, a restricted candidate list (some columns `0`, so the number of candidates is not the number of variables) always got one slot per candidate. With `GREMA_SLOT_DIV` set, such a list
is searched as a random subset of `round(candidates / GREMA_SLOT_DIV)` slots (the subset mode that the unknown-regulator case already uses). Without the variable nothing changes.

## Reproduce

```shell
# data (GEO, journal site); the xlsx is read with the Python standard library only (no external package)
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142864/matrix/GSE142864-GPL16244_series_matrix.txt.gz   # save as GPL16244.txt.gz
curl -o GPL16244_annot.txt "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16244&targ=self&form=text&view=full"
curl -o gold.xlsx "https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0106479.s020"
python3 -I evaluation/analysis/scripts/hackett_prep.py DATA_DIR gold.xlsx OUT --per-tf 10 --network 3 --step 15
GREMA_ROOT=$PWD GREMA_SLOT_DIV=5 evaluation/analysis/scripts/run_one.sh RUNDIR OUT/expression.txt OUT/knowledge.txt -g 500 -f 1 -t 3
python3 evaluation/analysis/scripts/eval_runs.py --gold OUT/gold.tsv "label=RUNDIR*/out/final_results.txt"
```

## Limits

- 3 runs per setting and 18 gold edges: only large differences can be seen. One data set, one design.
- Log ratios of a 2-colour array (not counts) were scaled to [0.01, 1]; 5 of the 7 time points are measured, the others interpolated.
- Each replicate is the induction of a different TF, not a repeat of the same experiment (like the 4 experiments of the SOS data, but with different perturbations); GREMA's model has no notion of this.
- The gold standard is inclusive (network 3); stricter networks were not run with GREMA (only counted).
- This is not a test on bulk RNA-seq; a real RNA-seq series with a gold standard (candidate: GSE97500, Arabidopsis nitrate response, 10 time points x 3 replicates) was not run.
