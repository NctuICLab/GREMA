# Can GREMA build a GRN from the hMSC time-series RNA-seq and ATAC-seq? A feasibility trial

Data: hMSC day 1, 4, 7, 10, 15 (internal): RNA-seq (nf-core, STAR + RSEM, one sample per day) and ATAC-seq (peaks, HINT-ATAC footprints with JASPAR motifs, per day).
Archive of inputs and results: `evaluation/analysis/results/hmsc_feasibility_results.tgz` (see `evaluation/analysis/results/README.md`). Scripts: `evaluation/analysis/scripts/hmsc_*.py`.
HFODE, `-f 1`, `-g 500`, real times (`-time 1,4,7,10,15 -tunit 3`), `GREMA_SLOT_DIV=5`, 4 threads.

**Answer: not with these data.** The data have 1 series of 5 time points (4 samples); no method that was tried recovers the regulations that the gold standard lists, except a plain correlation that mostly measures co-expression.
This is a limit of the data size and of the gold standard (11 true pairs, 7 of them of one TF), not a measured failure of GREMA in general; section 5 lists what would change it.

## 1. The data and the check (`check_data.py`)

- 78,932 genes; 12,113 with TPM > 5 on some day; 4,093 of them change more than 2-fold. 5 time points centred have at most 4 independent dimensions; the first two components of the changing genes hold 83% of the variance.
- `check_data.py` on the two designs below: only 1 series, 5 time points, 4 samples, mean |correlation| between the genes 0.65 (RUNX2 design) and 0.75 (multi design): all four warnings apply.
- ATAC: the footprint files have about 1.1-1.2 million footprints per day (days 1 and 7 counted; 579 JASPAR motifs on day 1) (one truncated line in the day-7 file is skipped); the differential analysis of the TF activities exists only against day 1 and the activity of RUNX2 is not significant (p 0.57, 0.57, 0.82 for day 4, 10, 15).
- `merged_peaks_annotated_header.txt` has only the header line, so the peaks were linked to genes by the distance to the TSS (GENCODE v35, protein coding): a footprint within +-10 kb.

## 2. Gold standard and prior (rules fixed before the GREMA runs; `hmsc_prep.py`)

- **Gold standard**: TRRUST v2 (human, 8,427 TF -> gene pairs, from the literature; downloaded, `external/trrust` of the run folder, not in the repository; the sign is Activation / Repression / Unknown). It is not specific to hMSC.
  There is no gold standard from hMSC data: the ChIP-seq of hMSC in `ifar/ChIP-seq` are histone marks, not TFs.
- **Design `runx2`** (the request): RUNX2 is the only regulator; candidates = genes with a RUNX2 footprint near the TSS and fold >= 2, the 60 with the largest fold. The gold standard has **1** true pair among the 60: it cannot be evaluated.
- **Design `multi`**: regulators = the 15 TRRUST TFs (TFs of the gold standard) with the largest fold that are expressed and have a motif; genes = the other genes with the largest fold (>= 2). With 60 genes: 6 true pairs; with 150 genes (decided after seeing the 6 and before any GREMA run, to have at least about 10): **2,460 candidate pairs, 11 true pairs, 1 of them with a footprint**; 5 of the 11 have a sign.
  Of the 11 true pairs, **7 are E2F1 -> MYBL2, DDIT3, TOP2A, AURKB, RRM2, CDK1, CCNB1** (mostly cell-cycle genes), the others EGR2 -> IGF2, SREBF1 -> SCD, DDIT3 -> TRIB3 and DDIT3 -> CDK2: effectively 4 regulators.
- Data: log2(TPM + 1), one global scaling to [0.01, 1]. Knowledge `atac`: TF -> gene is `?` only when the TF has a footprint near the gene (275 of 2,460 pairs); knowledge `all`: every TF -> gene pair is `?`.

## 3. Is the footprint prior informative? (`hmsc_prior_enrichment.py`, no GREMA)

204 TRRUST TFs with a motif and expression x 10,853 expressed genes = 2.2 million pairs, 2,690 of them in TRRUST: **32% of the TRRUST pairs have a footprint near the gene against 15% of all pairs** (odds ratio 2.7);
per-TF AUROC of the footprint 0.55 (106 TFs with >= 5 targets), 0.500 when the genes are shuffled within each TF (20 shuffles, range 0.491-0.508). A real but weak signal; not corrected for the accessibility and the expression level of the promoter, which also raise the number of footprints and the number of literature reports.
In the 15 TFs x 150 genes design the footprint is no use: 1 of the 11 true pairs has one (AUROC 0.49).

## 4. Results on the multi design (2,460 pairs, 11 true; random AUROC 0.50, sd about 0.09; `eval_runs.py`, `hmsc_baselines.py`, `hmsc_time_control.py`, `hmsc_diagnose.py`)

| Method | runs | AUROC | note |
|---|---|---|---|
| random | | 0.50 +- 0.09 | |
| footprint prior alone (0/1) | | 0.49 | |
| \|Pearson\| of the two 5-point profiles | | **0.83** | 8 of the 11 true pairs rank 5-196 of 2,460 (6 of them E2F1 targets) |
| \|Spearman\| | | 0.72 | |
| \|Pearson\| of the TF level and the change of the gene | | 0.38 | the time order does not help |
| dynGENIE3 (RF, score) | 5 seeds | 0.36 | same with the `atac` knowledge (it only selects regulators, not pairs); sd over seeds 0.017 |
| **GREMA, all pairs candidates** | 3 | **0.48 (binary), 0.48 (conf), 0.55 (vote)**; sd 0.02-0.03 | 1.3 true pairs called of 11, 382 false; all 4 signs of the called true pairs were `-`; the pair called by most runs, DDIT3 -> TRIB3, is an activation in TRRUST: wrong |
| GREMA, `atac` prior | 5 | 0.47, 0.47, 0.51 | 0 of 11 called: the prior removes 10 of the 11 true pairs; 127 pairs called (one regulator for every gene that has a candidate) |
| GREMA, `atac` prior, time-shuffled data | 3 | 0.47, 0.47, 0.47 | same as with the real data |

The ranks of the 11 true pairs by the GREMA vote are 153-2,294 of 2,460 (median about 1,000): random. Eight of the 11 true pairs are in the best 200 by the correlation; by the GREMA vote the same eight are at ranks 153-1,842.

**Does GREMA use the time order?** Yes: the vote scores of the 275 `atac` candidates correlate 0.97 between runs on the real data, 0.41 between a run on the real and a run on the time-shuffled data, 0.44 between two shufflings.
The GREMA vote does not follow the correlation either (Spearman with \|Pearson\| -0.04, 2,460 pairs, `all` runs; reproducibility between the `all` runs 0.69). The number of pairs called (127 = the genes with a candidate, `atac`) is set by the slots, not by the data.
So GREMA scores depend on the dynamics, but for these data they do not point to the pairs of the gold standard.

## 5. What this says, and what to change

**What the trial shows.** (a) With 1 series of 5 points `check_data.py` was right: nothing that uses the dynamics recovers the gold pairs; the only signal is co-expression of cell-cycle genes (E2F1 and its targets), which a correlation finds and the dynamics-based scores do not.
(b) The footprint prior has a real but weak association with TRRUST (odds ratio 2.7). (c) As a hard filter (`0` without a footprint) it removes most of the true pairs; it cannot be used in this way.
(d) No hMSC-specific gold standard exists; TRRUST has 11 pairs here, 7 of one TF, so differences of AUROC smaller than about 0.2 cannot be read.

**Changes that are likely to matter, in the order of the expected effect** (none of them is tested here except where said):

1. **More time points and replicates (data).** From the simulations of this repository, one series is weaker than four (AUROC 0.82 against 0.91, simulations with 20 points) and the data check asks for at least 8 points and 30 samples. The effect of 5 against 8-12 points on HFODE data was not simulated; it is the cheapest test to run before the experiment is planned (a simulation with `sim_gen.py`).
2. **The ATAC prior as a soft prior, not a filter** (needs a change of GREMA): the random slot subsets of the genetic algorithm are drawn with a higher probability for the candidates with a footprint, and all the other candidates stay possible. The knowledge file would need a third state for "supported" next to `?`.
   The gain is bounded by section 3 (odds ratio 2.7, AUROC 0.55 per TF), so it improves the search space, not the ranking, by itself.
3. **A regulator activity from ATAC**: the footprint score (protection, number of footprints in the promoter of the targets) per day as the measured variable of the TF in place of its mRNA. The earlier simulations (`PROTEIN_REGULATOR_REPORT.md`) show that the measured variable should be the activity of the regulator. Here the differential analysis (RUNX2 activity not significant) is too weak and the activity is only given against day 1.
4. **Modules instead of genes**: the changing genes follow 2-3 profiles; a GRN between 5-10 module profiles (clusters of the expressed genes) with HFODE has 5 points for few variables; the edges inside a module are then chosen by the ATAC prior. Needs a clustering step before GREMA and an edge expansion after; no change of the solver.
5. **Parameters shared inside a module** (a change of the solver: one degradation and one threshold for the genes of a module) so that the 4 samples constrain fewer parameters.
6. **An hMSC gold standard**: regulations from TF perturbation or ChIP of TFs in MSC or osteoblast (public data: for example ChIP-Atlas; the CRISPRi of hFOB in `ifar/CRISPRi/Conery2025_hFOB` is a perturbation of enhancers, not of TFs, and was not examined). With the gold standard of TRRUST only, any design needs more than 11 pairs (more TFs, not more genes of the same TF).

## Limits and deviations

- One series, 5 points; 11 true pairs, 7 of E2F1; TRRUST is a literature set from many cell types; every number with 11 true pairs has a large uncertainty (random sd of the AUROC about 0.09).
- Deviations from the plan: the number of genes was raised from 60 to 150 (rule above); the `all` knowledge has 3 runs instead of 5 (`all_r4` was cancelled after 4 minutes and `all_r5` not started, because of the run time, about 2 hours each); the time-shuffled control has 3 shufflings.
- dynGENIE3 ran with equally spaced times (0, 1, 2, 3, 4), not the real days.
- The footprint rule (+-10 kb of the TSS on any of the 5 days) and the genes (TPM >= 5, fold >= 2) were fixed before the runs; the window was not varied. The design choices (15 TFs, the largest fold) selected weak TFs (EGR2, BHLHE41, NPAS2, ...) because the largest fold among the TRRUST TFs was the rule.
- Not examined: the differences between the days in accessibility itself (the peak signal as a time series), the promoters of the genes with the longer distance regulation (enhancers), and the ATAC replicates.
