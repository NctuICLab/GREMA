# Regulator slots on IRMA (HFODE and S-system) and the direction of the S-system output

Archive: `evaluation/analysis/results/ssystem_slots_results.tgz` (see `evaluation/analysis/results/README.md`). IRMA, switch-on (16 points) and switch-off (21), one replicate each,
`IRMA_no_prior_knowledge.txt` (diagonal `0`, so 4 candidate regulators per gene), `-f 1`, `-g 500`, `-t 5`, 10 runs per setting; AUROC / AUPRC averaged over the runs
(random 0.5 / 0.4), gold standard `input/IRMA/IRMA_GoldStandard.tsv` (8 edges); differences below about 0.05 AUROC are noise.

## 1. What the default number of regulator slots is

A slot holds one regulator that is searched for a gene at the same time (a slot can still end up with "no regulation"). Default rule of the HFODE solver (and, since this PR, of the S-system solver):
`round(candidates / 3)` slots, at least 1, when every gene may regulate (all genes, or all but the diagonal); one slot for every candidate when the candidate list is restricted.
`GREMA_SLOT_DIV` replaces the divisor 3 (README, section "Regulator slots").

| Data | candidates per gene | default slots |
|---|---|---|
| IRMA (5 genes, diagonal `0`) | 4 | **1** |
| SOS (8 genes, all `?`) | 8 | 3 |
| DREAM4 10 genes (diagonal `0`) | 9 | 3 |
| Hackett small network (only 4 TFs may regulate) | 4 | 4 (restricted list) |

All earlier IRMA results of HFODE (also those of the other reports) were therefore made with **one slot per gene**.
The S-system solver used `round(N/3)` with `N` the number of genes (2 slots on IRMA) and could not be changed; it now follows the HFODE rule (this PR).

## 2. IRMA: 1, 2 and 4 slots (`GREMA_SLOT_DIV` unset / 2 / 1)

AUROC / AUPRC, conf = CONFIDENCE_LEVEL, vote = VOTE_SCORE; TP / FP of the voted network.

**Switch-on**

| Model, slots | conf | vote | TP / FP |
|---|---|---|---|
| HFODE, 1 (default) | 0.705 / 0.596 | 0.781 / 0.754 | 4.0 / 1.1 |
| HFODE, 2 | 0.698 / 0.546 | 0.828 / 0.824 | 5.9 / 4.1 |
| HFODE, 4 | 0.679 / 0.513 | 0.675 / 0.532 | 7.7 / 8.8 |
| S-system, 1 | 0.386 / 0.377 | 0.484 / 0.417 | 1.0 / 4.3 |
| S-system, 2 | 0.316 / 0.347 | 0.280 / 0.336 | 2.5 / 7.4 |
| S-system, 4 | 0.464 / 0.384 | 0.500 / 0.400 | 8.0 / 12.0 |
| S-system, old rule `round(N/3)` (= 2 slots) | 0.290 / 0.345 | 0.307 / 0.367 | 2.2 / 7.7 |

**Switch-off**

| Model, slots | conf | vote | TP / FP |
|---|---|---|---|
| HFODE, 1 (default) | 0.683 / 0.572 | 0.605 / 0.544 | 4.0 / 1.8 |
| HFODE, 2 | 0.459 / 0.386 | 0.520 / 0.465 | 3.7 / 6.5 |
| HFODE, 4 | 0.419 / 0.369 | 0.503 / 0.441 | 6.5 / 10.7 |
| S-system, 1 | 0.500 / 0.400 | 0.493 / 0.494 | 2.0 / 3.0 |
| S-system, 2 | 0.478 / 0.393 | 0.335 / 0.353 | 3.8 / 6.0 |
| S-system, 4 | 0.426 / 0.372 | 0.500 / 0.400 | 8.0 / 12.0 |
| S-system, old rule (= 2 slots) | 0.468 / 0.387 | 0.349 / 0.356 | 3.8 / 6.2 |

- **HFODE: one slot is not too few on IRMA.** conf is best with 1 slot; more slots find more of the 8 true edges but even more false ones (switch-off: FP 1.8 -> 10.7), and the conf AUROC falls to or below random with 4 slots.
  The vote of switch-on is highest with 2 slots (0.828 vs 0.781), within noise. On SOS (8 / 3 / 2 slots) and the Hackett network (10 / 3 / 2) fewer slots than "every candidate" cut the false positives; IRMA's default is already the smallest setting (1 slot), and on none of the three data sets did more slots than the default give a better ranking.
- **S-system on IRMA does not give a usable ranking**: AUROC 0.25-0.5 (random level or below) for every number of slots; the new rule with 2 slots gives the same as the old `round(N/3)` (both are 2 slots).
- The new S-system slot rule itself behaves as designed (slot numbers 1 / 2 / 4 for the three settings, checked by running the solver directly).

## 3. Is the direction of the S-system output reversed? (it is not)

On IRMA the S-system AUROC was below 0.5, and with the gold standard transposed (TF and GENE swapped) it was above 0.5 (switch-on, 1 slot: conf 0.716, vote 0.793; HFODE with the transposed gold: 0.40 / 0.36).
That pattern looked like a reversed output, so a test with a known truth was made: simulated data with the IRMA topology (`evaluation/analysis/scripts/sim_gen.py`, 3 replicates, noise 0.02, seeds 1-5, HFODE-generated; archive `orient/`), HFODE and S-system, `-f 1`, `-g 500`:

| | gold as is (conf / vote) | gold transposed (conf / vote) |
|---|---|---|
| HFODE | 0.825 / 0.817 | 0.383 / 0.389 |
| S-system (old rule) | 0.702 / 0.677 | 0.498 / 0.524 |

Both models are above random only with the true orientation, so **the S-system labels are not reversed**, and the code agrees: `cal_g[Var][j]` is the effect of gene `j` on the target gene `Var`, as in the knowledge file and in `GREMA_main.pl`.
The below-random AUROC on the real IRMA data is therefore a property of the model on these data (not understood; possible reasons that were not tested: the feedback loops of IRMA, `h_ij` fixed to 0, a single time series).
**The transposed evaluation is a diagnostic only**: choosing the orientation of the output by the gold standard would be using the answer, and on the simulated data the reversed output is worse.

## Limits

- One data set (IRMA), `-f 1` only, 10 runs per setting (the orientation test: 5 seeds); the S-system `-f 3` runs of the first batch (20 runs) were made but not evaluated and are not in the archive.
- 8 true edges. The slot comparison on SOS and the Hackett network (HFODE only) is in `evaluation/analysis/SOS_REPORT.md` and `evaluation/analysis/HACKETT_REPORT.md`.
- S-system was run with the settings of this repository (`h_ij` fixed to 0); other settings were not tried.
